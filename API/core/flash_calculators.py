"""
Flash calculation system for REFPROP API.

This module defines base classes and implementations for different types
of flash calculations (PT, PH, TS, etc.).
"""

from typing import Dict, List, Any, Tuple, Optional, Union, Callable, Iterator
import numpy as np
import logging
from collections import defaultdict

# Import from property system
from API.core.property_system import PropertyRegistry

logger = logging.getLogger(__name__)

class FlashCalculator:
    """
    Base class for all flash calculations.
    
    This class provides a common interface and shared functionality for all
    types of flash calculations. Specific flash types (PT, PH, TS, etc.) should
    implement the abstract methods defined here.
    """
    
    def __init__(self, rp, property_registry):
        """
        Initialize the flash calculator.
        
        Args:
            rp: REFPROP instance
            property_registry: PropertyRegistry instance
        """
        self.rp = rp
        self.property_registry = property_registry
    
    def calculate_flash_grid(self, composition: List[Dict[str, Any]], 
                            variables: Dict[str, Dict[str, Any]], 
                            properties: List[str], 
                            **options) -> Tuple[List[Dict[str, Any]], Dict[str, Any], Dict[str, np.ndarray]]:
        """
        Calculate properties across a grid of points.
        
        Args:
            composition: List of fluid components and fractions
            variables: Dictionary of variables and their ranges
            properties: List of properties to calculate
            options: Additional options (grid_type, etc.)
            
        Returns:
            results: List of calculated property points
            grid_info: Information about the calculation grid
            grids: Dictionary of grid arrays used for calculation
        """
        # Setup mixture
        z = self._setup_mixture(composition)
        molar_mass = self.rp.WMOLdll(z)
        
        # Generate grids based on flash type
        grids = self._generate_grids(z, variables, options)
        
        # Calculate at each grid point
        results = []
        
        # Store progress information
        total_points = self._get_total_grid_points(grids)
        completed = 0
        error_count = 0
        
        logger.info(f"Starting {self.__class__.__name__} calculation with {total_points} grid points")
        
        # Use a try-except block to ensure we return partial results on error
        try:
            for point_idx, grid_point in enumerate(self._grid_iterator(grids)):
                try:
                    # Calculate base properties at this point
                    base_props = self._calculate_base_properties(z, grid_point)
                    
                    # Add molar mass for convenience in property calculations
                    base_props['molar_mass'] = molar_mass
                    
                    # Calculate all requested properties
                    calculated_props = self._calculate_all_properties(
                        base_props, properties, z, molar_mass
                    )
                    
                    # Add grid indices
                    grid_indices = self._get_grid_indices(grid_point, grids)
                    
                    # Add to results
                    results.append({
                        'index': point_idx,
                        **grid_indices,
                        **calculated_props
                    })
                    
                    completed += 1
                    
                    # Log progress every 10% or 100 points
                    if completed % max(1, total_points // 10) == 0 or completed % 100 == 0:
                        logger.info(f"Completed {completed}/{total_points} points ({completed/total_points*100:.1f}%)")
                    
                except Exception as e:
                    # Log error and continue with next point
                    error_count += 1
                    logger.warning(f"Error at point {grid_point}: {str(e)}")
                    continue
            
        except Exception as e:
            # Handle unexpected errors in the iteration itself
            logger.error(f"Error during grid calculation: {str(e)}")
        
        # Log final statistics
        logger.info(f"Calculation complete: {completed}/{total_points} points calculated "
                    f"({completed/total_points*100:.1f}%), {error_count} errors")
        
        # Prepare grid info
        grid_info = self._prepare_grid_info(grids, options, len(results))
        
        return results, grid_info, grids
    
    def _setup_mixture(self, composition: List[Dict[str, Any]]) -> List[float]:
        """
        Setup REFPROP mixture.
        
        Args:
            composition: List of fluid components and fractions
            
        Returns:
            z: Composition array
            
        Raises:
            ValueError: If mixture setup fails
        """
        fluid_string = '|'.join(f"{comp['fluid']}.FLD" for comp in composition)
        z = [comp['fraction'] for comp in composition] + [0] * (20 - len(composition))
        
        ierr, herr = self.rp.SETUPdll(len(composition), fluid_string, 'HMX.BNC', 'DEF')
        if ierr > 0:
            raise ValueError(f"Error setting up mixture: {herr}")
        return z
    
    def _calculate_all_properties(self, base_props: Dict[str, Any], 
                                 requested_properties: List[str], 
                                 z: List[float], 
                                 molar_mass: float) -> Dict[str, Dict[str, Any]]:
        """
        Calculate all requested properties.
        
        Args:
            base_props: Dictionary of base properties
            requested_properties: List of properties to calculate
            z: Composition array
            molar_mass: Molecular weight
            
        Returns:
            Dictionary of calculated properties with values and units
        """
        results = {}
        for prop in requested_properties:
            try:
                value = self.property_registry.calculate_property(
                    prop, base_props, self.rp, z
                )
                
                # Format property with units
                if value is not None:
                    unit_system = 'SI'  # TODO: Make configurable
                    unit = self.property_registry.get_property_unit(prop, unit_system)
                    
                    results[prop] = {
                        "value": float(value) if isinstance(value, (int, float, np.number)) else value,
                        "unit": unit
                    }
            except Exception as e:
                logger.warning(f"Error calculating property {prop}: {str(e)}")
                
        return results
    
    def _get_total_grid_points(self, grids: Dict[str, np.ndarray]) -> int:
        """
        Get the total number of grid points.
        
        Args:
            grids: Dictionary of grid arrays
            
        Returns:
            Total number of grid points
        """
        return np.prod([len(grid) for grid in grids.values()])
    
    # Methods to be implemented by subclasses
    def _generate_grids(self, z: List[float], variables: Dict[str, Dict[str, Any]], 
                       options: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """
        Generate calculation grids.
        
        Args:
            z: Composition array
            variables: Dictionary of variables and their ranges
            options: Additional options
            
        Returns:
            Dictionary of grid arrays
        """
        raise NotImplementedError("Subclasses must implement")
    
    def _grid_iterator(self, grids: Dict[str, np.ndarray]) -> Iterator[Tuple[Any, ...]]:
        """
        Iterate over grid points.
        
        Args:
            grids: Dictionary of grid arrays
            
        Yields:
            Tuple with indices and values for each grid point
        """
        raise NotImplementedError("Subclasses must implement")
    
    def _calculate_base_properties(self, z: List[float], grid_point: Tuple[Any, ...]) -> Dict[str, Any]:
        """
        Calculate base properties at a grid point.
        
        Args:
            z: Composition array
            grid_point: Grid point information (from grid_iterator)
            
        Returns:
            Dictionary of base properties
        """
        raise NotImplementedError("Subclasses must implement")
    
    def _get_grid_indices(self, grid_point: Tuple[Any, ...], 
                         grids: Dict[str, np.ndarray]) -> Dict[str, int]:
        """
        Get grid indices for a point.
        
        Args:
            grid_point: Grid point information (from grid_iterator)
            grids: Dictionary of grid arrays
            
        Returns:
            Dictionary of grid indices
        """
        raise NotImplementedError("Subclasses must implement")
    
    def _prepare_grid_info(self, grids: Dict[str, np.ndarray], options: Dict[str, Any], 
                          result_count: int) -> Dict[str, Any]:
        """
        Prepare grid information.
        
        Args:
            grids: Dictionary of grid arrays
            options: Additional options
            result_count: Number of successful results
            
        Returns:
            Dictionary with grid information
        """
        raise NotImplementedError("Subclasses must implement")


class PTFlashCalculator(FlashCalculator):
    """
    PT-Flash specific implementation.
    
    This class implements the pressure-temperature flash calculation.
    """
    
    def _generate_grids(self, z: List[float], variables: Dict[str, Dict[str, Any]], 
                       options: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """
        Generate PT grids.
        
        Args:
            z: Composition array
            variables: Dictionary of variables and their ranges
            options: Additional options
            
        Returns:
            Dictionary with pressure and temperature grids
        """
        from API.utils.grid_generator import generate_grid, get_phase_boundaries_pt
        
        # Extract variables
        p_range = variables["pressure"]["range"]
        t_range = variables["temperature"]["range"]
        p_res = variables["pressure"]["resolution"]
        t_res = variables["temperature"]["resolution"]
        
        # Get grid options
        grid_type = options.get("grid_type", "equidistant")
        enhancement_factor = options.get("enhancement_factor", 5.0)
        boundary_zone_width = options.get("boundary_zone_width")
        
        # Get phase boundaries for adaptive grid
        t_boundaries, p_boundaries = [], []
        if grid_type == "adaptive":
            try:
                t_boundaries, p_boundaries = get_phase_boundaries_pt(
                    self.rp, z, t_range, p_range
                )
                logger.info(f"Found {len(t_boundaries)} temperature and {len(p_boundaries)} pressure phase boundaries")
            except Exception as e:
                logger.warning(f"Error getting phase boundaries: {str(e)}")
        
        # Generate pressure grid (in bar)
        P_grid = generate_grid(
            p_range["from"], p_range["to"], p_res,
            grid_type, p_boundaries,
            enhancement_factor, boundary_zone_width
        )
        
        # Generate temperature grid (convert from °C to K)
        T_grid = generate_grid(
            t_range["from"] + 273.15, t_range["to"] + 273.15, t_res,
            grid_type, t_boundaries,
            enhancement_factor, boundary_zone_width
        )
        
        logger.info(f"Generated PT grid with {len(P_grid)} pressure points and {len(T_grid)} temperature points")
        
        return {"pressure": P_grid, "temperature": T_grid}
    
    def _grid_iterator(self, grids: Dict[str, np.ndarray]) -> Iterator[Tuple[int, int, float, float]]:
        """
        Iterate over PT grid points.
        
        Args:
            grids: Dictionary with pressure and temperature grids
            
        Yields:
            Tuple with (p_idx, t_idx, P_bar, T_K)
        """
        P_grid = grids["pressure"]
        T_grid = grids["temperature"]
        
        for p_idx, P in enumerate(P_grid):
            for t_idx, T in enumerate(T_grid):
                yield (p_idx, t_idx, P, T)
    
    def _calculate_base_properties(self, z: List[float], grid_point: Tuple[int, int, float, float]) -> Dict[str, Any]:
        """
        Calculate base properties at PT point.
        
        Args:
            z: Composition array
            grid_point: Tuple with (p_idx, t_idx, P_bar, T_K)
            
        Returns:
            Dictionary of base properties
            
        Raises:
            ValueError: If calculation fails
        """
        p_idx, t_idx, P_bar, T_K = grid_point
        
        # Convert pressure from bar to kPa for REFPROP
        P_kpa = P_bar * 100
        
        # Call REFPROP TPFLSHdll
        result = self.rp.TPFLSHdll(T_K, P_kpa, z)
        
        if result.ierr > 0:
            raise ValueError(f"Error in TPFLSHdll: {result.herr}")
        
        # Extract results
        D = result.D      # Overall density [mol/L]
        Dl = result.Dl    # Liquid density [mol/L]
        Dv = result.Dv    # Vapor density [mol/L]
        q = result.q      # Vapor fraction (quality) [mol/mol]
        e = result.e      # Internal energy [J/mol]
        h = result.h      # Enthalpy [J/mol]
        s = result.s      # Entropy [J/mol-K]
        Cv = result.Cv    # Isochoric heat capacity [J/mol-K]
        Cp = result.Cp    # Isobaric heat capacity [J/mol-K]
        w = result.w      # Speed of sound [m/s]
        x = result.x      # Liquid composition [mol/mol]
        y = result.y      # Vapor composition [mol/mol]
        
        # Get transport properties
        try:
            eta, tcx, ierr_trn, herr_trn = self.rp.TRNPRPdll(T_K, D, z)
            if ierr_trn > 0:
                logger.warning(f"Transport properties warning: {herr_trn}")
                eta = tcx = None
        except Exception as e:
            logger.warning(f"Error calculating transport properties: {str(e)}")
            eta = tcx = None
            
        # Get surface tension for two-phase
        surface_tension = None
        if 0 < q < 1:
            try:
                sigma, ierr_st, herr_st = self.rp.SURTENdll(T_K, Dl, Dv, x, y)
                if ierr_st == 0:
                    surface_tension = sigma
            except Exception as e:
                logger.warning(f"Error calculating surface tension: {str(e)}")
                
        # Get critical properties
        try:
            Tc, Pc, Dc, ierr_crit, herr_crit = self.rp.CRITPdll(z)
            if ierr_crit > 0:
                logger.warning(f"Critical properties warning: {herr_crit}")
                Tc = Pc = Dc = None
        except Exception as e:
            logger.warning(f"Error calculating critical properties: {str(e)}")
            Tc = Pc = Dc = None
            
        # Get thermodynamic derivatives
        try:
            derivatives = self.rp.DERVPVTdll(T_K, D, z)
            dPdD = derivatives.dPdD
            dPdT = derivatives.dPdT
            dDdP = derivatives.dDdP
            dDdT = derivatives.dDdT
        except Exception as e:
            logger.warning(f"Error calculating derivatives: {str(e)}")
            dPdD = dPdT = dDdP = dDdT = None
            
        # Calculate phase-specific properties for two-phase regions
        liquid_cp = vapor_cp = None
        liquid_cv = vapor_cv = None
        liquid_enthalpy = vapor_enthalpy = None
        liquid_entropy = vapor_entropy = None
        liquid_viscosity = vapor_viscosity = None
        liquid_thermal_conductivity = vapor_thermal_conductivity = None
        dDdP_liquid = dDdP_vapor = None
        dDdT_liquid = dDdT_vapor = None
        
        if 0 < q < 1:
            # Two-phase region - calculate properties for each phase
            try:
                # Liquid phase properties
                P_l, e_l, liquid_enthalpy, liquid_entropy, liquid_cv, liquid_cp, w_l, hjt_l = self.rp.THERMdll(T_K, Dl, z)
                
                # Vapor phase properties
                P_v, e_v, vapor_enthalpy, vapor_entropy, vapor_cv, vapor_cp, w_v, hjt_v = self.rp.THERMdll(T_K, Dv, z)
                
                # Transport properties
                eta_l, tcx_l, ierr_l, herr_l = self.rp.TRNPRPdll(T_K, Dl, z)
                if ierr_l == 0:
                    liquid_viscosity = eta_l
                    liquid_thermal_conductivity = tcx_l
                    
                eta_v, tcx_v, ierr_v, herr_v = self.rp.TRNPRPdll(T_K, Dv, z)
                if ierr_v == 0:
                    vapor_viscosity = eta_v
                    vapor_thermal_conductivity = tcx_v
                    
                # Derivatives for each phase
                try:
                    derivs_l = self.rp.DERVPVTdll(T_K, Dl, z)
                    dDdP_liquid = derivs_l.dDdP
                    dDdT_liquid = derivs_l.dDdT
                    
                    derivs_v = self.rp.DERVPVTdll(T_K, Dv, z)
                    dDdP_vapor = derivs_v.dDdP
                    dDdT_vapor = derivs_v.dDdT
                except Exception as e:
                    logger.warning(f"Error calculating phase-specific derivatives: {str(e)}")
            except Exception as e:
                logger.warning(f"Error calculating phase-specific properties: {str(e)}")
        elif q == 0:
            # Single-phase liquid
            liquid_cp = Cp
            liquid_cv = Cv
            liquid_enthalpy = h
            liquid_entropy = s
            liquid_viscosity = eta
            liquid_thermal_conductivity = tcx
            dDdP_liquid = dDdP
            dDdT_liquid = dDdT
        elif q == 1:
            # Single-phase vapor
            vapor_cp = Cp
            vapor_cv = Cv
            vapor_enthalpy = h
            vapor_entropy = s
            vapor_viscosity = eta
            vapor_thermal_conductivity = tcx
            dDdP_vapor = dDdP
            dDdT_vapor = dDdT
            
        # Return base properties
        return {
            "temperature": T_K - 273.15,  # K to °C
            "pressure": P_bar,      # Already in bar
            "density": D,
            "liquid_density": Dl,
            "vapor_density": Dv,
            "vapor_fraction": q,
            "internal_energy": e,
            "enthalpy": h,
            "entropy": s,
            "cv": Cv,
            "cp": Cp,
            "sound_speed": w,
            "viscosity": eta,
            "thermal_conductivity": tcx,
            "surface_tension": surface_tension,
            "critical_temperature": Tc,
            "critical_pressure": Pc / 100 if Pc is not None else None,  # kPa to bar
            "critical_density": Dc,
            "dDdP": dDdP,
            "dDdT": dDdT,
            "x": list(x[:len(z)]) if x is not None else None,  # Liquid composition 
            "y": list(y[:len(z)]) if y is not None else None,  # Vapor composition
            # Phase-specific properties
            "liquid_cp": liquid_cp,
            "vapor_cp": vapor_cp,
            "liquid_cv": liquid_cv,
            "vapor_cv": vapor_cv,
            "liquid_enthalpy": liquid_enthalpy,
            "vapor_enthalpy": vapor_enthalpy,
            "liquid_entropy": liquid_entropy,
            "vapor_entropy": vapor_entropy,
            "liquid_viscosity": liquid_viscosity, 
            "vapor_viscosity": vapor_viscosity,
            "liquid_thermal_conductivity": liquid_thermal_conductivity,
            "vapor_thermal_conductivity": vapor_thermal_conductivity,
            "dDdP_liquid": dDdP_liquid,
            "dDdP_vapor": dDdP_vapor,
            "dDdT_liquid": dDdT_liquid,
            "dDdT_vapor": dDdT_vapor
        }
    
    def _get_grid_indices(self, grid_point: Tuple[int, int, float, float], 
                         grids: Dict[str, np.ndarray]) -> Dict[str, int]:
        """
        Get grid indices for PT point.
        
        Args:
            grid_point: Tuple with (p_idx, t_idx, P_bar, T_K)
            grids: Dictionary with pressure and temperature grids
            
        Returns:
            Dictionary with p_idx and t_idx
        """
        p_idx, t_idx, _, _ = grid_point
        return {"p_idx": int(p_idx), "t_idx": int(t_idx)}
    
    def _prepare_grid_info(self, grids: Dict[str, np.ndarray], options: Dict[str, Any], 
                          result_count: int) -> Dict[str, Any]:
        """
        Prepare grid information for PT flash.
        
        Args:
            grids: Dictionary with pressure and temperature grids
            options: Additional options
            result_count: Number of successful results
            
        Returns:
            Dictionary with grid information
        """
        return {
            "type": options.get("grid_type", "equidistant"),
            "pressure_points": len(grids["pressure"]),
            "temperature_points": len(grids["temperature"]),
            "total_points": result_count
        }


# TODO: Implement other flash calculators (PH, TS, VT, UV) with similar structure

class PHFlashCalculator(FlashCalculator):
    """
    PH-Flash specific implementation.
    
    This class implements the pressure-enthalpy flash calculation.
    """
    
    def _generate_grids(self, z: List[float], variables: Dict[str, Dict[str, Any]], 
                       options: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """
        Generate PH grids.
        
        Args:
            z: Composition array
            variables: Dictionary of variables and their ranges
            options: Additional options
            
        Returns:
            Dictionary with pressure and enthalpy grids
        """
        from API.utils.grid_generator import generate_grid, get_phase_boundaries_ph
        
        # Extract variables
        p_range = variables["pressure"]["range"]
        h_range = variables["enthalpy"]["range"]
        p_res = variables["pressure"]["resolution"]
        h_res = variables["enthalpy"]["resolution"]
        
        # Get grid options
        grid_type = options.get("grid_type", "equidistant")
        enhancement_factor = options.get("enhancement_factor", 5.0)
        boundary_zone_width = options.get("boundary_zone_width")
        
        # Get phase boundaries for adaptive grid
        p_boundaries, h_boundaries = [], []
        if grid_type == "adaptive":
            try:
                p_boundaries, h_boundaries = get_phase_boundaries_ph(
                    self.rp, z, p_range, h_range
                )
                logger.info(f"Found {len(p_boundaries)} pressure and {len(h_boundaries)} enthalpy phase boundaries")
            except Exception as e:
                logger.warning(f"Error getting PH phase boundaries: {str(e)}")
        
        # Generate pressure grid (in bar)
        P_grid = generate_grid(
            p_range["from"], p_range["to"], p_res,
            grid_type, p_boundaries,
            enhancement_factor, boundary_zone_width
        )
        
        # Generate enthalpy grid (in J/mol)
        h_grid = generate_grid(
            h_range["from"], h_range["to"], h_res,
            grid_type, h_boundaries,
            enhancement_factor, boundary_zone_width
        )
        
        logger.info(f"Generated PH grid with {len(P_grid)} pressure points and {len(h_grid)} enthalpy points")
        
        return {"pressure": P_grid, "enthalpy": h_grid}
    
    def _grid_iterator(self, grids: Dict[str, np.ndarray]) -> Iterator[Tuple[int, int, float, float]]:
        """
        Iterate over PH grid points.
        
        Args:
            grids: Dictionary with pressure and enthalpy grids
            
        Yields:
            Tuple with (p_idx, h_idx, P_bar, h_J_mol)
        """
        P_grid = grids["pressure"]
        h_grid = grids["enthalpy"]
        
        for p_idx, P in enumerate(P_grid):
            for h_idx, h in enumerate(h_grid):
                yield (p_idx, h_idx, P, h)
    
    def _calculate_base_properties(self, z: List[float], grid_point: Tuple[int, int, float, float]) -> Dict[str, Any]:
        """
        Calculate base properties at PH point.
        
        Args:
            z: Composition array
            grid_point: Tuple with (p_idx, h_idx, P_bar, h_J_mol)
            
        Returns:
            Dictionary of base properties
            
        Raises:
            ValueError: If calculation fails
        """
        p_idx, h_idx, P_bar, h_J_mol = grid_point
        
        # Convert pressure from bar to kPa for REFPROP
        P_kpa = P_bar * 100
        
        # Call REFPROP PHFLSHdll
        result = self.rp.PHFLSHdll(P_kpa, h_J_mol, z)
        
        if result.ierr > 0:
            raise ValueError(f"Error in PHFLSHdll: {result.herr}")
        
        # Extract results
        T = result.T      # Temperature [K]
        D = result.D      # Overall density [mol/L]
        Dl = result.Dl    # Liquid density [mol/L]
        Dv = result.Dv    # Vapor density [mol/L]
        q = result.q      # Vapor fraction (quality) [mol/mol]
        e = result.e      # Internal energy [J/mol]
        s = result.s      # Entropy [J/mol-K]
        Cv = result.Cv    # Isochoric heat capacity [J/mol-K]
        Cp = result.Cp    # Isobaric heat capacity [J/mol-K]
        w = result.w      # Speed of sound [m/s]
        x = result.x      # Liquid composition [mol/mol]
        y = result.y      # Vapor composition [mol/mol]
        
        # Get transport properties
        try:
            transport = self.rp.TRNPRPdll(T, D, z)
            eta = transport.eta  # Viscosity [μPa·s]
            tcx = transport.tcx  # Thermal conductivity [W/(m·K)]
        except Exception as e:
            logger.warning(f"Error calculating transport properties: {str(e)}")
            eta = tcx = None
            
        # Get surface tension for two-phase
        surface_tension = None
        if 0 < q < 1:
            try:
                sigma_result = self.rp.SURTENdll(T, Dl, Dv, x, y)
                if sigma_result.ierr == 0:
                    surface_tension = sigma_result.sigma
            except Exception as e:
                logger.warning(f"Error calculating surface tension: {str(e)}")
                
        # Get critical properties
        try:
            crit_result = self.rp.CRITPdll(z)
            if crit_result.ierr == 0:
                Tc = crit_result.Tc
                Pc = crit_result.Pc
                Dc = crit_result.Dc
            else:
                logger.warning(f"Critical properties warning: {crit_result.herr}")
                Tc = Pc = Dc = None
        except Exception as e:
            logger.warning(f"Error calculating critical properties: {str(e)}")
            Tc = Pc = Dc = None
            
        # Get thermodynamic derivatives
        try:
            derivatives = self.rp.DERVPVTdll(T, D, z)
            dPdD = derivatives.dPdD
            dPdT = derivatives.dPdT
            dDdP = derivatives.dDdP
            dDdT = derivatives.dDdT
        except Exception as e:
            logger.warning(f"Error calculating derivatives: {str(e)}")
            dPdD = dPdT = dDdP = dDdT = None
            
        # Calculate phase-specific properties for two-phase regions
        liquid_cp = vapor_cp = None
        liquid_cv = vapor_cv = None
        liquid_enthalpy = vapor_enthalpy = None
        liquid_entropy = vapor_entropy = None
        liquid_viscosity = vapor_viscosity = None
        liquid_thermal_conductivity = vapor_thermal_conductivity = None
        dDdP_liquid = dDdP_vapor = None
        dDdT_liquid = dDdT_vapor = None
        
        if 0 < q < 1:
            # Two-phase region - calculate properties for each phase
            try:
                # Liquid phase properties
                therm_l = self.rp.THERMdll(T, Dl, z)
                liquid_enthalpy = therm_l.h
                liquid_entropy = therm_l.s
                liquid_cv = therm_l.Cv
                liquid_cp = therm_l.Cp
                
                # Vapor phase properties
                therm_v = self.rp.THERMdll(T, Dv, z)
                vapor_enthalpy = therm_v.h
                vapor_entropy = therm_v.s
                vapor_cv = therm_v.Cv
                vapor_cp = therm_v.Cp
                
                # Transport properties
                trans_l = self.rp.TRNPRPdll(T, Dl, z)
                if trans_l.ierr == 0:
                    liquid_viscosity = trans_l.eta
                    liquid_thermal_conductivity = trans_l.tcx
                    
                trans_v = self.rp.TRNPRPdll(T, Dv, z)
                if trans_v.ierr == 0:
                    vapor_viscosity = trans_v.eta
                    vapor_thermal_conductivity = trans_v.tcx
                    
                # Derivatives for each phase
                try:
                    derivs_l = self.rp.DERVPVTdll(T, Dl, z)
                    dDdP_liquid = derivs_l.dDdP
                    dDdT_liquid = derivs_l.dDdT
                    
                    derivs_v = self.rp.DERVPVTdll(T, Dv, z)
                    dDdP_vapor = derivs_v.dDdP
                    dDdT_vapor = derivs_v.dDdT
                except Exception as e:
                    logger.warning(f"Error calculating phase-specific derivatives: {str(e)}")
            except Exception as e:
                logger.warning(f"Error calculating phase-specific properties: {str(e)}")
        elif q == 0:
            # Single-phase liquid
            liquid_cp = Cp
            liquid_cv = Cv
            liquid_enthalpy = h_J_mol  # Use the input enthalpy
            liquid_entropy = s
            liquid_viscosity = eta
            liquid_thermal_conductivity = tcx
            dDdP_liquid = dDdP
            dDdT_liquid = dDdT
        elif q == 1:
            # Single-phase vapor
            vapor_cp = Cp
            vapor_cv = Cv
            vapor_enthalpy = h_J_mol  # Use the input enthalpy
            vapor_entropy = s
            vapor_viscosity = eta
            vapor_thermal_conductivity = tcx
            dDdP_vapor = dDdP
            dDdT_vapor = dDdT
            
        # Return base properties
        return {
            "temperature": T - 273.15,  # K to °C
            "pressure": P_bar,      # Already in bar
            "enthalpy": h_J_mol,    # Already in J/mol
            "density": D,
            "liquid_density": Dl,
            "vapor_density": Dv,
            "vapor_fraction": q,
            "internal_energy": e,
            "entropy": s,
            "cv": Cv,
            "cp": Cp,
            "sound_speed": w,
            "viscosity": eta,
            "thermal_conductivity": tcx,
            "surface_tension": surface_tension,
            "critical_temperature": Tc,
            "critical_pressure": Pc / 100 if Pc is not None else None,  # kPa to bar
            "critical_density": Dc,
            "dDdP": dDdP,
            "dDdT": dDdT,
            "x": list(x[:len(z)]) if x is not None else None,  # Liquid composition 
            "y": list(y[:len(z)]) if y is not None else None,  # Vapor composition
            # Phase-specific properties
            "liquid_cp": liquid_cp,
            "vapor_cp": vapor_cp,
            "liquid_cv": liquid_cv,
            "vapor_cv": vapor_cv,
            "liquid_enthalpy": liquid_enthalpy,
            "vapor_enthalpy": vapor_enthalpy,
            "liquid_entropy": liquid_entropy,
            "vapor_entropy": vapor_entropy,
            "liquid_viscosity": liquid_viscosity, 
            "vapor_viscosity": vapor_viscosity,
            "liquid_thermal_conductivity": liquid_thermal_conductivity,
            "vapor_thermal_conductivity": vapor_thermal_conductivity,
            "dDdP_liquid": dDdP_liquid,
            "dDdP_vapor": dDdP_vapor,
            "dDdT_liquid": dDdT_liquid,
            "dDdT_vapor": dDdT_vapor
        }
    
    def _get_grid_indices(self, grid_point: Tuple[int, int, float, float], 
                         grids: Dict[str, np.ndarray]) -> Dict[str, int]:
        """
        Get grid indices for PH point.
        
        Args:
            grid_point: Tuple with (p_idx, h_idx, P_bar, h_J_mol)
            grids: Dictionary with pressure and enthalpy grids
            
        Returns:
            Dictionary with p_idx and h_idx
        """
        p_idx, h_idx, _, _ = grid_point
        return {"p_idx": int(p_idx), "h_idx": int(h_idx)}
    
    def _prepare_grid_info(self, grids: Dict[str, np.ndarray], options: Dict[str, Any], 
                          result_count: int) -> Dict[str, Any]:
        """
        Prepare grid information for PH flash.
        
        Args:
            grids: Dictionary with pressure and enthalpy grids
            options: Additional options
            result_count: Number of successful results
            
        Returns:
            Dictionary with grid information
        """
        return {
            "type": options.get("grid_type", "equidistant"),
            "pressure_points": len(grids["pressure"]),
            "enthalpy_points": len(grids["enthalpy"]),
            "total_points": result_count
        }


class TSFlashCalculator(FlashCalculator):
    """
    TS-Flash specific implementation.
    
    This class implements the temperature-entropy flash calculation.
    """
    
    def _generate_grids(self, z: List[float], variables: Dict[str, Dict[str, Any]], 
                       options: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """
        Generate TS grids.
        
        Args:
            z: Composition array
            variables: Dictionary of variables and their ranges
            options: Additional options
            
        Returns:
            Dictionary with temperature and entropy grids
        """
        from API.utils.grid_generator import generate_grid, get_phase_boundaries_ts
        
        # Extract variables
        t_range = variables["temperature"]["range"]
        s_range = variables["entropy"]["range"]
        t_res = variables["temperature"]["resolution"]
        s_res = variables["entropy"]["resolution"]
        
        # Get grid options
        grid_type = options.get("grid_type", "equidistant")
        enhancement_factor = options.get("enhancement_factor", 5.0)
        boundary_zone_width = options.get("boundary_zone_width")
        
        # Get phase boundaries for adaptive grid
        t_boundaries, s_boundaries = [], []
        if grid_type == "adaptive":
            try:
                t_boundaries, s_boundaries = get_phase_boundaries_ts(
                    self.rp, z, t_range, s_range
                )
                logger.info(f"Found {len(t_boundaries)} temperature and {len(s_boundaries)} entropy phase boundaries")
            except Exception as e:
                logger.warning(f"Error getting TS phase boundaries: {str(e)}")
        
        # Generate temperature grid (in K)
        T_grid = generate_grid(
            t_range["from"] + 273.15, t_range["to"] + 273.15, t_res,
            grid_type, t_boundaries,
            enhancement_factor, boundary_zone_width
        )
        
        # Generate entropy grid (in J/(mol·K))
        S_grid = generate_grid(
            s_range["from"], s_range["to"], s_res,
            grid_type, s_boundaries,
            enhancement_factor, boundary_zone_width
        )
        
        logger.info(f"Generated TS grid with {len(T_grid)} temperature points and {len(S_grid)} entropy points")
        
        return {"temperature": T_grid, "entropy": S_grid}
    
    def _grid_iterator(self, grids: Dict[str, np.ndarray]) -> Iterator[Tuple[int, int, float, float]]:
        """
        Iterate over TS grid points.
        
        Args:
            grids: Dictionary with temperature and entropy grids
            
        Yields:
            Tuple with (t_idx, s_idx, T_K, S_J_mol_K)
        """
        T_grid = grids["temperature"]
        S_grid = grids["entropy"]
        
        for t_idx, T in enumerate(T_grid):
            for s_idx, S in enumerate(S_grid):
                yield (t_idx, s_idx, T, S)
    
    def _calculate_base_properties(self, z: List[float], grid_point: Tuple[int, int, float, float]) -> Dict[str, Any]:
        """
        Calculate base properties at TS point.
        
        Args:
            z: Composition array
            grid_point: Tuple with (t_idx, s_idx, T_K, S_J_mol_K)
            
        Returns:
            Dictionary of base properties
            
        Raises:
            ValueError: If calculation fails
        """
        t_idx, s_idx, T_K, S_J_mol_K = grid_point
        
        # Call REFPROP TSFLSHdll - Note: kr=1 for initial phase guess (liquid)
        result = self.rp.TSFLSHdll(T_K, S_J_mol_K, z, 1)
        
        if result.ierr > 0:
            # Try again with vapor phase guess (kr=2)
            result = self.rp.TSFLSHdll(T_K, S_J_mol_K, z, 2)
            if result.ierr > 0:
                raise ValueError(f"Error in TSFLSHdll: {result.herr}")
        
        # Extract results
        P = result.P      # Pressure [kPa]
        D = result.D      # Overall density [mol/L]
        Dl = result.Dl    # Liquid density [mol/L]
        Dv = result.Dv    # Vapor density [mol/L]
        q = result.q      # Vapor fraction (quality) [mol/mol]
        e = result.e      # Internal energy [J/mol]
        h = result.h      # Enthalpy [J/mol]
        Cv = result.Cv    # Isochoric heat capacity [J/mol-K]
        Cp = result.Cp    # Isobaric heat capacity [J/mol-K]
        w = result.w      # Speed of sound [m/s]
        x = result.x      # Liquid composition [mol/mol]
        y = result.y      # Vapor composition [mol/mol]
        
        # Get transport properties
        try:
            transport = self.rp.TRNPRPdll(T_K, D, z)
            eta = transport.eta  # Viscosity [μPa·s]
            tcx = transport.tcx  # Thermal conductivity [W/(m·K)]
        except Exception as e:
            logger.warning(f"Error calculating transport properties: {str(e)}")
            eta = tcx = None
            
        # Get surface tension for two-phase
        surface_tension = None
        if 0 < q < 1:
            try:
                sigma_result = self.rp.SURTENdll(T_K, Dl, Dv, x, y)
                if sigma_result.ierr == 0:
                    surface_tension = sigma_result.sigma
            except Exception as e:
                logger.warning(f"Error calculating surface tension: {str(e)}")
                
        # Get critical properties
        try:
            crit_result = self.rp.CRITPdll(z)
            if crit_result.ierr == 0:
                Tc = crit_result.Tc
                Pc = crit_result.Pc
                Dc = crit_result.Dc
            else:
                logger.warning(f"Critical properties warning: {crit_result.herr}")
                Tc = Pc = Dc = None
        except Exception as e:
            logger.warning(f"Error calculating critical properties: {str(e)}")
            Tc = Pc = Dc = None
            
        # Get thermodynamic derivatives
        try:
            derivatives = self.rp.DERVPVTdll(T_K, D, z)
            dPdD = derivatives.dPdD
            dPdT = derivatives.dPdT
            dDdP = derivatives.dDdP
            dDdT = derivatives.dDdT
        except Exception as e:
            logger.warning(f"Error calculating derivatives: {str(e)}")
            dPdD = dPdT = dDdP = dDdT = None
            
        # Calculate phase-specific properties for two-phase regions
        liquid_cp = vapor_cp = None
        liquid_cv = vapor_cv = None
        liquid_enthalpy = vapor_enthalpy = None
        liquid_entropy = vapor_entropy = None
        liquid_viscosity = vapor_viscosity = None
        liquid_thermal_conductivity = vapor_thermal_conductivity = None
        dDdP_liquid = dDdP_vapor = None
        dDdT_liquid = dDdT_vapor = None
        
        if 0 < q < 1:
            # Two-phase region - calculate properties for each phase
            try:
                # Liquid phase properties
                therm_l = self.rp.THERMdll(T_K, Dl, z)
                liquid_enthalpy = therm_l.h
                liquid_entropy = S_J_mol_K  # Same as input entropy for single phase
                liquid_cv = therm_l.Cv
                liquid_cp = therm_l.Cp
                
                # Vapor phase properties
                therm_v = self.rp.THERMdll(T_K, Dv, z)
                vapor_enthalpy = therm_v.h
                vapor_entropy = S_J_mol_K  # Same as input entropy for single phase
                vapor_cv = therm_v.Cv
                vapor_cp = therm_v.Cp
                
                # Transport properties
                trans_l = self.rp.TRNPRPdll(T_K, Dl, z)
                if trans_l.ierr == 0:
                    liquid_viscosity = trans_l.eta
                    liquid_thermal_conductivity = trans_l.tcx
                    
                trans_v = self.rp.TRNPRPdll(T_K, Dv, z)
                if trans_v.ierr == 0:
                    vapor_viscosity = trans_v.eta
                    vapor_thermal_conductivity = trans_v.tcx
                    
                # Derivatives for each phase
                try:
                    derivs_l = self.rp.DERVPVTdll(T_K, Dl, z)
                    dDdP_liquid = derivs_l.dDdP
                    dDdT_liquid = derivs_l.dDdT
                    
                    derivs_v = self.rp.DERVPVTdll(T_K, Dv, z)
                    dDdP_vapor = derivs_v.dDdP
                    dDdT_vapor = derivs_v.dDdT
                except Exception as e:
                    logger.warning(f"Error calculating phase-specific derivatives: {str(e)}")
            except Exception as e:
                logger.warning(f"Error calculating phase-specific properties: {str(e)}")
        elif q == 0:
            # Single-phase liquid
            liquid_cp = Cp
            liquid_cv = Cv
            liquid_enthalpy = h
            liquid_entropy = S_J_mol_K  # Same as input entropy
            liquid_viscosity = eta
            liquid_thermal_conductivity = tcx
            dDdP_liquid = dDdP
            dDdT_liquid = dDdT
        elif q == 1:
            # Single-phase vapor
            vapor_cp = Cp
            vapor_cv = Cv
            vapor_enthalpy = h
            vapor_entropy = S_J_mol_K  # Same as input entropy
            vapor_viscosity = eta
            vapor_thermal_conductivity = tcx
            dDdP_vapor = dDdP
            dDdT_vapor = dDdT
            
        # Return base properties
        return {
            "temperature": T_K - 273.15,  # K to °C
            "pressure": P / 100,          # kPa to bar
            "entropy": S_J_mol_K,         # Already in J/(mol·K)
            "density": D,
            "liquid_density": Dl,
            "vapor_density": Dv,
            "vapor_fraction": q,
            "internal_energy": e,
            "enthalpy": h,
            "cv": Cv,
            "cp": Cp,
            "sound_speed": w,
            "viscosity": eta,
            "thermal_conductivity": tcx,
            "surface_tension": surface_tension,
            "critical_temperature": Tc,
            "critical_pressure": Pc / 100 if Pc is not None else None,  # kPa to bar
            "critical_density": Dc,
            "dDdP": dDdP,
            "dDdT": dDdT,
            "x": list(x[:len(z)]) if x is not None else None,  # Liquid composition 
            "y": list(y[:len(z)]) if y is not None else None,  # Vapor composition
            # Phase-specific properties
            "liquid_cp": liquid_cp,
            "vapor_cp": vapor_cp,
            "liquid_cv": liquid_cv,
            "vapor_cv": vapor_cv,
            "liquid_enthalpy": liquid_enthalpy,
            "vapor_enthalpy": vapor_enthalpy,
            "liquid_entropy": liquid_entropy,
            "vapor_entropy": vapor_entropy,
            "liquid_viscosity": liquid_viscosity, 
            "vapor_viscosity": vapor_viscosity,
            "liquid_thermal_conductivity": liquid_thermal_conductivity,
            "vapor_thermal_conductivity": vapor_thermal_conductivity,
            "dDdP_liquid": dDdP_liquid,
            "dDdP_vapor": dDdP_vapor,
            "dDdT_liquid": dDdT_liquid,
            "dDdT_vapor": dDdT_vapor
        }
    
    def _get_grid_indices(self, grid_point: Tuple[int, int, float, float], 
                         grids: Dict[str, np.ndarray]) -> Dict[str, int]:
        """
        Get grid indices for TS point.
        
        Args:
            grid_point: Tuple with (t_idx, s_idx, T_K, S_J_mol_K)
            grids: Dictionary with temperature and entropy grids
            
        Returns:
            Dictionary with t_idx and s_idx
        """
        t_idx, s_idx, _, _ = grid_point
        return {"t_idx": int(t_idx), "s_idx": int(s_idx)}
    
    def _prepare_grid_info(self, grids: Dict[str, np.ndarray], options: Dict[str, Any], 
                          result_count: int) -> Dict[str, Any]:
        """
        Prepare grid information for TS flash.
        
        Args:
            grids: Dictionary with temperature and entropy grids
            options: Additional options
            result_count: Number of successful results
            
        Returns:
            Dictionary with grid information
        """
        return {
            "type": options.get("grid_type", "equidistant"),
            "temperature_points": len(grids["temperature"]),
            "entropy_points": len(grids["entropy"]),
            "total_points": result_count
        }


class VTFlashCalculator(FlashCalculator):
    """
    VT-Flash specific implementation.
    
    This class implements the volume-temperature (or density-temperature) flash calculation.
    """
    
    def _generate_grids(self, z: List[float], variables: Dict[str, Dict[str, Any]], 
                       options: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """
        Generate VT grids.
        
        Args:
            z: Composition array
            variables: Dictionary of variables and their ranges
            options: Additional options
            
        Returns:
            Dictionary with specific_volume and temperature grids
        """
        from API.utils.grid_generator import generate_grid, get_phase_boundaries_tv
        
        # Extract variables
        v_range = variables["specific_volume"]["range"]
        t_range = variables["temperature"]["range"]
        v_res = variables["specific_volume"]["resolution"]
        t_res = variables["temperature"]["resolution"]
        
        # Get grid options
        grid_type = options.get("grid_type", "equidistant")
        enhancement_factor = options.get("enhancement_factor", 5.0)
        boundary_zone_width = options.get("boundary_zone_width")
        
        # Get phase boundaries for adaptive grid
        t_boundaries, v_boundaries = [], []
        if grid_type == "adaptive":
            try:
                t_boundaries, v_boundaries = get_phase_boundaries_tv(
                    self.rp, z, t_range, v_range
                )
                logger.info(f"Found {len(t_boundaries)} temperature and {len(v_boundaries)} volume phase boundaries")
            except Exception as e:
                logger.warning(f"Error getting TV phase boundaries: {str(e)}")
        
        # Generate temperature grid (in K)
        T_grid = generate_grid(
            t_range["from"] + 273.15, t_range["to"] + 273.15, t_res,
            grid_type, t_boundaries,
            enhancement_factor, boundary_zone_width
        )
        
        # Generate specific volume grid (in m³/mol)
        V_grid = generate_grid(
            v_range["from"], v_range["to"], v_res,
            grid_type, v_boundaries,
            enhancement_factor, boundary_zone_width
        )
        
        logger.info(f"Generated VT grid with {len(V_grid)} specific volume points and {len(T_grid)} temperature points")
        
        return {"specific_volume": V_grid, "temperature": T_grid}
    
    def _grid_iterator(self, grids: Dict[str, np.ndarray]) -> Iterator[Tuple[int, int, float, float]]:
        """
        Iterate over VT grid points.
        
        Args:
            grids: Dictionary with specific_volume and temperature grids
            
        Yields:
            Tuple with (v_idx, t_idx, V_m3_mol, T_K)
        """
        V_grid = grids["specific_volume"]
        T_grid = grids["temperature"]
        
        for v_idx, V in enumerate(V_grid):
            for t_idx, T in enumerate(T_grid):
                yield (v_idx, t_idx, V, T)
    
    def _calculate_base_properties(self, z: List[float], grid_point: Tuple[int, int, float, float]) -> Dict[str, Any]:
        """
        Calculate base properties at VT point.
        
        Args:
            z: Composition array
            grid_point: Tuple with (v_idx, t_idx, V_m3_mol, T_K)
            
        Returns:
            Dictionary of base properties
            
        Raises:
            ValueError: If calculation fails
        """
        v_idx, t_idx, V_m3_mol, T_K = grid_point
        
        # Convert specific volume to density (mol/L)
        # REFPROP uses mol/L = 1000 mol/m³, so 1/V * 1000
        D = 1.0 / V_m3_mol  # mol/L
        
        # Call REFPROP TDFLSHdll
        result = self.rp.TDFLSHdll(T_K, D, z)
        
        if result.ierr > 0:
            raise ValueError(f"Error in TDFLSHdll: {result.herr}")
        
        # Extract results
        P = result.P      # Pressure [kPa]
        Dl = result.Dl    # Liquid density [mol/L]
        Dv = result.Dv    # Vapor density [mol/L]
        q = result.q      # Vapor fraction (quality) [mol/mol]
        e = result.e      # Internal energy [J/mol]
        h = result.h      # Enthalpy [J/mol]
        s = result.s      # Entropy [J/(mol·K)]
        Cv = result.Cv    # Isochoric heat capacity [J/(mol·K)]
        Cp = result.Cp    # Isobaric heat capacity [J/(mol·K)]
        w = result.w      # Speed of sound [m/s]
        x = result.x      # Liquid composition [mol/mol]
        y = result.y      # Vapor composition [mol/mol]
        
        # Get transport properties
        try:
            transport = self.rp.TRNPRPdll(T_K, D, z)
            eta = transport.eta  # Viscosity [μPa·s]
            tcx = transport.tcx  # Thermal conductivity [W/(m·K)]
        except Exception as e:
            logger.warning(f"Error calculating transport properties: {str(e)}")
            eta = tcx = None
            
        # Get surface tension for two-phase
        surface_tension = None
        if 0 < q < 1:
            try:
                sigma_result = self.rp.SURTENdll(T_K, Dl, Dv, x, y)
                if sigma_result.ierr == 0:
                    surface_tension = sigma_result.sigma
            except Exception as e:
                logger.warning(f"Error calculating surface tension: {str(e)}")
                
        # Get critical properties
        try:
            crit_result = self.rp.CRITPdll(z)
            if crit_result.ierr == 0:
                Tc = crit_result.Tc
                Pc = crit_result.Pc
                Dc = crit_result.Dc
            else:
                logger.warning(f"Critical properties warning: {crit_result.herr}")
                Tc = Pc = Dc = None
        except Exception as e:
            logger.warning(f"Error calculating critical properties: {str(e)}")
            Tc = Pc = Dc = None
            
        # Get thermodynamic derivatives
        try:
            derivatives = self.rp.DERVPVTdll(T_K, D, z)
            dPdD = derivatives.dPdD
            dPdT = derivatives.dPdT
            dDdP = derivatives.dDdP
            dDdT = derivatives.dDdT
        except Exception as e:
            logger.warning(f"Error calculating derivatives: {str(e)}")
            dPdD = dPdT = dDdP = dDdT = None
            
        # Calculate phase-specific properties for two-phase regions
        liquid_cp = vapor_cp = None
        liquid_cv = vapor_cv = None
        liquid_enthalpy = vapor_enthalpy = None
        liquid_entropy = vapor_entropy = None
        liquid_viscosity = vapor_viscosity = None
        liquid_thermal_conductivity = vapor_thermal_conductivity = None
        dDdP_liquid = dDdP_vapor = None
        dDdT_liquid = dDdT_vapor = None
        
        if 0 < q < 1:
            # Two-phase region - calculate properties for each phase
            try:
                # Liquid phase properties
                therm_l = self.rp.THERMdll(T_K, Dl, z)
                liquid_enthalpy = therm_l.h
                liquid_entropy = therm_l.s
                liquid_cv = therm_l.Cv
                liquid_cp = therm_l.Cp
                
                # Vapor phase properties
                therm_v = self.rp.THERMdll(T_K, Dv, z)
                vapor_enthalpy = therm_v.h
                vapor_entropy = therm_v.s
                vapor_cv = therm_v.Cv
                vapor_cp = therm_v.Cp
                
                # Transport properties
                trans_l = self.rp.TRNPRPdll(T_K, Dl, z)
                if trans_l.ierr == 0:
                    liquid_viscosity = trans_l.eta
                    liquid_thermal_conductivity = trans_l.tcx
                    
                trans_v = self.rp.TRNPRPdll(T_K, Dv, z)
                if trans_v.ierr == 0:
                    vapor_viscosity = trans_v.eta
                    vapor_thermal_conductivity = trans_v.tcx
                    
                # Derivatives for each phase
                try:
                    derivs_l = self.rp.DERVPVTdll(T_K, Dl, z)
                    dDdP_liquid = derivs_l.dDdP
                    dDdT_liquid = derivs_l.dDdT
                    
                    derivs_v = self.rp.DERVPVTdll(T_K, Dv, z)
                    dDdP_vapor = derivs_v.dDdP
                    dDdT_vapor = derivs_v.dDdT
                except Exception as e:
                    logger.warning(f"Error calculating phase-specific derivatives: {str(e)}")
            except Exception as e:
                logger.warning(f"Error calculating phase-specific properties: {str(e)}")
        elif q == 0:
            # Single-phase liquid
            liquid_cp = Cp
            liquid_cv = Cv
            liquid_enthalpy = h
            liquid_entropy = s
            liquid_viscosity = eta
            liquid_thermal_conductivity = tcx
            dDdP_liquid = dDdP
            dDdT_liquid = dDdT
        elif q == 1:
            # Single-phase vapor
            vapor_cp = Cp
            vapor_cv = Cv
            vapor_enthalpy = h
            vapor_entropy = s
            vapor_viscosity = eta
            vapor_thermal_conductivity = tcx
            dDdP_vapor = dDdP
            dDdT_vapor = dDdT
            
        # Return base properties
        return {
            "temperature": T_K - 273.15,    # K to °C
            "pressure": P / 100,            # kPa to bar
            "specific_volume": V_m3_mol,    # Already in m³/mol
            "density": D,                   # Already in mol/L
            "liquid_density": Dl,
            "vapor_density": Dv,
            "vapor_fraction": q,
            "internal_energy": e,
            "enthalpy": h,
            "entropy": s,
            "cv": Cv,
            "cp": Cp,
            "sound_speed": w,
            "viscosity": eta,
            "thermal_conductivity": tcx,
            "surface_tension": surface_tension,
            "critical_temperature": Tc,
            "critical_pressure": Pc / 100 if Pc is not None else None,  # kPa to bar
            "critical_density": Dc,
            "dDdP": dDdP,
            "dDdT": dDdT,
            "x": list(x[:len(z)]) if x is not None else None,  # Liquid composition 
            "y": list(y[:len(z)]) if y is not None else None,  # Vapor composition
            # Phase-specific properties
            "liquid_cp": liquid_cp,
            "vapor_cp": vapor_cp,
            "liquid_cv": liquid_cv,
            "vapor_cv": vapor_cv,
            "liquid_enthalpy": liquid_enthalpy,
            "vapor_enthalpy": vapor_enthalpy,
            "liquid_entropy": liquid_entropy,
            "vapor_entropy": vapor_entropy,
            "liquid_viscosity": liquid_viscosity, 
            "vapor_viscosity": vapor_viscosity,
            "liquid_thermal_conductivity": liquid_thermal_conductivity,
            "vapor_thermal_conductivity": vapor_thermal_conductivity,
            "dDdP_liquid": dDdP_liquid,
            "dDdP_vapor": dDdP_vapor,
            "dDdT_liquid": dDdT_liquid,
            "dDdT_vapor": dDdT_vapor
        }
    
    def _get_grid_indices(self, grid_point: Tuple[int, int, float, float], 
                         grids: Dict[str, np.ndarray]) -> Dict[str, int]:
        """
        Get grid indices for VT point.
        
        Args:
            grid_point: Tuple with (v_idx, t_idx, V_m3_mol, T_K)
            grids: Dictionary with specific_volume and temperature grids
            
        Returns:
            Dictionary with v_idx and t_idx
        """
        v_idx, t_idx, _, _ = grid_point
        return {"v_idx": int(v_idx), "t_idx": int(t_idx)}
    
    def _prepare_grid_info(self, grids: Dict[str, np.ndarray], options: Dict[str, Any], 
                          result_count: int) -> Dict[str, Any]:
        """
        Prepare grid information for VT flash.
        
        Args:
            grids: Dictionary with specific_volume and temperature grids
            options: Additional options
            result_count: Number of successful results
            
        Returns:
            Dictionary with grid information
        """
        return {
            "type": options.get("grid_type", "equidistant"),
            "specific_volume_points": len(grids["specific_volume"]),
            "temperature_points": len(grids["temperature"]),
            "total_points": result_count
        }


class UVFlashCalculator(FlashCalculator):
    """
    UV-Flash specific implementation.
    
    This class implements the internal energy-volume (or internal energy-density) flash calculation.
    """
    
    def _generate_grids(self, z: List[float], variables: Dict[str, Dict[str, Any]], 
                       options: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """
        Generate UV grids.
        
        Args:
            z: Composition array
            variables: Dictionary of variables and their ranges
            options: Additional options
            
        Returns:
            Dictionary with internal_energy and specific_volume grids
        """
        from API.utils.grid_generator import generate_grid, get_phase_boundaries_uv
        
        # Extract variables
        u_range = variables["internal_energy"]["range"]
        v_range = variables["specific_volume"]["range"]
        u_res = variables["internal_energy"]["resolution"]
        v_res = variables["specific_volume"]["resolution"]
        
        # Get grid options
        grid_type = options.get("grid_type", "equidistant")
        enhancement_factor = options.get("enhancement_factor", 5.0)
        boundary_zone_width = options.get("boundary_zone_width")
        
        # Get phase boundaries for adaptive grid
        u_boundaries, v_boundaries = [], []
        if grid_type == "adaptive":
            try:
                u_boundaries, v_boundaries = get_phase_boundaries_uv(
                    self.rp, z, u_range, v_range
                )
                logger.info(f"Found {len(u_boundaries)} energy and {len(v_boundaries)} volume phase boundaries")
            except Exception as e:
                logger.warning(f"Error getting UV phase boundaries: {str(e)}")
        
        # Generate internal energy grid (in J/mol)
        U_grid = generate_grid(
            u_range["from"], u_range["to"], u_res,
            grid_type, u_boundaries,
            enhancement_factor, boundary_zone_width
        )
        
        # Generate specific volume grid (in m³/mol)
        V_grid = generate_grid(
            v_range["from"], v_range["to"], v_res,
            grid_type, v_boundaries,
            enhancement_factor, boundary_zone_width
        )
        
        logger.info(f"Generated UV grid with {len(U_grid)} internal energy points and {len(V_grid)} specific volume points")
        
        return {"internal_energy": U_grid, "specific_volume": V_grid}
    
    def _grid_iterator(self, grids: Dict[str, np.ndarray]) -> Iterator[Tuple[int, int, float, float]]:
        """
        Iterate over UV grid points.
        
        Args:
            grids: Dictionary with internal_energy and specific_volume grids
            
        Yields:
            Tuple with (u_idx, v_idx, U_J_mol, V_m3_mol)
        """
        U_grid = grids["internal_energy"]
        V_grid = grids["specific_volume"]
        
        for u_idx, U in enumerate(U_grid):
            for v_idx, V in enumerate(V_grid):
                yield (u_idx, v_idx, U, V)
    
    def _calculate_base_properties(self, z: List[float], grid_point: Tuple[int, int, float, float]) -> Dict[str, Any]:
        """
        Calculate base properties at UV point.
        
        Args:
            z: Composition array
            grid_point: Tuple with (u_idx, v_idx, U_J_mol, V_m3_mol)
            
        Returns:
            Dictionary of base properties
            
        Raises:
            ValueError: If calculation fails
        """
        u_idx, v_idx, U_J_mol, V_m3_mol = grid_point
        
        # Convert specific volume to density (mol/L)
        # REFPROP uses mol/L = 1000 mol/m³, so 1/V * 1000
        D = 1.0 / V_m3_mol  # mol/L
        
        # Call REFPROP DEFLSHdll
        result = self.rp.DEFLSHdll(D, U_J_mol, z)
        
        if result.ierr > 0:
            raise ValueError(f"Error in DEFLSHdll: {result.herr}")
        
        # Extract results
        T = result.T      # Temperature [K]
        P = result.P      # Pressure [kPa]
        Dl = result.Dl    # Liquid density [mol/L]
        Dv = result.Dv    # Vapor density [mol/L]
        q = result.q      # Vapor fraction (quality) [mol/mol]
        h = result.h      # Enthalpy [J/mol]
        s = result.s      # Entropy [J/(mol·K)]
        Cv = result.Cv    # Isochoric heat capacity [J/(mol·K)]
        Cp = result.Cp    # Isobaric heat capacity [J/(mol·K)]
        w = result.w      # Speed of sound [m/s]
        x = result.x      # Liquid composition [mol/mol]
        y = result.y      # Vapor composition [mol/mol]
        
        # Get transport properties
        try:
            transport = self.rp.TRNPRPdll(T, D, z)
            eta = transport.eta  # Viscosity [μPa·s]
            tcx = transport.tcx  # Thermal conductivity [W/(m·K)]
        except Exception as e:
            logger.warning(f"Error calculating transport properties: {str(e)}")
            eta = tcx = None
            
        # Get surface tension for two-phase
        surface_tension = None
        if 0 < q < 1:
            try:
                sigma_result = self.rp.SURTENdll(T, Dl, Dv, x, y)
                if sigma_result.ierr == 0:
                    surface_tension = sigma_result.sigma
            except Exception as e:
                logger.warning(f"Error calculating surface tension: {str(e)}")
                
        # Get critical properties
        try:
            crit_result = self.rp.CRITPdll(z)
            if crit_result.ierr == 0:
                Tc = crit_result.Tc
                Pc = crit_result.Pc
                Dc = crit_result.Dc
            else:
                logger.warning(f"Critical properties warning: {crit_result.herr}")
                Tc = Pc = Dc = None
        except Exception as e:
            logger.warning(f"Error calculating critical properties: {str(e)}")
            Tc = Pc = Dc = None
            
        # Get thermodynamic derivatives
        try:
            derivatives = self.rp.DERVPVTdll(T, D, z)
            dPdD = derivatives.dPdD
            dPdT = derivatives.dPdT
            dDdP = derivatives.dDdP
            dDdT = derivatives.dDdT
        except Exception as e:
            logger.warning(f"Error calculating derivatives: {str(e)}")
            dPdD = dPdT = dDdP = dDdT = None
            
        # Calculate phase-specific properties for two-phase regions
        liquid_cp = vapor_cp = None
        liquid_cv = vapor_cv = None
        liquid_enthalpy = vapor_enthalpy = None
        liquid_entropy = vapor_entropy = None
        liquid_viscosity = vapor_viscosity = None
        liquid_thermal_conductivity = vapor_thermal_conductivity = None
        dDdP_liquid = dDdP_vapor = None
        dDdT_liquid = dDdT_vapor = None
        
        if 0 < q < 1:
            # Two-phase region - calculate properties for each phase
            try:
                # Liquid phase properties
                therm_l = self.rp.THERMdll(T, Dl, z)
                liquid_enthalpy = therm_l.h
                liquid_entropy = therm_l.s
                liquid_cv = therm_l.Cv
                liquid_cp = therm_l.Cp
                
                # Vapor phase properties
                therm_v = self.rp.THERMdll(T, Dv, z)
                vapor_enthalpy = therm_v.h
                vapor_entropy = therm_v.s
                vapor_cv = therm_v.Cv
                vapor_cp = therm_v.Cp
                
                # Transport properties
                trans_l = self.rp.TRNPRPdll(T, Dl, z)
                if trans_l.ierr == 0:
                    liquid_viscosity = trans_l.eta
                    liquid_thermal_conductivity = trans_l.tcx
                    
                trans_v = self.rp.TRNPRPdll(T, Dv, z)
                if trans_v.ierr == 0:
                    vapor_viscosity = trans_v.eta
                    vapor_thermal_conductivity = trans_v.tcx
                    
                # Derivatives for each phase
                try:
                    derivs_l = self.rp.DERVPVTdll(T, Dl, z)
                    dDdP_liquid = derivs_l.dDdP
                    dDdT_liquid = derivs_l.dDdT
                    
                    derivs_v = self.rp.DERVPVTdll(T, Dv, z)
                    dDdP_vapor = derivs_v.dDdP
                    dDdT_vapor = derivs_v.dDdT
                except Exception as e:
                    logger.warning(f"Error calculating phase-specific derivatives: {str(e)}")
            except Exception as e:
                logger.warning(f"Error calculating phase-specific properties: {str(e)}")
        elif q == 0:
            # Single-phase liquid
            liquid_cp = Cp
            liquid_cv = Cv
            liquid_enthalpy = h
            liquid_entropy = s
            liquid_viscosity = eta
            liquid_thermal_conductivity = tcx
            dDdP_liquid = dDdP
            dDdT_liquid = dDdT
        elif q == 1:
            # Single-phase vapor
            vapor_cp = Cp
            vapor_cv = Cv
            vapor_enthalpy = h
            vapor_entropy = s
            vapor_viscosity = eta
            vapor_thermal_conductivity = tcx
            dDdP_vapor = dDdP
            dDdT_vapor = dDdT
            
        # Return base properties
        return {
            "temperature": T - 273.15,      # K to °C
            "pressure": P / 100,            # kPa to bar
            "internal_energy": U_J_mol,     # Already in J/mol
            "specific_volume": V_m3_mol,    # Already in m³/mol
            "density": D,                   # Already in mol/L
            "liquid_density": Dl,
            "vapor_density": Dv,
            "vapor_fraction": q,
            "enthalpy": h,
            "entropy": s,
            "cv": Cv,
            "cp": Cp,
            "sound_speed": w,
            "viscosity": eta,
            "thermal_conductivity": tcx,
            "surface_tension": surface_tension,
            "critical_temperature": Tc,
            "critical_pressure": Pc / 100 if Pc is not None else None,  # kPa to bar
            "critical_density": Dc,
            "dDdP": dDdP,
            "dDdT": dDdT,
            "x": list(x[:len(z)]) if x is not None else None,  # Liquid composition 
            "y": list(y[:len(z)]) if y is not None else None,  # Vapor composition
            # Phase-specific properties
            "liquid_cp": liquid_cp,
            "vapor_cp": vapor_cp,
            "liquid_cv": liquid_cv,
            "vapor_cv": vapor_cv,
            "liquid_enthalpy": liquid_enthalpy,
            "vapor_enthalpy": vapor_enthalpy,
            "liquid_entropy": liquid_entropy,
            "vapor_entropy": vapor_entropy,
            "liquid_viscosity": liquid_viscosity, 
            "vapor_viscosity": vapor_viscosity,
            "liquid_thermal_conductivity": liquid_thermal_conductivity,
            "vapor_thermal_conductivity": vapor_thermal_conductivity,
            "dDdP_liquid": dDdP_liquid,
            "dDdP_vapor": dDdP_vapor,
            "dDdT_liquid": dDdT_liquid,
            "dDdT_vapor": dDdT_vapor
        }
    
    def _get_grid_indices(self, grid_point: Tuple[int, int, float, float], 
                         grids: Dict[str, np.ndarray]) -> Dict[str, int]:
        """
        Get grid indices for UV point.
        
        Args:
            grid_point: Tuple with (u_idx, v_idx, U_J_mol, V_m3_mol)
            grids: Dictionary with internal_energy and specific_volume grids
            
        Returns:
            Dictionary with u_idx and v_idx
        """
        u_idx, v_idx, _, _ = grid_point
        return {"u_idx": int(u_idx), "v_idx": int(v_idx)}
    
    def _prepare_grid_info(self, grids: Dict[str, np.ndarray], options: Dict[str, Any], 
                          result_count: int) -> Dict[str, Any]:
        """
        Prepare grid information for UV flash.
        
        Args:
            grids: Dictionary with internal_energy and specific_volume grids
            options: Additional options
            result_count: Number of successful results
            
        Returns:
            Dictionary with grid information
        """
        return {
            "type": options.get("grid_type", "equidistant"),
            "internal_energy_points": len(grids["internal_energy"]),
            "specific_volume_points": len(grids["specific_volume"]),
            "total_points": result_count
        }