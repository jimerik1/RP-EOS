from flask import request, jsonify
import numpy as np
import sys
import traceback
from typing import List, Dict, Any, Tuple

from API.endpoints import extended_pt_flash_bp
from API.refprop_setup import RP
from API.unit_converter import UnitConverter
from API.utils.helpers import get_phase, validate_composition

def setup_mixture(composition: List[Dict[str, Any]]) -> List[float]:
    """Setup REFPROP mixture"""
    fluid_string = '|'.join(f"{comp['fluid']}.FLD" for comp in composition)
    z = [comp['fraction'] for comp in composition] + [0] * (20 - len(composition))
    
    ierr, herr = RP.SETUPdll(len(composition), fluid_string, 'HMX.BNC', 'DEF')
    if ierr > 0:
        raise ValueError(f"Error setting up mixture: {herr}")
    return z

def is_below_triple_point(T: float, P: float, z: List[float], pure_component_idx: int) -> bool:
    """
    Check if the state point is below the triple point temperature.
    
    Args:
        T: Temperature in K
        P: Pressure in kPa
        z: Composition array
        pure_component_idx: Index of the pure component
        
    Returns:
        True if T is below triple point temperature, False otherwise
    """
    try:
        if pure_component_idx > 0:
            # Get triple point temperature
            wmm, Ttrp, Tnbpt, Tc, Pc, Dc, Zc, acf, dip, Rgas = RP.INFOdll(pure_component_idx)
            return T < Ttrp
    except Exception:
        pass
    
    return False

def is_in_solid_phase(T: float, P: float, z: List[float], pure_component_idx: int) -> bool:
    """
    Check if the state point is in the solid phase.
    
    Args:
        T: Temperature in K
        P: Pressure in kPa
        z: Composition array
        pure_component_idx: Index of the pure component
        
    Returns:
        True if state is in solid phase, False otherwise
    """
    if pure_component_idx <= 0:
        return False
        
    try:
        # Get triple point temperature
        wmm, Ttrp, Tnbpt, Tc, Pc, Dc, Zc, acf, dip, Rgas = RP.INFOdll(pure_component_idx)
        
        if T >= Ttrp:
            # Above triple point, check if it's below melting line
            try:
                P_melt, ierr, herr = RP.MELTTdll(T, z)
                return ierr == 0 and P < P_melt
            except Exception:
                return False
        else:
            # Below triple point, check if it's below sublimation line
            try:
                P_subl, ierr, herr = RP.SUBLTdll(T, z)
                return ierr == 0 and P < P_subl
            except Exception:
                # If we're below triple point and can't determine sublimation,
                # conservatively assume it's solid
                return True
    except Exception:
        return False

def calculate_properties_extended(
    z: List[float], T: float, P: float, 
    units_system: str = 'SI', 
    is_pure: bool = False, 
    pure_component_idx: int = 0
) -> Dict[str, Any]:
    """
    Calculate fluid properties at given temperature and pressure with solid phase support.
    
    Args:
        z: Composition array
        T: Temperature in K
        P: Pressure in kPa
        units_system: Unit system to use ('SI' or 'CGS')
        is_pure: Whether this is a pure fluid
        pure_component_idx: Index of the pure component
        
    Returns:
        Dictionary of calculated properties with values and units
    """
    # Initialize unit converter
    converter = UnitConverter()
    
    # Get molecular weight for unit conversions
    wmm = RP.WMOLdll(z)
    
    # Special handling for solid phase (pure fluids only)
    is_solid = False
    at_melting = False
    at_sublimation = False
    
    if is_pure:
        is_solid = is_in_solid_phase(T, P, z, pure_component_idx)
    
    # Initialize properties dictionary
    raw_properties = {}
    
    if is_solid:
        # For solid phase, use estimated properties
        # Get estimated density for the solid
        try:
            # Try to use liquid density as an approximation (slightly higher)
            # Get triple point temperature
            wmm, Ttrp, Tnbpt, Tc, Pc, Dc, Zc, acf, dip, Rgas = RP.INFOdll(pure_component_idx)
            
            # Get density near triple point as reference
            P_triple, ierr, herr = RP.SUBLTdll(Ttrp, z)
            if ierr == 0:
                # Use higher density than liquid at triple point
                D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr2, herr2 = RP.TPFLSHdll(Ttrp, P_triple, z)
                
                if ierr2 == 0:
                    solid_D = Dl * 1.1  # Solid typically 5-15% denser than liquid
                else:
                    # Fallback - use critical density as reference
                    Tc, Pc, Dc, ierr3, herr3 = RP.CRITPdll(z)
                    solid_D = Dc * 3.0  # Rough estimate
            else:
                # Fallback - use critical density as reference
                Tc, Pc, Dc, ierr3, herr3 = RP.CRITPdll(z)
                solid_D = Dc * 3.0  # Rough estimate
        except Exception:
            # Another fallback method
            # Try to calculate a reference state and use that as an approximation
            try:
                D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr = RP.TPFLSHdll(T+50, P, z)
                solid_D = D * 1.2  # Estimation
            except Exception:
                # Final fallback - use a reasonable value
                solid_D = 20.0  # mol/L (reasonable for many substances)
        
        # For enthalpy, entropy and other properties, use liquid as approximation
        # but adjust for heat of fusion
        try:
            # Calculate properties near, but above, melting point
            try:
                P_melt, ierr, herr = RP.MELTTdll(T+1, z)
                if ierr != 0:
                    P_melt = P
            except Exception:
                P_melt = P
                
            D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr = RP.TPFLSHdll(T+1, P_melt, z)
            
            if ierr == 0:
                # Adjust enthalpy for heat of fusion (typically 5-15% of total)
                # Heat of fusion typically lowers the enthalpy of the solid
                solid_enthalpy_adjustment = 0.9  # Factor to adjust liquid enthalpy to solid
                solid_entropy_adjustment = 0.9   # Similar adjustment for entropy
                solid_cv_adjustment = 0.8        # Solid heat capacity typically lower
                
                # Set the properties with adjustments
                raw_properties = {
                    'density': solid_D,
                    'liquid_density': None,
                    'vapor_density': None,
                    'vapor_fraction': -999,  # Flag for solid
                    'internal_energy': e * solid_enthalpy_adjustment,
                    'enthalpy': h * solid_enthalpy_adjustment,
                    'entropy': s * solid_entropy_adjustment,
                    'cv': Cv * solid_cv_adjustment,
                    'cp': Cp * solid_cv_adjustment,
                    'sound_speed': w * 1.5,  # Sound travels faster in solids
                    'viscosity': None,  # No viscosity for solids
                    'thermal_conductivity': None,  # Could estimate but less reliable
                    'surface_tension': None,
                    'temperature': T - 273.15,  # Convert to Celsius
                    'pressure': P / 100  # Convert kPa to bar
                }
            else:
                # Fallback with minimal properties
                raw_properties = {
                    'density': solid_D,
                    'vapor_fraction': -999,  # Flag for solid
                    'temperature': T - 273.15,
                    'pressure': P / 100
                }
        except Exception as e:
            # Final fallback with minimal properties
            raw_properties = {
                'density': solid_D,
                'vapor_fraction': -999,  # Flag for solid
                'temperature': T - 273.15,
                'pressure': P / 100
            }
    else:
        # Not in solid phase, use normal REFPROP calculations
        try:
            # Get basic thermodynamic properties
            D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr = RP.TPFLSHdll(T, P, z)
            
            if ierr > 0:
                # Error in calculation, return minimal set with error
                return {
                    'temperature': converter.convert_property(
                        'temperature', float(T - 273.15), wmm, 'SI', units_system
                    ),
                    'pressure': converter.convert_property(
                        'pressure', float(P / 100), wmm, 'SI', units_system
                    ),
                    'calculation_error': {
                        'value': herr,
                        'unit': None
                    }
                }
            
            # Check if we're at a phase boundary
            if is_pure:
                try:
                    # Check if on melting line
                    P_melt, ierr_melt, herr_melt = RP.MELTTdll(T, z)
                    if ierr_melt == 0 and abs(P - P_melt) / P < 0.01:  # Within 1%
                        at_melting = True
                        q = -997  # Special flag for solid-liquid boundary
                except Exception:
                    pass
                    
                try:
                    # Check if on sublimation line
                    P_subl, ierr_subl, herr_subl = RP.SUBLTdll(T, z)
                    if ierr_subl == 0 and abs(P - P_subl) / P < 0.01:  # Within 1%
                        at_sublimation = True
                        q = 997  # Special flag for solid-vapor boundary
                except Exception:
                    pass
                    
                # Check for triple point (both melting and sublimation curves match)
                if at_melting and at_sublimation:
                    q = 996  # Special flag for triple point
            
            # Get transport properties
            try:
                eta, tcx, ierr_trn, herr_trn = RP.TRNPRPdll(T, D, z)
                if ierr_trn > 0:
                    eta = tcx = None
            except Exception:
                eta = tcx = None
                
            # Get surface tension if in two-phase region
            surface_tension = None
            if 0 < q < 1:
                try:
                    sigma, ierr_st, herr_st = RP.SURTENdll(T, Dl, Dv, x, y)
                    if ierr_st == 0:
                        surface_tension = float(sigma)
                except Exception:
                    pass
                    
            # Get critical properties
            try:
                Tc, Pc, Dc, ierr_crit, herr_crit = RP.CRITPdll(z)
                if ierr_crit > 0:
                    Tc = Pc = Dc = None
            except Exception:
                Tc = Pc = Dc = None
                
            # Get additional thermodynamic derivatives
            try:
                dPdD, dPdT, d2PdD2, d2PdT2, d2PdTD, dDdP, dDdT, d2DdP2, d2DdT2, d2DdPT, \
                dTdP, dTdD, d2TdP2, d2TdD2, d2TdPD = RP.DERVPVTdll(T, D, z)
            except Exception:
                dPdD = dPdT = dDdP = dDdT = None
            
            # Calculate compressibility factor
            Z = P / (D * 8.31446261815324 * T)
            
            # Build raw properties dictionary
            raw_properties = {
                'density': D,
                'liquid_density': Dl,
                'vapor_density': Dv,
                'vapor_fraction': q,
                'internal_energy': e,
                'enthalpy': h,
                'entropy': s,
                'cv': Cv,
                'cp': Cp,
                'sound_speed': w,
                'viscosity': eta,
                'thermal_conductivity': tcx,
                'surface_tension': surface_tension,
                'critical_temperature': Tc,
                'critical_pressure': Pc / 100,  # Convert from kPa to bar
                'critical_density': Dc,
                'compressibility_factor': Z,
                'isothermal_compressibility': -1/D * dDdP if dDdP is not None else None,
                'volume_expansivity': 1/D * dDdT if dDdT is not None else None,
                'dp_dt_saturation': dPdT,
                'joule_thomson_coefficient': (T*dDdT/dDdP - 1)/(Cp * 100) if (dDdT is not None and dDdP is not None and Cp is not None) else None,
                'temperature': T - 273.15,  # Convert to Celsius
                'pressure': P / 100  # Convert kPa to bar
            }
        except Exception as e:
            # If calculation fails completely, return error with minimal set
            return {
                'temperature': converter.convert_property(
                    'temperature', float(T - 273.15), wmm, 'SI', units_system
                ),
                'pressure': converter.convert_property(
                    'pressure', float(P / 100), wmm, 'SI', units_system
                ),
                'calculation_error': {
                    'value': str(e),
                    'unit': None
                }
            }
    
    # Convert properties to requested unit system
    properties = {}
    for prop_id, value in raw_properties.items():
        if value is not None:  # Skip undefined properties
            try:
                properties[prop_id] = converter.convert_property(
                    prop_id, float(value), wmm, 'SI', units_system
                )
            except Exception as e:
                print(f"Warning: Could not convert property {prop_id}: {str(e)}")
                properties[prop_id] = {'value': float(value) if value is not None else None, 'unit': 'unknown'}
    
    # Add phase information
    properties['phase'] = {
        'value': get_phase(raw_properties.get('vapor_fraction')),
        'unit': None
    }
    
    return properties

@extended_pt_flash_bp.route('/extended_pt_flash', methods=['POST'])
def extended_pt_flash():
    try:
        data = request.get_json(force=True)
        
        # Validate the request structure
        required_fields = ['composition', 'variables']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing field: {field}'}), 400
                
        # Check if required variables exist
        variables = data.get('variables', {})
        if 'pressure' not in variables or 'temperature' not in variables:
            return jsonify({'error': 'Missing pressure or temperature variables'}), 400

        # Extract calculation settings
        calculation = data.get('calculation', {})
        properties = calculation.get('properties', [])
        units_system = calculation.get('units_system', 'SI')  # Default to SI
        include_solid_phase = calculation.get('include_solid_phase', True)  # Default to include solid phase
        
        if not properties:
            return jsonify({'error': 'No properties specified for calculation'}), 400

        # Validate composition
        if not validate_composition(data['composition']):
            return jsonify({'error': 'Invalid composition - fractions must sum to 1'}), 400

        # Setup mixture
        z = setup_mixture(data['composition'])

        # Determine if it's a pure fluid
        is_pure = sum(1 for component in z if component > 0) == 1
        pure_component_idx = next((i+1 for i, component in enumerate(z) if component > 0), 0) if is_pure else 0

        # Extract range and resolution parameters
        pressure_range = variables['pressure'].get('range', {})
        temperature_range = variables['temperature'].get('range', {})
        pressure_resolution = variables['pressure'].get('resolution')
        temperature_resolution = variables['temperature'].get('resolution')
        
        # Validate required parameters exist
        if not all([pressure_range.get('from'), pressure_range.get('to'), 
                   temperature_range.get('from'), temperature_range.get('to'),
                   pressure_resolution, temperature_resolution]):
            return jsonify({'error': 'Missing range or resolution parameters'}), 400

        # Debug log
        print(f"Calculating extended PT flash for temperature range: {temperature_range['from']} to {temperature_range['to']} °C, "
              f"pressure range: {pressure_range['from']} to {pressure_range['to']} bar")

        # Create arrays for calculations
        T_range = np.arange(
            float(temperature_range['from']) + 273.15,
            float(temperature_range['to']) + 273.15 + float(temperature_resolution),
            float(temperature_resolution)
        )
        P_range = np.arange(
            float(pressure_range['from']) * 100,  # Convert bar to kPa
            float(pressure_range['to']) * 100 + float(pressure_resolution) * 100,
            float(pressure_resolution) * 100
        )

        # Calculate properties
        results = []
        idx = 0
        for T in T_range:
            for P in P_range:
                try:
                    props = calculate_properties_extended(
                        z, float(T), float(P), units_system, is_pure, pure_component_idx
                    )
                    
                    # Filter properties
                    filtered_props = {k: v for k, v in props.items() 
                                   if k in properties or k in ['temperature', 'pressure', 'phase']}
                    
                    # Skip points with calculation errors
                    if 'calculation_error' in props and len(filtered_props) <= 3:  # Only T, P, and error
                        continue
                        
                    # Add result with index
                    results.append({
                        'index': idx,
                        **filtered_props
                    })
                    idx += 1
                except Exception as fe:
                    print(f"Error processing T={T-273.15}°C, P={P/100} bar: {fe}", file=sys.stderr)
                    continue

        return jsonify({'results': results})
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500