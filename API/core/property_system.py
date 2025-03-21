"""
Property registry system for REFPROP API.

This module defines a central registry for all properties available in the API,
including their metadata, dependencies, and calculation methods.
"""

from typing import Dict, List, Any, Callable, Optional, Union, Tuple
import numpy as np
import logging

logger = logging.getLogger(__name__)

# Define property calculation function type
PropertyCalculator = Callable[[Dict[str, Any], Any, List[float]], Any]

class PropertyRegistry:
    """
    Central registry for all REFPROP properties with metadata.
    
    This class:
    1. Maintains a registry of all available properties
    2. Defines calculation methods for derived properties
    3. Handles property calculation with proper dependency management
    4. Provides unit information for different unit systems
    """
    
    def __init__(self):
        """Initialize the property registry and register all available properties."""
        self.properties = {}
        self.register_all_properties()
    
    def register(self, name: str, metadata: Dict[str, Any]) -> None:
        """
        Register a property with its metadata and calculation method.
        
        Args:
            name: Property identifier string
            metadata: Dictionary containing property metadata:
                - description: Human-readable description
                - si_unit: Unit in SI system
                - cgs_unit: Unit in CGS system
                - calculation_method: Function to calculate the property (None for base properties)
                - dependencies: List of other properties this property depends on
                - phase_specific: Whether this property has phase-specific variants
                - aliases: Alternative names for this property
        """
        self.properties[name] = metadata
        
        # Also register any aliases
        for alias in metadata.get('aliases', []):
            # Don't overwrite existing properties with aliases
            if alias not in self.properties:
                alias_metadata = metadata.copy()
                alias_metadata['is_alias'] = True
                alias_metadata['alias_of'] = name
                self.properties[alias] = alias_metadata
    
    def register_all_properties(self) -> None:
        """Register all available properties with their metadata."""
        # ========================================================
        # I. BASIC THERMODYNAMIC PROPERTIES
        # ========================================================
        
        # 1. Base properties (directly from flash calculations)
        self.register("temperature", {
            "description": "Temperature",
            "si_unit": "°C", 
            "cgs_unit": "°C",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["T"],
        })
        
        self.register("pressure", {
            "description": "Pressure",
            "si_unit": "bar", 
            "cgs_unit": "dyn/cm²",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["P"],
        })
        
        self.register("density", {
            "description": "Bulk density",
            "si_unit": "mol/L", 
            "cgs_unit": "g/cm³",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": True,
            "aliases": ["D", "rho"],
        })
        
        self.register("liquid_density", {
            "description": "Liquid phase density",
            "si_unit": "mol/L", 
            "cgs_unit": "g/cm³",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["Dl", "rhol"],
        })
        
        self.register("vapor_density", {
            "description": "Vapor phase density",
            "si_unit": "mol/L", 
            "cgs_unit": "g/cm³",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["Dv", "rhov"],
        })
        
        self.register("vapor_fraction", {
            "description": "Vapor quality/fraction",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["quality", "q"],
        })
        
        self.register("specific_volume", {
            "description": "Specific volume",
            "si_unit": "m³/mol", 
            "cgs_unit": "cm³/g",
            "calculation_method": self._calculate_specific_volume,
            "dependencies": ["density"],
            "phase_specific": True,
            "aliases": ["v"],
        })
        
        self.register("internal_energy", {
            "description": "Specific internal energy",
            "si_unit": "J/mol", 
            "cgs_unit": "erg/g",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": True,
            "aliases": ["energy", "e"],
        })
        
        self.register("enthalpy", {
            "description": "Specific enthalpy",
            "si_unit": "J/mol", 
            "cgs_unit": "erg/g",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": True,
            "aliases": ["h"],
        })
        
        self.register("entropy", {
            "description": "Specific entropy",
            "si_unit": "J/(mol·K)", 
            "cgs_unit": "erg/(g·K)",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": True,
            "aliases": ["s"],
        })
        
        self.register("cv", {
            "description": "Isochoric specific heat capacity",
            "si_unit": "J/(mol·K)", 
            "cgs_unit": "erg/(g·K)",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": True,
            "aliases": ["isochoric_heat_capacity"],
        })
        
        self.register("cp", {
            "description": "Isobaric specific heat capacity",
            "si_unit": "J/(mol·K)", 
            "cgs_unit": "erg/(g·K)",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": True,
            "aliases": ["isobaric_heat_capacity"],
        })
        
        self.register("sound_speed", {
            "description": "Speed of sound",
            "si_unit": "m/s", 
            "cgs_unit": "cm/s",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": True,
            "aliases": ["w"],
        })
        
        # 2. Energy-related properties
        self.register("helmholtz_energy", {
            "description": "Specific Helmholtz free energy",
            "si_unit": "J/mol", 
            "cgs_unit": "erg/g",
            "calculation_method": self._calculate_helmholtz_energy,
            "dependencies": ["internal_energy", "temperature", "entropy"],
            "phase_specific": True,
            "aliases": ["a"],
        })
        
        self.register("gibbs_energy", {
            "description": "Specific Gibbs free energy",
            "si_unit": "J/mol", 
            "cgs_unit": "erg/g",
            "calculation_method": self._calculate_gibbs_energy,
            "dependencies": ["enthalpy", "temperature", "entropy"],
            "phase_specific": True,
            "aliases": ["g"],
        })
        
        self.register("heat_of_vaporization", {
            "description": "Heat of vaporization",
            "si_unit": "J/mol", 
            "cgs_unit": "erg/g",
            "calculation_method": self._calculate_heat_of_vaporization,
            "dependencies": ["vapor_enthalpy", "liquid_enthalpy", "vapor_fraction"],
            "phase_specific": False,
            "aliases": ["hvap", "latent_heat"],
        })
        
        # 3. Compositional properties (for mixtures)
        self.register("chemical_potential", {
            "description": "Chemical potential",
            "si_unit": "J/mol", 
            "cgs_unit": "erg/g",
            "calculation_method": self._calculate_chemical_potential,
            "dependencies": ["temperature", "pressure", "density", "x", "y"],
            "phase_specific": True,
            "aliases": ["mu"],
        })
        
        self.register("fugacity", {
            "description": "Fugacity",
            "si_unit": "bar", 
            "cgs_unit": "dyn/cm²",
            "calculation_method": self._calculate_fugacity,
            "dependencies": ["temperature", "pressure", "density", "x", "y"],
            "phase_specific": True,
            "aliases": ["f"],
        })
        
        self.register("fugacity_coefficient", {
            "description": "Fugacity coefficient",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": self._calculate_fugacity_coefficient,
            "dependencies": ["fugacity", "pressure"],
            "phase_specific": True,
            "aliases": ["phi"],
        })
        
        self.register("k_value", {
            "description": "K value (vapor/liquid distribution ratio)",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": self._calculate_k_value,
            "dependencies": ["x", "y"],
            "phase_specific": False,
            "aliases": ["K"],
        })
        
        self.register("molar_mass", {
            "description": "Molar mass",
            "si_unit": "g/mol", 
            "cgs_unit": "g/mol",
            "calculation_method": self._calculate_molar_mass,
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["molecular_weight", "wmm"],
        })
        
        # ========================================================
        # II. TRANSPORT PROPERTIES
        # ========================================================
        
        # 1. Viscosity and related properties
        self.register("viscosity", {
            "description": "Dynamic viscosity",
            "si_unit": "μPa·s", 
            "cgs_unit": "poise",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": True,
            "aliases": ["dynamic_viscosity", "eta"],
        })
        
        self.register("kinematic_viscosity", {
            "description": "Kinematic viscosity",
            "si_unit": "cm²/s", 
            "cgs_unit": "cm²/s",
            "calculation_method": self._calculate_kinematic_viscosity,
            "dependencies": ["viscosity", "density"],
            "phase_specific": True,
            "aliases": ["nu"],
        })
        
        # 2. Thermal transport properties
        self.register("thermal_conductivity", {
            "description": "Thermal conductivity",
            "si_unit": "W/(m·K)", 
            "cgs_unit": "cal/(s·cm·K)",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": True,
            "aliases": ["tcx"],
        })
        
        self.register("thermal_diffusivity", {
            "description": "Thermal diffusivity",
            "si_unit": "cm²/s", 
            "cgs_unit": "cm²/s",
            "calculation_method": self._calculate_thermal_diffusivity,
            "dependencies": ["thermal_conductivity", "density", "cp"],
            "phase_specific": True,
            "aliases": ["alpha"],
        })
        
        # 3. Dimensionless transport groups
        self.register("prandtl_number", {
            "description": "Prandtl number",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": self._calculate_prandtl,
            "dependencies": ["cp", "viscosity", "thermal_conductivity"],
            "phase_specific": True,
            "aliases": ["Pr"],
        })
        
        # 4. Interface properties
        self.register("surface_tension", {
            "description": "Surface tension (two-phase only)",
            "si_unit": "N/m", 
            "cgs_unit": "dyn/cm",
            "calculation_method": None,  # Base property
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["sigma"],
        })
        
        # 5. Electrical properties
        self.register("dielectric_constant", {
            "description": "Dielectric constant",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": self._calculate_dielectric_constant,
            "dependencies": ["temperature", "density"],
            "phase_specific": True,
            "aliases": ["epsilon"],
        })
        
        # ========================================================
        # III. DERIVATIVE PROPERTIES
        # ========================================================
        
        # 1. Thermodynamic derivatives
        self.register("dDdP", {
            "description": "Pressure derivative of density",
            "si_unit": "(mol/L)/kPa", 
            "cgs_unit": "(g/cm³)/(dyn/cm²)",
            "calculation_method": None,  # Usually calculated by REFPROP directly
            "dependencies": [],
            "phase_specific": True,
            "aliases": [],
        })
        
        self.register("dDdT", {
            "description": "Temperature derivative of density",
            "si_unit": "(mol/L)/K", 
            "cgs_unit": "(g/cm³)/K",
            "calculation_method": None,  # Usually calculated by REFPROP directly
            "dependencies": [],
            "phase_specific": True,
            "aliases": [],
        })
        
        self.register("dPdT", {
            "description": "Temperature derivative of pressure",
            "si_unit": "kPa/K", 
            "cgs_unit": "(dyn/cm²)/K",
            "calculation_method": None,  # Usually calculated by REFPROP directly
            "dependencies": [],
            "phase_specific": True,
            "aliases": [],
        })
        
        self.register("dPdD", {
            "description": "Density derivative of pressure",
            "si_unit": "kPa/(mol/L)", 
            "cgs_unit": "(dyn/cm²)/(g/cm³)",
            "calculation_method": None,  # Usually calculated by REFPROP directly
            "dependencies": [],
            "phase_specific": True,
            "aliases": [],
        })
        
        self.register("d2PdD2", {
            "description": "Second density derivative of pressure",
            "si_unit": "kPa/(mol/L)²", 
            "cgs_unit": "(dyn/cm²)/(g/cm³)²",
            "calculation_method": None,  # Usually calculated by REFPROP directly
            "dependencies": [],
            "phase_specific": True,
            "aliases": [],
        })
        
        # 2. Derived thermodynamic properties
        self.register("isothermal_compressibility", {
            "description": "Isothermal compressibility",
            "si_unit": "1/kPa", 
            "cgs_unit": "cm²/dyn",
            "calculation_method": self._calculate_isothermal_compressibility,
            "dependencies": ["density", "dDdP"],
            "phase_specific": True,
            "aliases": ["kappa_T", "beta_T"],
        })
        
        self.register("volume_expansivity", {
            "description": "Volume expansivity",
            "si_unit": "1/K", 
            "cgs_unit": "1/K",
            "calculation_method": self._calculate_volume_expansivity,
            "dependencies": ["density", "dDdT"],
            "phase_specific": True,
            "aliases": ["alpha_V"],
        })
        
        self.register("isentropic_coefficient", {
            "description": "Isentropic coefficient (ratio of specific heats)",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": self._calculate_isentropic_coefficient,
            "dependencies": ["cp", "cv"],
            "phase_specific": True,
            "aliases": ["gamma", "kappa"],
        })
        
        self.register("adiabatic_compressibility", {
            "description": "Adiabatic compressibility",
            "si_unit": "1/kPa", 
            "cgs_unit": "cm²/dyn",
            "calculation_method": self._calculate_adiabatic_compressibility,
            "dependencies": ["isothermal_compressibility", "isentropic_coefficient"],
            "phase_specific": True,
            "aliases": ["kappa_S", "beta_S"],
        })
        
        self.register("joule_thomson_coefficient", {
            "description": "Joule-Thomson coefficient",
            "si_unit": "K/bar", 
            "cgs_unit": "K·cm²/dyn",
            "calculation_method": self._calculate_joule_thomson,
            "dependencies": ["temperature", "density", "cp", "dDdT", "dDdP"],
            "phase_specific": True,
            "aliases": ["mu_JT"],
        })
        
        # 3. Advanced thermodynamic properties
        self.register("specific_heat_input", {
            "description": "Specific heat input",
            "si_unit": "J/(mol·K)", 
            "cgs_unit": "erg/(g·K)",
            "calculation_method": self._calculate_specific_heat_input,
            "dependencies": ["cp", "temperature", "pressure", "volume_expansivity", "isothermal_compressibility"],
            "phase_specific": True,
            "aliases": ["cstar"],
        })
        
        self.register("gruneisen_parameter", {
            "description": "Grüneisen parameter",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": self._calculate_gruneisen,
            "dependencies": ["volume_expansivity", "isothermal_compressibility", "cp", "density"],
            "phase_specific": True,
            "aliases": ["gruneisen", "gamma_G"],
        })
        
        self.register("critical_flow_factor", {
            "description": "Critical flow factor",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": self._calculate_critical_flow_factor,
            "dependencies": ["cp", "cv", "temperature", "pressure", "density"],
            "phase_specific": True,
            "aliases": ["cff"],
        })
        
        # 4. Virial coefficients
        self.register("second_virial_coefficient", {
            "description": "Second virial coefficient",
            "si_unit": "m³/mol", 
            "cgs_unit": "cm³/g",
            "calculation_method": self._calculate_second_virial,
            "dependencies": ["temperature"],
            "phase_specific": False,
            "aliases": ["B", "B_virial"],
        })
        
        self.register("third_virial_coefficient", {
            "description": "Third virial coefficient",
            "si_unit": "(m³/mol)²", 
            "cgs_unit": "(cm³/g)²",
            "calculation_method": self._calculate_third_virial,
            "dependencies": ["temperature"],
            "phase_specific": False,
            "aliases": ["C", "C_virial"],
        })
        
        self.register("second_acoustic_virial_coefficient", {
            "description": "Second acoustic virial coefficient",
            "si_unit": "m³/mol", 
            "cgs_unit": "cm³/g",
            "calculation_method": self._calculate_second_acoustic_virial,
            "dependencies": ["temperature"],
            "phase_specific": False,
            "aliases": ["Ba", "B_acoustic"],
        })
        
        self.register("third_acoustic_virial_coefficient", {
            "description": "Third acoustic virial coefficient",
            "si_unit": "(m³/mol)²", 
            "cgs_unit": "(cm³/g)²",
            "calculation_method": self._calculate_third_acoustic_virial,
            "dependencies": ["temperature"],
            "phase_specific": False,
            "aliases": ["Ca", "C_acoustic"],
        })
        
        # Binary interaction parameters
        self.register("binary_interaction_parameter", {
            "description": "Binary interaction parameter between components",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": self._calculate_binary_interaction,
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["B12", "k12", "BIP"],
        })
        
        # ========================================================
        # IV. COMBUSTION PROPERTIES
        # ========================================================
        
        self.register("gross_heating_value", {
            "description": "Gross heating value (higher heating value)",
            "si_unit": "J/mol", 
            "cgs_unit": "erg/g",
            "calculation_method": self._calculate_gross_heating_value,
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["HHV", "GHV"],
        })
        
        self.register("net_heating_value", {
            "description": "Net heating value (lower heating value)",
            "si_unit": "J/mol", 
            "cgs_unit": "erg/g",
            "calculation_method": self._calculate_net_heating_value,
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["LHV", "NHV"],
        })
        
        # ========================================================
        # V. STATE INFORMATION
        # ========================================================
        
        # 1. Phase information
        self.register("phase", {
            "description": "Phase state (Liquid, Vapor, etc.)",
            "si_unit": "text", 
            "cgs_unit": "text",
            "calculation_method": self._calculate_phase,
            "dependencies": ["vapor_fraction"],
            "phase_specific": False,
            "aliases": [],
        })
        
        # 2. State indicators
        self.register("compressibility_factor", {
            "description": "Compressibility factor (Z)",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": self._calculate_z_factor,
            "dependencies": ["pressure", "density", "temperature"],
            "phase_specific": True,
            "aliases": ["Z"],
        })
        
        # 3. Excess properties (for mixtures)
        self.register("excess_volume", {
            "description": "Excess volume of mixing",
            "si_unit": "m³/mol", 
            "cgs_unit": "cm³/g",
            "calculation_method": self._calculate_excess_volume,
            "dependencies": ["temperature", "pressure", "density", "x", "y"],
            "phase_specific": True,
            "aliases": ["v_excess"],
        })
        
        self.register("excess_enthalpy", {
            "description": "Excess enthalpy of mixing",
            "si_unit": "J/mol", 
            "cgs_unit": "erg/g",
            "calculation_method": self._calculate_excess_enthalpy,
            "dependencies": ["temperature", "pressure", "density", "x", "y"],
            "phase_specific": True,
            "aliases": ["h_excess"],
        })
        
        # ========================================================
        # VI. CRITICAL POINT PROPERTIES
        # ========================================================
        
        self.register("critical_temperature", {
            "description": "Critical temperature",
            "si_unit": "K", 
            "cgs_unit": "K",
            "calculation_method": None,  # Usually provided by REFPROP directly
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["Tc"],
        })
        
        self.register("critical_pressure", {
            "description": "Critical pressure",
            "si_unit": "bar", 
            "cgs_unit": "dyn/cm²",
            "calculation_method": None,  # Usually provided by REFPROP directly
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["Pc"],
        })
        
        self.register("critical_density", {
            "description": "Critical density",
            "si_unit": "mol/L", 
            "cgs_unit": "g/cm³",
            "calculation_method": None,  # Usually provided by REFPROP directly
            "dependencies": [],
            "phase_specific": False,
            "aliases": ["Dc", "rhoc"],
        })
        
        self.register("critical_compressibility_factor", {
            "description": "Critical compressibility factor",
            "si_unit": "dimensionless", 
            "cgs_unit": "dimensionless",
            "calculation_method": self._calculate_critical_z,
            "dependencies": ["critical_pressure", "critical_temperature", "critical_density"],
            "phase_specific": False,
            "aliases": ["Zc"],
        })
        
        # ========================================================
        # VII. MISCELLANEOUS PROPERTIES
        # ========================================================
        
        self.register("exergy", {
            "description": "Specific exergy",
            "si_unit": "J/mol", 
            "cgs_unit": "erg/g",
            "calculation_method": self._calculate_exergy,
            "dependencies": ["enthalpy", "entropy", "temperature"],
            "phase_specific": True,
            "aliases": ["availability"],
        })
        
        # Phase-specific variants of base properties
        # These are created automatically based on base properties
        self._register_phase_specific_properties()
        
    def _register_phase_specific_properties(self) -> None:
        """Register liquid and vapor variants of phase-specific properties."""
        phase_specific_props = [
            name for name, metadata in self.properties.items()
            if metadata.get("phase_specific", False) and not name.startswith("liquid_") and not name.startswith("vapor_")
        ]
        
        for prop_name in phase_specific_props:
            # Skip already registered phase-specific properties
            if f"liquid_{prop_name}" in self.properties or f"vapor_{prop_name}" in self.properties:
                continue
                
            base_metadata = self.properties[prop_name]
            
            # Register liquid version
            liquid_metadata = base_metadata.copy()
            liquid_metadata["description"] = f"Liquid phase {base_metadata['description'].lower()}"
            if base_metadata["calculation_method"] is None:
                liquid_metadata["calculation_method"] = None  # Base property
            else:
                liquid_metadata["calculation_method"] = lambda base_props, rp, z, prop=prop_name: self._calculate_phase_specific(
                    base_props, rp, z, prop, "liquid"
                )
            self.register(f"liquid_{prop_name}", liquid_metadata)
            
            # Register vapor version
            vapor_metadata = base_metadata.copy()
            vapor_metadata["description"] = f"Vapor phase {base_metadata['description'].lower()}"
            if base_metadata["calculation_method"] is None:
                vapor_metadata["calculation_method"] = None  # Base property
            else:
                vapor_metadata["calculation_method"] = lambda base_props, rp, z, prop=prop_name: self._calculate_phase_specific(
                    base_props, rp, z, prop, "vapor"
                )
            self.register(f"vapor_{prop_name}", vapor_metadata)
    
    def _calculate_phase_specific(self, base_props: Dict[str, Any], rp: Any, z: List[float], 
                                 property_name: str, phase: str) -> Any:
        """
        Calculate a phase-specific property by delegating to the base property method.
        
        This is used for properties that can be calculated for specific phases.
        
        Args:
            base_props: Dictionary of base properties (from flash calculation)
            rp: REFPROP instance
            z: Composition array
            property_name: Name of the base property
            phase: Phase to calculate for ('liquid' or 'vapor')
        
        Returns:
            Calculated property value
        """
        # Create phase-specific property dict with the same keys but phase-specific values
        phase_specific_props = {}
        for k, v in base_props.items():
            prefix = f"{phase}_"
            if k.startswith(prefix):
                # This is a phase-specific property, add base name to props
                base_name = k[len(prefix):]
                phase_specific_props[base_name] = v
            elif k == "density":
                # Use phase density for base density
                phase_specific_props[k] = base_props.get(f"{phase}_density")
            else:
                # Keep other properties as is
                phase_specific_props[k] = v
        
        # Get the calculation method for the base property
        calculation_method = self.properties[property_name]["calculation_method"]
        if calculation_method is None:
            # If it's a base property, just return the corresponding phase-specific value
            return base_props.get(f"{phase}_{property_name}")
            
        # Calculate using the base property method but with phase-specific inputs
        return calculation_method(phase_specific_props, rp, z)
        
    def get_property_info(self, name: str) -> Dict[str, Any]:
        """
        Get information about a property.
        
        Args:
            name: Property name
            
        Returns:
            Dictionary with property information
        
        Raises:
            ValueError: If property not found
        """
        if name not in self.properties:
            raise ValueError(f"Unknown property: {name}")
            
        # If it's an alias, get the canonical name
        if self.properties[name].get('is_alias', False):
            canonical_name = self.properties[name]['alias_of']
            return self.properties[canonical_name]
            
        return self.properties[name]
        
    def get_property_unit(self, name: str, unit_system: str = 'SI') -> str:
        """
        Get the unit for a property in the specified unit system.
        
        Args:
            name: Property name
            unit_system: Unit system ('SI' or 'CGS')
            
        Returns:
            Unit string
        """
        try:
            prop_info = self.get_property_info(name)
            
            if unit_system.upper() == 'CGS':
                return prop_info.get('cgs_unit', 'dimensionless')
            else:
                return prop_info.get('si_unit', 'dimensionless')
        except ValueError:
            return 'unknown'
        
    def calculate_property(self, name: str, base_props: Dict[str, Any], 
                          rp: Any, z: List[float]) -> Any:
        """
        Calculate a property value.
        
        Args:
            name: Property name
            base_props: Dictionary of base properties
            rp: REFPROP instance
            z: Composition array
            
        Returns:
            Calculated property value
            
        Raises:
            ValueError: If property not found or calculation fails
        """
        try:
            # Get the property info
            if name not in self.properties:
                raise ValueError(f"Unknown property: {name}")
                
            prop_info = self.get_property_info(name)
            
            # If it's a base property, just return from base_props
            if prop_info["calculation_method"] is None:
                if name in base_props:
                    return base_props[name]
                return None
                
            # Check dependencies
            for dep in prop_info["dependencies"]:
                if dep not in base_props:
                    # Try to calculate the dependency
                    dep_value = self.calculate_property(dep, base_props, rp, z)
                    
                    if dep_value is None:
                        raise ValueError(f"Missing dependency {dep} for property {name}")
                    
                    # Add to base_props for future use
                    base_props[dep] = dep_value
                    
            # Calculate the property
            return prop_info["calculation_method"](base_props, rp, z)
            
        except Exception as e:
            logger.warning(f"Error calculating property {name}: {str(e)}")
            return None

    # ========================================================
    # BASIC PROPERTY CALCULATION METHODS
    # ========================================================
    
    def _calculate_specific_volume(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate specific volume (inverse of density)."""
        D = base_props["density"]  # mol/L
        if D <= 0:
            raise ValueError("Density must be positive to calculate specific volume")
        
        # Convert mol/L to m³/mol: 1 mol/L = 1 mol/0.001 m³ => v = 0.001/D m³/mol
        return 0.001 / D  # m³/mol
        
    def _calculate_phase(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> str:
        """Determine the phase state based on vapor fraction."""
        q = base_props.get("vapor_fraction")
        if q is None:
            return "Unknown"
            
        if q == 0:
            return "Liquid"
        elif q == 1:
            return "Vapor"
        elif q == -998:
            return "Liquid"  # Subcooled liquid
        elif q == 998:
            return "Vapor"   # Superheated vapor
        elif q == 999:
            return "Supercritical"
        elif q == -999:
            return "Solid"   # Added support for solid phase
        elif q == -997:
            return "Solid-Liquid"  # Solid-liquid equilibrium (melting)
        elif q == 997:
            return "Solid-Vapor"   # Solid-vapor equilibrium (sublimation)
        elif q == 996:
            return "Triple-Point"  # Triple point (solid-liquid-vapor)
        elif 0 < q < 1:
            return "Two-Phase"
        elif q < 0:
            return "Liquid"  # Subcooled liquid
        else:
            return "Vapor"   # Superheated vapor
            
    def _calculate_z_factor(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate compressibility factor Z = PV/RT."""
        P_kpa = base_props["pressure"] * 100  # bar to kPa
        T_K = base_props["temperature"] + 273.15  # °C to K
        D = base_props["density"]
        
        # Z = P/(ρRT) with R = 8.31446 J/(mol·K)
        return P_kpa / (D * 8.31446261815324 * T_K)
    
    # ========================================================
    # ENERGY PROPERTY CALCULATIONS
    # ========================================================
    
    def _calculate_helmholtz_energy(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate Helmholtz free energy A = U - TS."""
        U = base_props["internal_energy"]  # J/mol
        T_K = base_props["temperature"] + 273.15  # °C to K
        S = base_props["entropy"]  # J/(mol·K)
        
        return U - T_K * S  # J/mol
    
    def _calculate_gibbs_energy(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate Gibbs free energy G = H - TS."""
        H = base_props["enthalpy"]  # J/mol
        T_K = base_props["temperature"] + 273.15  # °C to K
        S = base_props["entropy"]  # J/(mol·K)
        
        return H - T_K * S  # J/mol
    
    def _calculate_heat_of_vaporization(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate heat of vaporization (latent heat)."""
        # Only applicable in two-phase region or at saturation
        q = base_props["vapor_fraction"]
        if q != 0 and q != 1 and not (0 < q < 1):
            # Not at saturation or in two-phase region
            return 0.0
            
        # Get liquid and vapor enthalpies
        H_v = base_props.get("vapor_enthalpy")
        H_l = base_props.get("liquid_enthalpy")
        
        if H_v is None or H_l is None:
            # Try to get from REFPROP directly if at saturation
            try:
                T_K = base_props["temperature"] + 273.15  # °C to K
                result = rp.SATTdll(T_K, z, 1)  # Get saturation properties at T
                if result.ierr == 0:
                    Dl = result.Dl
                    Dv = result.Dv
                    
                    # Get enthalpies at saturation
                    result_l = rp.THERMdll(T_K, Dl, z)
                    result_v = rp.THERMdll(T_K, Dv, z)
                    
                    H_l = result_l.h
                    H_v = result_v.h
                else:
                    return 0.0  # Could not calculate
            except Exception:
                return 0.0  # Could not calculate
                
        if H_v is not None and H_l is not None:
            return H_v - H_l  # J/mol
        
        return 0.0  # Could not calculate
    
    # ========================================================
    # COMPOSITIONAL PROPERTY CALCULATIONS
    # ========================================================
    
    def _calculate_molar_mass(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate the molar mass (molecular weight) of the mixture."""
        return rp.WMOLdll(z)  # g/mol
    
    def _calculate_chemical_potential(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> List[float]:
        """Calculate chemical potentials for each component."""
        T_K = base_props["temperature"] + 273.15  # °C to K
        P_kpa = base_props["pressure"] * 100  # bar to kPa
        
        try:
            # Try using THERM2 if available
            result = rp.THERM2dll(T_K, P_kpa, z)
            return list(result.mu)
        except Exception:
            # Use PRESSdll as fallback approach
            # Note: This is a simplified approximation
            G = self._calculate_gibbs_energy(base_props, rp, z)
            return [G] * len([c for c in z if c > 0])  # Same value for all components
    
    def _calculate_fugacity(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> List[float]:
        """Calculate fugacities for each component."""
        T_K = base_props["temperature"] + 273.15  # °C to K
        P_kpa = base_props["pressure"] * 100  # bar to kPa
        D = base_props["density"]
        
        try:
            # Try to calculate using REFPROP directly
            result = rp.FUGCOFdll(T_K, D, z)
            if result.ierr == 0:
                # Convert fugacity coefficients to fugacities
                # f_i = φ_i * x_i * P
                x_i = z  # For a homogeneous phase, use bulk composition
                return [phi * x * P_kpa / 100 for phi, x in zip(result.f, x_i) if x > 0]  # Convert kPa to bar
        except Exception:
            pass
            
        # Use a simplified approach as fallback
        return [base_props["pressure"] * z_i for z_i in z if z_i > 0]  # bar
    
    def _calculate_fugacity_coefficient(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> List[float]:
        """Calculate fugacity coefficients for each component."""
        # Try to get fugacities first
        if "fugacity" in base_props:
            fugacities = base_props["fugacity"]
            P = base_props["pressure"]  # bar
            
            # φ_i = f_i / (x_i * P)
            return [f / (z_i * P) if z_i > 0 else 0 for f, z_i in zip(fugacities, z) if z_i > 0]
            
        # Calculate directly using REFPROP
        T_K = base_props["temperature"] + 273.15  # °C to K
        D = base_props["density"]
        
        try:
            result = rp.FUGCOFdll(T_K, D, z)
            if result.ierr == 0:
                return [result.f[i] for i, z_i in enumerate(z) if z_i > 0]
        except Exception:
            pass
            
        # Default to ideal gas (φ = 1)
        return [1.0] * len([c for c in z if c > 0])
    
    def _calculate_k_value(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> List[float]:
        """Calculate K values (vapor/liquid distribution ratios)."""
        # Only relevant in two-phase region
        q = base_props["vapor_fraction"]
        if not (0 < q < 1):
            return [1.0] * len([c for c in z if c > 0])  # Default to 1.0 for single-phase
            
        # Get liquid and vapor compositions
        x = base_props.get("x")
        y = base_props.get("y")
        
        if x is None or y is None:
            return [1.0] * len([c for c in z if c > 0])
            
        # K = y/x for each component
        return [y_i / x_i if x_i > 0 else 0.0 for x_i, y_i in zip(x, y)]
    
    # ========================================================
    # TRANSPORT PROPERTY CALCULATIONS
    # ========================================================
    
    def _calculate_kinematic_viscosity(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate kinematic viscosity ν = η/ρ."""
        if "molar_mass" not in base_props:
            # Get molecular weight from REFPROP
            wmm = rp.WMOLdll(z)
        else:
            wmm = base_props["molar_mass"]
            
        eta = base_props["viscosity"]  # μPa·s
        D = base_props["density"]      # mol/L
        
        # Convert units: η[μPa·s]/ρ[kg/m³] = [m²/s] * 10⁴ = [cm²/s]
        # Where ρ[kg/m³] = ρ[mol/L] * M[g/mol] / 1000
        # Thus: ν = η * 1e-6 / (D * wmm / 1000) * 10⁴
        return (eta * 1e-6) / (D * wmm / 1000) * 10000
    
    def _calculate_thermal_diffusivity(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate thermal diffusivity α = k/(ρCp)."""
        if "molar_mass" not in base_props:
            # Get molecular weight from REFPROP
            wmm = rp.WMOLdll(z)
        else:
            wmm = base_props["molar_mass"]
            
        k = base_props["thermal_conductivity"]  # W/(m·K)
        D = base_props["density"]               # mol/L
        Cp = base_props["cp"]                   # J/(mol·K)
        
        # α = k/(ρCp) where:
        # k = [W/(m·K)] = [J/(s·m·K)]
        # ρ = [kg/m³] = D[mol/L] * M[g/mol] / 1000 [g/kg] = D * wmm / 1000
        # Cp' = [J/(kg·K)] = Cp[J/(mol·K)] * 1000 / wmm
        # Thus: α = k / (D * wmm / 1000 * Cp * 1000 / wmm) = k / (D * Cp) in m²/s
        # Convert to cm²/s: * 10⁴
        return k / (D * Cp) * 10000
    
    def _calculate_prandtl(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate Prandtl number Pr = Cpη/k."""
        if "molar_mass" not in base_props:
            # Get molecular weight from REFPROP
            wmm = rp.WMOLdll(z)
        else:
            wmm = base_props["molar_mass"]
            
        Cp = base_props["cp"]                   # J/(mol·K)
        eta = base_props["viscosity"]           # μPa·s
        k = base_props["thermal_conductivity"]  # W/(m·K)
        
        # Pr = Cp'η/k where:
        # Cp' = [J/(kg·K)] = Cp[J/(mol·K)] * 1000 / wmm
        # η = [Pa·s] = η[μPa·s] * 1e-6
        # k = [W/(m·K)] = [J/(s·m·K)]
        # Thus: Pr = Cp * 1000 / wmm * η * 1e-6 / k
        wmm_kg = wmm / 1000  # g/mol to kg/mol
        return (Cp / wmm_kg) * (eta * 1e-6) / k
    
    def _calculate_dielectric_constant(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate dielectric constant (relative permittivity)."""
        T_K = base_props["temperature"] + 273.15  # °C to K
        D = base_props["density"]
        
        try:
            # REFPROP function for dielectric constant
            result = rp.DIELECdll(T_K, D, z)
            if result.ierr == 0:
                return result.de
        except Exception:
            pass
            
        # Default to 1.0 (vacuum) if calculation fails
        return 1.0
    
    # ========================================================
    # DERIVATIVE PROPERTY CALCULATIONS
    # ========================================================
    
    def _calculate_isothermal_compressibility(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate isothermal compressibility κT = -1/V(∂V/∂P)T = -1/ρ(∂ρ/∂P)T."""
        D = base_props["density"]      # mol/L
        dDdP = base_props["dDdP"]      # (mol/L)/kPa
        
        # κT = -1/ρ(∂ρ/∂P)T = -1/D * dDdP in 1/kPa
        return -1 / D * dDdP
    
    def _calculate_volume_expansivity(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate volume expansivity αV = 1/V(∂V/∂T)P = -1/ρ(∂ρ/∂T)P."""
        D = base_props["density"]      # mol/L
        dDdT = base_props["dDdT"]      # (mol/L)/K
        
        # αV = -1/ρ(∂ρ/∂T)P = -1/D * dDdT in 1/K
        return -1 / D * dDdT
    
    def _calculate_isentropic_coefficient(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate isentropic coefficient (ratio of specific heats) γ = Cp/Cv."""
        Cp = base_props["cp"]  # J/(mol·K)
        Cv = base_props["cv"]  # J/(mol·K)
        
        if Cv <= 0:
            raise ValueError("Cv must be positive to calculate isentropic coefficient")
            
        return Cp / Cv  # dimensionless
    
    def _calculate_adiabatic_compressibility(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate adiabatic compressibility κS = κT/γ."""
        kappa_T = base_props["isothermal_compressibility"]  # 1/kPa
        gamma = base_props["isentropic_coefficient"]  # dimensionless
        
        if gamma <= 0:
            raise ValueError("Isentropic coefficient must be positive")
            
        return kappa_T / gamma  # 1/kPa
    
    def _calculate_joule_thomson(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate Joule-Thomson coefficient μ_JT = (∂T/∂P)_H = (T(∂ρ/∂T)_P/ρ²Cp - 1/ρCp)."""
        T_K = base_props["temperature"] + 273.15  # °C to K
        D = base_props["density"]
        Cp = base_props["cp"]
        dDdT = base_props["dDdT"]
        dDdP = base_props["dDdP"]
        
        # μ_JT = T * [(∂V/∂T)_P] / Cp
        # But we have dDdT = (∂ρ/∂T)_P and (∂V/∂T)_P = -1/ρ² * (∂ρ/∂T)_P
        # So: μ_JT = -T/(D²*Cp) * dDdT 
        # Or alternatively: μ_JT = (T*dDdT/dDdP - 1)/(Cp * 100) in K/bar
        return (T_K * dDdT / dDdP - 1) / (Cp * 100)
    
    def _calculate_specific_heat_input(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate specific heat input C* = Cp - T*αV²/κT."""
        Cp = base_props["cp"]  # J/(mol·K)
        T_K = base_props["temperature"] + 273.15  # °C to K
        alpha_V = base_props["volume_expansivity"]  # 1/K
        kappa_T = base_props["isothermal_compressibility"]  # 1/kPa
        
        # Unit conversion factor for kappa_T (1/kPa to 1/bar): 100
        return Cp - T_K * alpha_V**2 / (kappa_T * 100)  # J/(mol·K)
    
    def _calculate_gruneisen(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate Grüneisen parameter Γ = αV/(ρ*Cp*κT)."""
        alpha_V = base_props["volume_expansivity"]  # 1/K
        D = base_props["density"]  # mol/L
        Cp = base_props["cp"]  # J/(mol·K)
        kappa_T = base_props["isothermal_compressibility"]  # 1/kPa
        
        # Unit conversion factor for kappa_T (1/kPa to 1/bar): 100
        return alpha_V / (D * Cp * kappa_T * 100)  # dimensionless
    
    def _calculate_critical_flow_factor(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate critical flow factor."""
        gamma = self._calculate_isentropic_coefficient(base_props, rp, z)
        
        # Simple approximation for ideal gases
        # For real gases, a more complex calculation would be needed
        return (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)))**0.5
    
    # ========================================================
    # VIRIAL COEFFICIENT CALCULATIONS
    # ========================================================
    
    def _calculate_second_virial(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate second virial coefficient B."""
        T_K = base_props["temperature"] + 273.15  # °C to K
        
        try:
            # Use VIRdll for second virial coefficient
            result = rp.VIRdll(T_K, z)
            if result.ierr == 0:
                return result.b  # m³/mol
        except Exception:
            pass
            
        # Fallback: Try to estimate from equation of state
        try:
            # At very low density, Z ≈ 1 + B*ρ
            # So we can estimate B by calculating Z at two low densities
            D1 = 1e-6  # Very low density
            D2 = 2e-6  # Double the first density
            
            P1_kpa, e1, h1, s1, cv1, cp1, w1, hjt1 = rp.THERMdll(T_K, D1, z)
            P2_kpa, e2, h2, s2, cv2, cp2, w2, hjt2 = rp.THERMdll(T_K, D2, z)
            
            Z1 = P1_kpa / (D1 * 8.31446261815324 * T_K)
            Z2 = P2_kpa / (D2 * 8.31446261815324 * T_K)
            
            # B = (Z - 1)/ρ ≈ (Z1 - 1)/D1 ≈ (Z2 - 1)/D2
            B1 = (Z1 - 1) / D1 * 0.001  # Convert to m³/mol
            B2 = (Z2 - 1) / D2 * 0.001  # Convert to m³/mol
            
            # Use average of the two estimates
            return (B1 + B2) / 2
        except Exception:
            pass
            
        # Default fallback - very rough estimate
        # For most gases, B is negative and on the order of -0.1 L/mol at moderate temperatures
        return -0.0001  # m³/mol
    
    def _calculate_third_virial(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate third virial coefficient C."""
        T_K = base_props["temperature"] + 273.15  # °C to K
        
        try:
            # Use VIRdll for third virial coefficient
            result = rp.VIRdll(T_K, z)
            if result.ierr == 0:
                return result.c  # (m³/mol)²
        except Exception:
            pass
            
        # Default fallback - very rough estimate
        # For most gases, C is positive and much smaller than B²
        B = self._calculate_second_virial(base_props, rp, z)
        return abs(B * B) * 0.1  # (m³/mol)²
    
    def _calculate_second_acoustic_virial(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate second acoustic virial coefficient."""
        # In many cases, Ba ≈ B - T*dB/dT
        T_K = base_props["temperature"] + 273.15  # °C to K
        dT = 0.1  # Small temperature increment
        
        try:
            # Calculate B at T and T+dT
            B_T = self._calculate_second_virial(base_props, rp, z)
            
            # Create new base_props with slightly higher T
            new_props = base_props.copy()
            new_props["temperature"] = base_props["temperature"] + dT
            
            B_T_plus_dT = self._calculate_second_virial(new_props, rp, z)
            
            # Estimate dB/dT
            dBdT = (B_T_plus_dT - B_T) / dT
            
            # Ba ≈ B - T*dB/dT
            return B_T - T_K * dBdT
        except Exception:
            pass
            
        # Fallback: Just return B
        return self._calculate_second_virial(base_props, rp, z)
    
    def _calculate_third_acoustic_virial(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate third acoustic virial coefficient."""
        # Much more complex, using a simplified approach
        C = self._calculate_third_virial(base_props, rp, z)
        return C  # Simplified approximation
    
    def _calculate_binary_interaction(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> List[List[float]]:
        """Calculate binary interaction parameters between components."""
        # Only applicable for mixtures
        if sum(1 for comp in z if comp > 0) <= 1:
            return [[0.0]]  # No binary interactions for pure fluids
            
        try:
            # Try to get binary interaction parameters from REFPROP
            # Example implementation - would need actual REFPROP function
            # This is a placeholder for the actual implementation
            n_comps = sum(1 for comp in z if comp > 0)
            bips = [[0.0 for _ in range(n_comps)] for _ in range(n_comps)]
            
            # In a real implementation, would iterate through component pairs
            # and get their BIPs from REFPROP
            
            return bips
        except Exception:
            pass
            
        # Default: zero interaction parameters
        n_comps = sum(1 for comp in z if comp > 0)
        return [[0.0 for _ in range(n_comps)] for _ in range(n_comps)]
    
    # ========================================================
    # COMBUSTION PROPERTY CALCULATIONS
    # ========================================================
    
    def _calculate_gross_heating_value(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate gross heating value (higher heating value)."""
        try:
            # Check if REFPROP has a direct implementation
            # This is a placeholder - would need the actual REFPROP function
            pass
        except Exception:
            pass
            
        # Default approach: estimate based on composition
        # This is a very simplified approximation
        # For a better implementation, would need heating values for each component
        
        # Approximate values for common components in kJ/mol
        heating_values = {
            "METHANE": 890.0,
            "ETHANE": 1560.0,
            "PROPANE": 2220.0,
            "BUTANE": 2880.0,
            "PENTANE": 3510.0,
            "HEXANE": 4160.0,
            "HYDROGEN": 286.0,
            "OXYGEN": 0.0,
            "NITROGEN": 0.0,
            "CO2": 0.0,
            "CO": 283.0,
            "WATER": 0.0
        }
        
        # Get component names - this would depend on how composition is represented
        # In a real implementation, would extract names from the mixture
        
        # Default to a moderate value for hydrocarbons (methane)
        return 890.0  # J/mol
    
    def _calculate_net_heating_value(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate net heating value (lower heating value)."""
        # Typically, NHV = GHV - latent heat of water formed during combustion
        # For a proper implementation, would need detailed composition analysis
        
        # Simplified approach: NHV is roughly 10% less than GHV for hydrocarbons
        ghv = self._calculate_gross_heating_value(base_props, rp, z)
        return ghv * 0.9
    
    # ========================================================
    # EXCESS PROPERTY CALCULATIONS
    # ========================================================
    
    def _calculate_excess_volume(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate excess volume of mixing."""
        # Only applicable for mixtures
        if sum(1 for comp in z if comp > 0) <= 1:
            return 0.0
            
        # This would require calculating the volume of the mixture
        # and subtracting the sum of the pure component volumes
        # weighted by their mole fractions
        
        # Simplified approach: use a small value
        # For a proper implementation, would need to calculate volumes
        # of pure components at the same T and P
        
        return 0.0001  # m³/mol
    
    def _calculate_excess_enthalpy(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate excess enthalpy of mixing."""
        # Only applicable for mixtures
        if sum(1 for comp in z if comp > 0) <= 1:
            return 0.0
            
        # This would require calculating the enthalpy of the mixture
        # and subtracting the sum of the pure component enthalpies
        # weighted by their mole fractions
        
        # Simplified approach: use a small value
        # For a proper implementation, would need to calculate enthalpies
        # of pure components at the same T and P
        
        return 0.0  # J/mol
    
    # ========================================================
    # CRITICAL PROPERTY CALCULATIONS
    # ========================================================
    
    def _calculate_critical_z(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate critical compressibility factor Zc."""
        Pc = base_props["critical_pressure"] * 100  # bar to kPa
        Tc = base_props["critical_temperature"]  # K
        Dc = base_props["critical_density"]  # mol/L
        
        # Zc = Pc/(ρc*R*Tc)
        return Pc / (Dc * 8.31446261815324 * Tc)
    
    # ========================================================
    # MISCELLANEOUS PROPERTY CALCULATIONS
    # ========================================================
    
    def _calculate_exergy(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate specific exergy (flow availability)."""
        h = base_props["enthalpy"]  # J/mol
        s = base_props["entropy"]  # J/(mol·K)
        T = base_props["temperature"] + 273.15  # °C to K
        
        # Reference state (ambient) properties
        T0 = 298.15  # K (25°C)
        P0 = 101.325  # kPa (1 atm)
        
        # A comprehensive exergy calculation would determine h0 and s0 at the reference state
        # but for simplicity, we use a physical exergy approximation
        
        # Physical exergy approximation: ex_ph ≈ h - h0 - T0*(s - s0) ≈ h - T0*s
        return h - T0 * s  # J/mol
    
    def _calculate_exergy_loss(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate exergy loss or destruction based on entropy generation."""
        # Reference temperature
        T0 = 298.15  # K (25°C)
        
        # If entropy generation data is available, use it
        if "entropy_generation" in base_props:
            return T0 * base_props["entropy_generation"]  # J/mol
        
        # Estimate entropy generation from current state vs. reference
        try:
            # Try to compute reference state properties
            P0 = 101.325  # kPa
            D0, Dl0, Dv0, x0, y0, q0, e0, h0, s0, Cv0, Cp0, w0, ierr, herr = rp.TPFLSHdll(T0, P0, z)
            
            if ierr == 0:
                # Entropy generation is increase from reference state
                s_generation = base_props["entropy"] - s0
                return T0 * max(0, s_generation)  # J/mol
        except Exception:
            pass
            
        # Without sufficient data, return zero
        return 0.0  # J/mol
    
    def _calculate_thermodynamic_efficiency(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate theoretical thermodynamic efficiency."""
        # For power cycles, efficiency relates to temperature
        T_hot = base_props["temperature"] + 273.15  # °C to K
        T_cold = 298.15  # K (25°C reference)
        
        # Carnot efficiency is the theoretical maximum 
        if T_hot > T_cold:
            return 1.0 - T_cold / T_hot  # dimensionless
        return 0.0
    
    def _calculate_available_energy(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate available energy for work extraction."""
        # Similar to exergy but focused on energy transformation potential
        exergy = self._calculate_exergy(base_props, rp, z)
        
        # Account for any process inefficiencies (simplified)
        process_efficiency = 0.75  # Typical process efficiency factor
        
        return exergy * process_efficiency  # J/mol
    
    def _calculate_second_law_efficiency(self, base_props: Dict[str, Any], rp: Any, z: List[float]) -> float:
        """Calculate second law efficiency."""
        # Second law efficiency = actual work / maximum theoretical work
        # Without actual work data, estimate from state properties
        
        exergy = self._calculate_exergy(base_props, rp, z)
        
        # Theoretical maximum work from Carnot cycle
        T = base_props["temperature"] + 273.15  # °C to K
        T0 = 298.15  # K (25°C)
        
        if "enthalpy" in base_props and T > T0:
            carnot_work = base_props["enthalpy"] * (1 - T0/T)
            
            if abs(carnot_work) > 1e-10:
                return min(exergy / carnot_work, 1.0)
        
        # Default fallback
        return 0.8  # Typical value for well-designed systems
    
    # ========================================================
    # ADDITIONAL UTILITY METHODS
    # ========================================================
    
    def get_available_properties(self) -> List[str]:
        """Get list of all available properties."""
        return list(self.properties.keys())
    
    def get_property_groups(self) -> Dict[str, List[str]]:
        """Get properties organized by category groups."""
        groups = {
            "basic": [],
            "energy": [],
            "transport": [],
            "derivative": [],
            "critical": [],
            "phase": [],
            "virial": [],
            "combustion": [],
            "miscellaneous": []
        }
        
        # Categorize properties
        for prop in self.properties:
            # Skip aliases
            if self.properties[prop].get('is_alias', False):
                continue
                
            if prop in ["temperature", "pressure", "density", "vapor_fraction", "enthalpy", "entropy"]:
                groups["basic"].append(prop)
            elif prop in ["internal_energy", "helmholtz_energy", "gibbs_energy", "heat_of_vaporization", "cp", "cv"]:
                groups["energy"].append(prop)
            elif prop in ["viscosity", "thermal_conductivity", "surface_tension", "kinematic_viscosity"]:
                groups["transport"].append(prop)
            elif prop in ["dDdP", "dDdT", "isothermal_compressibility", "volume_expansivity"]:
                groups["derivative"].append(prop)
            elif prop.startswith("critical_"):
                groups["critical"].append(prop)
            elif prop in ["phase", "compressibility_factor", "vapor_fraction"]:
                groups["phase"].append(prop)
            elif "virial" in prop:
                groups["virial"].append(prop)
            elif "heating_value" in prop:
                groups["combustion"].append(prop)
            elif prop in ["exergy", "exergy_loss", "thermodynamic_efficiency", "second_law_efficiency"]:
                groups["miscellaneous"].append(prop)
            else:
                # Categorize based on dependency patterns
                deps = self.properties[prop].get("dependencies", [])
                if "viscosity" in deps or "thermal_conductivity" in deps:
                    groups["transport"].append(prop)
                elif "enthalpy" in deps or "entropy" in deps:
                    groups["energy"].append(prop)
                elif "dDdP" in deps or "dDdT" in deps:
                    groups["derivative"].append(prop)
                else:
                    groups["miscellaneous"].append(prop)
        
        return groups
    
    def get_property_dependencies(self, property_name: str) -> List[str]:
        """Get all dependencies for a property, including nested dependencies."""
        if property_name not in self.properties:
            raise ValueError(f"Unknown property: {property_name}")
            
        # Get direct dependencies
        direct_deps = self.properties[property_name].get("dependencies", [])
        
        # Handle aliases
        if self.properties[property_name].get('is_alias', False):
            canonical_name = self.properties[property_name]['alias_of']
            direct_deps = self.properties[canonical_name].get("dependencies", [])
        
        # No dependencies or base property
        if not direct_deps or self.properties[property_name].get("calculation_method") is None:
            return []
            
        # Get all nested dependencies recursively
        all_deps = set(direct_deps)
        for dep in direct_deps:
            nested_deps = self.get_property_dependencies(dep)
            all_deps.update(nested_deps)
            
        return list(all_deps)
    
    def get_base_property_requirements(self, properties: List[str]) -> List[str]:
        """
        Determine which base properties need to be calculated for a set of requested properties.
        This helps optimize calculations by only computing necessary base properties.
        
        Args:
            properties: List of property names to calculate
            
        Returns:
            List of base property names needed
        """
        base_props = set()
        
        for prop in properties:
            # Add the property if it's a base property
            if prop in self.properties and self.properties[prop].get("calculation_method") is None:
                base_props.add(prop)
                
            # Add all dependencies recursively
            deps = self.get_property_dependencies(prop)
            for dep in deps:
                if dep in self.properties and self.properties[dep].get("calculation_method") is None:
                    base_props.add(dep)
        
        # Always include critical base properties
        essential_props = ["temperature", "pressure", "density", "vapor_fraction"]
        for prop in essential_props:
            base_props.add(prop)
            
        return list(base_props)