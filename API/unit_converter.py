"""
Unit conversion module for REFPROP wrapper.
Handles conversions between different unit systems with proper handling of molecular weights
and compound units. Supports both forward and reverse conversions.
"""

from typing import Dict, Any, Union, List, Optional
import numpy as np

class UnitConverter:
    def __init__(self):
        # Basic conversion factors
        self.PRESSURE_FACTORS = {
            'SI_TO_CGS': 10000,  # kPa to dyn/cm²
            'CGS_TO_SI': 1/10000,  # dyn/cm² to kPa
            'BAR_TO_KPA': 100,   # bar to kPa
            'KPA_TO_BAR': 0.01   # kPa to bar
        }
        
        self.ENERGY_FACTORS = {
            'J_TO_ERG': 1e7,     # Joules to ergs
            'ERG_TO_J': 1e-7,    # ergs to Joules
            'CAL_TO_J': 4.184,   # calories to Joules
            'J_TO_CAL': 1/4.184  # Joules to calories
        }
        
        self.LENGTH_FACTORS = {
            'M_TO_CM': 100,      # meters to centimeters
            'CM_TO_M': 0.01      # centimeters to meters
        }
        
        self.THERMAL_CONDUCTIVITY_FACTORS = {
            'W_M_K_TO_CAL_CM_S_K': 0.002388458966275,  # W/(m·K) to cal/(cm·s·K)
            'CAL_CM_S_K_TO_W_M_K': 1/0.002388458966275  # cal/(cm·s·K) to W/(m·K)
        }
        
        self.VISCOSITY_FACTORS = {
            'UPA_S_TO_POISE': 1e-7,  # μPa·s to poise
            'POISE_TO_UPA_S': 1e7    # poise to μPa·s
        }
        
        # Define unit mappings for each property in each unit system
        self.UNITS = {
            'SI': {
                'density': 'mol/L',
                'liquid_density': 'mol/L',
                'vapor_density': 'mol/L',
                'critical_density': 'mol/L',
                'pressure': 'bar',
                'critical_pressure': 'bar',
                'internal_energy': 'J/mol',
                'enthalpy': 'J/mol',
                'entropy': 'J/(mol·K)',
                'cv': 'J/(mol·K)',
                'cp': 'J/(mol·K)',
                'sound_speed': 'm/s',
                'viscosity': 'μPa·s',
                'thermal_conductivity': 'W/(m·K)',
                'surface_tension': 'N/m',
                'temperature': '°C',
                'critical_temperature': 'K',
                'isothermal_compressibility': '1/kPa',
                'volume_expansivity': '1/K',
                'dp_dt_saturation': 'kPa/K',
                'joule_thomson_coefficient': 'K/bar',
                'kinematic_viscosity': 'cm²/s',
                'thermal_diffusivity': 'cm²/s',
                'compressibility_factor': 'dimensionless',
                'vapor_fraction': 'dimensionless',
                'prandtl_number': 'dimensionless'
            },
            'CGS': {
                'density': 'g/cm³',
                'liquid_density': 'g/cm³',
                'vapor_density': 'g/cm³',
                'critical_density': 'g/cm³',
                'pressure': 'dyn/cm²',
                'critical_pressure': 'dyn/cm²',
                'internal_energy': 'erg/g',
                'enthalpy': 'erg/g',
                'entropy': 'erg/(g·K)',
                'cv': 'erg/(g·K)',
                'cp': 'erg/(g·K)',
                'sound_speed': 'cm/s',
                'viscosity': 'poise',
                'thermal_conductivity': 'cal/(s·cm·K)',
                'surface_tension': 'dyn/cm',
                'temperature': '°C',
                'critical_temperature': 'K',
                'isothermal_compressibility': 'cm²/dyn',
                'volume_expansivity': '1/K',
                'dp_dt_saturation': 'dyn/(cm²·K)',
                'joule_thomson_coefficient': 'K·cm²/dyn',
                'kinematic_viscosity': 'cm²/s',
                'thermal_diffusivity': 'cm²/s',
                'compressibility_factor': 'dimensionless',
                'vapor_fraction': 'dimensionless',
                'prandtl_number': 'dimensionless'
            }
        }

    def get_unit(self, property_id: str, unit_system: str) -> str:
        """Get the appropriate unit string for a property in the specified unit system"""
        return self.UNITS.get(unit_system, {}).get(property_id, 'dimensionless')

    def convert_density(self, value: float, wmm: float, from_unit: str, to_unit: str) -> float:
        """
        Convert density between different units considering molecular weight
        
        Args:
            value: Density value
            wmm: Molecular weight [g/mol]
            from_unit: Original unit ('mol/L' or 'kg/m³')
            to_unit: Target unit ('g/cm³')
            
        Returns:
            Converted density value
        """
        if from_unit == 'mol/L' and to_unit == 'g/cm³':
            return value * wmm / 1000  # mol/L * (g/mol) / (1000 cm³/L)
        elif from_unit == 'g/cm³' and to_unit == 'mol/L':
            return value * 1000 / wmm  # g/cm³ * (1000 cm³/L) / (g/mol)
        elif from_unit == 'kg/m³' and to_unit == 'g/cm³':
            return value / 1000
        elif from_unit == 'g/cm³' and to_unit == 'kg/m³':
            return value * 1000
        return value

    def convert_energy_per_mass(self, value: float, wmm: float, from_unit: str, to_unit: str) -> float:
        """
        Convert energy per mass units (e.g., specific energy, enthalpy, etc.)
        
        Args:
            value: Energy value
            wmm: Molecular weight [g/mol]
            from_unit: Original unit ('J/mol' or 'erg/g')
            to_unit: Target unit ('erg/g' or 'J/mol')
            
        Returns:
            Converted energy value
        """
        if from_unit == 'J/mol' and to_unit == 'erg/g':
            return value * self.ENERGY_FACTORS['J_TO_ERG'] / wmm
        elif from_unit == 'erg/g' and to_unit == 'J/mol':
            return value * self.ENERGY_FACTORS['ERG_TO_J'] * wmm
        return value

    def convert_pressure(self, value: float, from_unit: str, to_unit: str) -> float:
        """
        Convert pressure between different units
        
        Args:
            value: Pressure value
            from_unit: Original unit ('kPa', 'bar', or 'dyn/cm²')
            to_unit: Target unit ('kPa', 'bar', or 'dyn/cm²')
            
        Returns:
            Converted pressure value
        """
        if from_unit == 'kPa' and to_unit == 'dyn/cm²':
            return value * self.PRESSURE_FACTORS['SI_TO_CGS']
        elif from_unit == 'dyn/cm²' and to_unit == 'kPa':
            return value * self.PRESSURE_FACTORS['CGS_TO_SI']
        elif from_unit == 'bar' and to_unit == 'dyn/cm²':
            return value * self.PRESSURE_FACTORS['BAR_TO_KPA'] * self.PRESSURE_FACTORS['SI_TO_CGS']
        elif from_unit == 'dyn/cm²' and to_unit == 'bar':
            return value * self.PRESSURE_FACTORS['CGS_TO_SI'] * self.PRESSURE_FACTORS['KPA_TO_BAR']
        elif from_unit == 'bar' and to_unit == 'kPa':
            return value * self.PRESSURE_FACTORS['BAR_TO_KPA']
        elif from_unit == 'kPa' and to_unit == 'bar':
            return value * self.PRESSURE_FACTORS['KPA_TO_BAR']
        return value

    def convert_thermal_conductivity(self, value: float, from_unit: str, to_unit: str) -> float:
        """
        Convert thermal conductivity between different units
        
        Args:
            value: Thermal conductivity value
            from_unit: Original unit ('W/(m·K)' or 'cal/(s·cm·K)')
            to_unit: Target unit ('cal/(s·cm·K)' or 'W/(m·K)')
            
        Returns:
            Converted thermal conductivity value
        """
        if from_unit == 'W/(m·K)' and to_unit == 'cal/(s·cm·K)':
            return value * self.THERMAL_CONDUCTIVITY_FACTORS['W_M_K_TO_CAL_CM_S_K']
        elif from_unit == 'cal/(s·cm·K)' and to_unit == 'W/(m·K)':
            return value * self.THERMAL_CONDUCTIVITY_FACTORS['CAL_CM_S_K_TO_W_M_K']
        return value

    def convert_viscosity(self, value: float, from_unit: str, to_unit: str) -> float:
        """
        Convert viscosity between different units
        
        Args:
            value: Viscosity value
            from_unit: Original unit ('μPa·s' or 'poise')
            to_unit: Target unit ('poise' or 'μPa·s')
            
        Returns:
            Converted viscosity value
        """
        if from_unit == 'μPa·s' and to_unit == 'poise':
            return value * self.VISCOSITY_FACTORS['UPA_S_TO_POISE']
        elif from_unit == 'poise' and to_unit == 'μPa·s':
            return value * self.VISCOSITY_FACTORS['POISE_TO_UPA_S']
        return value

    def convert_property(self, 
                          property_id: str, 
                          value: float, 
                          wmm: float, 
                          from_system: str = 'SI', 
                          to_system: str = 'SI') -> Dict[str, Any]:
        """
        Convert a property value between unit systems
        
        Args:
            property_id: Property identifier (e.g., 'density', 'pressure')
            value: Property value
            wmm: Molecular weight [g/mol]
            from_system: Original unit system ('SI' or 'CGS')
            to_system: Target unit system ('SI' or 'CGS')
        
        Returns:
            Dictionary with converted value and unit
        """
        # If both systems are the same, no conversion needed
        if from_system.upper() == to_system.upper():
            return {
                'value': value,
                'unit': self.get_unit(property_id, from_system)
            }

        # Get the units for this property in both systems
        from_unit = self.get_unit(property_id, from_system)
        to_unit = self.get_unit(property_id, to_system)
        
        # Handle dimensionless properties (no conversion needed)
        dimensionless_properties = [
            'vapor_fraction', 'compressibility_factor', 'prandtl_number'
        ]
        if property_id in dimensionless_properties:
            return {'value': value, 'unit': 'dimensionless'}
            
        # Properties with the same unit in both systems
        same_unit_properties = ['temperature', 'critical_temperature']
        if property_id in same_unit_properties:
            return {'value': value, 'unit': to_unit}
        
        # Apply conversions based on property type
        try:
            if property_id in ['density', 'liquid_density', 'vapor_density', 'critical_density']:
                converted = self.convert_density(value, wmm, from_unit, to_unit)
            elif property_id in ['internal_energy', 'enthalpy', 'entropy', 'cv', 'cp']:
                converted = self.convert_energy_per_mass(value, wmm, from_unit, to_unit)
            elif property_id in ['pressure', 'critical_pressure']:
                converted = self.convert_pressure(value, from_unit, to_unit)
            elif property_id == 'thermal_conductivity':
                converted = self.convert_thermal_conductivity(value, from_unit, to_unit)
            elif property_id == 'viscosity':
                converted = self.convert_viscosity(value, from_unit, to_unit)
            elif property_id == 'sound_speed':
                # m/s to cm/s or vice versa
                if from_unit == 'm/s' and to_unit == 'cm/s':
                    converted = value * self.LENGTH_FACTORS['M_TO_CM']
                elif from_unit == 'cm/s' and to_unit == 'm/s':
                    converted = value * self.LENGTH_FACTORS['CM_TO_M']
                else:
                    converted = value
            elif property_id == 'surface_tension':
                # N/m to dyn/cm or vice versa
                if from_unit == 'N/m' and to_unit == 'dyn/cm':
                    converted = value * 1000  # 1 N/m = 1000 dyn/cm
                elif from_unit == 'dyn/cm' and to_unit == 'N/m':
                    converted = value / 1000  # 1 dyn/cm = 0.001 N/m
                else:
                    converted = value
            elif property_id == 'isothermal_compressibility':
                # 1/kPa to cm²/dyn or vice versa
                if from_unit == '1/kPa' and to_unit == 'cm²/dyn':
                    converted = value / self.PRESSURE_FACTORS['SI_TO_CGS']
                elif from_unit == 'cm²/dyn' and to_unit == '1/kPa':
                    converted = value * self.PRESSURE_FACTORS['SI_TO_CGS']
                else:
                    converted = value
            elif property_id == 'dp_dt_saturation':
                # kPa/K to dyn/(cm²·K) or vice versa
                if from_unit == 'kPa/K' and to_unit == 'dyn/(cm²·K)':
                    converted = value * self.PRESSURE_FACTORS['SI_TO_CGS']
                elif from_unit == 'dyn/(cm²·K)' and to_unit == 'kPa/K':
                    converted = value * self.PRESSURE_FACTORS['CGS_TO_SI']
                else:
                    converted = value
            elif property_id == 'joule_thomson_coefficient':
                # K/bar to K·cm²/dyn or vice versa
                if from_unit == 'K/bar' and to_unit == 'K·cm²/dyn':
                    converted = value / (self.PRESSURE_FACTORS['BAR_TO_KPA'] * self.PRESSURE_FACTORS['SI_TO_CGS'])
                elif from_unit == 'K·cm²/dyn' and to_unit == 'K/bar':
                    converted = value * (self.PRESSURE_FACTORS['BAR_TO_KPA'] * self.PRESSURE_FACTORS['SI_TO_CGS'])
                else:
                    converted = value
            else:
                # No specific conversion defined, return as is
                converted = value
                
            return {'value': converted, 'unit': to_unit}
        except Exception as e:
            # If conversion fails, return with original value and unit
            print(f"Conversion error for {property_id}: {str(e)}")
            return {'value': value, 'unit': from_unit}

    def convert_property_reverse(self, 
                                 property_id: str, 
                                 value: float, 
                                 wmm: float, 
                                 from_system: str, 
                                 to_system: str = 'SI') -> Dict[str, Any]:
        """
        Convert a property value from user units back to SI units.
        
        Args:
            property_id: Property identifier (e.g., 'enthalpy', 'pressure')
            value: Property value in user units
            wmm: Molecular weight [g/mol]
            from_system: Source unit system ('SI' or 'CGS')
            to_system: Target unit system ('SI' or 'CGS')
        
        Returns:
            Dictionary with converted value and unit
        """
        # This is just calling convert_property with reversed from/to systems
        return self.convert_property(property_id, value, wmm, from_system, to_system)