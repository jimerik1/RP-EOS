from flask import Blueprint, jsonify, request
import sys
import traceback
from typing import List, Dict, Any

# Create a blueprint for the available_properties endpoint
available_properties_bp = Blueprint('available_properties', __name__)

# Define the properties that can be requested in flash calculations
AVAILABLE_PROPERTIES = [
    {
        'id': 'temperature',
        'name': 'Temperature',
        'description': 'Temperature',
        'si_unit': '°C',
        'cgs_unit': '°C',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': True,
        'is_output': True
    },
    {
        'id': 'pressure',
        'name': 'Pressure',
        'description': 'Pressure',
        'si_unit': 'bar',
        'cgs_unit': 'dyn/cm²',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': True,
        'is_output': True
    },
    {
        'id': 'density',
        'name': 'Bulk density',
        'description': 'Overall mixture density',
        'si_unit': 'mol/L',
        'cgs_unit': 'g/cm³',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'liquid_density',
        'name': 'Liquid density',
        'description': 'Density of the liquid phase',
        'si_unit': 'mol/L',
        'cgs_unit': 'g/cm³',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'vapor_density',
        'name': 'Vapor density',
        'description': 'Density of the vapor phase',
        'si_unit': 'mol/L',
        'cgs_unit': 'g/cm³',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'critical_density',
        'name': 'Critical density',
        'description': 'Density at the critical point',
        'si_unit': 'mol/L',
        'cgs_unit': 'g/cm³',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'critical_pressure',
        'name': 'Critical pressure',
        'description': 'Pressure at the critical point',
        'si_unit': 'bar',
        'cgs_unit': 'dyn/cm²',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'critical_temperature',
        'name': 'Critical temperature',
        'description': 'Temperature at the critical point',
        'si_unit': 'K',
        'cgs_unit': 'K',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'internal_energy',
        'name': 'Internal energy',
        'description': 'Specific internal energy',
        'si_unit': 'J/mol',
        'cgs_unit': 'erg/g',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': True,
        'is_output': True
    },
    {
        'id': 'enthalpy',
        'name': 'Enthalpy',
        'description': 'Specific enthalpy',
        'si_unit': 'J/mol',
        'cgs_unit': 'erg/g',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': True,
        'is_output': True
    },
    {
        'id': 'entropy',
        'name': 'Entropy',
        'description': 'Specific entropy',
        'si_unit': 'J/(mol·K)',
        'cgs_unit': 'erg/(g·K)',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': True,
        'is_output': True
    },
    {
        'id': 'cv',
        'name': 'Isochoric specific heat',
        'description': 'Isochoric (constant volume) specific heat capacity',
        'si_unit': 'J/(mol·K)',
        'cgs_unit': 'erg/(g·K)',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'cp',
        'name': 'Isobaric specific heat',
        'description': 'Isobaric (constant pressure) specific heat capacity',
        'si_unit': 'J/(mol·K)',
        'cgs_unit': 'erg/(g·K)',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'sound_speed',
        'name': 'Speed of sound',
        'description': 'Speed of sound in the fluid',
        'si_unit': 'm/s',
        'cgs_unit': 'cm/s',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'viscosity',
        'name': 'Dynamic viscosity',
        'description': 'Dynamic viscosity of the fluid',
        'si_unit': 'μPa·s',
        'cgs_unit': 'poise',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'thermal_conductivity',
        'name': 'Thermal conductivity',
        'description': 'Thermal conductivity of the fluid',
        'si_unit': 'W/(m·K)',
        'cgs_unit': 'cal/(s·cm·K)',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'surface_tension',
        'name': 'Surface tension',
        'description': 'Surface tension between liquid and vapor phases',
        'si_unit': 'N/m',
        'cgs_unit': 'dyn/cm',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'vapor_fraction',
        'name': 'Vapor quality',
        'description': 'Fraction of vapor phase (quality)',
        'si_unit': 'dimensionless',
        'cgs_unit': 'dimensionless',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': True,
        'is_output': True
    },
    {
        'id': 'compressibility_factor',
        'name': 'Compressibility factor',
        'description': 'Compressibility factor (Z)',
        'si_unit': 'dimensionless',
        'cgs_unit': 'dimensionless',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'prandtl_number',
        'name': 'Prandtl number',
        'description': 'Prandtl number',
        'si_unit': 'dimensionless',
        'cgs_unit': 'dimensionless',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'isothermal_compressibility',
        'name': 'Isothermal compressibility',
        'description': 'Isothermal compressibility',
        'si_unit': '1/kPa',
        'cgs_unit': 'cm²/dyn',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'volume_expansivity',
        'name': 'Volume expansivity',
        'description': 'Volume expansivity (thermal expansion coefficient)',
        'si_unit': '1/K',
        'cgs_unit': '1/K',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'dp_dt_saturation',
        'name': 'dP/dT saturation',
        'description': 'Derivative of pressure with respect to temperature along saturation line',
        'si_unit': 'kPa/K',
        'cgs_unit': 'dyn/(cm²·K)',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'joule_thomson_coefficient',
        'name': 'Joule-Thomson coefficient',
        'description': 'Joule-Thomson (Joule-Kelvin) coefficient',
        'si_unit': 'K/bar',
        'cgs_unit': 'K·cm²/dyn',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'kinematic_viscosity',
        'name': 'Kinematic viscosity',
        'description': 'Kinematic viscosity (dynamic viscosity divided by density)',
        'si_unit': 'cm²/s',
        'cgs_unit': 'cm²/s',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'thermal_diffusivity',
        'name': 'Thermal diffusivity',
        'description': 'Thermal diffusivity',
        'si_unit': 'cm²/s',
        'cgs_unit': 'cm²/s',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'phase',
        'name': 'Phase state',
        'description': 'Phase state of the fluid (Liquid, Vapor, Two-Phase, etc.)',
        'si_unit': 'text',
        'cgs_unit': 'text',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'specific_volume',
        'name': 'Specific volume',
        'description': 'Specific volume (1/density)',
        'si_unit': 'm³/mol',
        'cgs_unit': 'cm³/g',
        'flash_types': ['vt_flash', 'uv_flash'],
        'is_input': True,
        'is_output': True
    },
    {
        'id': 'dDdP',
        'name': 'dD/dP (density derivative)',
        'description': 'Pressure derivative of density',
        'si_unit': '(mol/L)/kPa',
        'cgs_unit': '(g/cm³)/(dyn/cm²)',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    },
    {
        'id': 'dDdT',
        'name': 'dD/dT (density derivative)',
        'description': 'Temperature derivative of density',
        'si_unit': '(mol/L)/K',
        'cgs_unit': '(g/cm³)/K',
        'flash_types': ['pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash'],
        'is_input': False,
        'is_output': True
    }
]

# Define flash calculation types
FLASH_TYPES = [
    {
        'id': 'pt_flash',
        'name': 'Pressure-Temperature Flash',
        'description': 'Calculate properties based on pressure and temperature inputs',
        'input_properties': ['pressure', 'temperature'],
        'endpoint': '/pt_flash'
    },
    {
        'id': 'ph_flash',
        'name': 'Pressure-Enthalpy Flash',
        'description': 'Calculate properties based on pressure and enthalpy inputs',
        'input_properties': ['pressure', 'enthalpy'],
        'endpoint': '/ph_flash'
    },
    {
        'id': 'ts_flash',
        'name': 'Temperature-Entropy Flash',
        'description': 'Calculate properties based on temperature and entropy inputs',
        'input_properties': ['temperature', 'entropy'],
        'endpoint': '/ts_flash'
    },
    {
        'id': 'vt_flash',
        'name': 'Specific Volume-Temperature Flash',
        'description': 'Calculate properties based on specific volume and temperature inputs',
        'input_properties': ['specific_volume', 'temperature'],
        'endpoint': '/vt_flash'
    },
    {
        'id': 'uv_flash',
        'name': 'Internal Energy-Specific Volume Flash',
        'description': 'Calculate properties based on internal energy and specific volume inputs',
        'input_properties': ['internal_energy', 'specific_volume'],
        'endpoint': '/uv_flash'
    }
]

@available_properties_bp.route('/available_properties', methods=['GET'])
def available_properties():
    """
    Get a list of all available properties that can be used in the API.
    
    Query parameters:
    - flash_type: Filter properties by compatible flash type
    - input_only: If 'true', return only input properties
    - output_only: If 'true', return only output properties
    """
    try:
        # Parse query parameters
        flash_type = request.args.get('flash_type')
        input_only = request.args.get('input_only', 'false').lower() == 'true'
        output_only = request.args.get('output_only', 'false').lower() == 'true'
        
        # Filter properties if flash_type is provided
        filtered_properties = AVAILABLE_PROPERTIES
        if flash_type:
            filtered_properties = [
                prop for prop in AVAILABLE_PROPERTIES 
                if flash_type in prop['flash_types']
            ]
        
        # Filter by input/output if requested
        if input_only:
            filtered_properties = [
                prop for prop in filtered_properties 
                if prop['is_input']
            ]
        elif output_only:
            filtered_properties = [
                prop for prop in filtered_properties 
                if prop['is_output']
            ]
        
        # Filter flash types if flash_type is provided
        filtered_flash_types = FLASH_TYPES
        if flash_type:
            filtered_flash_types = [
                ft for ft in FLASH_TYPES 
                if ft['id'] == flash_type
            ]
        
        return jsonify({
            'properties': filtered_properties,
            'flash_types': filtered_flash_types
        })
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500