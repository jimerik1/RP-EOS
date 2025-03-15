from flask import request, jsonify
import numpy as np
import sys
import traceback
from typing import List, Dict, Any

from API.endpoints import pt_flash_bp
from API.refprop_setup import RP
from API.unit_converter import UnitConverter
from API.utils.helpers import get_phase

def validate_composition(composition: List[Dict[str, Any]]) -> bool:
    """Validate composition data"""
    total = sum(comp['fraction'] for comp in composition)
    return abs(total - 1.0) < 1e-6

def setup_mixture(composition: List[Dict[str, Any]]) -> List[float]:
    """Setup REFPROP mixture"""
    fluid_string = '|'.join(f"{comp['fluid']}.FLD" for comp in composition)
    z = [comp['fraction'] for comp in composition] + [0] * (20 - len(composition))
    
    ierr, herr = RP.SETUPdll(len(composition), fluid_string, 'HMX.BNC', 'DEF')
    if ierr > 0:
        raise ValueError(f"Error setting up mixture: {herr}")
    return z

def calculate_properties(z: List[float], T: float, P: float, units_system: str = 'SI') -> Dict[str, Any]:
    """
    Calculate fluid properties at given temperature and pressure with specified unit system.
    
    Args:
        z: Composition array
        T: Temperature in K
        P: Pressure in bar
        units_system: Unit system to use ('SI' or 'CGS')
    
    Returns:
        Dictionary of calculated properties with values and units
    """
    # Initialize unit converter
    converter = UnitConverter()
    
    # Get molecular weight for unit conversions
    wmm = RP.WMOLdll(z)
    
    # Convert pressure to kPa for REFPROP
    P_kpa = P * 100  # bar to kPa
    
    # Get basic thermodynamic properties
    D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr = RP.TPFLSHdll(T, P_kpa, z)
    if ierr > 0:
        raise ValueError(f"Error in TPFLSHdll: {herr}")
    
    # Get transport properties
    eta, tcx, ierr, herr = RP.TRNPRPdll(T, D, z)
    if ierr > 0:
        raise ValueError(f"Error in TRNPRPdll: {herr}")
        
    # Get surface tension if in two-phase region
    if 0 < q < 1:
        sigma, ierr, herr = RP.SURTENdll(T, Dl, Dv, x, y)
        surface_tension = float(sigma) if ierr == 0 else None
    else:
        surface_tension = None
        
    # Get critical properties
    Tc, Pc, Dc, ierr, herr = RP.CRITPdll(z)
    if ierr > 0:
        Tc = Pc = Dc = None
        
    # Get additional thermodynamic derivatives
    dPdD, dPdT, d2PdD2, d2PdT2, d2PdTD, dDdP, dDdT, d2DdP2, d2DdT2, d2DdPT, \
    dTdP, dTdD, d2TdP2, d2TdD2, d2TdPD = RP.DERVPVTdll(T, D, z)
    
    # Calculate compressibility factor
    Z = P_kpa / (D * 8.31446261815324 * T)
    
    # Get molecular weight for conversions
    wmm_kg = wmm / 1000  # Convert g/mol to kg/mol

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
        'isothermal_compressibility': -1/D * dDdP,
        'volume_expansivity': 1/D * dDdT,
        'dp_dt_saturation': dPdT,
        'joule_thomson_coefficient': (T*dDdT/dDdP - 1)/(Cp * 100),  # Convert to K/bar
        'kinematic_viscosity': (eta * 1e-6) / (D * wmm / 1000) * 10000,
        'thermal_diffusivity': tcx / (D * Cp) * 10000,  # Convert to cm²/s
        'prandtl_number': (Cp / wmm_kg) * (eta * 1e-6) / tcx,
        'temperature': T - 273.15,  # Convert to Celsius
        'pressure': P
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
                properties[prop_id] = {'value': float(value), 'unit': 'unknown'}
    
    # Add phase information
    properties['phase'] = {
        'value': get_phase(q),
        'unit': None
    }
    
    return properties

@pt_flash_bp.route('/pt_flash', methods=['POST'])
def pt_flash():
    try:
        data = request.get_json(force=True)
        
        # Validate the new structure
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
        
        if not properties:
            return jsonify({'error': 'No properties specified for calculation'}), 400

        # Validate composition
        if not validate_composition(data['composition']):
            return jsonify({'error': 'Invalid composition - fractions must sum to 1'}), 400

        # Setup mixture
        z = setup_mixture(data['composition'])

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
        print("Received request with temperature range:", temperature_range,
              "and pressure range:", pressure_range)

        # Create arrays for calculations
        T_range = np.arange(
            float(temperature_range['from']) + 273.15,
            float(temperature_range['to']) + 273.15 + float(temperature_resolution),
            float(temperature_resolution)
        )
        P_range = np.arange(
            float(pressure_range['from']),
            float(pressure_range['to']) + float(pressure_resolution),
            float(pressure_resolution)
        )

        # Calculate properties (rest of the function remains the same)
        results = []
        idx = 0
        for T in T_range:
            for P in P_range:
                try:
                    props = calculate_properties(z, float(T), float(P), units_system)
                    filtered_props = {k: v for k, v in props.items() 
                                   if k in properties or k in ['temperature', 'pressure']}
                    results.append({
                        'index': idx,
                        **filtered_props
                    })
                    idx += 1
                except Exception as fe:
                    print(f"Error processing T={T}, P={P}: {fe}", file=sys.stderr)
                    continue

        return jsonify({'results': results})
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500