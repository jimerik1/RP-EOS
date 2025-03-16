from flask import request, jsonify, Response
import numpy as np
import sys
import traceback
from typing import List, Dict, Any

from API.endpoints import ph_flash_bp
from API.refprop_setup import RP
from API.unit_converter import UnitConverter
from API.utils.helpers import get_phase, convert_for_json

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

def calculate_properties_ph(z: List[float], P: float, h: float, units_system: str = 'SI') -> Dict[str, Any]:
    """
    Calculate fluid properties at given pressure and enthalpy with specified unit system.
    
    Args:
        z: Composition array
        P: Pressure in bar
        h: Enthalpy in J/mol
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
    
    # Get basic thermodynamic properties using PH flash
    T, D, Dl, Dv, x, y, q, e, s, Cv, Cp, w, ierr, herr = RP.PHFLSHdll(P_kpa, h, z)
    if ierr > 0:
        raise ValueError(f"Error in PHFLSHdll: {herr}")
    
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
        'critical_pressure': Pc / 100 if Pc is not None else None,  # Convert from kPa to bar
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
        'pressure': P,
        'x': list(x[:len(z)]),  # Liquid composition
        'y': list(y[:len(z)]),  # Vapor composition
        'dDdP': dDdP,          # Add pressure derivative of density
        'dDdT': dDdT           # Add temperature derivative of density
    }
    
    # Convert properties to requested unit system
    properties = {}
    for prop_id, value in raw_properties.items():
        if value is not None:  # Skip undefined properties
            try:
                if prop_id in ['x', 'y']:
                    # Composition vectors
                    properties[prop_id] = {'value': value, 'unit': 'mole fraction'}
                else:
                    properties[prop_id] = converter.convert_property(
                        prop_id, float(value), wmm, 'SI', units_system
                    )
            except Exception as e:
                print(f"Warning: Could not convert property {prop_id}: {str(e)}")
                properties[prop_id] = {'value': float(value) if value is not None else None, 'unit': 'unknown'}
    
    # Add phase information
    properties['phase'] = {
        'value': get_phase(q),
        'unit': None
    }
    
    return properties

@ph_flash_bp.route('/ph_flash', methods=['POST'])
def ph_flash():
    try:
        data = request.get_json(force=True)
        
        # Validate the new structure
        required_fields = ['composition', 'variables']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing field: {field}'}), 400
                
        # Check if required variables exist
        variables = data.get('variables', {})
        if 'pressure' not in variables or 'enthalpy' not in variables:
            return jsonify({'error': 'Missing pressure or enthalpy variables'}), 400

        # Extract calculation settings
        calculation = data.get('calculation', {})
        properties = calculation.get('properties', [])
        units_system = calculation.get('units_system', 'SI')  # Default to SI
        response_format = calculation.get('response_format', 'json')  # Default to JSON
        
        if not properties:
            return jsonify({'error': 'No properties specified for calculation'}), 400

        # Validate composition
        if not validate_composition(data['composition']):
            return jsonify({'error': 'Invalid composition - fractions must sum to 1'}), 400

        # Setup mixture
        z = setup_mixture(data['composition'])
        
        # Get molecular weight for unit conversions
        wmm = RP.WMOLdll(z)

        # Extract range and resolution parameters with robust error handling
        try:
            pressure_range = variables['pressure'].get('range', {})
            enthalpy_range = variables['enthalpy'].get('range', {})
            
            # Ensure all necessary values exist with defaults if not
            p_from = float(pressure_range.get('from', 1.0))
            p_to = float(pressure_range.get('to', 100.0))
            h_from = float(enthalpy_range.get('from', 100.0))
            h_to = float(enthalpy_range.get('to', 1000.0))
            
            pressure_resolution = float(variables['pressure'].get('resolution', 10.0))
            enthalpy_resolution = float(variables['enthalpy'].get('resolution', 100.0))
            
            # Ensure to > from
            if p_to <= p_from:
                p_to = p_from + pressure_resolution
            if h_to <= h_from:
                h_to = h_from + enthalpy_resolution
                
            # Update the ranges for later use
            pressure_range = {'from': p_from, 'to': p_to}
            enthalpy_range = {'from': h_from, 'to': h_to}
            
        except (ValueError, TypeError) as ve:
            return jsonify({'error': f'Invalid range or resolution parameters: {str(ve)}'}), 400

        # Debug log
        print("Received request with enthalpy range:", enthalpy_range,
              "and pressure range:", pressure_range)

        # Create arrays for calculations
        P_range = np.arange(p_from, p_to + pressure_resolution, pressure_resolution)
        h_range = np.arange(h_from, h_to + enthalpy_resolution, enthalpy_resolution)

        # Calculate properties
        results = []
        idx = 0
        for p_idx, P in enumerate(P_range):
            for h_idx, h in enumerate(h_range):
                try:
                    props = calculate_properties_ph(z, float(P), float(h), units_system)
                    filtered_props = {k: v for k, v in props.items() 
                                   if k in properties or k in ['temperature', 'pressure', 'enthalpy']}
                    
                    # Add grid indices for OLGA TAB formatting
                    results.append({
                        'index': idx,
                        'p_idx': p_idx,
                        'h_idx': h_idx,
                        **filtered_props
                    })
                    idx += 1
                except Exception as fe:
                    print(f"Error processing P={P}, h={h}: {fe}", file=sys.stderr)
                    continue

        # Return response in the requested format
        if response_format.lower() == 'olga_tab':
            from API.utils.olga_formatter import format_olga_tab
            try:
                # Create structured variable dictionaries for formatter
                pressure_vars = {
                    'range': pressure_range,
                    'resolution': pressure_resolution
                }
                
                enthalpy_vars = {
                    'range': enthalpy_range,
                    'resolution': enthalpy_resolution
                }
                
                response = format_olga_tab(
                    pressure_vars,     # x-axis (pressure) 
                    enthalpy_vars,     # y-axis (enthalpy)
                    results,
                    data['composition'],
                    wmm,
                    endpoint_type='ph_flash'  # Specify endpoint type for correct grid variables
                )
                return response  # Return the Response object directly
            except Exception as e:
                print(f"Error formatting OLGA TAB: {e}", file=sys.stderr)
                traceback.print_exc()
                return jsonify({'error': f'Error formatting OLGA TAB response: {str(e)}'}), 500
        else:
            return jsonify({'results': results})
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500