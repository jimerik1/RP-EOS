from flask import request, jsonify, Response
import numpy as np
import sys
import traceback
from typing import List, Dict, Any

from API.endpoints import vt_flash_bp
from API.refprop_setup import RP
from API.unit_converter import UnitConverter
from API.utils.helpers import get_phase
from API.utils.grid_generator import generate_grid, get_phase_boundaries_tv

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

def calculate_properties_vt(z: List[float], T: float, v: float, units_system: str = 'SI') -> Dict[str, Any]:
    """
    Calculate fluid properties at given temperature and specific volume with specified unit system.
    
    Args:
        z: Composition array
        T: Temperature in K
        v: Specific volume in m³/mol
        units_system: Unit system to use ('SI' or 'CGS')
    
    Returns:
        Dictionary of calculated properties with values and units
    """
    # Initialize unit converter
    converter = UnitConverter()
    
    # Get molecular weight for unit conversions
    wmm = RP.WMOLdll(z)
    
    # Convert specific volume to density
    D = 1.0 / v  # mol/L = mol/(0.001 m³) = 1000 / (m³/mol)
    
    # Use TDFLSHdll (Temperature-Density) for VT flash
    P, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr = RP.TDFLSHdll(T, D, z)
    if ierr > 0:
        raise ValueError(f"Error in TDFLSHdll: {herr}")
    
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
    Z = P / (D * 8.31446261815324 * T)
    
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
        'pressure': P / 100,  # Convert kPa to bar
        'specific_volume': v,
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

@vt_flash_bp.route('/vt_flash', methods=['POST'])
def vt_flash():
    try:
        data = request.get_json(force=True)
        
        # Validate the request structure
        required_fields = ['composition', 'variables']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing field: {field}'}), 400
                
        # Check if required variables exist
        variables = data.get('variables', {})
        if 'temperature' not in variables or 'specific_volume' not in variables:
            return jsonify({'error': 'Missing temperature or specific_volume variables'}), 400

        # Extract calculation settings
        calculation = data.get('calculation', {})
        properties = calculation.get('properties', [])
        units_system = calculation.get('units_system', 'SI')  # Default to SI
        response_format = calculation.get('response_format', 'json')  # Default to JSON
        
        # Extract grid_type parameter and related options
        grid_type = calculation.get('grid_type', 'equidistant')  # Default to equidistant grid
        enhancement_factor = calculation.get('enhancement_factor', 5.0)  
        boundary_zone_width = calculation.get('boundary_zone_width', None)  
        
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
            temperature_range = variables['temperature'].get('range', {})
            volume_range = variables['specific_volume'].get('range', {})
            
            # Ensure all necessary values exist with defaults if not
            t_from = float(temperature_range.get('from', 0.0))
            t_to = float(temperature_range.get('to', 100.0))
            v_from = float(volume_range.get('from', 0.01))
            v_to = float(volume_range.get('to', 1.0))
            
            temperature_resolution = float(variables['temperature'].get('resolution', 5.0))
            volume_resolution = float(variables['specific_volume'].get('resolution', 0.05))
            
            # Ensure to > from
            if t_to <= t_from:
                t_to = t_from + temperature_resolution
            if v_to <= v_from:
                v_to = v_from + volume_resolution
                
            # Update the ranges for later use
            temperature_range = {'from': t_from, 'to': t_to}
            volume_range = {'from': v_from, 'to': v_to}
            
        except (ValueError, TypeError) as ve:
            return jsonify({'error': f'Invalid range or resolution parameters: {str(ve)}'}), 400

        # Debug log
        print(f"Calculating VT flash for temperature range: {temperature_range['from']} to {temperature_range['to']} °C, "
              f"specific volume range: {volume_range['from']} to {volume_range['to']} m³/mol, grid_type: {grid_type}")

        # Generate grids based on grid_type
        if grid_type.lower() != 'equidistant':
            # If we're using an adaptive grid, determine phase boundaries first
            if grid_type.lower() == 'adaptive':
                try:
                    # Get phase boundaries in T-V space
                    t_boundaries, v_boundaries = get_phase_boundaries_tv(
                        RP, z, temperature_range, volume_range
                    )
                    
                    print(f"Identified phase boundaries: {len(t_boundaries)} temperature points, "
                          f"{len(v_boundaries)} specific volume points")
                except Exception as e:
                    print(f"Error determining phase boundaries: {e}")
                    t_boundaries, v_boundaries = [], []
            else:
                t_boundaries, v_boundaries = [], []
            
            # Generate grids using the utility function
            T_range = generate_grid(
                t_from + 273.15,  # Convert from °C to K
                t_to + 273.15, 
                temperature_resolution,
                grid_type, 
                t_boundaries, 
                enhancement_factor, 
                boundary_zone_width
            )
            
            v_range = generate_grid(
                v_from, 
                v_to, 
                volume_resolution,
                grid_type, 
                v_boundaries, 
                enhancement_factor, 
                boundary_zone_width
            )
            
            # Debug info about the grid
            print(f"Generated {len(T_range)} temperature points and {len(v_range)} specific volume points")
            
        else:
            # Create regular (equidistant) grids as before
            T_range = np.arange(
                t_from + 273.15,  # Convert from °C to K
                t_to + 273.15 + temperature_resolution,
                temperature_resolution
            )
            
            v_range = np.arange(
                v_from,
                v_to + volume_resolution,
                volume_resolution
            )

        # Calculate properties
        results = []
        idx = 0
        for t_idx, T in enumerate(T_range):
            for v_idx, v in enumerate(v_range):
                try:
                    props = calculate_properties_vt(z, float(T), float(v), units_system)
                    filtered_props = {k: v for k, v in props.items() 
                                   if k in properties or k in ['temperature', 'pressure', 'specific_volume', 'phase']}
                    
                    # Add grid indices for formatting
                    results.append({
                        'index': idx,
                        't_idx': t_idx,
                        'v_idx': v_idx,
                        **filtered_props
                    })
                    idx += 1
                except Exception as fe:
                    print(f"Error processing T={T-273.15}°C, v={v} m³/mol: {fe}", file=sys.stderr)
                    continue

        # Return response in the requested format
        if response_format.lower() == 'olga_tab':
            from API.utils.olga_formatter import format_olga_tab
            try:
                # Create structured variable dictionaries for formatter
                temperature_vars = {
                    'range': {'from': (T_range.min() - 273.15), 'to': (T_range.max() - 273.15)},
                    'resolution': temperature_resolution,
                    'values': T_range - 273.15  # Pass the actual grid values in °C
                }
                
                volume_vars = {
                    'range': {'from': v_range.min(), 'to': v_range.max()},
                    'resolution': volume_resolution,
                    'values': v_range  # Pass the actual grid values
                }
                
                response = format_olga_tab(
                    temperature_vars,  # x-axis (temperature)
                    volume_vars,       # y-axis (specific volume)
                    results,
                    data['composition'],
                    wmm,
                    endpoint_type='vt_flash'
                )
                return response  # Return the Response object directly
            except Exception as e:
                print(f"Error formatting OLGA TAB: {e}", file=sys.stderr)
                traceback.print_exc()
                return jsonify({'error': f'Error formatting OLGA TAB response: {str(e)}'}), 500
        else:
            # For JSON responses, include grid information
            return jsonify({
                'results': results,
                'grid_info': {
                    'type': grid_type,
                    'temperature_points': len(T_range),
                    'specific_volume_points': len(v_range),
                    'total_points': len(results)
                }
            })
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500