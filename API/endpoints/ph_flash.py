from flask import request, jsonify
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
        'critical_pressure': Pc,
        'critical_density': Dc,
        'compressibility_factor': Z,
        'isothermal_compressibility': -1/D * dDdP,
        'volume_expansivity': 1/D * dDdT,
        'dp_dt_saturation': dPdT,
        'joule_thomson_coefficient': (T*dDdT/dDdP - 1)/Cp,
        'kinematic_viscosity': eta/D,
        'thermal_diffusivity': tcx/(D*Cp),
        'prandtl_number': Cp*eta/tcx,
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

@ph_flash_bp.route('/ph_flash', methods=['POST'])
def ph_flash():
    try:
        data = request.get_json(force=True)
        # Validate that required fields exist
        required_fields = ['composition', 'pressure_range', 'enthalpy_range', 
                          'pressure_resolution', 'enthalpy_resolution', 'properties']
        units_system = data.get('units_system', 'SI')  # Default to SI if not specified

        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing field: {field}'}), 400

        if not validate_composition(data['composition']):
            return jsonify({'error': 'Invalid composition - fractions must sum to 1'}), 400

        z = setup_mixture(data['composition'])

        # Debug log the inputs before performing Fortran calls
        print("Received request with enthalpy range:", data["enthalpy_range"],
              "and pressure range:", data["pressure_range"])

        # Create arrays for pressure and enthalpy values
        P_range = np.arange(
            float(data['pressure_range']['from']),
            float(data['pressure_range']['to']) + float(data['pressure_resolution']),
            float(data['pressure_resolution'])
        )
        
        h_range = np.arange(
            float(data['enthalpy_range']['from']),
            float(data['enthalpy_range']['to']) + float(data['enthalpy_resolution']),
            float(data['enthalpy_resolution'])
        )

        results = []
        idx = 0
        for P in P_range:
            for h in h_range:
                try:
                    props = calculate_properties_ph(z, float(P), float(h), units_system)
                    filtered_props = {k: v for k, v in props.items() 
                                   if k in data['properties'] or k in ['temperature', 'pressure', 'enthalpy']}
                    results.append({
                        'index': idx,
                        **filtered_props
                    })
                    idx += 1
                except Exception as fe:
                    print(f"Error processing P={P}, h={h}: {fe}", file=sys.stderr)
                    continue

        return jsonify({'results': results})
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500