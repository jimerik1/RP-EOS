from flask import request, jsonify
import numpy as np
import sys
import traceback
from typing import List, Dict, Any, Tuple

from API.endpoints import phase_boundaries_bp
from API.refprop_setup import RP
from API.unit_converter import UnitConverter
from API.utils.helpers import validate_composition

def setup_mixture(composition: List[Dict[str, Any]]) -> List[float]:
    """Setup REFPROP mixture"""
    fluid_string = '|'.join(f"{comp['fluid']}.FLD" for comp in composition)
    z = [comp['fraction'] for comp in composition] + [0] * (20 - len(composition))
    
    ierr, herr = RP.SETUPdll(len(composition), fluid_string, 'HMX.BNC', 'DEF')
    if ierr > 0:
        raise ValueError(f"Error setting up mixture: {herr}")
    return z

def calculate_phase_boundaries(z: List[float], temp_range: Dict[str, float], 
                              temp_resolution: float, boundary_types: List[str],
                              units_system: str = 'SI') -> Dict[str, List[Dict[str, Any]]]:
    """
    Calculate phase boundaries (melting, sublimation, vapor-liquid) for the given composition.
    
    Args:
        z: Composition array
        temp_range: Dictionary with 'from' and 'to' temperature values (°C)
        temp_resolution: Temperature step size (°C)
        boundary_types: List of boundary types to calculate ('melting', 'sublimation', 'vaporization')
        units_system: Unit system to use ('SI' or 'CGS')
    
    Returns:
        Dictionary with phase boundary data
    """
    # Initialize unit converter
    converter = UnitConverter()
    
    # Get molecular weight for unit conversions
    wmm = RP.WMOLdll(z)
    
    # Generate temperature range in K
    temps = np.arange(
        float(temp_range['from']) + 273.15,
        float(temp_range['to']) + 273.15 + float(temp_resolution),
        float(temp_resolution)
    )
    
    # Initialize result arrays
    melting_curve = []
    sublimation_curve = []
    vaporization_curve = []
    triple_point = None
    critical_point = None
    
    # Get critical point
    try:
        Tc, Pc, Dc, ierr, herr = RP.CRITPdll(z)
        if ierr == 0:
            critical_temp = Tc - 273.15  # Convert to Celsius
            critical_press = Pc / 100    # Convert kPa to bar
            critical_point = {
                "temperature": converter.convert_property("temperature", float(critical_temp), wmm, 'SI', units_system),
                "pressure": converter.convert_property("pressure", float(critical_press), wmm, 'SI', units_system)
            }
    except Exception as e:
        print(f"Error getting critical point: {str(e)}")
    
    # Check for supported substances for solid calculations
    # For pure fluids only
    is_pure = sum(1 for component in z if component > 0) == 1
    pure_component_idx = next((i+1 for i, component in enumerate(z) if component > 0), 0) if is_pure else 0
    
    # Get component name for logging
    comp_info = RP.NAMEdll(pure_component_idx) if is_pure else None
    comp_name = comp_info[0] if comp_info else "Mixture"
    
    print(f"Calculating phase boundaries for {comp_name}")
    
    # Calculate melting curve if requested
    if 'melting' in boundary_types and is_pure:
        for T in temps:
            try:
                # Use MELTTdll to calculate melting pressure at given temperature
                P, ierr, herr = RP.MELTTdll(T, z)
                if ierr == 0:
                    # Pressures less than 1e-6 are essentially zero, which can create issues in log plots
                    if P < 1e-6:
                        continue
                        
                    # Convert units as needed
                    T_unit = converter.convert_property("temperature", float(T - 273.15), wmm, 'SI', units_system)
                    P_unit = converter.convert_property("pressure", float(P / 100), wmm, 'SI', units_system)
                    
                    # Get enthalpy at the melting point using TPFLSHdll
                    try:
                        D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr2, herr2 = RP.TPFLSHdll(T, P, z)
                        if ierr2 == 0:
                            # Convert enthalpy
                            h_unit = converter.convert_property("enthalpy", float(h), wmm, 'SI', units_system)
                            
                            # Add to melting curve
                            melting_curve.append({
                                "temperature": T_unit,
                                "pressure": P_unit,
                                "enthalpy": h_unit
                            })
                    except Exception as e:
                        # If enthalpy calculation fails, just record T and P
                        melting_curve.append({
                            "temperature": T_unit,
                            "pressure": P_unit
                        })
            except Exception as e:
                # Skip if MELTTdll fails for this temperature
                continue
        
        print(f"Generated {len(melting_curve)} points for melting curve")
    
    # Calculate sublimation curve if requested
    if 'sublimation' in boundary_types and is_pure:
        for T in temps:
            try:
                # Use SUBLTdll to calculate sublimation pressure at given temperature
                P, ierr, herr = RP.SUBLTdll(T, z)
                if ierr == 0:
                    # Pressures less than 1e-6 are essentially zero
                    if P < 1e-6:
                        continue
                        
                    # Convert units as needed
                    T_unit = converter.convert_property("temperature", float(T - 273.15), wmm, 'SI', units_system)
                    P_unit = converter.convert_property("pressure", float(P / 100), wmm, 'SI', units_system)
                    
                    # Get enthalpy at the sublimation point using TPFLSHdll
                    try:
                        D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr2, herr2 = RP.TPFLSHdll(T, P, z)
                        if ierr2 == 0:
                            # Convert enthalpy
                            h_unit = converter.convert_property("enthalpy", float(h), wmm, 'SI', units_system)
                            
                            # Add to sublimation curve
                            sublimation_curve.append({
                                "temperature": T_unit,
                                "pressure": P_unit,
                                "enthalpy": h_unit
                            })
                    except Exception as e:
                        # If enthalpy calculation fails, just record T and P
                        sublimation_curve.append({
                            "temperature": T_unit,
                            "pressure": P_unit
                        })
            except Exception as e:
                # Skip if SUBLTdll fails for this temperature
                continue
        
        print(f"Generated {len(sublimation_curve)} points for sublimation curve")
    
    # Calculate vapor-liquid (saturation) curve if requested
    if 'vaporization' in boundary_types:
        for T in temps:
            try:
                # Use SATTdll to calculate saturation properties at given temperature
                # kph=1 gets the bubble point, kph=2 gets the dew point
                P, Dl, Dv, x, y, ierr, herr = RP.SATTdll(T, z, 1)  # Bubble point
                
                if ierr == 0:
                    # Convert units
                    T_unit = converter.convert_property("temperature", float(T - 273.15), wmm, 'SI', units_system)
                    P_unit = converter.convert_property("pressure", float(P / 100), wmm, 'SI', units_system)
                    
                    # Get enthalpy at bubble point
                    try:
                        D, Dl, Dv, x_vle, y_vle, q, e, h, s, Cv, Cp, w, ierr2, herr2 = RP.TPFLSHdll(T, P, z)
                        if ierr2 == 0:
                            h_unit = converter.convert_property("enthalpy", float(h), wmm, 'SI', units_system)
                            
                            vaporization_curve.append({
                                "temperature": T_unit,
                                "pressure": P_unit,
                                "enthalpy": h_unit,
                                "density_liquid": converter.convert_property("density", float(Dl), wmm, 'SI', units_system),
                                "density_vapor": converter.convert_property("density", float(Dv), wmm, 'SI', units_system),
                                "type": "bubble"  # Indicate this is a bubble point
                            })
                    except Exception as e:
                        # If enthalpy calculation fails, add without enthalpy
                        vaporization_curve.append({
                            "temperature": T_unit,
                            "pressure": P_unit,
                            "type": "bubble"
                        })
            except Exception as e:
                # Skip if SATTdll fails for this temperature
                continue
        
        print(f"Generated {len(vaporization_curve)} points for vapor-liquid curve")
    
    # Try to find the triple point if we're dealing with a pure fluid
    if is_pure and ('sublimation' in boundary_types or 'melting' in boundary_types):
        try:
            # Get triple point temperature from INFO
            wmm, Ttrp, Tnbpt, Tc, Pc, Dc, Zc, acf, dip, Rgas = RP.INFOdll(pure_component_idx)
            
            if Ttrp > 0:  # Valid triple point found
                # Get pressure at triple point using SUBLTdll
                P_triple, ierr, herr = RP.SUBLTdll(Ttrp, z)
                
                if ierr == 0:
                    # Convert units
                    Ttrp_unit = converter.convert_property("temperature", float(Ttrp - 273.15), wmm, 'SI', units_system)
                    P_triple_unit = converter.convert_property("pressure", float(P_triple / 100), wmm, 'SI', units_system)
                    
                    # Get enthalpy at triple point
                    try:
                        D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr2, herr2 = RP.TPFLSHdll(Ttrp, P_triple, z)
                        if ierr2 == 0:
                            h_unit = converter.convert_property("enthalpy", float(h), wmm, 'SI', units_system)
                            
                            triple_point = {
                                "temperature": Ttrp_unit,
                                "pressure": P_triple_unit,
                                "enthalpy": h_unit
                            }
                    except Exception as e:
                        # If enthalpy calculation fails
                        triple_point = {
                            "temperature": Ttrp_unit,
                            "pressure": P_triple_unit
                        }
        except Exception as e:
            print(f"Error finding triple point: {str(e)}")
    
    # Build final response
    result = {}
    
    if melting_curve:
        result["melting_curve"] = melting_curve
        
    if sublimation_curve:
        result["sublimation_curve"] = sublimation_curve
        
    if vaporization_curve:
        result["vaporization_curve"] = vaporization_curve
        
    if triple_point:
        result["triple_point"] = triple_point
        
    if critical_point:
        result["critical_point"] = critical_point
    
    return result

@phase_boundaries_bp.route('/phase_boundaries', methods=['POST'])
def phase_boundaries():
    try:
        data = request.get_json(force=True)
        
        # Validate the request structure
        required_fields = ['composition', 'variables']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing field: {field}'}), 400
                
        # Extract variables
        variables = data.get('variables', {})
        if 'temperature' not in variables:
            return jsonify({'error': 'Missing temperature variable'}), 400
            
        # Extract temperature range and resolution
        temp_range = variables['temperature'].get('range', {})
        temp_resolution = variables['temperature'].get('resolution')
        
        # Validate required parameters exist
        if not all([temp_range.get('from'), temp_range.get('to'), temp_resolution]):
            return jsonify({'error': 'Missing temperature range or resolution parameters'}), 400
        
        # Extract calculation settings
        calculation = data.get('calculation', {})
        boundary_types = calculation.get('boundary_types', ['melting', 'sublimation', 'vaporization'])
        units_system = calculation.get('units_system', 'SI')  # Default to SI
        
        # Validate composition
        if not validate_composition(data['composition']):
            return jsonify({'error': 'Invalid composition - fractions must sum to 1'}), 400

        # Setup mixture
        z = setup_mixture(data['composition'])
        
        # Debug log
        print(f"Calculating phase boundaries for temperature range: {temp_range['from']} to {temp_range['to']} °C")
        
        # Calculate phase boundaries
        boundaries = calculate_phase_boundaries(
            z, temp_range, temp_resolution, boundary_types, units_system
        )
        
        # Return results
        return jsonify(boundaries)
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500