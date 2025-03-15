# file: API/endpoints/phase_envelope_pt.py
from flask import request, jsonify
import numpy as np
import sys
import traceback
from ctypes import c_double, c_int, byref, create_string_buffer

from API.endpoints import phase_envelope_pt_bp  # <-- Define this blueprint in __init__.py
from API.refprop_setup import RP
from API.unit_converter import UnitConverter   # or wherever you handle your unit conversions
from API.utils.helpers import validate_composition, get_phase, convert_for_json

@phase_envelope_pt_bp.route('/phase_envelope_pt', methods=['POST'])
def phase_envelope_pt():
    try:
        data = request.get_json(force=True)
        
        # Validate the new structure
        required_fields = ['composition', 'variables']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing field: {field}'}), 400
                
        # Check if required variables exist
        variables = data.get('variables', {})
        if 'temperature' not in variables:
            return jsonify({'error': 'Missing temperature variable'}), 400

        # Extract calculation settings
        calculation = data.get('calculation', {})
        desired_curve = calculation.get('curve_type', 'both')  # Default to both curves
        
        # Validate composition
        if not validate_composition(data['composition']):
            return jsonify({'error': 'Invalid composition - fractions must sum to 1'}), 400

        # Extract range and resolution parameters
        temperature_range = variables['temperature'].get('range', {})
        temperature_resolution = variables['temperature'].get('resolution')
        
        # Validate required parameters exist
        if not all([temperature_range.get('from'), temperature_range.get('to'), 
                    temperature_resolution]):
            return jsonify({'error': 'Missing temperature range or resolution parameters'}), 400

        # Setup mixture:
        fluid_string = '|'.join(f"{c['fluid']}.FLD" for c in data['composition'])
        z = [c['fraction'] for c in data['composition']] + [0]*(20-len(data['composition']))
        
        # Convert z to c_double array for REFPROP
        z_array = (len(z)*c_double)(*z)

        ierr, herr = RP.SETUPdll(len(data['composition']),
                                 fluid_string,
                                 'HMX.BNC',
                                 'DEF')
        if ierr > 0:
            return jsonify({"error": f"REFPROP SETUPdll error: {herr}"}), 400

        results_bubble = []
        results_dew = []
        
        # Loop over temperature range:
        Ts = np.arange(
            float(temperature_range['from']) + 273.15, 
            float(temperature_range['to']) + 273.15, 
            float(temperature_resolution)
        )
        
        # Rest of the function remains the same
        for T in Ts:
            # Bubble point calculation (kph=1)
            if desired_curve in ["both", "bubble"]:
                try:
                    result = RP.SATTdll(T, z_array, 1)
                    if result.ierr == 0:
                        # Convert to our units and store
                        T_C = T - 273.15
                        results_bubble.append({
                            "temperature": T_C,
                            "pressure": result.P / 100,  # Convert kPa to bar
                            "liquid_density": result.Dl,
                            "vapor_density": result.Dv,
                            "liquid_composition": list(result.x[:len(data['composition'])]),
                            "vapor_composition": list(result.y[:len(data['composition'])])
                        })
                except Exception as e:
                    print(f"Bubble point error at T={T-273.15}°C: {str(e)}")
            
            # Dew point calculation (kph=2)
            if desired_curve in ["both", "dew"]:
                try:
                    result = RP.SATTdll(T, z_array, 2)
                    if result.ierr == 0:
                        T_C = T - 273.15
                        results_dew.append({
                            "temperature": T_C,
                            "pressure": result.P / 100,  # Convert kPa to bar
                            "liquid_density": result.Dl,
                            "vapor_density": result.Dv,
                            "liquid_composition": list(result.x[:len(data['composition'])]),
                            "vapor_composition": list(result.y[:len(data['composition'])])
                        })
                except Exception as e:
                    print(f"Dew point error at T={T-273.15}°C: {str(e)}")

        # Return final results in JSON
        response = {
            "bubble_curve": results_bubble,
            "dew_curve": results_dew
        }
        return jsonify(response)

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500