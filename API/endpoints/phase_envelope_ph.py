# file: API/endpoints/phase_envelope_ph.py
from flask import request, jsonify
import numpy as np
import sys
import traceback
from ctypes import c_double, c_int, byref

from API.endpoints import phase_envelope_ph_bp
from API.refprop_setup import RP
from API.utils.helpers import validate_composition, get_phase, convert_for_json

@phase_envelope_ph_bp.route('/phase_envelope_ph', methods=['POST'])
def phase_envelope_ph():
    try:
        data = request.get_json(force=True)
        
        # Validate the new structure
        required_fields = ['composition', 'variables']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing field: {field}'}), 400
                
        # Check if required variables exist
        variables = data.get('variables', {})
        if 'pressure' not in variables:
            return jsonify({'error': 'Missing pressure variable'}), 400

        # Extract calculation settings
        calculation = data.get('calculation', {})
        desired_curve = calculation.get('curve_type', 'both')  # Default to both curves
        
        # Validate composition
        if not validate_composition(data['composition']):
            return jsonify({'error': 'Invalid composition - fractions must sum to 1'}), 400

        # Extract range and resolution parameters
        pressure_range = variables['pressure'].get('range', {})
        pressure_resolution = variables['pressure'].get('resolution')
        
        # Validate required parameters exist
        if not all([pressure_range.get('from'), pressure_range.get('to'), 
                    pressure_resolution]):
            return jsonify({'error': 'Missing pressure range or resolution parameters'}), 400

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

        # Convert pressures from bar to kPa for REFPROP
        pressures = np.arange(
            float(pressure_range['from'])*100, 
            (float(pressure_range['to']) + float(pressure_resolution))*100, 
            float(pressure_resolution)*100
        )
        
        # Rest of the function remains the same
        for p_kpa in pressures:
            # Calculate bubble point (q=0)
            if desired_curve in ["both", "bubble"]:
                try:
                    # Using the ctREFPROP wrapper correctly
                    result = RP.PQFLSHdll(p_kpa, 0.0, z_array, 1)
                    
                    if result.ierr == 0:
                        results_bubble.append({
                            "pressure": p_kpa / 100,  # Convert kPa to bar
                            "enthalpy": result.h,
                            "temperature": result.T - 273.15,  # Convert K to °C
                            "density": result.D,
                            "vapor_fraction": 0.0  # By definition, bubble point has q=0
                        })
                except Exception as e:
                    print(f"Bubble point error at P={p_kpa/100} bar: {str(e)}")
            
            # Calculate dew point (q=1)
            if desired_curve in ["both", "dew"]:
                try:
                    # Using the ctREFPROP wrapper correctly
                    result = RP.PQFLSHdll(p_kpa, 1.0, z_array, 1)
                    
                    if result.ierr == 0:
                        results_dew.append({
                            "pressure": p_kpa / 100,  # Convert kPa to bar
                            "enthalpy": result.h,
                            "temperature": result.T - 273.15,  # Convert K to °C
                            "density": result.D,
                            "vapor_fraction": 1.0  # By definition, dew point has q=1
                        })
                except Exception as e:
                    print(f"Dew point error at P={p_kpa/100} bar: {str(e)}")

        response = {
            "bubble_curve": results_bubble,
            "dew_curve": results_dew
        }
        return jsonify(response)

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500