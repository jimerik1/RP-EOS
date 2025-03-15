# file: API/endpoints/phase_envelope_pt.py
from flask import request, jsonify
import numpy as np
import sys
import traceback
from ctypes import c_double, c_int, byref, create_string_buffer

from API.endpoints import phase_envelope_pt_bp  # <-- Define this blueprint in __init__.py
from API.refprop_setup import RP
from API.unit_converter import UnitConverter   # or wherever you handle your unit conversions
from API.utils.helpers import convert_for_json

@phase_envelope_pt_bp.route('/phase_envelope_pt', methods=['POST'])
def phase_envelope_pt():
    """
    Calculate bubble/dew curves in the PT plane (Phase Envelope).
    Example JSON payload might look like:
    {
      "composition": [ {"fluid": "CO2", "fraction": 0.4},
                       {"fluid": "N2",  "fraction": 0.6} ],
      "temperature_range": { "from": -50, "to": 100 },
      "temperature_resolution": 2,
      "desired_curve": "both"      # or "bubble" or "dew"
    }
    """
    try:
        data = request.get_json(force=True)
        composition = data["composition"]
        T_min = float(data["temperature_range"]["from"])
        T_max = float(data["temperature_range"]["to"])
        dT = float(data["temperature_resolution"])
        desired_curve = data.get("desired_curve", "both")  # bubble, dew, or both

        # Setup mixture:
        fluid_string = '|'.join(f"{c['fluid']}.FLD" for c in composition)
        z = [c['fraction'] for c in composition] + [0]*(20-len(composition))
        
        # Convert z to c_double array for REFPROP
        z_array = (len(z)*c_double)(*z)
        
        ierr, herr = RP.SETUPdll(len(composition),
                                 fluid_string,  # fluid string
                                 'HMX.BNC',     # mixture file
                                 'DEF')
        if ierr > 0:
            return jsonify({"error": f"REFPROP SETUPdll error: {herr}"}), 400

        results_bubble = []
        results_dew = []
        
        # Loop over temperature range:
        Ts = np.arange(T_min+273.15, T_max+273.15, dT)  # in K if your input is in °C
        
        # Note: The ctREFPROP wrapper provides a different interface than calling the
        # DLL directly. We need to use the wrapper's methods correctly.
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
                            "liquid_composition": list(result.x[:len(composition)]),
                            "vapor_composition": list(result.y[:len(composition)])
                        })
                except Exception as e:
                    print(f"Bubble point error at T={T-273.15}°C: {str(e)}")
                    # Continue to next point
            
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
                            "liquid_composition": list(result.x[:len(composition)]),
                            "vapor_composition": list(result.y[:len(composition)])
                        })
                except Exception as e:
                    print(f"Dew point error at T={T-273.15}°C: {str(e)}")
                    # Continue to next point

        # Return final results in JSON
        response = {
            "bubble_curve": results_bubble,
            "dew_curve": results_dew
        }
        return jsonify(response)

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500