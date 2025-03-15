# file: API/endpoints/phase_envelope_ph.py
from flask import request, jsonify
import numpy as np
import sys
import traceback
from ctypes import c_double, c_int, byref

from API.endpoints import phase_envelope_ph_bp
from API.refprop_setup import RP
from API.utils.helpers import convert_for_json

@phase_envelope_ph_bp.route('/phase_envelope_ph', methods=['POST'])
def phase_envelope_ph():
    """
    Calculate bubble/dew curves in the p-h plane (Phase Envelope).
    Example JSON payload:
    {
      "composition": [ {"fluid": "CO2", "fraction": 1.0} ],
      "pressure_range": { "from": 10, "to": 200 },
      "pressure_resolution": 5,
      "desired_curve": "both"      # or "bubble"/"dew"
    }
    """
    try:
        data = request.get_json(force=True)
        composition = data["composition"]
        p_min = float(data["pressure_range"]["from"])
        p_max = float(data["pressure_range"]["to"])
        dp = float(data["pressure_resolution"])
        desired_curve = data.get("desired_curve", "both")

        # Setup mixture:
        fluid_string = '|'.join(f"{c['fluid']}.FLD" for c in composition)
        z = [c['fraction'] for c in composition] + [0]*(20-len(composition))
        
        # Convert z to c_double array for REFPROP
        z_array = (len(z)*c_double)(*z)

        ierr, herr = RP.SETUPdll(len(composition),
                                 fluid_string,
                                 'HMX.BNC',
                                 'DEF')
        if ierr > 0:
            return jsonify({"error": f"REFPROP SETUPdll error: {herr}"}), 400
        
        results_bubble = []
        results_dew = []

        # Convert pressures from bar to kPa for REFPROP
        pressures = np.arange(p_min*100, (p_max+dp)*100, dp*100)
        
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
                    # Continue to next pressure
            
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
                    # Continue to next pressure

        response = {
            "bubble_curve": results_bubble,
            "dew_curve": results_dew
        }
        return jsonify(response)

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500