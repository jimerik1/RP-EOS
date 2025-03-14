from flask import request, jsonify
import numpy as np
import sys
from typing import List, Dict, Any

from API.refprop_setup import RP
from API.unit_converter import UnitConverter
from API.utils.helpers import get_phase, convert_for_json

# Create a blueprint for this endpoint
# This line assumes you've added vle_bp to endpoints/__init__.py
from flask import Blueprint
vle_bp = Blueprint('vle', __name__)

@vle_bp.route('/vle', methods=['POST'])
def calculate_vle():
    """
    Calculate vapor-liquid equilibrium for a specified mixture.
    
    Expects JSON with:
    {
        "composition": [{"fluid": "name", "fraction": value}],
        "temperature_range": {"from": value, "to": value},
        "temperature_resolution": value,
        "units_system": "SI" or "CGS"
    }
    """
    try:
        data = request.get_json(force=True)
        
        # Validate request
        if 'composition' not in data or 'temperature_range' not in data:
            return jsonify({'error': 'Missing required fields'}), 400
            
        # Setup the mixture
        fluid_string = '|'.join(f"{comp['fluid']}.FLD" for comp in data['composition'])
        z = [comp['fraction'] for comp in data['composition']] + [0] * (20 - len(data['composition']))
        
        ierr, herr = RP.SETUPdll(len(data['composition']), fluid_string, 'HMX.BNC', 'DEF')
        if ierr > 0:
            return jsonify({'error': f'Error setting up mixture: {herr}'}), 400
            
        # Get temperature range
        t_start = float(data['temperature_range']['from']) + 273.15  # Convert to K
        t_end = float(data['temperature_range']['to']) + 273.15
        t_step = float(data.get('temperature_resolution', 1.0))
        
        # Unit system
        units_system = data.get('units_system', 'SI')
        converter = UnitConverter()
        
        # Calculate vapor pressure curve
        temperatures = np.arange(t_start, t_end + t_step, t_step)
        results = []
        
        for T in temperatures:
            try:
                # Calculate saturation state at temperature T
                P, Dl, Dv, x, y, ierr, herr = RP.SATTdll(T, z, 1)
                
                if ierr == 0:
                    # Convert to requested units
                    wmm = RP.WMOLdll(z)
                    
                    # Get properties for liquid and vapor phases
                    e_liq, h_liq, s_liq, cv_liq, cp_liq, w_liq = RP.THERMdll(T, Dl, z)
                    e_vap, h_vap, s_vap, cv_vap, cp_vap, w_vap = RP.THERMdll(T, Dv, z)
                    
                    # Transport properties
                    eta_liq, tcx_liq, ierr_liq, herr_liq = RP.TRNPRPdll(T, Dl, z)
                    eta_vap, tcx_vap, ierr_vap, herr_vap = RP.TRNPRPdll(T, Dv, z)
                    
                    # Surface tension
                    sigma, ierr_surf, herr_surf = RP.SURTENdll(T, Dl, Dv, x, y)
                    
                    # Build response
                    result = {
                        'temperature': converter.convert_property('temperature', T-273.15, wmm, 'SI', units_system),
                        'pressure': converter.convert_property('pressure', P/100, wmm, 'SI', units_system),  # kPa to bar
                        'liquid_density': converter.convert_property('density', Dl, wmm, 'SI', units_system),
                        'vapor_density': converter.convert_property('density', Dv, wmm, 'SI', units_system),
                        'liquid_enthalpy': converter.convert_property('enthalpy', h_liq, wmm, 'SI', units_system),
                        'vapor_enthalpy': converter.convert_property('enthalpy', h_vap, wmm, 'SI', units_system),
                        'surface_tension': converter.convert_property('surface_tension', sigma, wmm, 'SI', units_system),
                        'liquid_composition': [float(val) for val in x[:len(data['composition'])]],
                        'vapor_composition': [float(val) for val in y[:len(data['composition'])]]
                    }
                    
                    results.append(result)
            except Exception as e:
                print(f"Error at T={T}: {str(e)}", file=sys.stderr)
                continue
                
        return jsonify({'vle_data': results})
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500