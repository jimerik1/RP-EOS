# API/endpoints/enthalpy_range_from_temp.py
from flask import request, jsonify
import numpy as np
import sys
import traceback
from typing import Dict, List, Any, Tuple, Optional
import logging

# Import the blueprint 
from API.endpoints import enthalpy_range_from_temp_bp
from API.refprop_setup import RP
from API.unit_converter import UnitConverter
from API.utils.helpers import validate_composition

logger = logging.getLogger(__name__)

def setup_mixture(composition: List[Dict[str, Any]]) -> List[float]:
    """Setup REFPROP mixture (Helper function)"""
    fluid_string = '|'.join(f"{comp['fluid']}.FLD" for comp in composition)
    z = [comp['fraction'] for comp in composition] + [0] * (20 - len(composition))
    ierr, herr = RP.SETUPdll(len(composition), fluid_string, 'HMX.BNC', 'DEF')
    if ierr > 0:
        raise ValueError(f"Error setting up mixture: {herr}")
    return z

def _calculate_saturation_enthalpy(T_K: float, z: List[float]) -> Tuple[Optional[float], Optional[float]]:
    """
    Calculate saturated liquid and vapor enthalpy at a given temperature.
    """
    h_l, h_v = None, None
    logger.debug(f"Calculating saturation enthalpy at T = {T_K:.2f} K") # Log input T

    # --- Liquid Phase ---
    try:
        logger.debug(f"Calling SATTdll(T={T_K}, z, kph=1)")
        P_sat_l, Dl, Dv_l, x_l, y_l, ierr_sat_l, herr_sat_l = RP.SATTdll(T_K, z, 1)
        herr_sat_l_decoded = herr_sat_l.decode().strip() if isinstance(herr_sat_l, bytes) else str(herr_sat_l).strip() # Decode error message
        logger.debug(f"SATTdll(kph=1) Result: P={P_sat_l}, Dl={Dl}, Dv={Dv_l}, ierr={ierr_sat_l}, herr='{herr_sat_l_decoded}'") # Log SATT result

        if ierr_sat_l == 0 and Dl > 0:
            logger.debug(f"Calling THERMdll(T={T_K}, D={Dl}, z) for liquid")
            P_therm_l, e_l, h_l_calc, s_l, Cv_l, Cp_l, w_l, hjt_l, ierr_therm_l, herr_therm_l = RP.THERMdll(T_K, Dl, z)
            herr_therm_l_decoded = herr_therm_l.decode().strip() if isinstance(herr_therm_l, bytes) else str(herr_therm_l).strip() # Decode error message
            logger.debug(f"THERMdll(liquid) Result: h={h_l_calc}, ierr={ierr_therm_l}, herr='{herr_therm_l_decoded}'") # Log THERM result

            if ierr_therm_l == 0:
                h_l = h_l_calc
                logger.debug(f"Successfully calculated h_l = {h_l}")
            else:
                logger.warning(f"THERMdll(liquid) failed at T={T_K} K, D={Dl}. Error: {herr_therm_l_decoded}")
        elif ierr_sat_l != 0:
             logger.warning(f"SATTdll(bubble) failed at T={T_K} K. Error: {herr_sat_l_decoded}")
        elif Dl <= 0:
             logger.warning(f"SATTdll(bubble) returned non-positive liquid density Dl={Dl} at T={T_K} K.")

    except Exception as e:
        logger.error(f"Exception calculating liquid saturation enthalpy at T={T_K} K: {e}", exc_info=True)

    # --- Vapor Phase ---
    try:
        logger.debug(f"Calling SATTdll(T={T_K}, z, kph=2)")
        P_sat_v, Dl_v, Dv, x_v, y_v, ierr_sat_v, herr_sat_v = RP.SATTdll(T_K, z, 2)
        herr_sat_v_decoded = herr_sat_v.decode().strip() if isinstance(herr_sat_v, bytes) else str(herr_sat_v).strip() # Decode error message
        logger.debug(f"SATTdll(kph=2) Result: P={P_sat_v}, Dl={Dl_v}, Dv={Dv}, ierr={ierr_sat_v}, herr='{herr_sat_v_decoded}'") # Log SATT result

        if ierr_sat_v == 0 and Dv > 0:
            logger.debug(f"Calling THERMdll(T={T_K}, D={Dv}, z) for vapor")
            P_therm_v, e_v, h_v_calc, s_v, Cv_v, Cp_v, w_v, hjt_v, ierr_therm_v, herr_therm_v = RP.THERMdll(T_K, Dv, z)
            herr_therm_v_decoded = herr_therm_v.decode().strip() if isinstance(herr_therm_v, bytes) else str(herr_therm_v).strip() # Decode error message
            logger.debug(f"THERMdll(vapor) Result: h={h_v_calc}, ierr={ierr_therm_v}, herr='{herr_therm_v_decoded}'") # Log THERM result

            if ierr_therm_v == 0:
                h_v = h_v_calc
                logger.debug(f"Successfully calculated h_v = {h_v}")
            else:
                logger.warning(f"THERMdll(vapor) failed at T={T_K} K, D={Dv}. Error: {herr_therm_v_decoded}")
        elif ierr_sat_v != 0:
             logger.warning(f"SATTdll(dew) failed at T={T_K} K. Error: {herr_sat_v_decoded}")
        elif Dv <= 0:
             logger.warning(f"SATTdll(dew) returned non-positive vapor density Dv={Dv} at T={T_K} K.")

    except Exception as e:
        logger.error(f"Exception calculating vapor saturation enthalpy at T={T_K} K: {e}", exc_info=True)

    logger.debug(f"Returning enthalpies: h_l={h_l}, h_v={h_v}")
    return h_l, h_v


@enthalpy_range_from_temp_bp.route('/enthalpy_range_from_temperature', methods=['POST'])
def enthalpy_range_from_temperature():
    """
    Calculates the enthalpy values corresponding to a given temperature range
    at a specified pressure (default: 101.325 kPa / 1 atm).
    """
    try:
        data = request.get_json(force=True)

        # --- Validation ---
        if 'composition' not in data:
            return jsonify({'error': 'Missing composition field'}), 400
        if 'temperature_range' not in data:
            return jsonify({'error': 'Missing temperature_range field'}), 400
        if 'from' not in data['temperature_range'] or 'to' not in data['temperature_range']:
             return jsonify({'error': 'temperature_range must include "from" and "to"'}), 400

        if not validate_composition(data['composition']):
            return jsonify({'error': 'Invalid composition - fractions must sum to 1'}), 400

        units_system = data.get('units_system', 'SI').upper()
        if units_system not in ['SI', 'CGS']:
            return jsonify({'error': 'Invalid units_system. Use "SI" or "CGS".'}), 400

        # --- Setup ---
        composition = data['composition']
        temp_range = data['temperature_range']
        T_min_C = float(temp_range['from'])
        T_max_C = float(temp_range['to'])
        
        # Default pressure (1 atm / 101.325 kPa)
        pressure = data.get('pressure', 101.325)  # kPa

        if T_min_C >= T_max_C:
             return jsonify({'error': '"from" temperature must be less than "to" temperature'}), 400

        z = setup_mixture(composition)
        wmm = RP.WMOLdll(z) # Needed for unit conversion
        converter = UnitConverter()

        # --- Get Critical and Triple Point Info ---
        notes = []
        Tc, Pc, Dc, ierr_crit, herr_crit = RP.CRITPdll(z)
        T_crit_K = Tc if ierr_crit == 0 else None
        
        # Get triple point info if it's a pure fluid
        is_pure = sum(1 for frac in z if frac > 0) == 1
        T_triple_K = None
        if is_pure:
            try:
                comp_idx = next(i + 1 for i, frac in enumerate(z) if frac > 0)
                wmm_info, Ttrp, Tnbpt, Tc_info, Pc_info, Dc_info, Zc, acf, dip, Rgas = RP.INFOdll(comp_idx)
                T_triple_K = Ttrp if Ttrp > 0 else None
            except Exception as e:
                logger.warning(f"Could not get triple point info: {e}")

        # --- Calculation ---
        T_min_K = T_min_C + 273.15
        T_max_K = T_max_C + 273.15
        
        # Add informational notes about phase region
        if T_crit_K:
            if T_min_K > T_crit_K:
                notes.append(f"Temperature range is entirely in supercritical region (above critical temperature of {T_crit_K - 273.15:.2f}°C)")
            elif T_max_K > T_crit_K:
                notes.append(f"Temperature range spans across critical point ({T_crit_K - 273.15:.2f}°C)")
        
        if T_triple_K and T_min_K < T_triple_K:
            notes.append(f"Min temperature ({T_min_C}°C) is below triple point temperature ({T_triple_K - 273.15:.2f}°C). Some calculations may be extrapolations.")
        
        # Generate a reasonable number of temperature points to evaluate
        num_points = 10
        temp_step = (T_max_K - T_min_K) / (num_points - 1)
        temperatures_K = [T_min_K + i * temp_step for i in range(num_points)]
        
        # Calculate enthalpies at each temperature
        enthalpies = []
        
        for T_K in temperatures_K:
            try:
                # Use TPFLSHdll to calculate thermodynamic properties at T and P
                D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr = RP.TPFLSHdll(T_K, pressure, z)
                
                if ierr == 0:
                    enthalpies.append(h)
                else:
                    logger.warning(f"TPFLSHdll failed at T={T_K} K, P={pressure} kPa. Error: {herr}")
            except Exception as e:
                logger.error(f"Exception calculating properties at T={T_K} K: {e}", exc_info=True)
        
        # --- Format Response ---
        result = {
            "temperature_range_input": {
                "from": T_min_C,
                "to": T_max_C,
                "unit": "°C"
            },
            "pressure": {
                "value": pressure,
                "unit": "kPa"
            },
            "notes": notes
        }
        
        enthalpy_unit = converter.get_unit("enthalpy", units_system)
        
        if enthalpies:
            h_min = min(enthalpies)
            h_max = max(enthalpies)
            result["enthalpy_range"] = {
                "min": converter.convert_property("enthalpy", h_min, wmm, 'SI', units_system)['value'],
                "max": converter.convert_property("enthalpy", h_max, wmm, 'SI', units_system)['value'],
                "unit": enthalpy_unit
            }
        else:
            notes.append("Could not determine enthalpy range (calculation failed)")
        
        return jsonify(result)

    except ValueError as ve:
        logger.error(f"Validation error: {str(ve)}")
        return jsonify({'error': str(ve)}), 400
    except Exception as e:
        logger.error("Error processing request:", exc_info=True)
        traceback.print_exc()
        return jsonify({'error': f"An unexpected error occurred: {str(e)}"}), 500
