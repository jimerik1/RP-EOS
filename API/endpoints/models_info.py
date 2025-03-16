# API/endpoints/models_info.py
from flask import request, jsonify
import traceback
import sys
from typing import Dict, Any, List

from API.endpoints import models_info_bp  # You'll need to create this blueprint in __init__.py
from API.refprop_setup import RP
from API.utils.helpers import validate_composition, trim

def get_model_info(z: List[float]) -> Dict[str, Any]:
    """
    Get detailed information about thermodynamic models used for a given composition
    
    Args:
        z: Composition array (mole fractions)
    
    Returns:
        Dictionary with model information
    """
    models = {}
    
    try:
        # Get the primary EOS model
        hcode, hcite = RP.GETMODdll(0, 'EOS')
        models['primary_eos'] = {
            'code': trim(hcode.raw),
            'reference': trim(hcite.raw)
        }
        
        # Get transport property models
        property_codes = {
            'viscosity': 'VIS',
            'thermal_conductivity': 'TCX',
            'surface_tension': 'STN',
            'dielectric_constant': 'DIE'
        }
        
        models['property_models'] = {}
        for prop_name, code in property_codes.items():
            try:
                hcode, hcite = RP.GETMODdll(0, code)
                models['property_models'][prop_name] = {
                    'code': trim(hcode.raw),
                    'reference': trim(hcite.raw)
                }
            except:
                # Skip properties that don't have specific models
                pass
        
        # Check if GERG models are being used
        try:
            # Check for GERG-2004
            ierr, herr = RP.GERG04dll(0, 0)  # Just check status, don't modify
            if ierr == 0:
                models['gerg_2004_available'] = True
            
            # Check for GERG-2008
            ierr, herr = RP.GERG08dll(0, 0)  # Just check status, don't modify
            if ierr == 0:
                models['gerg_2008_available'] = True
        except:
            pass
        
        # Get mixing rules if it's a mixture
        if sum(1 for comp in z if comp > 0) > 1:
            try:
                # Get first two components with non-zero fractions
                i, j = [idx+1 for idx, val in enumerate(z) if val > 0][:2]
                hmixrule = RP.GETKTVdll(i, j).hmxrul
                models['mixing_rule'] = trim(hmixrule)
                
                # Get binary interaction parameters
                hmodij = RP.GETKTVdll(i, j).hmodij
                fij = RP.GETKTVdll(i, j).fij
                
                models['binary_interaction'] = {
                    'model': trim(hmodij),
                    'parameters': [float(f) for f in fij]
                }
            except Exception as e:
                models['mixing_rule_error'] = str(e)
        
        # Get information about components
        models['components'] = []
        for i, comp_fraction in enumerate(z):
            if comp_fraction > 0:
                comp_idx = i + 1  # REFPROP uses 1-based indexing
                try:
                    # Get component name
                    hnam, hn80, hcasn = RP.NAMEdll(comp_idx)
                    
                    # Get component info
                    wmm, Ttrp, Tnbpt, Tc, Pc, Dc, Zc, acf, dip, Rgas = RP.INFOdll(comp_idx)
                    
                    # Get component EOS
                    hcode, hcite = RP.GETMODdll(comp_idx, 'EOS')
                    
                    # Store component information
                    models['components'].append({
                        'name': trim(hnam),
                        'cas_number': trim(hcasn),
                        'molar_mass': wmm,
                        'critical_temperature': Tc,
                        'critical_pressure': Pc,
                        'critical_density': Dc,
                        'acentric_factor': acf,
                        'eos_model': {
                            'code': trim(hcode),
                            'reference': trim(hcite)
                        }
                    })
                except Exception as e:
                    models['components'].append({
                        'index': comp_idx,
                        'error': str(e)
                    })
                
    except Exception as e:
        models['error'] = str(e)
        
    return models

@models_info_bp.route('/models_info', methods=['POST'])
def models_info():
    """
    Endpoint to get information about thermodynamic models used for a composition
    """
    try:
        data = request.get_json(force=True)
        
        # Validate the composition field exists
        if 'composition' not in data:
            return jsonify({'error': 'Missing composition field'}), 400
                
        # Validate composition
        if not validate_composition(data['composition']):
            return jsonify({'error': 'Invalid composition - fractions must sum to 1'}), 400

        # Setup the mixture in REFPROP
        fluid_string = '|'.join(f"{comp['fluid']}.FLD" for comp in data['composition'])
        z = [comp['fraction'] for comp in data['composition']] + [0] * (20 - len(data['composition']))
        
        ierr, herr = RP.SETUPdll(len(data['composition']), fluid_string, 'HMX.BNC', 'DEF')
        if ierr > 0:
            return jsonify({'error': f"Error setting up mixture: {herr}"}), 400
        
        # Get model information
        models = get_model_info(z)
        
        # Return the model information
        return jsonify({
            'composition': data['composition'],
            'models': models
        })
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500