from flask import request, jsonify
import sys
import traceback
from typing import List, Dict, Any, Tuple

from API.endpoints import critical_point_bp
from API.refprop_setup import RP
from API.unit_converter import UnitConverter

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

def find_critical_point(z: List[float], units_system: str = 'SI') -> Dict[str, Any]:
    """
    Find the critical point for the given composition.
    
    Args:
        z: Composition array
        units_system: Unit system to use ('SI' or 'CGS')
    
    Returns:
        Dictionary of critical properties with values and units
    """
    # Initialize unit converter
    converter = UnitConverter()
    
    # Get molecular weight for unit conversions
    wmm = RP.WMOLdll(z)
    
    # Get critical properties
    Tc, Pc, Dc, ierr, herr = RP.CRITPdll(z)
    if ierr > 0:
        raise ValueError(f"Error in CRITPdll: {herr}")
        
    # For pure fluids, get more detailed critical information
    if sum(1 for zz in z if zz > 0) == 1:
        # Get the index of the single component
        icomp = next(i for i, zz in enumerate(z) if zz > 0) + 1
        
        # Get more detailed information from INFO function
        wmm, Ttrp, Tnbpt, Tc, Pc_info, Dc_info, Zc, acf, dip, Rgas = RP.INFOdll(icomp)
    
    # Calculate compressibility factor at critical point
    Zc = Pc * 100 / (Dc * 8.31446261815324 * Tc)
    
    # Get additional properties at critical point if possible
    try:
        # Try to get transport properties at critical point
        eta, tcx, ierr_trn, herr_trn = RP.TRNPRPdll(Tc, Dc, z)
        
        # Get thermodynamic properties at critical point
        P, e, h, s, Cv, Cp, w, hjt = RP.THERMdll(Tc, Dc, z)
        
        # Additional properties including acentric factor, etc. could be added here
    except Exception as e:
        # Critical point is often a difficult region for property calculation
        # so we'll handle exceptions gracefully
        eta = tcx = None
        P = Pc
        e = h = s = Cv = Cp = w = hjt = None
        
    # Build raw properties dictionary
    raw_properties = {
        'critical_temperature': Tc,
        'critical_pressure': Pc / 100,  # Convert from kPa to bar for SI
        'critical_density': Dc,
        'critical_compressibility': Zc,
        'critical_enthalpy': h,
        'critical_entropy': s,
        'critical_viscosity': eta,
        'critical_thermal_conductivity': tcx,
        'critical_sound_speed': w
    }
    
    # Convert properties to requested unit system
    properties = {}
    for prop_id, value in raw_properties.items():
        if value is not None:  # Skip undefined properties
            try:
                # Special handling for properties that aren't in the converter
                if prop_id == 'critical_compressibility':
                    properties[prop_id] = {'value': float(value), 'unit': 'dimensionless'}
                else:
                    # Try to derive the standard property name by removing 'critical_' prefix
                    std_prop_id = prop_id.replace('critical_', '')
                    properties[prop_id] = converter.convert_property(
                        std_prop_id, float(value), wmm, 'SI', units_system
                    )
            except Exception as e:
                print(f"Warning: Could not convert property {prop_id}: {str(e)}")
                properties[prop_id] = {'value': float(value) if value is not None else None, 'unit': 'unknown'}
    
    return properties

@critical_point_bp.route('/critical_point', methods=['POST'])
def critical_point():
    try:
        data = request.get_json(force=True)
        
        # Validate required fields
        if 'composition' not in data:
            return jsonify({'error': 'Missing composition field'}), 400
                
        # Extract calculation settings
        units_system = data.get('units_system', 'SI')  # Default to SI
        
        # Validate composition
        if not validate_composition(data['composition']):
            return jsonify({'error': 'Invalid composition - fractions must sum to 1'}), 400

        # Setup mixture
        z = setup_mixture(data['composition'])
        
        # Debug log
        print(f"Calculating critical point for composition: {data['composition']}")
        
        # Calculate critical point
        critical_props = find_critical_point(z, units_system)
            
        # Return results
        return jsonify({'critical_point': critical_props})
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500