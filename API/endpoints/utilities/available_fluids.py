from flask import Blueprint, jsonify, request
import os
import sys
import traceback
from pathlib import Path
from typing import List, Dict, Any

# Create a blueprint for the available_fluids endpoint
available_fluids_bp = Blueprint('available_fluids', __name__)

def parse_fluid_file(file_path: str) -> Dict[str, Any]:
    """
    Parse a REFPROP fluid file (.FLD) to extract information.
    
    Args:
        file_path: Path to the fluid file
        
    Returns:
        Dictionary with parsed fluid information
    """
    fluid_info = {}
    
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
            
            # Extract basic information from the first lines
            if len(lines) >= 14:
                # Parse short name (remove any comments)
                short_name_parts = lines[0].split('!')
                fluid_info['short_name'] = short_name_parts[0].strip()
                if len(short_name_parts) > 1:
                    fluid_info['short_name_comment'] = short_name_parts[1].strip()
                
                # Parse CAS number
                cas_parts = lines[1].split('!')
                fluid_info['cas_number'] = cas_parts[0].strip()
                
                # Parse full name
                full_name_parts = lines[2].split('!')
                fluid_info['full_name'] = full_name_parts[0].strip()
                
                # Parse chemical formula
                formula_parts = lines[3].split('!')
                fluid_info['chemical_formula'] = formula_parts[0].strip()
                if len(formula_parts) > 1 and '{' in formula_parts[1]:
                    # Extract formula in braces like {C2H6O2}
                    import re
                    braces_match = re.search(r'\{([^}]+)\}', formula_parts[1])
                    if braces_match:
                        fluid_info['formula_alt'] = braces_match.group(1)
                
                # Parse synonym
                synonym_parts = lines[4].split('!')
                fluid_info['synonym'] = synonym_parts[0].strip()
                
                # Extract numeric properties with improved error handling
                # Molar mass
                try:
                    molar_mass_parts = lines[5].split('!')
                    fluid_info['molar_mass'] = float(molar_mass_parts[0].strip())
                    if len(molar_mass_parts) > 1 and '[g/mol]' in molar_mass_parts[1]:
                        fluid_info['molar_mass_unit'] = 'g/mol'
                except (ValueError, IndexError):
                    fluid_info['molar_mass'] = None
                
                # Triple point temperature
                try:
                    tp_parts = lines[6].split('!')
                    fluid_info['triple_point_temperature'] = float(tp_parts[0].strip())
                    if len(tp_parts) > 1 and '[K]' in tp_parts[1]:
                        fluid_info['triple_point_temperature_unit'] = 'K'
                except (ValueError, IndexError):
                    fluid_info['triple_point_temperature'] = None
                
                # Normal boiling point
                try:
                    bp_parts = lines[7].split('!')
                    fluid_info['normal_boiling_point'] = float(bp_parts[0].strip())
                    if len(bp_parts) > 1 and '[K]' in bp_parts[1]:
                        fluid_info['normal_boiling_point_unit'] = 'K'
                except (ValueError, IndexError):
                    fluid_info['normal_boiling_point'] = None
                
                # Critical temperature
                try:
                    ct_parts = lines[8].split('!')
                    fluid_info['critical_temperature'] = float(ct_parts[0].strip())
                    if len(ct_parts) > 1 and '[K]' in ct_parts[1]:
                        fluid_info['critical_temperature_unit'] = 'K'
                except (ValueError, IndexError):
                    fluid_info['critical_temperature'] = None
                
                # Critical pressure
                try:
                    cp_parts = lines[9].split('!')
                    fluid_info['critical_pressure'] = float(cp_parts[0].strip())
                    if len(cp_parts) > 1 and '[kPa]' in cp_parts[1]:
                        fluid_info['critical_pressure_unit'] = 'kPa'
                except (ValueError, IndexError):
                    fluid_info['critical_pressure'] = None
                
                # Critical density
                try:
                    cd_parts = lines[10].split('!')
                    fluid_info['critical_density'] = float(cd_parts[0].strip())
                    if len(cd_parts) > 1 and '[mol/L]' in cd_parts[1]:
                        fluid_info['critical_density_unit'] = 'mol/L'
                except (ValueError, IndexError):
                    fluid_info['critical_density'] = None
                
                # Acentric factor
                try:
                    af_parts = lines[11].split('!')
                    fluid_info['acentric_factor'] = float(af_parts[0].strip())
                except (ValueError, IndexError):
                    fluid_info['acentric_factor'] = None
                
                # Try to extract dipole moment
                try:
                    dipole_line = lines[12].strip()
                    dipole_parts = dipole_line.split('!')
                    dipole_value = dipole_parts[0].strip()
                    fluid_info['dipole_moment'] = float(dipole_value)
                    
                    # Try to extract dipole moment reference if present
                    if len(dipole_parts) > 1:
                        fluid_info['dipole_moment_reference'] = dipole_parts[1].strip()
                        if 'Debye' in dipole_parts[1]:
                            fluid_info['dipole_moment_unit'] = 'Debye'
                except (ValueError, IndexError):
                    fluid_info['dipole_moment'] = None
                
                # Reference state
                ref_parts = lines[13].split('!')
                fluid_info['reference_state'] = ref_parts[0].strip()
                
                # Try to extract version number
                if len(lines) > 14:
                    try:
                        ver_parts = lines[14].split('!')
                        fluid_info['version'] = float(ver_parts[0].strip())
                    except (ValueError, IndexError):
                        fluid_info['version'] = None
            
            # Extract additional metadata from tagged lines
            for line in lines:
                line = line.strip()
                if ':UN:' in line:
                    un_parts = line.split(':UN:')
                    if len(un_parts) > 0:
                        un_value = un_parts[0].strip()
                        if un_value != '????':
                            fluid_info['un_number'] = un_value
                
                elif ':Family:' in line:
                    family_parts = line.split(':Family:')
                    if len(family_parts) > 0:
                        family_value = family_parts[0].strip()
                        if family_value != '????':
                            fluid_info['family'] = family_value
                
                elif ':Heat:' in line:
                    heat_parts = line.split(':Heat:')
                    if len(heat_parts) > 0:
                        heat_value = heat_parts[0].strip()
                        if heat_value != '????':
                            try:
                                fluid_info['heating_value'] = float(heat_value)
                                if '[kJ/mol]' in line:
                                    fluid_info['heating_value_unit'] = 'kJ/mol'
                            except ValueError:
                                fluid_info['heating_value'] = heat_value
                
                elif ':InChi:' in line:
                    inchi_parts = line.split(':InChi:')
                    if len(inchi_parts) > 0:
                        inchi_value = inchi_parts[0].strip()
                        if inchi_value != '????':
                            fluid_info['inchi'] = inchi_value
                
                elif ':InChiKey:' in line:
                    inchikey_parts = line.split(':InChiKey:')
                    if len(inchikey_parts) > 0:
                        inchikey_value = inchikey_parts[0].strip()
                        if inchikey_value != '????':
                            fluid_info['inchikey'] = inchikey_value
                
                elif ':AltID:' in line:
                    altid_parts = line.split(':AltID:')
                    if len(altid_parts) > 0:
                        altid_value = altid_parts[0].strip()
                        if altid_value != '????':
                            fluid_info['alt_id'] = altid_value
                
                elif ':Hash:' in line:
                    hash_parts = line.split(':Hash:')
                    if len(hash_parts) > 0:
                        hash_value = hash_parts[0].strip()
                        if hash_value != '????':
                            fluid_info['hash'] = hash_value
                            
                elif ':Safety:' in line:
                    safety_parts = line.split(':Safety:')
                    if len(safety_parts) > 0:
                        safety_value = safety_parts[0].strip()
                        if safety_value != '????':
                            fluid_info['safety_group'] = safety_value
        
        # Add file name without extension for API reference
        fluid_info['fluid_id'] = os.path.splitext(os.path.basename(file_path))[0]
        
        return fluid_info
    
    except Exception as e:
        print(f"Error parsing fluid file {file_path}: {str(e)}", file=sys.stderr)
        traceback.print_exc()
        return {'error': str(e), 'file_name': os.path.basename(file_path)}

def get_fluids_directory():
    """Get the FLUIDS directory path"""
    base_dir = Path(__file__).resolve().parent.parent.parent.parent
    return base_dir / 'FLUIDS'

def get_fluid_path(fluid_id):
    """Get path to a specific fluid file by ID"""
    fluids_dir = get_fluids_directory()
    potential_paths = [
        fluids_dir / f"{fluid_id}.FLD",
        fluids_dir / f"{fluid_id.upper()}.FLD",
        fluids_dir / f"{fluid_id.lower()}.FLD"
    ]
    
    for path in potential_paths:
        if path.exists():
            return path
    
    return None

@available_fluids_bp.route('/available_fluids', methods=['GET'])
def available_fluids():
    """
    Get a list of all available fluids and their information.
    
    Query parameters:
    - short_only: If 'true', returns only basic information (default: false)
    - search: Filter fluids by name, formula, or CAS number (case-insensitive)
    """
    try:
        # Parse query parameters
        short_only = request.args.get('short_only', 'false').lower() == 'true'
        search_term = request.args.get('search', '').lower()
        
        # Get FLUIDS directory
        fluids_dir = get_fluids_directory()
        
        if not fluids_dir.exists():
            return jsonify({'error': f'FLUIDS directory not found at {str(fluids_dir)}'}), 404
        
        # Get all .FLD files in the FLUIDS directory
        fld_files = list(fluids_dir.glob('*.FLD'))
        
        # Parse each fluid file
        fluids = []
        for fld_file in fld_files:
            fluid_info = parse_fluid_file(str(fld_file))
            
            # Filter by search term if provided
            if search_term:
                searchable_fields = [
                    fluid_info.get('short_name', '').lower(),
                    fluid_info.get('full_name', '').lower(),
                    fluid_info.get('chemical_formula', '').lower(),
                    fluid_info.get('cas_number', '').lower(),
                    fluid_info.get('synonym', '').lower()
                ]
                if not any(search_term in field for field in searchable_fields):
                    continue
            
            # If short_only is true, only include basic information
            if short_only:
                fluids.append({
                    'fluid_id': fluid_info.get('fluid_id'),
                    'short_name': fluid_info.get('short_name'),
                    'full_name': fluid_info.get('full_name'),
                    'chemical_formula': fluid_info.get('chemical_formula'),
                    'cas_number': fluid_info.get('cas_number')
                })
            else:
                fluids.append(fluid_info)
        
        # Sort fluids by short name
        fluids.sort(key=lambda x: x.get('short_name', '').lower())
        
        return jsonify({
            'count': len(fluids),
            'fluids': fluids
        })
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500
        
@available_fluids_bp.route('/available_fluids/<fluid_id>', methods=['GET'])
def get_fluid_by_id(fluid_id):
    """
    Get detailed information about a specific fluid by ID.
    
    Parameters:
    - fluid_id: The ID of the fluid to retrieve
    """
    try:
        # Find the fluid file
        fluid_path = get_fluid_path(fluid_id)
        
        if not fluid_path:
            return jsonify({
                'error': f'Fluid not found: {fluid_id}',
                'suggestion': 'Use /available_fluids to get a list of all available fluids'
            }), 404
        
        # Parse the fluid file
        fluid_info = parse_fluid_file(str(fluid_path))
        
        # Return the fluid information
        return jsonify(fluid_info)
        
    except Exception as e:
        print(f"Error getting fluid {fluid_id}:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500