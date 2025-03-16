from flask import Response
import numpy as np
import re
import sys
import traceback

def format_olga_tab(x_vars, y_vars, results, composition, molar_mass, endpoint_type='pt_flash'):
    """
    Format calculation results in OLGA TAB format using generic grid variables.
    
    Args:
        x_vars (dict): Dictionary with 'range', 'resolution', and optional 'values' for x-axis variable
        y_vars (dict): Dictionary with 'range', 'resolution', and optional 'values' for y-axis variable
        results (list): List of result dictionaries from the calculation
        composition (list): List of fluid compositions (dict with 'fluid' and 'fraction')
        molar_mass (float): Molar mass of the mixture in g/mol
        endpoint_type (str): Type of endpoint ('pt_flash', 'ph_flash', 'ts_flash')
        
    Returns:
        Response: Flask response object with OLGA TAB formatted content
    """
    try:
        # Determine the x and y variable names based on endpoint type
        if endpoint_type == 'pt_flash':
            x_name = 'pressure'
            y_name = 'temperature'
            x_header = 'Pressure (Pa)'
            y_header = 'Temperature (C)'
            x_idx_name = 'p_idx'
            y_idx_name = 't_idx'
            # For PT flash, convert bar to Pa for OLGA
            x_multiplier = 1e5  # bar to Pa
            y_multiplier = 1.0  # already in C
        elif endpoint_type == 'ph_flash':
            x_name = 'pressure'
            y_name = 'enthalpy'
            x_header = 'Pressure (Pa)'
            y_header = 'Enthalpy (J/mol)'
            x_idx_name = 'p_idx'
            y_idx_name = 'h_idx'
            # For PH flash, convert bar to Pa for OLGA
            x_multiplier = 1e5  # bar to Pa
            y_multiplier = 1.0  # J/mol
        elif endpoint_type == 'ts_flash':
            x_name = 'temperature'
            y_name = 'entropy'
            x_header = 'Temperature (C)'
            y_header = 'Entropy (J/mol-K)'
            x_idx_name = 't_idx'
            y_idx_name = 's_idx'
            x_multiplier = 1.0  # already in C
            y_multiplier = 1.0  # J/mol-K
        else:
            # Default to PT flash
            x_name = 'pressure'
            y_name = 'temperature'
            x_header = 'Pressure (Pa)'
            y_header = 'Temperature (C)'
            x_idx_name = 'p_idx'
            y_idx_name = 't_idx'
            x_multiplier = 1e5  # bar to Pa
            y_multiplier = 1.0  # C

        # Use provided grid values if available, otherwise generate from range
        if 'values' in x_vars:
            # Use the provided grid values directly
            x_grid = x_vars['values']
            nx_grid = len(x_grid)
        else:
            # Extract range and resolution with robust error handling
            try:
                # Get x-axis range with defaults
                x_range = x_vars.get('range', {})
                if not x_range:
                    x_range = {'from': 1.0, 'to': 100.0}
                    
                x_from = float(x_range.get('from', 1.0))
                x_to = float(x_range.get('to', 100.0))
                
                # Ensure valid range
                if x_to <= x_from:
                    x_to = x_from + 10.0
                    
                # Get resolution with default
                x_resolution = float(x_vars.get('resolution', 10.0))
                if x_resolution <= 0:
                    x_resolution = 10.0
                    
                # Generate regular grid
                x_grid = np.arange(x_from, x_to + x_resolution, x_resolution)
                nx_grid = len(x_grid)
                    
            except (TypeError, ValueError) as e:
                # If there's any issue parsing the ranges, use defaults
                print(f"Error parsing grid ranges: {e}. Using defaults.")
                x_from, x_to, x_resolution = 1.0, 100.0, 10.0
                x_grid = np.linspace(x_from, x_to, 10)
                nx_grid = len(x_grid)

        # Same for y-axis
        if 'values' in y_vars:
            # Use the provided grid values directly
            y_grid = y_vars['values']
            ny_grid = len(y_grid)
        else:
            try:
                y_range = y_vars.get('range', {})
                if not y_range:
                    y_range = {'from': 0.0, 'to': 100.0}
                    
                y_from = float(y_range.get('from', 0.0))
                y_to = float(y_range.get('to', 100.0))
                
                # Ensure valid range
                if y_to <= y_from:
                    y_to = y_from + 10.0
                    
                # Get resolution with default
                y_resolution = float(y_vars.get('resolution', 5.0))
                if y_resolution <= 0:
                    y_resolution = 5.0
                    
                # Generate regular grid
                y_grid = np.arange(y_from, y_to + y_resolution, y_resolution)
                ny_grid = len(y_grid)
                
            except (TypeError, ValueError) as e:
                # If there's any issue parsing the ranges, use defaults
                print(f"Error parsing grid ranges: {e}. Using defaults.")
                y_from, y_to, y_resolution = 0.0, 100.0, 5.0
                y_grid = np.linspace(y_from, y_to, 20)
                ny_grid = len(y_grid)
            
        # Debug info
        print(f"Generating OLGA TAB with {nx_grid} {x_name} points and {ny_grid} {y_name} points")
        
        # Compose fluid description
        fluid_desc = " ".join([f"{comp['fluid']}-{comp['fraction']:.4f}" for comp in composition])
        
        # Create OLGA TAB header - customized for the grid variables
        olga_tab = f"'Span-Wagner EOS {fluid_desc} ({x_name.upper()}-{y_name.upper()} grid)'\n"
        olga_tab += f"   {nx_grid}  {ny_grid}    .276731E-08\n"
        
        # X-axis grid with appropriate multiplier
        x_grid_formatted = [x * x_multiplier for x in x_grid]
        olga_tab += format_grid_values(x_grid_formatted)
        
        # Y-axis grid with appropriate multiplier
        y_grid_formatted = [y * y_multiplier for y in y_grid]
        olga_tab += format_grid_values(y_grid_formatted)
        
        # Define OLGA TAB property headers and corresponding API properties
        olga_properties = [
            {'name': 'LIQUID DENSITY (KG/M3)', 'key': 'liquid_density', 'converter': lambda x: x * molar_mass / 1000.0},
            {'name': 'GAS DENSITY (KG/M3)', 'key': 'vapor_density', 'converter': lambda x: x * molar_mass / 1000.0},
            {'name': 'PRES. DERIV. OF LIQUID DENS.', 'key': 'dDdP_liquid', 'converter': lambda x: x * molar_mass / 1000.0},
            {'name': 'PRES. DERIV. OF GAS DENS.', 'key': 'dDdP_vapor', 'converter': lambda x: x * molar_mass / 1000.0},
            {'name': 'TEMP. DERIV. OF LIQUID DENS.', 'key': 'dDdT_liquid', 'converter': lambda x: x * molar_mass / 1000.0},
            {'name': 'TEMP. DERIV. OF GAS DENS.', 'key': 'dDdT_vapor', 'converter': lambda x: x * molar_mass / 1000.0},
            {'name': 'GAS MASS FRACTION OF GAS +LIQUID', 'key': 'vapor_fraction', 'converter': convert_vapor_fraction_to_mass},
            {'name': 'WATER MASS FRACTION OF LIQUID', 'key': 'water_liquid', 'converter': lambda x: x},
            {'name': 'WATER MASS FRACTION OF GAS', 'key': 'water_vapor', 'converter': lambda x: x},
            {'name': 'LIQUID VISCOSITY (N S/M2)', 'key': 'liquid_viscosity', 'converter': lambda x: x * 1e-6},
            {'name': 'GAS VISCOSITY (N S/M2)', 'key': 'vapor_viscosity', 'converter': lambda x: x * 1e-6},
            {'name': 'LIQUID SPECIFIC HEAT (J/KG K)', 'key': 'liquid_cp', 'converter': lambda x: x * 1000.0 / molar_mass},
            {'name': 'GAS SPECIFIC HEAT (J/KG K)', 'key': 'vapor_cp', 'converter': lambda x: x * 1000.0 / molar_mass},
            {'name': 'LIQUID ENTHALPY (J/KG)', 'key': 'liquid_enthalpy', 'converter': lambda x: x * 1000.0 / molar_mass},
            {'name': 'GAS ENTHALPY (J/KG)', 'key': 'vapor_enthalpy', 'converter': lambda x: x * 1000.0 / molar_mass},
            {'name': 'LIQUID THERMAL COND. (W/M K)', 'key': 'liquid_thermal_conductivity', 'converter': lambda x: x},
            {'name': 'GAS THERMAL COND. (W/M K)', 'key': 'vapor_thermal_conductivity', 'converter': lambda x: x},
            {'name': 'SURFACE TENSION GAS/LIQUID (N/M)', 'key': 'surface_tension', 'converter': lambda x: x},
            {'name': 'LIQUID ENTROPY (J/KG/C)', 'key': 'liquid_entropy', 'converter': lambda x: x * 1000.0 / molar_mass},
            {'name': 'GAS ENTROPY (J/KG/C)', 'key': 'vapor_entropy', 'converter': lambda x: x * 1000.0 / molar_mass}
        ]
        
        # Extract and format properties
        property_arrays = {}
        for prop in olga_properties:
            # Initialize with zeros for all properties
            property_arrays[prop['name']] = np.zeros((nx_grid, ny_grid))
        
        # For the specific flash type, use the correct indices
        mapped_points = 0
        
        # Populate property arrays
        for result in results:
            try:
                # Get indices based on the endpoint type
                x_idx = result.get(x_idx_name)
                y_idx = result.get(y_idx_name)
                
                # If specific indices are not found, try to derive them
                if x_idx is None or y_idx is None:
                    if 'index' in result:
                        # Try to derive indices from index (assuming row-major order)
                        idx = int(result['index'])
                        x_idx = idx // ny_grid
                        y_idx = idx % ny_grid
                    else:
                        # Try to map by values
                        if x_name in result and y_name in result:
                            x_val = get_value_from_field(result[x_name])
                            y_val = get_value_from_field(result[y_name])
                            x_idx = find_nearest_index(x_grid, x_val)
                            y_idx = find_nearest_index(y_grid, y_val)
                        else:
                            # Can't determine position, skip this result
                            continue
                
                # Convert to integers if needed
                try:
                    x_idx = int(x_idx)
                    y_idx = int(y_idx)
                except (TypeError, ValueError):
                    continue
                    
                # Ensure indices are within bounds
                if 0 <= x_idx < nx_grid and 0 <= y_idx < ny_grid:
                    mapped_points += 1
                    
                    # Extract phase information
                    phase = get_phase_from_result(result)
                    is_two_phase = phase == 'two-phase'
                    
                    # Process each property
                    for prop in olga_properties:
                        try:
                            value = extract_property_value(result, prop['key'], phase, composition, molar_mass)
                            if value is not None:
                                converted_value = prop['converter'](value)
                                property_arrays[prop['name']][x_idx, y_idx] = converted_value
                            else:
                                # Default for missing values based on property
                                if prop['name'] == 'GAS MASS FRACTION OF GAS +LIQUID':
                                    if phase == 'vapor':
                                        property_arrays[prop['name']][x_idx, y_idx] = 1.0
                                    elif phase == 'liquid':
                                        property_arrays[prop['name']][x_idx, y_idx] = 0.0
                        except Exception as e:
                            print(f"Error processing property {prop['key']} for point ({x_idx}, {y_idx}): {e}")
            except Exception as e:
                print(f"Error processing result: {e}")
                continue
                
        print(f"Successfully mapped {mapped_points} points to the grid")
        
        # Add properties to OLGA TAB
        for prop in olga_properties:
            olga_tab += f" {prop['name']}                \n"
            olga_tab += format_grid_values(property_arrays[prop['name']].flatten())
        
        # Return as a Response object with proper MIME type
        filename = f"{endpoint_type.replace('_', '-')}.tab"
        return Response(olga_tab, mimetype='text/plain', headers={
            'Content-Disposition': f'attachment; filename={filename}'
        })
        
    except Exception as e:
        traceback.print_exc(file=sys.stderr)
        print(f"Error in format_olga_tab: {str(e)}")
        # Return an error Response to prevent socket hang-up
        error_msg = f"Error generating OLGA TAB format: {str(e)}\n\n{traceback.format_exc()}"
        return Response(error_msg, status=500, mimetype='text/plain')

def find_nearest_index(array, value):
    """Find the index of the value closest to the target in the array"""
    try:
        array = np.asarray(array)
        idx = (np.abs(array - float(value))).argmin()
        return idx
    except Exception as e:
        print(f"Error in find_nearest_index: {e}")
        return 0  # Return 0 as a fallback

def extract_property_value(result, key, phase, composition, molar_mass):
    """
    Extract property value from result, handling various property naming schemes.
    
    Args:
        result (dict): Result dictionary from calculation
        key (str): Property key to extract
        phase (str): Phase state ('liquid', 'vapor', 'two-phase', etc.)
        composition (list): List of fluid compositions
        molar_mass (float): Molar mass of the mixture in g/mol
        
    Returns:
        float or None: Extracted property value or None if not found
    """
    try:
        # Direct property match
        if key in result:
            return get_value_from_field(result[key])
        
        # Handle phase-specific properties
        if key.startswith('liquid_') and (phase == 'liquid' or phase == 'two-phase'):
            base_key = key[7:]  # Remove 'liquid_' prefix
            if base_key in result:
                return get_value_from_field(result[base_key])
        
        if key.startswith('vapor_') and (phase == 'vapor' or phase == 'two-phase'):
            base_key = key[6:]  # Remove 'vapor_' prefix
            if base_key in result:
                return get_value_from_field(result[base_key])
        
        # Special handling for specific properties
        if key == 'liquid_density' and 'Dl' in result:
            return get_value_from_field(result['Dl'])
        elif key == 'liquid_density' and 'liquid_density' in result:
            return get_value_from_field(result['liquid_density'])
        
        if key == 'vapor_density' and 'Dv' in result:
            return get_value_from_field(result['Dv'])
        elif key == 'vapor_density' and 'vapor_density' in result:
            return get_value_from_field(result['vapor_density'])
        
        if key == 'liquid_cp' and 'cp' in result and phase == 'liquid':
            return get_value_from_field(result['cp'])
        
        if key == 'vapor_cp' and 'cp' in result and phase == 'vapor':
            return get_value_from_field(result['cp'])
        
        if key == 'liquid_enthalpy' and 'enthalpy' in result and phase == 'liquid':
            return get_value_from_field(result['enthalpy'])
        
        if key == 'vapor_enthalpy' and 'enthalpy' in result and phase == 'vapor':
            return get_value_from_field(result['enthalpy'])
        
        if key == 'liquid_entropy' and 'entropy' in result and phase == 'liquid':
            return get_value_from_field(result['entropy'])
        
        if key == 'vapor_entropy' and 'entropy' in result and phase == 'vapor':
            return get_value_from_field(result['entropy'])
        
        if key == 'liquid_viscosity' and 'viscosity' in result and phase == 'liquid':
            return get_value_from_field(result['viscosity'])
        
        if key == 'vapor_viscosity' and 'viscosity' in result and phase == 'vapor':
            return get_value_from_field(result['viscosity'])
        
        if key == 'liquid_thermal_conductivity' and 'thermal_conductivity' in result and phase == 'liquid':
            return get_value_from_field(result['thermal_conductivity'])
        
        if key == 'vapor_thermal_conductivity' and 'thermal_conductivity' in result and phase == 'vapor':
            return get_value_from_field(result['thermal_conductivity'])
        
        if key == 'dDdP_liquid' and 'dDdP' in result and phase == 'liquid':
            return get_value_from_field(result['dDdP'])
        
        if key == 'dDdP_vapor' and 'dDdP' in result and phase == 'vapor':
            return get_value_from_field(result['dDdP'])
        
        if key == 'dDdT_liquid' and 'dDdT' in result and phase == 'liquid':
            return get_value_from_field(result['dDdT'])
        
        if key == 'dDdT_vapor' and 'dDdT' in result and phase == 'vapor':
            return get_value_from_field(result['dDdT'])
        
        if key == 'water_liquid' and 'x' in result:
            x = get_value_from_field(result['x'])
            if isinstance(x, list) and len(x) > 0:
                # Find water component by name
                water_idx = next((i for i, comp in enumerate(composition) if "WATER" in comp['fluid'].upper()), -1)
                return x[water_idx] if water_idx >= 0 and water_idx < len(x) else 0.0
            return 0.0
        
        if key == 'water_vapor' and 'y' in result:
            y = get_value_from_field(result['y'])
            if isinstance(y, list) and len(y) > 0:
                # Find water component by name
                water_idx = next((i for i, comp in enumerate(composition) if "WATER" in comp['fluid'].upper()), -1)
                return y[water_idx] if water_idx >= 0 and water_idx < len(y) else 0.0
            return 0.0
        
        # Not found
        return None
    except Exception as e:
        print(f"Error extracting property {key}: {e}")
        return None

def get_value_from_field(field):
    """Extract value from field, handling both direct values and dictionaries with 'value' key."""
    try:
        if isinstance(field, dict) and 'value' in field:
            return field['value']
        return field
    except Exception as e:
        print(f"Error getting value from field: {e}")
        return None

def get_phase_from_result(result):
    """Determine the phase state from the result dictionary."""
    try:
        # Try to get phase from 'phase' field
        if 'phase' in result:
            phase_value = get_value_from_field(result['phase'])
            if isinstance(phase_value, str):
                phase_lower = phase_value.lower()
                if 'liquid' in phase_lower:
                    return 'liquid'
                elif 'vapor' in phase_lower or 'gas' in phase_lower:
                    return 'vapor'
                elif 'two' in phase_lower or 'mixed' in phase_lower:
                    return 'two-phase'
        
        # Try to determine from vapor fraction
        vapor_fraction_key = next((k for k in ['q', 'vapor_fraction'] if k in result), None)
        if vapor_fraction_key:
            q = get_value_from_field(result[vapor_fraction_key])
            if q == 0:
                return 'liquid'
            elif q == 1:
                return 'vapor'
            elif 0 < q < 1:
                return 'two-phase'
        
        # Default to unknown
        return 'unknown'
    except Exception as e:
        print(f"Error determining phase: {e}")
        return 'unknown'

def convert_vapor_fraction_to_mass(q, result=None, composition=None, molar_mass=None):
    """
    Convert vapor fraction from mole basis to mass basis.
    
    In a simple case, this would need molar masses of each component,
    but for simplicity we'll return the mole fraction directly.
    A proper implementation would use liquid and vapor compositions
    along with component molar masses to do the conversion.
    """
    try:
        # Special cases for quality
        if q is None:
            return 0.0
        
        # For flagging regions outside phase envelope
        if q == 998:  # Special flag for vapor
            return 998.0
        elif q == -998:  # Special flag for liquid
            return -998.0
        elif q == 999:  # Special flag for supercritical
            return 999.0
        elif q < 0 or q > 1:  # Special flags or errors
            return 998.0  # Return as vapor by default
        
        # For a proper implementation:
        # q_mass = (q * M_vapor) / ((1-q) * M_liquid + q * M_vapor)
        # where M_vapor and M_liquid are the molar masses of vapor and liquid phases
        
        return q  # Simplified implementation
    except Exception as e:
        print(f"Error converting vapor fraction: {e}")
        return 0.0  # Default to liquid

def format_grid_values(values, values_per_line=5):
    """
    Format grid values for OLGA TAB format.
    
    Args:
        values (list): List of values to format
        values_per_line (int): Number of values per line
        
    Returns:
        str: Formatted grid values
    """
    try:
        formatted = ""
        values_array = np.asarray(values).flatten()
        
        for i in range(0, len(values_array), values_per_line):
            line_values = values_array[i:i + values_per_line]
            formatted += "    " + "    ".join([f".{value:.6E}" for value in line_values]) + "\n"
        
        return formatted
    except Exception as e:
        print(f"Error formatting grid values: {e}")
        # Return a minimal valid grid as fallback
        return "    .0.000000E+00\n"