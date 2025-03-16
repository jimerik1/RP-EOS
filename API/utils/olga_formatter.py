from flask import Response
import numpy as np
import re

def format_olga_tab(pressure_vars, temperature_vars, results, composition, molar_mass):
    """
    Format calculation results in OLGA TAB format.
    
    Args:
        pressure_vars (dict): Dictionary with 'range' and 'resolution' for pressure
        temperature_vars (dict): Dictionary with 'range' and 'resolution' for temperature
        results (list): List of result dictionaries from the calculation
        composition (list): List of fluid compositions (dict with 'fluid' and 'fraction')
        molar_mass (float): Molar mass of the mixture in g/mol
        
    Returns:
        Response: Flask response object with OLGA TAB formatted content
    """
    try:
        # Extract range and resolution
        pressure_range = pressure_vars.get('range', {})
        pressure_resolution = pressure_vars.get('resolution')
        
        temperature_range = temperature_vars.get('range', {})
        temperature_resolution = temperature_vars.get('resolution')
        
        # Generate pressure and temperature grids
        pressures = np.arange(
            float(pressure_range.get('from', 0)),
            float(pressure_range.get('to', 100)) + float(pressure_resolution),
            float(pressure_resolution)
        )
        temperatures = np.arange(
            float(temperature_range.get('from', 0)),
            float(temperature_range.get('to', 100)) + float(temperature_resolution),
            float(temperature_resolution)
        )
        
        np_grid = len(pressures)
        nt_grid = len(temperatures)
        
        # Compose fluid description
        fluid_desc = " ".join([f"{comp['fluid']}-{comp['fraction']:.4f}" for comp in composition])
        
        # Create OLGA TAB header
        olga_tab = f"'Span-Wagner EOS {fluid_desc}'\n"
        olga_tab += f"   {np_grid}  {nt_grid}    .276731E-08\n"
        
        # Pressure grid in Pa (convert from bar)
        pressure_pa = [p * 1e5 for p in pressures]
        olga_tab += format_grid_values(pressure_pa)
        
        # Temperature grid in °C
        olga_tab += format_grid_values(temperatures)
        
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
            property_arrays[prop['name']] = np.zeros((np_grid, nt_grid))
        
        # Populate property arrays
        for result in results:
            # Get indices
            p_idx = result.get('p_idx', result['index'] // nt_grid)
            t_idx = result.get('t_idx', result['index'] % nt_grid)
            
            # Ensure indices are within bounds
            if 0 <= p_idx < np_grid and 0 <= t_idx < nt_grid:
                # Extract phase information
                phase = get_phase_from_result(result)
                is_two_phase = phase == 'two-phase'
                
                # Process each property
                for prop in olga_properties:
                    value = extract_property_value(result, prop['key'], phase, composition, molar_mass)
                    if value is not None:
                        try:
                            converted_value = prop['converter'](value)
                            property_arrays[prop['name']][p_idx, t_idx] = converted_value
                        except Exception as e:
                            print(f"Error converting {prop['key']}: {e}")
                            # Use a sentinel value to indicate error or undefined
                            if prop['name'] == 'GAS MASS FRACTION OF GAS +LIQUID':
                                # Special case for quality - use value that indicates outside phase envelope
                                property_arrays[prop['name']][p_idx, t_idx] = 998.0 if phase == 'vapor' else -998.0 if phase == 'liquid' else 0.0
        
        # Add properties to OLGA TAB
        for prop in olga_properties:
            olga_tab += f" {prop['name']}                \n"
            olga_tab += format_grid_values(property_arrays[prop['name']].flatten())
        
        return Response(olga_tab, mimetype='text/plain')
    except Exception as e:
        import traceback
        print(f"Error in format_olga_tab: {str(e)}")
        print(traceback.format_exc())
        return Response(f"Error generating OLGA TAB: {str(e)}", status=500, mimetype='text/plain')

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

def get_value_from_field(field):
    """Extract value from field, handling both direct values and dictionaries with 'value' key."""
    if isinstance(field, dict) and 'value' in field:
        return field['value']
    return field

def get_phase_from_result(result):
    """Determine the phase state from the result dictionary."""
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

def convert_vapor_fraction_to_mass(q, result=None, composition=None, molar_mass=None):
    """
    Convert vapor fraction from mole basis to mass basis.
    
    In a simple case, this would need molar masses of each component,
    but for simplicity we'll return the mole fraction directly.
    A proper implementation would use liquid and vapor compositions
    along with component molar masses to do the conversion.
    """
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

def format_grid_values(values, values_per_line=5):
    """
    Format grid values for OLGA TAB format.
    
    Args:
        values (list): List of values to format
        values_per_line (int): Number of values per line
        
    Returns:
        str: Formatted grid values
    """
    formatted = ""
    values_array = np.asarray(values).flatten()
    
    for i in range(0, len(values_array), values_per_line):
        line_values = values_array[i:i + values_per_line]
        formatted += "    " + "    ".join([f".{value:.6E}" for value in line_values]) + "\n"
    
    return formatted