"""
Enhanced OLGA TAB formatter for REFPROP API.
"""

from flask import Response
import numpy as np
import sys
import traceback
import logging
from typing import Dict, List, Any, Optional, Tuple, Union

# Import the configuration
from API.utils.olga_config import (
    OLGA_PROPERTY_MAPPINGS, 
    OLGA_REQUIRED_PROPERTIES, 
    PHASE_MAPPING,
    SPECIAL_VALUES,
    DEFAULT_OLGA_OPTIONS
)

def format_olga_tab(
    x_vars: Dict, 
    y_vars: Dict, 
    results: List[Dict], 
    composition: List[Dict], 
    molar_mass: float, 
    endpoint_type: str = 'pt_flash',
    requested_properties: Optional[List[str]] = None,
    options: Optional[Dict] = None
) -> Response:
    """
    Format calculation results in OLGA TAB format using generic grid variables.
    
    Args:
        x_vars (dict): Dictionary with 'range', 'resolution', and optional 'values' for x-axis variable
        y_vars (dict): Dictionary with 'range', 'resolution', and optional 'values' for y-axis variable
        results (list): List of result dictionaries from the calculation
        composition (list): List of fluid compositions (dict with 'fluid' and 'fraction')
        molar_mass (float): Molar mass of the mixture in g/mol
        endpoint_type (str): Type of endpoint ('pt_flash', 'ph_flash', 'ts_flash', 'vt_flash', 'uv_flash')
        requested_properties (list, optional): List of properties specifically requested by the user
        options (dict, optional): Dictionary of formatting options
        
    Returns:
        Response: Flask response object with OLGA TAB formatted content
    """
    # Initialize logging
    logger = logging.getLogger('olga_formatter')
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    
    # Merge options with defaults
    if options is None:
        options = {}
    for key, value in DEFAULT_OLGA_OPTIONS.items():
        if key not in options:
            options[key] = value
            
    # Set logger level based on debug_level option
    if options['debug_level'] == 0:
        logger.setLevel(logging.ERROR)
    elif options['debug_level'] == 1:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)

    try:
        # Determine grid variable names and units based on endpoint type
        grid_config = get_grid_config(endpoint_type)
        x_name = grid_config['x_name']
        y_name = grid_config['y_name']
        x_header = grid_config['x_header']
        y_header = grid_config['y_header']
        x_idx_name = grid_config['x_idx_name']
        y_idx_name = grid_config['y_idx_name']
        x_multiplier = grid_config['x_multiplier']
        y_multiplier = grid_config['y_multiplier']

        # Use provided grid values if available, otherwise generate from range
        x_grid, nx_grid = get_grid_values(x_vars, logger)
        y_grid, ny_grid = get_grid_values(y_vars, logger)
        
        logger.info(f"Generating OLGA TAB with {nx_grid} {x_name} points and {ny_grid} {y_name} points")
        
        # Compose fluid description
        fluid_desc = " ".join([f"{comp['fluid']}-{comp['fraction']:.4f}" for comp in composition])
        
        # Create OLGA TAB header
        olga_tab = f"'Span-Wagner EOS {fluid_desc}'\n"
        # Use exactly 4 spaces before the dimensions, and include the constant
        olga_tab += f"    {nx_grid}  {ny_grid}    .802294E-08\n"
        
        # X-axis grid with appropriate multiplier
        x_grid_formatted = [x * x_multiplier for x in x_grid]
        olga_tab += format_grid_values(x_grid_formatted, options['values_per_line'], options['indent_spaces'])
        
        # Y-axis grid with appropriate multiplier
        y_grid_formatted = [y * y_multiplier for y in y_grid]
        olga_tab += format_grid_values(y_grid_formatted, options['values_per_line'], options['indent_spaces'])
        
        # Add phase envelope (bubble and dew curves) as special blocks
        # First the bubble point block (filled with special not_applicable values)
        bubble_values = np.ones(ny_grid) * 1.0E+10
        olga_tab += format_grid_values(bubble_values, options['values_per_line'], options['indent_spaces'])
        
        # Then the dew point block (filled with 0.0 values)
        dew_values = np.zeros(ny_grid)
        olga_tab += format_grid_values(dew_values, options['values_per_line'], options['indent_spaces'])
        
        # Prioritize and filter property mappings
        olga_properties = prioritize_properties(OLGA_PROPERTY_MAPPINGS, requested_properties)
        
        # Extract and format properties
        property_arrays = {}
        for prop in olga_properties:
            # Initialize with zeros for all properties
            property_arrays[prop['name']] = np.zeros((nx_grid, ny_grid))
        
        # Create a 2D grid for phase detection
        phase_grid = np.zeros((nx_grid, ny_grid))
        
        # Map results to property arrays
        mapped_points = map_results_to_grid(
            results, property_arrays, phase_grid, 
            x_grid, y_grid, nx_grid, ny_grid,
            x_idx_name, y_idx_name, x_name, y_name,
            olga_properties, molar_mass, composition,
            options, logger
        )
        
        logger.info(f"Successfully mapped {mapped_points} points to the grid")
        
        # Apply any post-processing or fallback calculations
        if options['use_fallbacks']:
            apply_fallback_calculations(
                property_arrays, phase_grid, 
                x_grid_formatted, y_grid_formatted,
                nx_grid, ny_grid, logger
            )
        
        # Add property blocks to OLGA TAB
        for prop in olga_properties:
            olga_tab += f" {prop['name']}                \n"
            olga_tab += format_grid_values(
                property_arrays[prop['name']].flatten(), 
                options['values_per_line'], 
                options['indent_spaces']
            )
        
        # Log statistics if debug enabled
        if options['debug_level'] >= 1:
            log_property_statistics(property_arrays, olga_properties, nx_grid, ny_grid, logger)
        
        # Return as a Response object with proper MIME type
        filename = f"{endpoint_type.replace('_', '-')}.tab"
        return Response(olga_tab, mimetype='text/plain', headers={
            'Content-Disposition': f'attachment; filename={filename}'
        })
        
    except Exception as e:
        logger.error(f"Error in format_olga_tab: {str(e)}")
        traceback.print_exc(file=sys.stderr)
        
        # Return an error Response to prevent socket hang-up
        error_msg = f"Error generating OLGA TAB format: {str(e)}\n\n{traceback.format_exc()}"
        return Response(error_msg, status=500, mimetype='text/plain')

def get_grid_config(endpoint_type: str) -> Dict[str, Any]:
    """
    Get the grid configuration for the specified endpoint type.
    """
    configs = {
        'pt_flash': {
            'x_name': 'pressure',
            'y_name': 'temperature',
            'x_header': 'Pressure (Pa)',
            'y_header': 'Temperature (C)',
            'x_idx_name': 'p_idx',
            'y_idx_name': 't_idx',
            'x_multiplier': 1e5,  # bar to Pa
            'y_multiplier': 1.0   # already in C
        },
        'ph_flash': {
            'x_name': 'pressure',
            'y_name': 'enthalpy',
            'x_header': 'Pressure (Pa)',
            'y_header': 'Enthalpy (J/mol)',
            'x_idx_name': 'p_idx',
            'y_idx_name': 'h_idx',
            'x_multiplier': 1e5,  # bar to Pa
            'y_multiplier': 1.0   # J/mol
        },
        'ts_flash': {
            'x_name': 'temperature',
            'y_name': 'entropy',
            'x_header': 'Temperature (C)',
            'y_header': 'Entropy (J/mol-K)',
            'x_idx_name': 't_idx',
            'y_idx_name': 's_idx',
            'x_multiplier': 1.0,  # already in C
            'y_multiplier': 1.0   # J/mol-K
        },
        'vt_flash': {
            'x_name': 'temperature',
            'y_name': 'specific_volume',
            'x_header': 'Temperature (C)',
            'y_header': 'Specific Volume (m3/mol)',
            'x_idx_name': 't_idx',
            'y_idx_name': 'v_idx',
            'x_multiplier': 1.0,  # already in C
            'y_multiplier': 1.0   # m3/mol
        },
        'uv_flash': {
            'x_name': 'internal_energy',
            'y_name': 'specific_volume',
            'x_header': 'Internal Energy (J/mol)',
            'y_header': 'Specific Volume (m3/mol)',
            'x_idx_name': 'u_idx',
            'y_idx_name': 'v_idx',
            'x_multiplier': 1.0,  # J/mol
            'y_multiplier': 1.0   # m3/mol
        }
    }
    
    # Default to PT flash if endpoint type is not recognized
    return configs.get(endpoint_type, configs['pt_flash'])

def get_grid_values(vars_dict: Dict, logger: logging.Logger) -> Tuple[np.ndarray, int]:
    """
    Extract grid values from the variables dictionary.
    """
    # Use provided grid values if available
    if 'values' in vars_dict:
        grid = vars_dict['values']
        return grid, len(grid)
    
    # Extract range and resolution with robust error handling
    try:
        # Get range with defaults
        range_dict = vars_dict.get('range', {})
        if not range_dict:
            range_dict = {'from': 1.0, 'to': 100.0}
            
        range_from = float(range_dict.get('from', 1.0))
        range_to = float(range_dict.get('to', 100.0))
        
        # Ensure valid range
        if range_to <= range_from:
            range_to = range_from + 10.0
            logger.warning(f"Invalid range (to <= from), adjusted to [{range_from}, {range_to}]")
            
        # Get resolution with default
        resolution = float(vars_dict.get('resolution', 10.0))
        if resolution <= 0:
            resolution = 10.0
            logger.warning(f"Invalid resolution (<=0), set to default {resolution}")
            
        # Generate regular grid
        num_points = int(round((range_to - range_from) / resolution)) + 1
        grid = np.linspace(range_from, range_to, num_points)
        
        return grid, len(grid)
    except (TypeError, ValueError) as e:
        # If there's any issue parsing the ranges, use defaults
        logger.error(f"Error parsing grid ranges: {e}. Using defaults.")
        grid = np.linspace(1.0, 100.0, 10)
        return grid, len(grid)

def prioritize_properties(
    olga_properties: List[Dict], 
    requested_properties: Optional[List[str]] = None
) -> List[Dict]:
    """
    Prioritize and filter OLGA properties based on user requests.
    """
    if not requested_properties:
        # If no specific properties requested, use all defined OLGA properties
        return olga_properties
    
    # Create a prioritized list based on user requests
    requested_olga_props = []
    standard_olga_props = []
    
    for prop in olga_properties:
        # Check if this property relates to any user-requested property
        is_requested = False
        for req_prop in requested_properties:
            if req_prop in prop['key'] or any(req_prop in fb for fb in prop['fallbacks']):
                is_requested = True
                break
                
        if is_requested:
            requested_olga_props.append(prop)
        else:
            standard_olga_props.append(prop)
    
    # Return the prioritized list
    return requested_olga_props + standard_olga_props

def map_results_to_grid(
    results: List[Dict],
    property_arrays: Dict[str, np.ndarray],
    phase_grid: np.ndarray,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    nx_grid: int,
    ny_grid: int,
    x_idx_name: str,
    y_idx_name: str,
    x_name: str,
    y_name: str,
    olga_properties: List[Dict],
    molar_mass: float,
    composition: List[Dict],
    options: Dict,
    logger: logging.Logger
) -> int:
    """
    Map calculation results to the 2D grid for OLGA TAB format.
    """
    mapped_points = 0
    
    for result in results:
        try:
            # Get grid indices from the result
            x_idx, y_idx = get_grid_indices(
                result, x_idx_name, y_idx_name, 
                x_name, y_name, x_grid, y_grid, 
                nx_grid, ny_grid
            )
            
            # Skip if indices are out of bounds
            if x_idx is None or y_idx is None or x_idx < 0 or x_idx >= nx_grid or y_idx < 0 or y_idx >= ny_grid:
                continue
                
            mapped_points += 1
            
            # Extract phase information
            phase = get_phase_from_result(result)
            
            # Store phase information for grid point
            if phase == 'liquid':
                phase_grid[x_idx, y_idx] = 0.0
            elif phase == 'vapor':
                phase_grid[x_idx, y_idx] = 1.0
                
                    # Debug log all vapor properties for verification
                if options['debug_level'] >= 2:
                    for prop in olga_properties:
                        if prop['name'] == 'GAS DENSITY (KG/M3)':
                            value = extract_property_value(result, prop['key'], prop['fallbacks'], phase, composition, molar_mass)
                            converted_value = prop['converter'](value, molar_mass) if value is not None else 0.0
                            logger.info(f"Gas density at point ({x_idx}, {y_idx}): {value} mol/L → {converted_value} kg/m³")

            elif phase == 'two-phase':
                # Use vapor fraction if available
                vf = get_value_from_field(result.get('vapor_fraction', {'value': 0.5}))
                phase_grid[x_idx, y_idx] = vf if 0.0 <= vf <= 1.0 else 0.5
            else:
                # For other phases (supercritical, etc.), use standard mapping
                phase_grid[x_idx, y_idx] = PHASE_MAPPING.get(phase, 0.5)
            
            # Process each property
            for prop in olga_properties:
                try:
                    # Check if property is applicable for this phase
                    if prop['condition'](phase):
                        # Extract the property value
                        value = extract_property_value(
                            result, prop['key'], prop['fallbacks'], 
                            phase, composition, molar_mass
                        )
                        
                        if value is not None:
                            # Apply conversion function
                            converted_value = prop['converter'](value, molar_mass)
                            property_arrays[prop['name']][x_idx, y_idx] = converted_value
                        elif options['debug_level'] >= 2:
                            logger.info(f"No value found for {prop['key']} at point ({x_idx}, {y_idx})")
                except Exception as prop_error:
                    if options['debug_level'] >= 1:
                        logger.warning(f"Error processing property {prop['key']} for point ({x_idx}, {y_idx}): {prop_error}")
        except Exception as e:
            if options['debug_level'] >= 1:
                logger.warning(f"Error processing result: {e}")
            continue
            
    return mapped_points

def get_grid_indices(
    result: Dict,
    x_idx_name: str,
    y_idx_name: str,
    x_name: str,
    y_name: str,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    nx_grid: int,
    ny_grid: int
) -> Tuple[Optional[int], Optional[int]]:
    """
    Get grid indices from result dictionary.
    """
    # First, try to get indices directly from the result
    x_idx = result.get(x_idx_name)
    y_idx = result.get(y_idx_name)
    
    # If specific indices are found, convert them to integers and return
    if x_idx is not None and y_idx is not None:
        try:
            return int(x_idx), int(y_idx)
        except (TypeError, ValueError):
            pass
    
    # If not found, try to derive from 'index' assuming row-major order
    if 'index' in result:
        try:
            idx = int(result['index'])
            x_idx = idx // ny_grid
            y_idx = idx % ny_grid
            
            # Make sure indices are within bounds
            if 0 <= x_idx < nx_grid and 0 <= y_idx < ny_grid:
                return x_idx, y_idx
        except (TypeError, ValueError):
            pass
    
    # Final attempt: map by coordinate values
    if x_name in result and y_name in result:
        try:
            x_val = get_value_from_field(result[x_name])
            y_val = get_value_from_field(result[y_name])
            
            if x_val is not None and y_val is not None:
                x_idx = find_nearest_index(x_grid, x_val)
                y_idx = find_nearest_index(y_grid, y_val)
                return x_idx, y_idx
        except Exception:
            pass
    
    # Could not determine position
    return None, None

def extract_property_value(
    result: Dict,
    key: str,
    fallbacks: List[str],
    phase: str,
    composition: List[Dict],
    molar_mass: float
) -> Optional[float]:
    """
    Extract property value from result, with fallback mechanisms.
    """
    try:
        # Direct property match
        if key in result:
            return get_value_from_field(result[key])
        
        # Try fallback keys
        for fallback in fallbacks:
            if fallback in result:
                value = get_value_from_field(result[fallback])
                if value is not None:
                    # For phase-specific properties, check if context is appropriate
                    if key.startswith('liquid_') and phase != 'liquid' and phase != 'two-phase':
                        continue
                    if key.startswith('vapor_') and phase != 'vapor' and phase != 'two-phase':
                        continue
                    return value
        
        # Special handling for specific property types
        
        # Phase-specific properties in two-phase region
        if phase == 'two-phase':
            q = get_value_from_field(result.get('vapor_fraction', {'value': 0.5}))
            
            # For liquid properties, try to derive from base and vapor properties
            if key.startswith('liquid_'):
                base_key = key[7:]  # Remove 'liquid_' prefix
                if base_key in result and 'vapor_' + base_key in result:
                    base_val = get_value_from_field(result[base_key])
                    vapor_val = get_value_from_field(result['vapor_' + base_key])
                    if base_val is not None and vapor_val is not None and q != 1.0:
                        # Calculate liquid value based on mixing equation
                        # base = q * vapor + (1-q) * liquid
                        # => liquid = (base - q * vapor) / (1-q)
                        try:
                            if q < 1.0:  # Avoid division by zero
                                return (base_val - q * vapor_val) / (1-q)
                        except Exception:
                            pass
            
            # For vapor properties, try to derive from base and liquid properties
            if key.startswith('vapor_'):
                base_key = key[6:]  # Remove 'vapor_' prefix
                if base_key in result and 'liquid_' + base_key in result:
                    base_val = get_value_from_field(result[base_key])
                    liquid_val = get_value_from_field(result['liquid_' + base_key])
                    if base_val is not None and liquid_val is not None and q != 0.0:
                        # Calculate vapor value based on mixing equation
                        # base = q * vapor + (1-q) * liquid
                        # => vapor = (base - (1-q) * liquid) / q
                        try:
                            if q > 0.0:  # Avoid division by zero
                                return (base_val - (1-q) * liquid_val) / q
                        except Exception:
                            pass
        
        # Water component properties
        if key == 'water_density' or key.startswith('water_'):
            # Find water component by name
            water_idx = next((i for i, comp in enumerate(composition) 
                             if "WATER" in comp['fluid'].upper()), -1)
            
            if water_idx >= 0:
                if phase == 'liquid' or phase == 'two-phase':
                    # Get water fraction in liquid
                    if 'x' in result:
                        x = get_value_from_field(result['x'])
                        if isinstance(x, list) and len(x) > water_idx:
                            water_fraction = x[water_idx]
                            # Get liquid density
                            liquid_density = get_value_from_field(result.get('liquid_density', 
                                                                 result.get('density')))
                            if liquid_density is not None and water_fraction > 0:
                                # Water density in liquid phase
                                return liquid_density * water_fraction
                
                if phase == 'vapor' or phase == 'two-phase':
                    # Get water fraction in vapor
                    if 'y' in result:
                        y = get_value_from_field(result['y'])
                        if isinstance(y, list) and len(y) > water_idx:
                            water_fraction = y[water_idx]
                            # Get vapor density
                            vapor_density = get_value_from_field(result.get('vapor_density', 
                                                                result.get('density')))
                            if vapor_density is not None and water_fraction > 0:
                                # Water density in vapor phase
                                return vapor_density * water_fraction
            
            # Default to zero for water properties if water not found
            return 0.0
        
        # Density derivatives for specific phases
        if key == 'dDdP_liquid' and 'dDdP' in result and phase == 'liquid':
            return get_value_from_field(result['dDdP'])
            
        if key == 'dDdP_vapor' and 'dDdP' in result and phase == 'vapor':
            return get_value_from_field(result['dDdP'])
            
        if key == 'dDdT_liquid' and 'dDdT' in result and phase == 'liquid':
            return get_value_from_field(result['dDdT'])
            
        if key == 'dDdT_vapor' and 'dDdT' in result and phase == 'vapor':
            return get_value_from_field(result['dDdT'])
        
        # Not found after all attempts
        return None
    except Exception as e:
        print(f"Error extracting property {key}: {e}")
        return None

def apply_fallback_calculations(
    property_arrays: Dict[str, np.ndarray],
    phase_grid: np.ndarray,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    nx_grid: int,
    ny_grid: int,
    logger: logging.Logger
) -> None:
    """
    Apply fallback calculations for missing or invalid property values.
    """
    # Identify zero values that need fallback calculations
    for i in range(nx_grid):
        for j in range(ny_grid):
            # Current phase for this grid point
            phase_value = phase_grid[i, j]
            
            # For density derivatives, use finite differences if direct values not available
            if 'DRHOG/DP (S2/M2)' in property_arrays and property_arrays['DRHOG/DP (S2/M2)'][i, j] == 0.0:
                if phase_value > 0.0:  # Has vapor component
                    # Use finite differences with neighboring points if available
                    if i > 0 and i < nx_grid - 1 and 'GAS DENSITY (KG/M3)' in property_arrays:
                        try:
                            dp = (x_grid[i+1] - x_grid[i-1]) / 2
                            if dp > 0:
                                drho_g = property_arrays['GAS DENSITY (KG/M3)'][i+1, j] - property_arrays['GAS DENSITY (KG/M3)'][i-1, j]
                                property_arrays['DRHOG/DP (S2/M2)'][i, j] = drho_g / dp
                        except Exception as e:
                            logger.info(f"Finite difference fallback failed for DRHOG/DP at ({i},{j}): {e}")
            
            if 'DRHOL/DP (S2/M2)' in property_arrays and property_arrays['DRHOL/DP (S2/M2)'][i, j] == 0.0:
                if phase_value < 1.0:  # Has liquid component
                    # Similar approach for liquid density derivatives
                    if i > 0 and i < nx_grid - 1 and 'LIQUID DENSITY (KG/M3)' in property_arrays:
                        try:
                            dp = (x_grid[i+1] - x_grid[i-1]) / 2
                            if dp > 0:
                                drho_l = property_arrays['LIQUID DENSITY (KG/M3)'][i+1, j] - property_arrays['LIQUID DENSITY (KG/M3)'][i-1, j]
                                property_arrays['DRHOL/DP (S2/M2)'][i, j] = drho_l / dp
                        except Exception as e:
                            logger.info(f"Finite difference fallback failed for DRHOL/DP at ({i},{j}): {e}")
                            
            # Similar approach for temperature derivatives
            if 'DRHOG/DT (KG/M3/K)' in property_arrays and property_arrays['DRHOG/DT (KG/M3/K)'][i, j] == 0.0:
                if phase_value > 0.0:  # Has vapor component
                    if j > 0 and j < ny_grid - 1 and 'GAS DENSITY (KG/M3)' in property_arrays:
                        try:
                            dt = (y_grid[j+1] - y_grid[j-1]) / 2
                            if dt > 0:
                                drho_g = property_arrays['GAS DENSITY (KG/M3)'][i, j+1] - property_arrays['GAS DENSITY (KG/M3)'][i, j-1]
                                property_arrays['DRHOG/DT (KG/M3/K)'][i, j] = drho_g / dt
                        except Exception as e:
                            logger.info(f"Finite difference fallback failed for DRHOG/DT at ({i},{j}): {e}")

def log_property_statistics(
    property_arrays: Dict[str, np.ndarray],
    olga_properties: List[Dict],
    nx_grid: int,
    ny_grid: int,
    logger: logging.Logger
) -> None:
    """
    Log statistics about property arrays.
    """
    for prop in olga_properties:
        if prop['name'] in property_arrays:
            array = property_arrays[prop['name']]
            zero_count = np.count_nonzero(array == 0)
            total_count = nx_grid * ny_grid
            if zero_count == total_count:
                logger.warning(f"Property {prop['name']} is all zeros")
            elif zero_count > 0:
                logger.info(f"Property {prop['name']} has {zero_count}/{total_count} zero values")

def get_value_from_field(field: Union[Dict, Any]) -> Any:
    """
    Extract value from field, handling both direct values and dictionaries with 'value' key.
    """
    try:
        if isinstance(field, dict) and 'value' in field:
            return field['value']
        return field
    except Exception as e:
        print(f"Error getting value from field: {e}")
        return None

def get_phase_from_result(result: Dict) -> str:
    """
    Determine the phase state from the result dictionary.
    """
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
                elif 'super' in phase_lower and 'critic' in phase_lower:
                    return 'supercritical'
                elif 'solid' in phase_lower:
                    return 'solid'
        
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
            elif q > 998:  # Special codes used in some implementations
                return 'vapor'
            elif q < -998:
                return 'liquid'
            elif q == 999:
                return 'supercritical'
            elif q == -999:
                return 'solid'
        
        # Default to unknown
        return 'unknown'
    except Exception as e:
        print(f"Error determining phase: {e}")
        return 'unknown'

def find_nearest_index(array: np.ndarray, value: float) -> int:
    """
    Find the index of the value closest to the target in the array.
    """
    try:
        array = np.asarray(array)
        idx = (np.abs(array - float(value))).argmin()
        return idx
    except Exception as e:
        print(f"Error in find_nearest_index: {e}")
        return 0  # Return 0 as a fallback

# In olga_formatter.py, find the format_grid_values function and update it:

def format_grid_values(values: np.ndarray, values_per_line: int = 5, indent_spaces: int = 5) -> str:
    """
    Format grid values for OLGA TAB format with the correct notation.
    """
    try:
        formatted = ""
        values_array = np.asarray(values).flatten()
        
        for i in range(0, len(values_array), values_per_line):
            line_values = values_array[i:i + values_per_line]
            line_strings = []
            
            for value in line_values:
                if value == 0:
                    line_strings.append(".000000E+00")
                    continue
                
                # Check if value is negative
                is_negative = value < 0
                abs_value = abs(value)
                
                # Get standard scientific notation for absolute value
                sci_notation = f"{abs_value:.6E}"
                
                # Extract mantissa and exponent parts
                mantissa_str, exponent_str = sci_notation.split('E')
                exponent = int(exponent_str)
                
                # CRITICAL FIX: Increment exponent by 1 to compensate for decimal point shift
                # This ensures: 123.456 → .123456E+03 (not .123456E+02)
                mantissa_without_point = mantissa_str.replace('.', '')
                exponent += 1
                exponent_str = f"{exponent:+03d}"
                
                # Add negative sign before decimal if needed
                if is_negative:
                    olga_formatted = f"-.{mantissa_without_point}E{exponent_str}"
                else:
                    olga_formatted = f".{mantissa_without_point}E{exponent_str}"
                
                line_strings.append(olga_formatted)
            
            formatted += " " * indent_spaces + "    ".join(line_strings) + "\n"
        
        return formatted
    except Exception as e:
        print(f"Error formatting grid values: {e}")
        return " " * indent_spaces + ".000000E+00\n"
    
