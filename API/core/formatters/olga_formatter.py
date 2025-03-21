"""
OLGA TAB response formatter wrapper for REFPROP API.

This module provides a wrapper around the existing OLGA TAB formatter
to adapt it to the new architecture.
"""

from flask import Response
from typing import Dict, List, Any, Optional
import numpy as np
import logging

logger = logging.getLogger(__name__)

def format_olga_response(results: List[Dict[str, Any]], 
                        grids: Dict[str, np.ndarray],
                        composition: List[Dict[str, Any]], 
                        endpoint_type: str = 'pt_flash',
                        requested_properties: Optional[List[str]] = None,
                        options: Optional[Dict[str, Any]] = None) -> Response:
    """
    Format calculation results as an OLGA TAB response.
    
    Args:
        results: List of calculated property points
        grids: Dictionary of grid arrays used for calculation
        composition: List of fluid components and fractions
        endpoint_type: Type of flash calculation ('pt_flash', 'ph_flash', etc.)
        requested_properties: List of property names explicitly requested
        options: Additional formatting options
            
    Returns:
        Flask response with OLGA TAB data
    """
    from API.utils.olga_formatter import format_olga_tab
    from API.utils.olga_config import OLGA_REQUIRED_PROPERTIES, DEFAULT_OLGA_OPTIONS
    from API.refprop_setup import RP
    
    # Merge options with defaults
    if options is None:
        options = {}
    merged_options = DEFAULT_OLGA_OPTIONS.copy()
    merged_options.update(options)
    
    # Get molecular weight
    z = [comp['fraction'] for comp in composition] + [0] * (20 - len(composition))
    molar_mass = RP.WMOLdll(z)
    
    # Configure default file name based on endpoint type
    filename = f"{endpoint_type.replace('_', '-')}.tab"
    
    # Extract grid variables based on endpoint type
    if endpoint_type == 'pt_flash':
        x_vars = _prepare_grid_vars(grids, "pressure", 100.0)  # Convert bar to kPa
        y_vars = _prepare_grid_vars(grids, "temperature", -273.15, offset=True)  # Convert K to °C
    elif endpoint_type == 'ph_flash':
        x_vars = _prepare_grid_vars(grids, "pressure", 100.0)  # Convert bar to kPa
        y_vars = _prepare_grid_vars(grids, "enthalpy")
    elif endpoint_type == 'ts_flash':
        x_vars = _prepare_grid_vars(grids, "temperature", -273.15, offset=True)  # Convert K to °C
        y_vars = _prepare_grid_vars(grids, "entropy")
    elif endpoint_type == 'vt_flash':
        x_vars = _prepare_grid_vars(grids, "temperature", -273.15, offset=True)  # Convert K to °C
        y_vars = _prepare_grid_vars(grids, "specific_volume")
    elif endpoint_type == 'uv_flash':
        x_vars = _prepare_grid_vars(grids, "internal_energy")
        y_vars = _prepare_grid_vars(grids, "specific_volume")
    else:
        # Default (fallback) for unknown types
        logger.warning(f"Unknown endpoint type: {endpoint_type}. Using first two grid variables.")
        keys = list(grids.keys())
        if len(keys) < 2:
            return Response(f"Error: Not enough grid variables for {endpoint_type}", 
                           status=400, mimetype='text/plain')
        x_vars = _prepare_grid_vars(grids, keys[0])
        y_vars = _prepare_grid_vars(grids, keys[1])
    
    # Call the existing OLGA formatter
    try:
        olga_response = format_olga_tab(
            x_vars, y_vars, results, composition, molar_mass, 
            endpoint_type=endpoint_type,
            requested_properties=requested_properties,
            options=merged_options
        )
        
        # If format_olga_tab returns a Response, use it directly
        if isinstance(olga_response, Response):
            return olga_response
            
        # Otherwise, create a new Response
        return Response(
            olga_response,
            mimetype='text/plain',
            headers={'Content-Disposition': f'attachment; filename={filename}'}
        )
        
    except Exception as e:
        logger.error(f"Error formatting OLGA TAB: {str(e)}")
        error_msg = f"Error generating OLGA TAB format: {str(e)}"
        return Response(error_msg, status=500, mimetype='text/plain')


def _prepare_grid_vars(grids: Dict[str, np.ndarray], var_name: str, 
                      conversion_factor: float = 1.0, offset: bool = False) -> Dict[str, Any]:
    """
    Prepare grid variable information for OLGA formatter.
    
    Args:
        grids: Dictionary of grid arrays
        var_name: Name of the variable
        conversion_factor: Factor to multiply values by (for unit conversion)
        offset: If True, apply conversion as an offset instead of a factor
            
    Returns:
        Dictionary with grid variable information
    """
    if var_name not in grids:
        logger.warning(f"Grid variable {var_name} not found. Using first available grid.")
        var_name = next(iter(grids))
    
    grid = grids[var_name]
    
    # Apply conversion (multiply or add)
    if offset:
        values = grid + conversion_factor
    else:
        values = grid * conversion_factor
    
    return {
        'range': {'from': float(values.min()), 'to': float(values.max())},
        'resolution': _estimate_resolution(values),
        'values': values
    }


def _estimate_resolution(values: np.ndarray) -> float:
    """
    Estimate the grid resolution from values.
    
    Args:
        values: Array of grid values
            
    Returns:
        Estimated resolution
    """
    if len(values) <= 1:
        return 1.0
    
    # Calculate all step sizes
    steps = np.diff(values)
    
    # Use the most common step size
    unique_steps, counts = np.unique(np.round(steps, 10), return_counts=True)
    
    if len(unique_steps) == 0:
        # Fallback if no steps found
        return (values[-1] - values[0]) / (len(values) - 1)
    
    # Return the most common step size
    return float(unique_steps[np.argmax(counts)])