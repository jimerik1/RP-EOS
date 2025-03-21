"""
UV-Flash OLGA TAB endpoint for the REFPROP API.

This endpoint calculates fluid properties at given internal energy and volume
conditions and returns the results in OLGA TAB format.
"""

from flask import request, jsonify
import sys
import traceback
import logging
from typing import Dict, List, Any, Optional

# Import the endpoint blueprint
from API.endpoints import uv_flash_olga_bp

# Import from the new core modules
from API.core.property_system import PropertyRegistry
from API.core.flash_calculators import UVFlashCalculator
from API.core.formatters.olga_formatter import format_olga_response

# Import REFPROP instance
from API.refprop_setup import RP

# Import utility functions
from API.utils.helpers import validate_composition
from API.utils.olga_config import OLGA_REQUIRED_PROPERTIES

# Configure logging
logger = logging.getLogger(__name__)

@uv_flash_olga_bp.route('/uv_flash_olga', methods=['POST'])
def uv_flash_olga():
    """
    Calculate properties using UV flash and return in OLGA TAB format.
    
    Request Format:
    {
        "composition": [
            {"fluid": "FLUID_NAME1", "fraction": X1},
            {"fluid": "FLUID_NAME2", "fraction": X2}
        ],
        "variables": {
            "internal_energy": {
                "range": {"from": U_MIN, "to": U_MAX},
                "resolution": U_STEP
            },
            "specific_volume": {
                "range": {"from": V_MIN, "to": V_MAX},
                "resolution": V_STEP
            }
        },
        "calculation": {
            "properties": ["property1", "property2", ...],
            "grid_type": "equidistant" | "adaptive" | "logarithmic" | "exponential",
            "enhancement_factor": 5.0,
            "boundary_zone_width": null
        }
    }
    
    Returns:
        OLGA TAB formatted text response
    """
    try:
        # Parse the request
        data = request.get_json(force=True)
        
        # Validate the request structure
        _validate_request(data)
        
        # Extract calculation options
        composition = data['composition']
        variables = data['variables']
        calculation = data.get('calculation', {})
        
        # Get requested properties, ensuring we include all OLGA required properties
        requested_properties = calculation.get('properties', [])
        all_properties = list(set(requested_properties).union(set(OLGA_REQUIRED_PROPERTIES)))
        
        # Get grid options
        grid_options = {
            'grid_type': calculation.get('grid_type', 'equidistant'),
            'enhancement_factor': calculation.get('enhancement_factor', 5.0),
            'boundary_zone_width': calculation.get('boundary_zone_width')
        }
        
        # Initialize the property registry and calculator
        registry = PropertyRegistry()
        calculator = UVFlashCalculator(RP, registry)
        
        # Calculate the flash
        logger.info(f"Starting UV flash calculation for OLGA TAB format with {len(all_properties)} properties")
        results, grid_info, grids = calculator.calculate_flash_grid(
            composition, variables, all_properties, **grid_options
        )
        
        logger.info(f"UV flash calculation complete: {len(results)} points calculated")
        
        # Format as OLGA TAB
        return format_olga_response(
            results, grids, composition, 
            endpoint_type='uv_flash',
            requested_properties=requested_properties
        )
        
    except ValueError as ve:
        # Handle validation errors
        logger.error(f"Validation error: {str(ve)}")
        return jsonify({'error': str(ve)}), 400
        
    except Exception as e:
        # Handle unexpected errors
        logger.error("Error processing request:", exc_info=True)
        return jsonify({'error': str(e)}), 500


def _validate_request(data: Dict[str, Any]) -> None:
    """
    Validate the request data.
    
    Args:
        data: Request data dictionary
        
    Raises:
        ValueError: If validation fails
    """
    # Check for required fields
    required_fields = ['composition', 'variables']
    for field in required_fields:
        if field not in data:
            raise ValueError(f'Missing field: {field}')
    
    # Check composition
    if not validate_composition(data['composition']):
        raise ValueError('Invalid composition - fractions must sum to 1')
    
    # Check variables
    variables = data.get('variables', {})
    if 'internal_energy' not in variables or 'specific_volume' not in variables:
        raise ValueError('Missing internal_energy or specific_volume variables')
    
    # Check ranges
    energy_range = variables['internal_energy'].get('range', {})
    volume_range = variables['specific_volume'].get('range', {})
    
    if not all([
        'from' in energy_range, 'to' in energy_range,
        'from' in volume_range, 'to' in volume_range
    ]):
        raise ValueError('Missing range parameters (from/to)')
    
    # Check resolutions
    if not all([
        'resolution' in variables['internal_energy'],
        'resolution' in variables['specific_volume']
    ]):
        raise ValueError('Missing resolution parameters')