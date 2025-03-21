"""
TS-Flash OLGA TAB endpoint for the REFPROP API.

This endpoint calculates fluid properties at given temperature and entropy
conditions and returns the results in OLGA TAB format.
"""

from flask import request, jsonify
import sys
import traceback
import logging
from typing import Dict, List, Any, Optional

# Import the endpoint blueprint
from API.endpoints import ts_flash_olga_bp

# Import from the new core modules
from API.core.property_system import PropertyRegistry
from API.core.flash_calculators import TSFlashCalculator
from API.core.formatters.olga_formatter import format_olga_response

# Import REFPROP instance
from API.refprop_setup import RP

# Import utility functions
from API.utils.helpers import validate_composition
from API.utils.olga_config import OLGA_REQUIRED_PROPERTIES

# Configure logging
logger = logging.getLogger(__name__)

@ts_flash_olga_bp.route('/ts_flash_olga', methods=['POST'])
def ts_flash_olga():
    """
    Calculate properties using TS flash and return in OLGA TAB format.
    
    Request Format:
    {
        "composition": [
            {"fluid": "FLUID_NAME1", "fraction": X1},
            {"fluid": "FLUID_NAME2", "fraction": X2}
        ],
        "variables": {
            "temperature": {
                "range": {"from": T_MIN, "to": T_MAX},
                "resolution": T_STEP
            },
            "entropy": {
                "range": {"from": S_MIN, "to": S_MAX},
                "resolution": S_STEP
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
        calculator = TSFlashCalculator(RP, registry)
        
        # Calculate the flash
        logger.info(f"Starting TS flash calculation for OLGA TAB format with {len(all_properties)} properties")
        results, grid_info, grids = calculator.calculate_flash_grid(
            composition, variables, all_properties, **grid_options
        )
        
        logger.info(f"TS flash calculation complete: {len(results)} points calculated")
        
        # Format as OLGA TAB
        return format_olga_response(
            results, grids, composition, 
            endpoint_type='ts_flash',
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
    if 'temperature' not in variables or 'entropy' not in variables:
        raise ValueError('Missing temperature or entropy variables')
    
    # Check ranges
    temperature_range = variables['temperature'].get('range', {})
    entropy_range = variables['entropy'].get('range', {})
    
    if not all([
        'from' in temperature_range, 'to' in temperature_range,
        'from' in entropy_range, 'to' in entropy_range
    ]):
        raise ValueError('Missing range parameters (from/to)')
    
    # Check resolutions
    if not all([
        'resolution' in variables['temperature'],
        'resolution' in variables['entropy']
    ]):
        raise ValueError('Missing resolution parameters')