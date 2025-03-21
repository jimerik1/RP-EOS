"""
JSON response formatter for REFPROP API.

This module provides functions to format calculation results as JSON responses.
"""

from flask import jsonify
from typing import Dict, List, Any, Optional
import logging

logger = logging.getLogger(__name__)

def format_json_response(results: List[Dict[str, Any]], 
                        grid_info: Optional[Dict[str, Any]] = None) -> Any:
    """
    Format calculation results as a JSON response.
    
    Args:
        results: List of calculated property points
        grid_info: Optional dictionary with grid information
            
    Returns:
        Flask response with JSON data
    """
    response = {'results': results}
    
    if grid_info:
        response['grid_info'] = grid_info
        
    return jsonify(response)


def filter_properties(results: List[Dict[str, Any]], 
                     requested_properties: List[str]) -> List[Dict[str, Any]]:
    """
    Filter results to include only requested properties.
    
    Args:
        results: List of calculated property points
        requested_properties: List of property names to include
            
    Returns:
        Filtered results
    """
    # Always include basic attributes even if not explicitly requested
    core_attributes = ["index", "p_idx", "t_idx", "temperature", "pressure", "phase"]
    
    # Build the final property list
    props_to_keep = list(set(requested_properties).union(set(core_attributes)))
    
    # Filter each result
    filtered_results = []
    for result in results:
        filtered_result = {
            prop: value for prop, value in result.items()
            if prop in props_to_keep
        }
        filtered_results.append(filtered_result)
        
    return filtered_results