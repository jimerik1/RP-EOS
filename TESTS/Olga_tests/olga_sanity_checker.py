"""
Test script for the enhanced OLGA TAB formatter.

This script makes an API call to the PT-Flash endpoint with OLGA TAB format 
and verifies the output contains valid property data.
"""

import requests
import argparse
import os
import re
import numpy as np
from typing import Dict, List, Any

def test_pt_flash_olga_tab(base_url: str, composition: List[Dict], output_file: str = None):
    """
    Test the PT-Flash endpoint with OLGA TAB format.
    
    Args:
        base_url: Base URL for the API
        composition: List of fluid compositions
        output_file: Optional output file to save the TAB file
    """
    print(f"Testing PT-Flash with OLGA TAB format at: {base_url}")
    
    endpoint = f"{base_url}/pt_flash"
    
    # Create test payload
    payload = {
        "composition": composition,
        "variables": {
            "pressure": {
                "range": {"from": 10, "to": 80},
                "resolution": 10
            },
            "temperature": {
                "range": {"from": -20, "to": 30},
                "resolution": 10
            }
        },
        "calculation": {
            "properties": [
                "temperature", "pressure", "density", "liquid_density", "vapor_density",
                "enthalpy", "entropy", "cp", "viscosity", "thermal_conductivity", 
                "surface_tension", "vapor_fraction", "phase", "dDdP", "dDdT"
            ],
            "units_system": "SI",
            "response_format": "olga_tab",
            "grid_type": "adaptive",
            "enhancement_factor": 5.0
        }
    }
    
    try:
        # Make the API call
        response = requests.post(endpoint, json=payload)
        
        if response.status_code != 200:
            print(f"Error: Response status code {response.status_code}")
            try:
                error_data = response.json()
                print(f"Error message: {error_data.get('error', 'Unknown error')}")
            except:
                print(f"Error content: {response.content}")
            return False
            
        # Get the content as text
        olga_tab_content = response.text
        
        # Save to file if requested
        if output_file:
            with open(output_file, 'w') as f:
                f.write(olga_tab_content)
            print(f"OLGA TAB content saved to: {output_file}")
            
        # Analyze the content
        analyze_olga_tab(olga_tab_content)
        
        return True
    
    except Exception as e:
        print(f"Error during test: {str(e)}")
        return False

def analyze_olga_tab(content: str) -> None:
    """
    Analyze the OLGA TAB content.
    
    Args:
        content: The OLGA TAB content as a string
    """
    lines = content.split('\n')
    
    # Check header
    if not lines[0].startswith("'Span-Wagner EOS"):
        print("Warning: Header doesn't match expected format")
    
    # Extract grid dimensions
    grid_dimensions = re.match(r'\s+(\d+)\s+(\d+)\s+', lines[1])
    if grid_dimensions:
        nx = int(grid_dimensions.group(1))
        ny = int(grid_dimensions.group(2))
        print(f"Grid dimensions: {nx}x{ny}")
    else:
        print("Warning: Could not extract grid dimensions")
        return
    
    # Count property blocks
    property_count = 0
    zero_blocks = []
    
    current_property = None
    property_lines = []
    
    for i, line in enumerate(lines[2:]):
        # Check if line starts a new property block
        if line.strip().endswith('                '):
            # If we were collecting a property, analyze it
            if current_property:
                # Check if all values are zero
                values = extract_values_from_lines(property_lines)
                if np.all(np.isclose(values, 0.0)):
                    zero_blocks.append(current_property)
                    
                # Reset for next property
                property_lines = []
                
            # Extract property name
            current_property = line.strip()
            property_count += 1
        elif current_property and line.strip():
            # Collect lines for current property
            property_lines.append(line)
    
    # Check if any blocks are all zeros
    print(f"Found {property_count} property blocks")
    if zero_blocks:
        print(f"Warning: {len(zero_blocks)} properties have all zero values:")
        for prop in zero_blocks:
            print(f"  - {prop}")
    else:
        print("All property blocks contain non-zero values.")

def extract_values_from_lines(lines: List[str]) -> np.ndarray:
    """
    Extract floating-point values from OLGA TAB formatted lines.
    
    Args:
        lines: List of OLGA TAB lines containing formatted values
        
    Returns:
        Array of extracted values
    """
    values = []
    
    for line in lines:
        # Extract values using regular expression
        # Pattern matches scientific notation with leading decimal point: .123456E+01
        matches = re.findall(r'[+-]?\.\d+E[+-]\d+', line)
        
        for match in matches:
            # Process the number:
            # 1. Find sign (- or implied +)
            sign = -1.0 if match.startswith('-') else 1.0
            
            # 2. Extract mantissa and exponent
            mantissa_exp = match[1:] if match.startswith('-') else match
            mantissa_str, exponent_str = mantissa_exp.split('E')
            
            # 3. Add decimal point to mantissa
            mantissa = float('0.' + mantissa_str[1:])
            
            # 4. Parse exponent
            exponent = int(exponent_str)
            
            # 5. Calculate final value
            value = sign * mantissa * (10 ** exponent)
            values.append(value)
    
    return np.array(values)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Test OLGA TAB formatter')
    parser.add_argument('--url', type=str, default='http://localhost:5051',
                       help='Base URL for the API (default: http://localhost:5051)')
    parser.add_argument('--output', type=str, default='test_output.tab',
                       help='Output file for TAB content (default: test_output.tab)')
    
    args = parser.parse_args()
    
    # Test composition
    composition = [
        {"fluid": "CO2", "fraction": 0.9},
        {"fluid": "METHANE", "fraction": 0.1}
    ]
    
    test_pt_flash_olga_tab(args.url, composition, args.output)