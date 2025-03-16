#!/usr/bin/env python3
"""
Example script to generate an OLGA TAB formatted file from the Span-Wagner EOS API.

This script demonstrates how to:
1. Make a request to the PT-Flash endpoint
2. Request output in OLGA TAB format
3. Save the result to a file

OLGA TAB files are commonly used in multiphase flow simulators like OLGA
to represent fluid property tables across pressure and temperature ranges.
"""

import requests
import os
import json
import sys

# API endpoint URL (modify as needed)
API_URL = "http://localhost:5051/pt_flash"

def generate_olga_tab(output_file="output.tab"):
    """
    Generate an OLGA TAB file for water using the Span-Wagner EOS API.
    
    Args:
        output_file (str): Path to the output TAB file
    """
    # Define input parameters for API request
    request_data = {
        "composition": [
            {"fluid": "WATER", "fraction": 1.0}
        ],
        "variables": {
            "pressure": {
                "range": {"from": 1.0, "to": 200.0},
                "resolution": 10.0
            },
            "temperature": {
                "range": {"from": 0.0, "to": 200.0},
                "resolution": 10.0
            }
        },
        "calculation": {
            "properties": [
                "density", "liquid_density", "vapor_density",
                "enthalpy", "entropy", "cp", "cv",
                "viscosity", "thermal_conductivity", "surface_tension",
                "dDdP", "dDdT", "vapor_fraction"
            ],
            "units_system": "SI",
            "response_format": "olga_tab"
        }
    }

    # Make API request
    print(f"Sending request to {API_URL}...")
    response = requests.post(API_URL, json=request_data)
    
    # Check response status
    if response.status_code != 200:
        print(f"Error: Received status code {response.status_code}")
        try:
            error_data = response.json()
            print(f"Error message: {error_data.get('error', 'Unknown error')}")
        except:
            print(f"Raw response: {response.text}")
        sys.exit(1)
    
    # Save response to file
    with open(output_file, "w") as f:
        f.write(response.text)
    
    print(f"Successfully generated OLGA TAB file: {output_file}")
    print(f"File size: {os.path.getsize(output_file)} bytes")

def main():
    # Parse command-line arguments
    import argparse
    parser = argparse.ArgumentParser(description="Generate OLGA TAB file from Span-Wagner EOS API")
    parser.add_argument("--output", "-o", default="output.tab", help="Output file path")
    args = parser.parse_args()
    
    # Generate the TAB file
    generate_olga_tab(args.output)

if __name__ == "__main__":
    main()