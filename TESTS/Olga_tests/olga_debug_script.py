#!/usr/bin/env python3
"""
Simplified OLGA TAB debug script to examine the raw format.
This will print the raw OLGA TAB content to help diagnose parsing issues.
"""

import requests
import argparse

def query_pt_flash_olga(
    base_url: str, 
    pressure_range: tuple,
    temperature: float,
    num_points: int = 25,
    fluid: str = "METHANE"
):
    """Query the PT flash endpoint with OLGA TAB response."""
    print(f"Querying PT flash OLGA at {temperature}°C, {pressure_range[0]}-{pressure_range[1]} bar")
    
    # Create request payload
    payload = {
        "composition": [
            {"fluid": fluid, "fraction": 1.0}
        ],
        "variables": {
            "pressure": {
                "range": {"from": pressure_range[0], "to": pressure_range[1]},
                "resolution": (pressure_range[1] - pressure_range[0]) / (num_points - 1)
            },
            "temperature": {
                "range": {"from": temperature, "to": temperature},
                "resolution": 1
            }
        },
        "calculation": {
            "properties": [
                "temperature", "pressure", "density", "vapor_density", 
                "liquid_density", "phase", "vapor_fraction", "enthalpy",
                "entropy", "cp", "viscosity", "thermal_conductivity"
            ],
            "units_system": "SI",
            "response_format": "olga_tab"
        }
    }
    
    try:
        # Make the API call
        response = requests.post(f"{base_url}/pt_flash", json=payload)
        
        if response.status_code != 200:
            print(f"Error: Response status code {response.status_code}")
            print(f"Response content: {response.content}")
            return ""
            
        # Return the OLGA TAB content
        return response.text
        
    except Exception as e:
        print(f"Error querying OLGA endpoint: {str(e)}")
        return ""

def examine_olga_tab(content: str):
    """Examine the OLGA TAB content to help diagnose parsing issues."""
    if not content:
        print("No OLGA TAB content received")
        return
        
    lines = content.split('\n')
    
    print("\n===== OLGA TAB CONTENT =====")
    print(f"Total number of lines: {len(lines)}")
    
    # Print the first few lines
    print("\nFirst 15 lines:")
    for i in range(min(15, len(lines))):
        print(f"Line {i+1}: {lines[i]}")
        
    # Print some lines from the middle
    mid = len(lines) // 2
    print(f"\nMiddle lines (around line {mid}):")
    for i in range(max(0, mid-5), min(len(lines), mid+5)):
        print(f"Line {i+1}: {lines[i]}")
    
    # Look for property blocks
    print("\nSearching for property blocks...")
    property_blocks = []
    for i, line in enumerate(lines):
        # Look for lines that might be property headers
        line_stripped = line.strip()
        if line_stripped and line.startswith(' ') and not line_stripped.startswith('.'):
            # This might be a property header
            print(f"Potential property block at line {i+1}: '{line_stripped}'")
            property_blocks.append((i+1, line_stripped))
    
    print(f"\nFound {len(property_blocks)} potential property blocks")
    
    # Check for grid dimensions line
    grid_line = None
    for i, line in enumerate(lines):
        if i > 0 and i < 5:  # Grid dimensions usually in the first few lines
            if len(line.strip().split()) == 3:  # Should have 3 numbers
                grid_line = line
                print(f"\nPotential grid dimensions at line {i+1}: '{grid_line}'")
                break
    
    # Check for OLGA scientific notation
    sci_notation_count = 0
    for line in lines[:100]:  # Check first 100 lines
        if '.' in line and 'E' in line:
            sci_notation_count += 1
    
    print(f"\nLines with potential scientific notation in first 100 lines: {sci_notation_count}")
    
    # Save the content to a file for manual inspection
    with open("olga_raw_output.tab", "w") as f:
        f.write(content)
    print("\nRaw OLGA output saved to 'olga_raw_output.tab'")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Examine raw OLGA TAB format')
    parser.add_argument('--url', type=str, default='http://localhost:5051',
                       help='Base URL for the API (default: http://localhost:5051)')
    parser.add_argument('--pmin', type=float, default=1.0,
                       help='Minimum pressure in bar (default: 1.0)')
    parser.add_argument('--pmax', type=float, default=1500.0,
                       help='Maximum pressure in bar (default: 1500.0)')
    parser.add_argument('--temp', type=float, default=50.0,
                       help='Temperature in °C (default: 50.0)')
    parser.add_argument('--points', type=int, default=25,
                       help='Number of pressure points (default: 25)')
    parser.add_argument('--fluid', type=str, default='METHANE',
                       help='Fluid name (default: METHANE)')
    
    args = parser.parse_args()
    
    # Query OLGA endpoint
    olga_content = query_pt_flash_olga(
        args.url, 
        (args.pmin, args.pmax), 
        args.temp,
        args.points,
        args.fluid
    )
    
    # Examine the content
    examine_olga_tab(olga_content)
    
    print("\nExamination complete!")

if __name__ == "__main__":
    main()