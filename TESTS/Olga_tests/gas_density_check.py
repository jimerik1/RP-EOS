"""
Validation script to compare gas density values between JSON and OLGA TAB formats.
This script queries the PT flash endpoint twice - once with JSON response and once with OLGA TAB,
then plots and compares the gas density values.
"""

import requests
import numpy as np
import matplotlib.pyplot as plt
import re
import argparse
import io
import pandas as pd
from typing import Dict, List, Tuple, Any

def query_pt_flash_json(
    base_url: str, 
    pressure_range: Tuple[float, float], 
    temperature: float,
    num_points: int = 20
) -> List[Dict[str, Any]]:
    """
    Query the PT flash endpoint for methane with JSON response.
    
    Args:
        base_url: Base URL for the API
        pressure_range: (min_pressure, max_pressure) in bar
        temperature: Fixed temperature in °C
        num_points: Number of pressure points
        
    Returns:
        List of result points
    """
    print(f"Querying PT flash JSON at {temperature}°C, {pressure_range[0]}-{pressure_range[1]} bar")
    
    # Create pressure points
    pressures = np.linspace(pressure_range[0], pressure_range[1], num_points)
    
    # Create request payload
    payload = {
        "composition": [
            {"fluid": "METHANE", "fraction": 1.0}
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
                "liquid_density", "phase", "vapor_fraction"
            ],
            "units_system": "SI",
            "response_format": "json"
        }
    }
    
    try:
        # Make the API call
        response = requests.post(f"{base_url}/pt_flash", json=payload)
        
        if response.status_code != 200:
            print(f"Error: Response status code {response.status_code}")
            print(f"Response content: {response.content}")
            return []
            
        # Parse the JSON response
        data = response.json()
        return data.get("results", [])
        
    except Exception as e:
        print(f"Error querying JSON endpoint: {str(e)}")
        return []

def query_pt_flash_olga(
    base_url: str, 
    pressure_range: Tuple[float, float], 
    temperature: float,
    num_points: int = 20
) -> str:
    """
    Query the PT flash endpoint for methane with OLGA TAB response.
    
    Args:
        base_url: Base URL for the API
        pressure_range: (min_pressure, max_pressure) in bar
        temperature: Fixed temperature in °C
        num_points: Number of pressure points
        
    Returns:
        OLGA TAB content as string
    """
    print(f"Querying PT flash OLGA at {temperature}°C, {pressure_range[0]}-{pressure_range[1]} bar")
    
    # Create pressure points with slight adjustment to match JSON points
    # OLGA grid generation might differ slightly, so we use the exact same grid
    pressures = np.linspace(pressure_range[0], pressure_range[1], num_points)
    
    # Create request payload
    payload = {
        "composition": [
            {"fluid": "METHANE", "fraction": 1.0}
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
                "liquid_density", "phase", "vapor_fraction"
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

def parse_olga_tab(content: str, property_name: str = "GAS DENSITY (KG/M3)") -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse OLGA TAB file to extract pressure and gas density values.
    
    Args:
        content: OLGA TAB content as string
        property_name: Name of the property to extract
        
    Returns:
        (pressures, densities) arrays
    """
    lines = content.split('\n')
    
    # Parse header to get grid dimensions
    grid_line = re.search(r'\s+(\d+)\s+(\d+)', lines[1])
    if not grid_line:
        print("Error: Could not parse OLGA TAB grid dimensions")
        return np.array([]), np.array([])
    
    n_pressure = int(grid_line.group(1))
    n_temp = int(grid_line.group(2))
    
    # Parse pressure values from the file (x-axis)
    pressure_values = []
    line_index = 2
    while len(pressure_values) < n_pressure:
        line = lines[line_index].strip()
        if line:
            # Extract values using regex for OLGA scientific notation
            values = re.findall(r'[+-]?\.\d+E[+-]\d+', line)
            for value_str in values:
                # Convert OLGA scientific notation to float
                if value_str.startswith('-'):
                    sign = -1.0
                    value_str = value_str[1:]
                else:
                    sign = 1.0
                    
                mantissa_str, exponent_str = value_str.split('E')
                # In OLGA format, the decimal point is at the start, so we need to adjust
                mantissa = float("0" + mantissa_str)
                exponent = int(exponent_str)
                
                # Convert to float
                value = sign * mantissa * (10 ** exponent)
                
                # Store the pressure in bar (convert from Pa)
                pressure_values.append(value / 1e5)  # Convert Pa to bar
                
                if len(pressure_values) >= n_pressure:
                    break
        line_index += 1
    
    # Find the property block
    property_start = None
    for i, line in enumerate(lines):
        if property_name in line:
            property_start = i + 1
            break
    
    if property_start is None:
        print(f"Error: Property '{property_name}' not found in OLGA TAB file")
        return np.array([]), np.array([])
    
    # Parse the property values - we want the first row (temperature = first value)
    property_values = []
    line_index = property_start
    
    while len(property_values) < n_pressure and line_index < len(lines):
        line = lines[line_index].strip()
        if line:
            # Extract values using regex for OLGA scientific notation
            values = re.findall(r'[+-]?\.\d+E[+-]\d+', line)
            for value_str in values:
                # Convert OLGA scientific notation to float
                if value_str.startswith('-'):
                    sign = -1.0
                    value_str = value_str[1:]
                else:
                    sign = 1.0
                    
                mantissa_str, exponent_str = value_str.split('E')
                # In OLGA format, the decimal point is at the start, so we need to adjust
                mantissa = float("0" + mantissa_str)
                exponent = int(exponent_str)
                
                # Convert to float
                value = sign * mantissa * (10 ** exponent)
                
                property_values.append(value)
                
                if len(property_values) >= n_pressure:
                    break
        line_index += 1
    
    # We only want the first row (first temperature)
    if len(property_values) > n_pressure:
        property_values = property_values[:n_pressure]
    
    return np.array(pressure_values), np.array(property_values)

def extract_json_values(json_results: List[Dict[str, Any]]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract pressure and gas density values from JSON results.
    
    Args:
        json_results: List of result points from JSON endpoint
        
    Returns:
        (pressures, densities) arrays
    """
    pressures = []
    densities = []
    
    for result in json_results:
        # Extract pressure and gas density
        pressure = result.get('pressure', {}).get('value')
        
        # Try vapor_density first, then density (for supercritical)
        if 'vapor_density' in result:
            density = result.get('vapor_density', {}).get('value')
        else:
            density = result.get('density', {}).get('value')
            
            # Also check the phase to ensure we're getting gas density
            phase = result.get('phase', {}).get('value', '')
            if 'liquid' in phase.lower() and 'vapor' not in phase.lower():
                # Skip liquid points
                continue
        
        if pressure is not None and density is not None:
            # Convert mol/L to kg/m³ using methane molecular weight (16.04 g/mol)
            density_kg_m3 = density * 16.04
            
            pressures.append(pressure)
            densities.append(density_kg_m3)
    
    # Sort by pressure
    sorted_indices = np.argsort(pressures)
    sorted_pressures = np.array(pressures)[sorted_indices]
    sorted_densities = np.array(densities)[sorted_indices]
    
    return sorted_pressures, sorted_densities

def plot_comparison(
    json_pressures: np.ndarray, 
    json_densities: np.ndarray,
    olga_pressures: np.ndarray, 
    olga_densities: np.ndarray,
    temperature: float,
    output_file: str = None
):
    """
    Plot comparison of gas density values from JSON and OLGA.
    
    Args:
        json_pressures: Pressure values from JSON response
        json_densities: Gas density values from JSON response
        olga_pressures: Pressure values from OLGA TAB response
        olga_densities: Gas density values from OLGA TAB response
        temperature: Temperature in °C
        output_file: Optional output file path
    """
    plt.figure(figsize=(12, 8))
    
    # Plot JSON data
    plt.plot(json_pressures, json_densities, 'bo-', label='JSON Response', markersize=6)
    
    # Plot OLGA data
    plt.plot(olga_pressures, olga_densities, 'rx--', label='OLGA TAB Response', markersize=6)
    
    plt.xlabel('Pressure (bar)', fontsize=12)
    plt.ylabel('Gas Density (kg/m³)', fontsize=12)
    plt.title(f'Methane Gas Density Comparison at {temperature}°C', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=12)
    
    # Calculate density difference statistics
    if len(json_pressures) > 0 and len(olga_pressures) > 0:
        # Find common pressure values (within tolerance)
        common_points = []
        for jp, jd in zip(json_pressures, json_densities):
            for op, od in zip(olga_pressures, olga_densities):
                if abs(jp - op) < 1e-5:
                    common_points.append((jp, jd, od))
                    break
        
        if common_points:
            # Calculate differences
            diffs = [abs(jd - od) for _, jd, od in common_points]
            max_diff = max(diffs)
            avg_diff = sum(diffs) / len(diffs)
            
            # Add statistics to plot
            plt.figtext(0.15, 0.15, 
                      f'Maximum Difference: {max_diff:.6f} kg/m³\n'
                      f'Average Difference: {avg_diff:.6f} kg/m³',
                      bbox=dict(facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")
        
    plt.show()

def save_data_to_csv(
    json_pressures: np.ndarray, 
    json_densities: np.ndarray,
    olga_pressures: np.ndarray, 
    olga_densities: np.ndarray,
    temperature: float,
    output_file: str
):
    """
    Save comparison data to CSV file.
    
    Args:
        json_pressures: Pressure values from JSON response
        json_densities: Gas density values from JSON response
        olga_pressures: Pressure values from OLGA TAB response
        olga_densities: Gas density values from OLGA TAB response
        temperature: Temperature in °C
        output_file: Output CSV file path
    """
    # Create a combined DataFrame
    # Since the pressure values might not match exactly, we'll create separate entries
    data = {
        'Source': [],
        'Pressure (bar)': [],
        'Gas Density (kg/m³)': [],
        'Temperature (°C)': []
    }
    
    # Add JSON data
    for p, d in zip(json_pressures, json_densities):
        data['Source'].append('JSON')
        data['Pressure (bar)'].append(p)
        data['Gas Density (kg/m³)'].append(d)
        data['Temperature (°C)'].append(temperature)
    
    # Add OLGA data
    for p, d in zip(olga_pressures, olga_densities):
        data['Source'].append('OLGA')
        data['Pressure (bar)'].append(p)
        data['Gas Density (kg/m³)'].append(d)
        data['Temperature (°C)'].append(temperature)
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    print(f"Data saved to {output_file}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Compare gas density values from JSON and OLGA TAB formats')
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
    parser.add_argument('--plot', type=str, default='density_comparison.png',
                       help='Output plot file (default: density_comparison.png)')
    parser.add_argument('--csv', type=str, default='density_comparison.csv',
                       help='Output CSV file (default: density_comparison.csv)')
    
    args = parser.parse_args()
    
    # Query both endpoints
    json_results = query_pt_flash_json(
        args.url, 
        (args.pmin, args.pmax), 
        args.temp,
        args.points
    )
    
    olga_content = query_pt_flash_olga(
        args.url, 
        (args.pmin, args.pmax), 
        args.temp,
        args.points
    )
    
    if not json_results:
        print("Error: No JSON results received")
        return
        
    if not olga_content:
        print("Error: No OLGA TAB content received")
        return
    
    # Extract values from both responses
    json_pressures, json_densities = extract_json_values(json_results)
    olga_pressures, olga_densities = parse_olga_tab(olga_content)
    
    if len(json_pressures) == 0 or len(olga_pressures) == 0:
        print("Error: Failed to extract valid data from responses")
        return
    
    # Save the comparison data to CSV
    save_data_to_csv(
        json_pressures, 
        json_densities,
        olga_pressures, 
        olga_densities,
        args.temp,
        args.csv
    )
    
    # Plot the comparison
    plot_comparison(
        json_pressures, 
        json_densities,
        olga_pressures, 
        olga_densities,
        args.temp,
        args.plot
    )
    
    print("Validation complete!")

if __name__ == "__main__":
    main()