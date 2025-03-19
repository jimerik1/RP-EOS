#!/usr/bin/env python3
"""
Enhanced debugging script for OLGA vs JSON format conversion issues.
This script helps identify exactly where and how the unit conversion mismatch occurs.
"""

import requests
import numpy as np
import matplotlib.pyplot as plt
import re
import sys
from typing import Dict, List, Tuple, Any, Optional

def query_pt_flash_json(
    base_url: str, 
    pressure_range: Tuple[float, float], 
    temperature: float,
    num_points: int = 20,
    verbose: bool = True
) -> Dict[str, Any]:
    """Query the PT flash endpoint with JSON response."""
    if verbose:
        print(f"Querying PT flash JSON at {temperature}°C, {pressure_range[0]}-{pressure_range[1]} bar")
    
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
                "liquid_density", "phase", "vapor_fraction", "critical_density"
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
            return {}
            
        # Parse the JSON response
        return response.json()
        
    except Exception as e:
        print(f"Error querying JSON endpoint: {str(e)}")
        return {}

def query_pt_flash_olga(
    base_url: str, 
    pressure_range: Tuple[float, float], 
    temperature: float,
    num_points: int = 20,
    verbose: bool = True
) -> str:
    """Query the PT flash endpoint with OLGA TAB response."""
    if verbose:
        print(f"Querying PT flash OLGA at {temperature}°C, {pressure_range[0]}-{pressure_range[1]} bar")
    
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
            "response_format": "olga_tab",
            "debug_level": 2  # Enable detailed debug output
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

def query_models_info(base_url: str) -> Dict[str, Any]:
    """Query the models_info endpoint to get fluid information including molecular weight."""
    payload = {
        "composition": [
            {"fluid": "METHANE", "fraction": 1.0}
        ]
    }
    
    try:
        response = requests.post(f"{base_url}/models_info", json=payload)
        if response.status_code != 200:
            print(f"Error: Response status code {response.status_code}")
            return {}
            
        return response.json()
    except Exception as e:
        print(f"Error querying models_info endpoint: {str(e)}")
        return {}

def parse_olga_scientific(value_str: str) -> float:
    """Parse a value in OLGA scientific notation."""
    if not value_str:
        return 0.0
        
    # Check for negative sign
    if value_str.startswith('-'):
        sign = -1.0
        value_str = value_str[1:]
    else:
        sign = 1.0
        
    try:
        mantissa_str, exponent_str = value_str.split('E')
        mantissa = float("0" + mantissa_str)
        exponent = int(exponent_str)
        
        # CRITICAL FIX: Decrement exponent to convert from OLGA format
        exponent -= 1  # This line is essential
        
        return sign * mantissa * (10 ** exponent)
    except Exception as e:
        print(f"Error parsing OLGA value '{value_str}': {e}")
        return 0.0

        
def parse_olga_tab(
    content: str, 
    property_name: str = "GAS DENSITY (KG/M3)",
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse OLGA TAB file to extract pressure and gas density values.
    
    Args:
        content: OLGA TAB content as string
        property_name: Name of the property to extract
        verbose: Whether to print verbose debug info
        
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
    
    if verbose:
        print(f"OLGA grid dimensions: {n_pressure} pressure points, {n_temp} temperature points")
    
    # Parse pressure values from the file (x-axis)
    pressure_values = []
    line_index = 2
    while len(pressure_values) < n_pressure and line_index < len(lines):
        line = lines[line_index].strip()
        if line:
            # Extract values using regex for OLGA scientific notation
            values = re.findall(r'[+-]?\.\d+E[+-]\d+', line)
            
            if verbose and line_index < 5:
                print(f"Line {line_index}: Found {len(values)} values: {values[:3]}...")
                
            for value_str in values:
                if not value_str:
                    continue
                    
                try:
                    value = parse_olga_scientific(value_str)
                    
                    # Store the pressure in bar (convert from Pa)
                    pressure_values.append(value / 1e5)  # Convert Pa to bar
                    
                    if verbose and len(pressure_values) <= 3:
                        print(f"Parsed pressure: {value_str} → {value} Pa → {value/1e5} bar")
                        
                    if len(pressure_values) >= n_pressure:
                        break
                except Exception as e:
                    print(f"Error processing value '{value_str}': {e}")
                    
        line_index += 1
    
    if len(pressure_values) < n_pressure:
        print(f"Warning: Only found {len(pressure_values)} pressure values, expected {n_pressure}")
    
    # Find the property block
    property_start = None
    for i, line in enumerate(lines):
        if property_name in line and line.strip() == " " + property_name:
            property_start = i + 1
            if verbose:
                print(f"Found property block '{property_name}' at line {i+1}")
            break
    
    if property_start is None:
        # Try alternative search
        for i, line in enumerate(lines):
            if property_name in line:
                property_start = i + 1
                if verbose:
                    print(f"Found property block '{property_name}' at line {i+1} using fallback search")
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
            
            if verbose and len(property_values) == 0:
                print(f"Property line {line_index}: Found {len(values)} values: {values[:3]}...")
                
            for value_str in values:
                if not value_str:
                    continue
                
                try:
                    value = parse_olga_scientific(value_str)
                    property_values.append(value)
                    
                    if verbose and len(property_values) <= 3:
                        print(f"Parsed property: {value_str} → {value}")
                        
                    if len(property_values) >= n_pressure:
                        break
                except Exception as e:
                    print(f"Error processing property value '{value_str}': {e}")
                    
        line_index += 1
    
    # We only want the first row (first temperature)
    if len(property_values) > n_pressure:
        property_values = property_values[:n_pressure]
    
    if verbose:
        if len(property_values) > 0:
            print(f"First few property values: {property_values[:5]}")
        else:
            print("No property values found")
    
    if len(property_values) < n_pressure:
        print(f"Warning: Only found {len(property_values)} property values, expected {n_pressure}")
        # Pad with zeros if necessary
        property_values.extend([0.0] * (n_pressure - len(property_values)))
    
    return np.array(pressure_values), np.array(property_values)

def parse_olga_tab(
    content: str, 
    property_name: str = "GAS DENSITY (KG/M3)",
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """Parse OLGA TAB file to extract pressure and property values."""
    lines = content.split('\n')
    
    # Parse header to get grid dimensions
    grid_line = re.search(r'\s+(\d+)\s+(\d+)', lines[1])
    if not grid_line:
        print("Error: Could not parse OLGA TAB grid dimensions")
        return np.array([]), np.array([])
    
    n_pressure = int(grid_line.group(1))
    n_temp = int(grid_line.group(2))
    
    if verbose:
        print(f"OLGA grid dimensions: {n_pressure} pressure points, {n_temp} temperature points")
    
    # Parse pressure values from the file (x-axis)
    pressure_values = []
    line_index = 2
    while len(pressure_values) < n_pressure and line_index < len(lines):
        line = lines[line_index].strip()
        if line:
            # Extract values using regex for OLGA scientific notation
            values = re.findall(r'[+-]?\.\d+E[+-]\d+', line)
            for value_str in values:
                value = parse_olga_scientific(value_str)
                
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
            if verbose:
                print(f"Found property block '{property_name}' at line {i+1}")
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
                value = parse_olga_scientific(value_str)
                property_values.append(value)
                
                if len(property_values) >= n_pressure:
                    break
        line_index += 1
    
    # We only want the first row (first temperature)
    if len(property_values) > n_pressure:
        property_values = property_values[:n_pressure]
    
    if verbose and len(property_values) > 0:
        print(f"First few property values: {property_values[:5]}")
    
    return np.array(pressure_values), np.array(property_values)

def extract_json_values(
    json_data: Dict[str, Any], 
    molecular_weight: float = 16.04, 
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """Extract pressure and gas density values from JSON results."""
    pressures = []
    densities_mol_l = []  # Store original density values
    densities_kg_m3 = []  # Store converted density values
    
    if not json_data or 'results' not in json_data:
        return np.array([]), np.array([])
    
    for result in json_data['results']:
        # Extract pressure
        if 'pressure' in result:
            pressure = result.get('pressure', {}).get('value')
        else:
            continue
            
        # Try vapor_density first, then density (for supercritical)
        phase = result.get('phase', {}).get('value', '')
        
        if 'vapor_density' in result:
            density = result.get('vapor_density', {}).get('value')
            source = 'vapor_density'
        else:
            density = result.get('density', {}).get('value')
            source = 'density'
            
            # Skip if this is explicitly liquid phase
            if 'liquid' in phase.lower() and 'vapor' not in phase.lower() and 'super' not in phase.lower():
                if verbose:
                    print(f"Skipping liquid phase point at {pressure} bar")
                continue
        
        if pressure is not None and density is not None:
            # Store original and converted values
            density_kg_m3 = density * molecular_weight
            
            pressures.append(pressure)
            densities_mol_l.append(density)
            densities_kg_m3.append(density_kg_m3)
            
            if verbose and len(pressures) <= 5:
                print(f"JSON: P={pressure:.2f} bar, phase={phase}, source={source}, "
                      f"density={density:.6f} mol/L, converted={density_kg_m3:.6f} kg/m³")
    
    # Sort by pressure
    if pressures:
        sorted_indices = np.argsort(pressures)
        sorted_pressures = np.array(pressures)[sorted_indices]
        sorted_densities_mol_l = np.array(densities_mol_l)[sorted_indices]
        sorted_densities_kg_m3 = np.array(densities_kg_m3)[sorted_indices]
        return sorted_pressures, sorted_densities_kg_m3, sorted_densities_mol_l
    
    return np.array([]), np.array([]), np.array([])

def compare_specific_points(
    json_pressures: np.ndarray, 
    json_densities: np.ndarray,
    json_densities_mol_l: np.ndarray,
    olga_pressures: np.ndarray, 
    olga_densities: np.ndarray,
    num_points: int = 5
) -> None:
    """Compare specific points between JSON and OLGA outputs."""
    print("\n===== POINT-BY-POINT COMPARISON =====")
    print("Pressure (bar) | JSON (kg/m³) | OLGA (kg/m³) | Difference | Ratio | JSON (mol/L)")
    print("---------------|--------------|--------------|------------|-------|------------")
    
    # Find points with matching pressures
    comparison_points = []
    for i, jp in enumerate(json_pressures):
        for j, op in enumerate(olga_pressures):
            if abs(jp - op) < 1e-3:  # Match within 0.001 bar
                comparison_points.append((jp, json_densities[i], olga_densities[j], json_densities_mol_l[i]))
                break
    
    # Sort by pressure and limit to num_points
    comparison_points.sort(key=lambda x: x[0])
    for p, jd, od, jd_mol_l in comparison_points[:num_points]:
        diff = jd - od
        ratio = jd / od if od != 0 else float('inf')
        print(f"{p:13.2f} | {jd:12.6f} | {od:12.6f} | {diff:10.6f} | {ratio:5.3f} | {jd_mol_l:12.6f}")
    
    # Calculate average ratio if possible
    if comparison_points:
        ratios = [jd/od for _, jd, od, _ in comparison_points if od != 0]
        if ratios:
            avg_ratio = sum(ratios) / len(ratios)
            print(f"\nAverage ratio (JSON/OLGA): {avg_ratio:.4f}")
            print(f"This suggests OLGA values might need to be multiplied by {avg_ratio:.4f}")

def debug_potential_conversion_issue(
    json_pressures: np.ndarray, 
    json_densities: np.ndarray,
    json_densities_mol_l: np.ndarray,
    olga_pressures: np.ndarray, 
    olga_densities: np.ndarray,
    mw: float
) -> None:
    """Debug potential conversion issues by trying different scaling factors."""
    print("\n===== DEBUGGING CONVERSION ISSUES =====")
    
    # Calculate statistics on ratios
    ratios = []
    for i, jp in enumerate(json_pressures):
        for j, op in enumerate(olga_pressures):
            if abs(jp - op) < 1e-3:  # Match within 0.001 bar
                if olga_densities[j] != 0:
                    ratios.append(json_densities[i] / olga_densities[j])
                break
    
    if not ratios:
        print("No matching points found for ratio analysis.")
        return
    
    avg_ratio = sum(ratios) / len(ratios)
    
    print(f"Average ratio of JSON/OLGA values: {avg_ratio:.4f}")
    
    # Test different scenarios
    print("\nTesting possible conversion issues:")
    
    # Scenario 1: OLGA values are already in mol/L despite kg/m³ label
    print("\nScenario 1: OLGA values are actually in mol/L (need conversion):")
    olga_converted = olga_densities * mw
    print_comparison(json_densities, olga_converted, title="JSON (kg/m³) vs OLGA * MW (kg/m³)")
    
    # Scenario 2: OLGA converter function didn't multiply by MW
    print("\nScenario 2: OLGA converter didn't apply MW conversion:")
    print_comparison(json_densities_mol_l, olga_densities/mw, title="JSON (mol/L) vs OLGA / MW (mol/L)")
    
    # Scenario 3: OLGA has incorrect handling of unit prefixes
    print("\nScenario 3: OLGA has 1000x scale error:")
    print_comparison(json_densities, olga_densities * 1000, title="JSON (kg/m³) vs OLGA * 1000 (kg/m³)")
    
    # Scenario 4: OLGA is using wrong molecular weight
    test_mw = 1.0
    while True:
        new_ratio = avg_ratio / test_mw
        if abs(new_ratio - 1.0) < 0.05:  # Within 5% of 1.0
            break
        test_mw += 0.1
        if test_mw > 50:  # Safety limit
            test_mw = None
            break
    
    if test_mw is not None:
        print(f"\nScenario 4: OLGA might be using MW of {test_mw:.1f} instead of {mw}:")
        olga_fixed = olga_densities * (mw/test_mw)
        print_comparison(json_densities, olga_fixed, title=f"JSON (kg/m³) vs OLGA adjusted with MW={test_mw:.1f} (kg/m³)")
    
    # Scenario 5: OLGA converter is being applied twice
    print("\nScenario 5: OLGA conversion is being applied twice:")
    olga_fixed = olga_densities * avg_ratio
    print_comparison(json_densities, olga_fixed, title="JSON (kg/m³) vs OLGA * correction factor (kg/m³)")

def print_comparison(
    values1: np.ndarray, 
    values2: np.ndarray, 
    title: str = "Comparison",
    num_points: int = 5
) -> None:
    """Print side-by-side comparison of two value arrays."""
    print(f"\n{title}")
    print("Value 1    | Value 2    | Ratio    | Difference")
    print("-----------|------------|----------|------------")
    
    # Take samples from different regions (start, middle, end)
    if len(values1) > 0 and len(values2) > 0:
        min_len = min(len(values1), len(values2))
        
        # Generate indices to sample from different regions
        if min_len <= num_points:
            indices = list(range(min_len))
        else:
            step = min_len // (num_points-1)
            indices = list(range(0, min_len, step))
            if indices[-1] < min_len - 1:
                indices.append(min_len - 1)
        
        for i in indices:
            if i < min_len:
                v1 = values1[i]
                v2 = values2[i]
                ratio = v1 / v2 if v2 != 0 else float('inf')
                diff = v1 - v2
                print(f"{v1:9.6f} | {v2:10.6f} | {ratio:8.4f} | {diff:10.6f}")
    else:
        print("Insufficient data for comparison")

def debug_molar_mass(base_url: str) -> float:
    """Debug molar mass calculation by querying the models_info endpoint."""
    print("\n===== DEBUGGING MOLAR MASS =====")
    
    # Get molecular weight from models_info
    models_info = query_models_info(base_url)
    api_mw = None
    
    if models_info and 'models' in models_info and 'components' in models_info['models']:
        for comp in models_info['models']['components']:
            if 'molar_mass' in comp:
                api_mw = comp['molar_mass']
                print(f"API reports molar mass: {api_mw} g/mol")
                break
    
    # Standard value for methane
    standard_mw = 16.04  # g/mol
    print(f"Standard methane molar mass: {standard_mw} g/mol")
    
    if api_mw and abs(api_mw - standard_mw) > 0.01:
        print(f"WARNING: API molar mass ({api_mw}) differs from standard value ({standard_mw})")
    
    return api_mw if api_mw else standard_mw

def main() -> None:
    # Parse command-line arguments
    if len(sys.argv) > 1:
        base_url = sys.argv[1]
    else:
        base_url = "http://localhost:5051"
    
    print(f"Using API URL: {base_url}")
    
    # Debug molar mass first
    mw = debug_molar_mass(base_url)
    
    # Query both endpoints
    pressure_range = (1.0, 1500.0)
    temperature = 50.0
    num_points = 25
    
    json_data = query_pt_flash_json(base_url, pressure_range, temperature, num_points)
    olga_content = query_pt_flash_olga(base_url, pressure_range, temperature, num_points)
    
    if not json_data or 'results' not in json_data:
        print("Error: No JSON results received")
        return
        
    if not olga_content:
        print("Error: No OLGA TAB content received")
        return
    
    # Print summary of the API response
    print(f"\n===== JSON RESPONSE SUMMARY =====")
    print(f"Number of result points: {len(json_data.get('results', []))}")
    
    # Extract values from both responses
    json_pressures, json_densities, json_densities_mol_l = extract_json_values(
        json_data, molecular_weight=mw, verbose=True
    )
    
    olga_pressures, olga_densities = parse_olga_tab(olga_content, verbose=True)
    
    if len(json_pressures) == 0 or len(olga_pressures) == 0:
        print("Error: Failed to extract valid data from responses")
        return
    
    # Compare specific points
    compare_specific_points(
        json_pressures, json_densities, json_densities_mol_l,
        olga_pressures, olga_densities
    )
    
    # Debug potential conversion issues
    debug_potential_conversion_issue(
        json_pressures, json_densities, json_densities_mol_l,
        olga_pressures, olga_densities, mw
    )
    
    # Plot comparison with correction
    plt.figure(figsize=(12, 8))
    
    # Original values
    plt.plot(json_pressures, json_densities, 'bo-', label='JSON Response', markersize=6)
    plt.plot(olga_pressures, olga_densities, 'rx--', label='OLGA TAB Response', markersize=6)
    
    # Corrected values based on average ratio
    ratios = []
    for i, jp in enumerate(json_pressures):
        for j, op in enumerate(olga_pressures):
            if abs(jp - op) < 1e-3 and olga_densities[j] != 0:
                ratios.append(json_densities[i] / olga_densities[j])
    
    if ratios:
        avg_ratio = sum(ratios) / len(ratios)
        plt.plot(olga_pressures, olga_densities * avg_ratio, 'g^--', 
                 label=f'OLGA TAB (corrected ×{avg_ratio:.4f})', markersize=6)
    
    plt.xlabel('Pressure (bar)', fontsize=12)
    plt.ylabel('Gas Density (kg/m³)', fontsize=12)
    plt.title(f'Methane Gas Density Comparison at {temperature}°C', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('density_comparison_debug.png')
    print(f"\nPlot saved to density_comparison_debug.png")
    
    print("\nDebugging complete!")

if __name__ == "__main__":
    main()