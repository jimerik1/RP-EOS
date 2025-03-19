#!/usr/bin/env python3
"""
OLGA vs JSON Format Test Script - Revised Version

This script compares properties between JSON and OLGA TAB format responses from the PT flash endpoint.
It handles unit conversions properly and ensures accurate comparison of values.
"""

import requests
import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import json
from typing import Dict, List, Tuple, Any, Optional
from dataclasses import dataclass

# Configuration
DEFAULT_API_URL = "http://localhost:5051"
TEMP_CELSIUS = 0.0  # Temperature for tests (fixed to 0°C for consistency)
PRESSURE_RANGE = (1.0, 150.0)  # Pressure range for tests (bar)
NUM_POINTS = 15  # Number of pressure points
DEFAULT_MOLECULAR_WEIGHT = 16.04  # g/mol for methane

@dataclass
class PropertyConfig:
    """Configuration for property testing."""
    json_name: str
    olga_name: str
    unit: str
    json_unit: Optional[str] = None  # Unit in JSON response
    olga_unit: Optional[str] = None  # Unit in OLGA response
    converter: Optional[callable] = None  # Function to convert JSON value to OLGA unit
    description: Optional[str] = None

# Define properties to test with their conversion functions
PROPERTIES_TO_TEST = [
    PropertyConfig(
        json_name="vapor_density", 
        olga_name="GAS DENSITY (KG/M3)", 
        unit="kg/m³",
        json_unit="mol/L",
        converter=lambda x, mw: x * mw,  # mol/L → kg/m³
        description="Gas density"
    ),
    PropertyConfig(
        json_name="liquid_density", 
        olga_name="LIQUID DENSITY (KG/M3)", 
        unit="kg/m³",
        json_unit="mol/L",
        converter=lambda x, mw: x * mw,  # mol/L → kg/m³
        description="Liquid density"
    ),
    PropertyConfig(
        json_name="vapor_viscosity", 
        olga_name="GAS VISCOSITY (NS/M2)", 
        unit="mPa·s",
        json_unit="μPa·s",
        converter=lambda x, mw: x * 1e-3,  # μPa·s → mPa·s (NS/M2)
        description="Gas viscosity"
    ),
    PropertyConfig(
        json_name="viscosity", 
        olga_name="GAS VISCOSITY (NS/M2)", 
        unit="mPa·s",
        json_unit="μPa·s",
        converter=lambda x, mw: x * 1e-3,  # μPa·s → mPa·s (NS/M2)
        description="Viscosity (fallback for gas)"
    ),
    PropertyConfig(
        json_name="liquid_viscosity", 
        olga_name="LIQUID VISCOSITY (NS/M2)", 
        unit="mPa·s",
        json_unit="μPa·s",
        converter=lambda x, mw: x * 1e-3,  # μPa·s → mPa·s (NS/M2)
        description="Liquid viscosity"
    ),
    PropertyConfig(
        json_name="vapor_cp", 
        olga_name="GAS HEAT CAPACITY (J/KG/K)", 
        unit="J/(kg·K)",
        json_unit="J/(mol·K)",
        converter=lambda x, mw: x * 1000.0 / mw,  # J/(mol·K) → J/(kg·K)
        description="Gas heat capacity"
    ),
    PropertyConfig(
        json_name="cp", 
        olga_name="GAS HEAT CAPACITY (J/KG/K)", 
        unit="J/(kg·K)",
        json_unit="J/(mol·K)",
        converter=lambda x, mw: x * 1000.0 / mw,  # J/(mol·K) → J/(kg·K)
        description="Heat capacity (fallback for gas)"
    ),
    PropertyConfig(
        json_name="liquid_cp", 
        olga_name="LIQUID HEAT CAPACITY (J/KG/K)", 
        unit="J/(kg·K)",
        json_unit="J/(mol·K)",
        converter=lambda x, mw: x * 1000.0 / mw,  # J/(mol·K) → J/(kg·K)
        description="Liquid heat capacity"
    ),
    PropertyConfig(
        json_name="vapor_enthalpy", 
        olga_name="GAS ENTHALPY (J/KG)", 
        unit="J/kg",
        json_unit="J/mol",
        converter=lambda x, mw: x * 1000.0 / mw,  # J/mol → J/kg
        description="Gas enthalpy"
    ),
    PropertyConfig(
        json_name="enthalpy", 
        olga_name="GAS ENTHALPY (J/KG)", 
        unit="J/kg",
        json_unit="J/mol",
        converter=lambda x, mw: x * 1000.0 / mw,  # J/mol → J/kg
        description="Enthalpy (fallback for gas)"
    ),
    PropertyConfig(
        json_name="liquid_enthalpy", 
        olga_name="LIQUID ENTHALPY (J/KG)", 
        unit="J/kg",
        json_unit="J/mol",
        converter=lambda x, mw: x * 1000.0 / mw,  # J/mol → J/kg
        description="Liquid enthalpy"
    ),
    PropertyConfig(
        json_name="vapor_thermal_conductivity", 
        olga_name="GAS THERMAL CONDUCTIVITY (W/M/K)", 
        unit="W/(m·K)",
        description="Gas thermal conductivity"
    ),
    PropertyConfig(
        json_name="thermal_conductivity", 
        olga_name="GAS THERMAL CONDUCTIVITY (W/M/K)", 
        unit="W/(m·K)",
        description="Thermal conductivity (fallback for gas)"
    ),
    PropertyConfig(
        json_name="liquid_thermal_conductivity", 
        olga_name="LIQUID THERMAL CONDUCTIVITY (W/M/K)", 
        unit="W/(m·K)",
        description="Liquid thermal conductivity"
    ),
    PropertyConfig(
        json_name="surface_tension", 
        olga_name="SURFACE TENSION (N/M)", 
        unit="N/m",
        description="Surface tension"
    ),
    PropertyConfig(
        json_name="vapor_entropy", 
        olga_name="GAS ENTROPY (J/KG/C)", 
        unit="J/(kg·K)",
        json_unit="J/(mol·K)",
        converter=lambda x, mw: x * 1000.0 / mw,  # J/(mol·K) → J/(kg·K)
        description="Gas entropy"
    ),
    PropertyConfig(
        json_name="entropy", 
        olga_name="GAS ENTROPY (J/KG/C)", 
        unit="J/(kg·K)",
        json_unit="J/(mol·K)",
        converter=lambda x, mw: x * 1000.0 / mw,  # J/(mol·K) → J/(kg·K)
        description="Entropy (fallback for gas)"
    ),
    PropertyConfig(
        json_name="liquid_entropy", 
        olga_name="LIQUID ENTROPY (J/KG/C)", 
        unit="J/(kg·K)",
        json_unit="J/(mol·K)",
        converter=lambda x, mw: x * 1000.0 / mw,  # J/(mol·K) → J/(kg·K)
        description="Liquid entropy"
    )
]

def get_models_info(base_url: str) -> Dict[str, Any]:
    """Query the models_info endpoint to get fluid properties."""
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

def get_molecular_weight(base_url: str) -> float:
    """Get the molecular weight of the fluid from the API."""
    models_info = get_models_info(base_url)
    mw = DEFAULT_MOLECULAR_WEIGHT
    
    if models_info and 'models' in models_info and 'components' in models_info['models']:
        for comp in models_info['models']['components']:
            if 'molar_mass' in comp:
                mw = comp['molar_mass']
                print(f"API reports molar mass: {mw} g/mol")
                break
    
    print(f"Using molecular weight: {mw} g/mol")
    return mw

def query_pt_flash_json(
    base_url: str, 
    pressure_range: Tuple[float, float], 
    temperature: float,
    num_points: int = 20,
    property_list: Optional[List[str]] = None,
    verbose: bool = True
) -> Dict[str, Any]:
    """Query the PT flash endpoint with JSON response."""
    if verbose:
        print(f"Querying PT flash JSON at {temperature}°C, {pressure_range[0]}-{pressure_range[1]} bar")
    
    # If no property list specified, use all available
    if property_list is None:
        property_list = [
            "temperature", "pressure", "density", "vapor_density", "liquid_density",
            "phase", "vapor_fraction", "viscosity", "vapor_viscosity", "liquid_viscosity",
            "cp", "vapor_cp", "liquid_cp", "enthalpy", "vapor_enthalpy", "liquid_enthalpy",
            "thermal_conductivity", "vapor_thermal_conductivity", "liquid_thermal_conductivity",
            "surface_tension", "entropy", "vapor_entropy", "liquid_entropy"
        ]
    
    # Create request payload
    payload = {
        "composition": [
            {"fluid": "METHANE", "fraction": 1.0}
        ],
        "variables": {
            "pressure": {
                "range": {"from": pressure_range[0], "to": pressure_range[1]},
                "resolution": (pressure_range[1] - pressure_range[0]) / (num_points - 1) if num_points > 1 else 1
            },
            "temperature": {
                "range": {"from": temperature, "to": temperature},
                "resolution": 1
            }
        },
        "calculation": {
            "properties": property_list,
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
    property_list: Optional[List[str]] = None,
    verbose: bool = True
) -> str:
    """Query the PT flash endpoint with OLGA TAB response."""
    if verbose:
        print(f"Querying PT flash OLGA at {temperature}°C, {pressure_range[0]}-{pressure_range[1]} bar")
    
    # If no property list specified, use all available
    if property_list is None:
        property_list = [
            "temperature", "pressure", "density", "vapor_density", "liquid_density",
            "phase", "vapor_fraction", "viscosity", "vapor_viscosity", "liquid_viscosity",
            "cp", "vapor_cp", "liquid_cp", "enthalpy", "vapor_enthalpy", "liquid_enthalpy",
            "thermal_conductivity", "vapor_thermal_conductivity", "liquid_thermal_conductivity",
            "surface_tension", "entropy", "vapor_entropy", "liquid_entropy"
        ]
    
    # Create request payload
    payload = {
        "composition": [
            {"fluid": "METHANE", "fraction": 1.0}
        ],
        "variables": {
            "pressure": {
                "range": {"from": pressure_range[0], "to": pressure_range[1]},
                "resolution": (pressure_range[1] - pressure_range[0]) / (num_points - 1) if num_points > 1 else 1
            },
            "temperature": {
                "range": {"from": temperature, "to": temperature},
                "resolution": 1
            }
        },
        "calculation": {
            "properties": property_list,
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

def parse_olga_scientific(value_str: str) -> float:
    """
    Parse a value in OLGA scientific notation.
    
    Args:
        value_str: Value string in OLGA notation (.123456E+01)
        
    Returns:
        Float value correctly parsed
    """
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
        # In OLGA format, the decimal point is at the start, add a leading zero
        mantissa = float("0" + mantissa_str)
        exponent = int(exponent_str)
        
        # No exponent adjustment needed - OLGA format is correctly parsed this way
        # Calculate the final value
        return sign * mantissa * (10 ** exponent)
    except Exception as e:
        print(f"Error parsing OLGA value '{value_str}': {e}")
        return 0.0  # Return zero instead of None on error

def parse_olga_tab(
    content: str, 
    property_name: str,
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse OLGA TAB file to extract pressure and property values.
    
    Args:
        content: OLGA TAB content as string
        property_name: Name of the property to extract
        verbose: Whether to print verbose debug info
        
    Returns:
        (pressures, property_values) arrays
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
            
            if verbose and line_index <= 4:
                print(f"Line {line_index} (Pressures): Found {len(values)} values: {values[:3] if len(values) >= 3 else values}")
                
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
    
    # Show the raw property block data for debugging
    if verbose:
        print("\nRAW OLGA VALUES (first few):")
        for i in range(min(3, len(lines) - property_start)):
            if property_start + i < len(lines):
                print(f"Line {property_start + i}: {lines[property_start + i]}")
    
    # Parse the property values - we want the first row (temperature = first value)
    property_values = []
    line_index = property_start
    
    while len(property_values) < n_pressure and line_index < len(lines):
        line = lines[line_index].strip()
        if line:
            # Extract values using regex for OLGA scientific notation
            values = re.findall(r'[+-]?\.\d+E[+-]\d+', line)
            
            if verbose and len(property_values) == 0:
                print(f"Property line {line_index}: Found {len(values)} values: {values[:3] if len(values) >= 3 else values}")
                
            for value_str in values:
                if not value_str:
                    continue
                
                try:
                    value = parse_olga_scientific(value_str)
                    property_values.append(value)
                    
                    if verbose and len(property_values) <= 3:
                        print(f"Parsed property value: {value_str} → {value}")
                        
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
            print(f"First few property values: {property_values[:min(5, len(property_values))]}")
        else:
            print("No property values found")
    
    if len(property_values) < n_pressure:
        print(f"Warning: Only found {len(property_values)} property values, expected {n_pressure}")
        # Pad with zeros if necessary
        property_values.extend([0.0] * (n_pressure - len(property_values)))
    
    return np.array(pressure_values), np.array(property_values)

def extract_json_values(
    json_data: Dict[str, Any], 
    property_config: PropertyConfig,
    molecular_weight: float = 16.04, 
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract pressure and property values from JSON results.
    
    Args:
        json_data: JSON response from PT flash endpoint
        property_config: Property configuration
        molecular_weight: Molecular weight for unit conversions
        verbose: Whether to print verbose debug info
        
    Returns:
        (pressures, converted_values, original_values) arrays
    """
    pressures = []
    original_values = []  # Values in original units
    converted_values = []  # Values in OLGA-compatible units
    
    if not json_data or 'results' not in json_data:
        return np.array([]), np.array([]), np.array([])
    
    property_name = property_config.json_name
    
    for result in json_data['results']:
        # Extract pressure
        if 'pressure' in result:
            if isinstance(result['pressure'], dict):
                pressure = result['pressure'].get('value')
            else:
                pressure = result['pressure']
        else:
            continue
            
        # Get phase information if available
        if 'phase' in result:
            if isinstance(result['phase'], dict):
                phase = result['phase'].get('value', '')
            else:
                phase = result['phase']
        else:
            phase = ""
            
        phase = phase.lower() if isinstance(phase, str) else ''
            
        # Try to get the property value
        original_value = None
        
        # Try the primary property name
        if property_name in result:
            prop_field = result[property_name]
            if isinstance(prop_field, dict):
                original_value = prop_field.get('value')
                prop_unit = prop_field.get('unit', '')
            else:
                original_value = prop_field
        
        # Apply phase filtering if needed
        if original_value is not None:
            # Skip if this is a phase-specific property that doesn't apply to the current phase
            if property_name.startswith('vapor_') and 'liquid' in phase and 'two' not in phase and 'super' not in phase:
                continue
            if property_name.startswith('liquid_') and 'vapor' in phase and 'two' not in phase and 'super' not in phase:
                continue
                
            # Convert to OLGA units if needed
            if property_config.converter:
                converted_value = property_config.converter(original_value, molecular_weight)
            else:
                converted_value = original_value
                
            # Store values
            pressures.append(pressure)
            original_values.append(original_value)
            converted_values.append(converted_value)
            
            if verbose and len(pressures) <= 5:
                print(f"JSON: P={pressure:.2f} bar, phase={phase}, {property_name}={original_value:.6f}")
                if property_config.converter:
                    print(f"  - Converted to {property_config.unit}: {converted_value:.6f}")
    
    # Sort by pressure
    if pressures:
        sorted_indices = np.argsort(pressures)
        sorted_pressures = np.array(pressures)[sorted_indices]
        sorted_original = np.array(original_values)[sorted_indices]
        sorted_converted = np.array(converted_values)[sorted_indices]
        return sorted_pressures, sorted_converted, sorted_original
    
    return np.array([]), np.array([]), np.array([])

def find_closest_pressure_index(olga_pressures, json_pressure, tolerance=1e-2):
    """Find the index of the closest pressure value within tolerance."""
    differences = np.abs(olga_pressures - json_pressure)
    min_diff_idx = np.argmin(differences)
    
    if differences[min_diff_idx] <= tolerance:
        return min_diff_idx
    return None

def compare_specific_points(
    json_pressures: np.ndarray, 
    json_values: np.ndarray,
    olga_pressures: np.ndarray, 
    olga_values: np.ndarray,
    property_unit: str,
    num_points: int = 5,
    tolerance: float = 1e-2  # Tolerance for matching pressure points (in bar)
) -> None:
    """
    Compare specific points between JSON and OLGA outputs.
    
    Args:
        json_pressures: Pressure values from JSON response
        json_values: Property values from JSON response
        olga_pressures: Pressure values from OLGA response
        olga_values: Property values from OLGA response
        property_unit: Unit string for the property
        num_points: Number of points to show in the comparison
        tolerance: Tolerance for matching pressure points
    """
    print("\n===== POINT-BY-POINT COMPARISON =====")
    print(f"Pressure (bar) | JSON ({property_unit}) | OLGA ({property_unit}) | Difference | Ratio")
    print("---------------|-----------------|-----------------|------------|-------")
    
    # Find points with matching pressures using tolerance
    comparison_points = []
    
    for i, jp in enumerate(json_pressures):
        olga_idx = find_closest_pressure_index(olga_pressures, jp, tolerance)
        if olga_idx is not None:
            comparison_points.append((jp, json_values[i], olga_values[olga_idx], olga_pressures[olga_idx]))
    
    # Sort by pressure and limit to num_points
    comparison_points.sort(key=lambda x: x[0])
    for p_json, jv, ov, p_olga in comparison_points[:num_points]:
        # Skip if olga value is zero (to avoid division by zero)
        if ov == 0 and jv == 0:
            ratio = 1.0
        elif ov == 0:
            ratio = float('inf')
        else:
            ratio = jv / ov
            
        diff = jv - ov
        print(f"{p_json:13.2f} | {jv:15.6f} | {ov:15.6f} | {diff:10.6f} | {ratio:7.4f}")
    
    # Calculate average ratio if possible
    if comparison_points:
        valid_ratios = []
        for _, jv, ov, _ in comparison_points:
            if ov != 0 and jv != 0:  # Only include non-zero values
                valid_ratios.append(jv/ov)
                
        if valid_ratios:
            avg_ratio = sum(valid_ratios) / len(valid_ratios)
            print(f"\nAverage ratio (JSON/OLGA): {avg_ratio:.4f}")
            if abs(avg_ratio - 1.0) > 0.05:  # 5% threshold
                print(f"WARNING: Values differ by a factor of approximately {avg_ratio:.4f}")
            else:
                print("Values match within 5% - conversion is working correctly")
        else:
            print("\nNo valid ratios found for comparison (all values are zero)")
    else:
        print("No matching pressure points found for comparison")

def plot_comparison(
    json_pressures: np.ndarray, 
    json_values: np.ndarray,
    olga_pressures: np.ndarray, 
    olga_values: np.ndarray,
    property_config: PropertyConfig,
    temperature: float,
    output_file: Optional[str] = None
) -> None:
    """
    Plot comparison of property values between JSON and OLGA outputs.
    
    Args:
        json_pressures: Pressure values from JSON response
        json_values: Property values from JSON response
        olga_pressures: Pressure values from OLGA response
        olga_values: Property values from OLGA response
        property_config: Property configuration
        temperature: Temperature for the plot title
        output_file: Filename for saving the plot (optional)
    """
    plt.figure(figsize=(12, 8))
    
    # Plot both datasets
    plt.plot(json_pressures, json_values, 'bo-', label='JSON Response', markersize=6)
    plt.plot(olga_pressures, olga_values, 'rx--', label='OLGA TAB Response', markersize=6)
    
    # Labels and title
    plt.xlabel('Pressure (bar)', fontsize=12)
    plt.ylabel(f'{property_config.description} ({property_config.unit})', fontsize=12)
    plt.title(f'Methane {property_config.description} Comparison at {temperature}°C', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=12)
    plt.tight_layout()
    
    # Save or show
    if output_file:
        plt.savefig(output_file)
        print(f"\nPlot saved to {output_file}")
    else:
        plt.show()

def test_property(
    base_url: str,
    property_config: PropertyConfig,
    molecular_weight: float,
    pressure_range: Tuple[float, float] = PRESSURE_RANGE,
    temperature: float = TEMP_CELSIUS,
    num_points: int = NUM_POINTS,
    verbose: bool = True
) -> None:
    """
    Test a specific property by comparing JSON and OLGA TAB responses.
    
    Args:
        base_url: API base URL
        property_config: Property configuration
        molecular_weight: Molecular weight for unit conversions
        pressure_range: Pressure range for testing
        temperature: Temperature for testing
        num_points: Number of pressure points
        verbose: Whether to print verbose debug info
    """
    print(f"\n\n{'='*80}")
    print(f"TESTING PROPERTY: {property_config.json_name} ({property_config.olga_name})")
    print(f"{'='*80}")
    
    # Define property list with phase information
    property_list = [
        "temperature", "pressure", "phase", "vapor_fraction",
        property_config.json_name
    ]
    
    # Add fallback properties if needed
    if property_config.json_name.startswith('vapor_'):
        base_property = property_config.json_name[6:]  # Remove 'vapor_' prefix
        if base_property not in property_list:
            property_list.append(base_property)
            
    if property_config.json_name.startswith('liquid_'):
        base_property = property_config.json_name[7:]  # Remove 'liquid_' prefix
        if base_property not in property_list:
            property_list.append(base_property)
    
    # Query both endpoints
    json_data = query_pt_flash_json(base_url, pressure_range, temperature, num_points, property_list, verbose)
    olga_content = query_pt_flash_olga(base_url, pressure_range, temperature, num_points, property_list, verbose)
    
    if not json_data or 'results' not in json_data:
        print("Error: No JSON results received")
        return
        
    if not olga_content:
        print("Error: No OLGA TAB content received")
        return
    
    # Extract values from both responses
    json_pressures, json_values, json_original = extract_json_values(
        json_data, property_config, molecular_weight, verbose
    )
    
    olga_pressures, olga_values = parse_olga_tab(olga_content, property_config.olga_name, verbose)
    
    if len(json_pressures) == 0:
        print(f"Error: No data found for {property_config.json_name} in JSON response")
        return
        
    if len(olga_pressures) == 0:
        print(f"Error: No data found for {property_config.olga_name} in OLGA response")
        return
    
    # Compare specific points
    compare_specific_points(
        json_pressures, json_values,
        olga_pressures, olga_values,
        property_config.unit
    )
    
    # Plot comparison
    output_file = f"{property_config.json_name}_comparison.png"
    plot_comparison(
        json_pressures, json_values,
        olga_pressures, olga_values,
        property_config, temperature,
        output_file
    )
    
    print(f"\nTest complete for {property_config.json_name}")

def main() -> None:
    """Main function to run property tests."""
    # Parse command-line arguments
    if len(sys.argv) > 1:
        base_url = sys.argv[1]
    else:
        base_url = DEFAULT_API_URL
    
    # Additional arguments (optional)
    property_to_test = sys.argv[2] if len(sys.argv) > 2 else None
    
    print(f"Using API URL: {base_url}")
    
    # Get molecular weight for unit conversions
    molecular_weight = get_molecular_weight(base_url)
    
    # Run tests for specified properties
    if property_to_test:
        # Test a single property if specified
        for prop_config in PROPERTIES_TO_TEST:
            if prop_config.json_name == property_to_test:
                test_property(base_url, prop_config, molecular_weight)
                break
        else:
            print(f"Property {property_to_test} not found in test configuration")
            print("Available properties:")
            for prop_config in PROPERTIES_TO_TEST:
                print(f"  - {prop_config.json_name}")
    else:
        # Test all primary properties (skipping fallbacks)
        tested_properties = set()
        for prop_config in PROPERTIES_TO_TEST:
            # Skip fallback properties (those that don't have vapor_ or liquid_ prefix)
            if (not prop_config.json_name.startswith('vapor_') and 
                not prop_config.json_name.startswith('liquid_') and
                prop_config.json_name in tested_properties):
                continue
                
            test_property(base_url, prop_config, molecular_weight)
            tested_properties.add(prop_config.json_name)
    
    print("\nAll property tests complete.")

if __name__ == "__main__":
    main()