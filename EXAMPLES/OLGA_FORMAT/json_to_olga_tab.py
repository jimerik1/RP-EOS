#!/usr/bin/env python3
"""
Convert an existing Span-Wagner EOS API JSON response to OLGA TAB format.

This script can be used to convert previously calculated PT-flash, PH-flash, 
or other calculation results from JSON format to OLGA TAB format without 
recalculating the values.
"""

import json
import sys
import numpy as np
import argparse
import re

def format_grid_values(values, values_per_line=5):
    """
    Format grid values for OLGA TAB format.
    
    Args:
        values (list): List of values to format
        values_per_line (int): Number of values per line
        
    Returns:
        str: Formatted grid values
    """
    formatted = ""
    values_array = np.asarray(values).flatten()
    
    for i in range(0, len(values_array), values_per_line):
        line_values = values_array[i:i + values_per_line]
        formatted += "    " + "    ".join([f".{value:.6E}" for value in line_values]) + "\n"
    
    return formatted

def extract_value(item, key):
    """Extract value from a nested dictionary"""
    if key in item:
        value = item[key]
        if isinstance(value, dict) and 'value' in value:
            return value['value']
        return value
    return None

def convert_json_to_olga_tab(input_file, output_file, composition_str=None, molar_mass=None):
    """
    Convert a JSON API response to OLGA TAB format.
    
    Args:
        input_file (str): Path to input JSON file
        output_file (str): Path to output TAB file
        composition_str (str): Optional composition string (e.g., "WATER-1.0")
        molar_mass (float): Optional molar mass in g/mol for unit conversions
    """
    # Load JSON data
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    # Check if the data has the expected format
    if 'results' not in data:
        print("Error: Input JSON does not contain 'results' array. Is this a valid API response?")
        sys.exit(1)
    
    results = data['results']
    if not results:
        print("Error: Results array is empty.")
        sys.exit(1)
    
    # Try to determine grid dimensions by finding max indices
    max_index = max(item.get('index', i) for i, item in enumerate(results))
    
    # Try to extract pressure and temperature points from results
    pressures = sorted(list(set([extract_value(item, 'pressure') for item in results if extract_value(item, 'pressure') is not None])))
    temperatures = sorted(list(set([extract_value(item, 'temperature') for item in results if extract_value(item, 'temperature') is not None])))
    
    np_grid = len(pressures)
    nt_grid = len(temperatures)
    
    if np_grid * nt_grid != len(results):
        print(f"Warning: Grid dimensions ({np_grid}x{nt_grid}) do not match the number of results ({len(results)}).")
        print("Attempting to continue, but the output file may be incorrect.")
    
    # Default composition string if not provided
    if not composition_str:
        composition_str = "UNKNOWN FLUID"
    
    # Create OLGA TAB header
    olga_tab = f"'Span-Wagner EOS {composition_str}'\n"
    olga_tab += f"   {np_grid}  {nt_grid}    .276731E-08\n"
    
    # Pressure grid in Pa (convert from bar)
    pressure_pa = [p * 1e5 for p in pressures]
    olga_tab += format_grid_values(pressure_pa)
    
    # Temperature grid in °C
    olga_tab += format_grid_values(temperatures)
    
    # Define OLGA TAB property headers to look for in the results
    olga_properties = [
        'LIQUID DENSITY (KG/M3)',
        'GAS DENSITY (KG/M3)',
        'PRES. DERIV. OF LIQUID DENS.',
        'PRES. DERIV. OF GAS DENS.',
        'TEMP. DERIV. OF LIQUID DENS.',
        'TEMP. DERIV. OF GAS DENS.',
        'GAS MASS FRACTION OF GAS +LIQUID',
        'WATER MASS FRACTION OF LIQUID',
        'WATER MASS FRACTION OF GAS',
        'LIQUID VISCOSITY (N S/M2)',
        'GAS VISCOSITY (N S/M2)',
        'LIQUID SPECIFIC HEAT (J/KG K)',
        'GAS SPECIFIC HEAT (J/KG K)',
        'LIQUID ENTHALPY (J/KG)',
        'GAS ENTHALPY (J/KG)',
        'LIQUID THERMAL COND. (W/M K)',
        'GAS THERMAL COND. (W/M K)',
        'SURFACE TENSION GAS/LIQUID (N/M)',
        'LIQUID ENTROPY (J/KG/C)',
        'GAS ENTROPY (J/KG/C)'
    ]
    
    # Mapping of property names to field names in the JSON
    property_mapping = {
        'LIQUID DENSITY (KG/M3)': 'liquid_density',
        'GAS DENSITY (KG/M3)': 'vapor_density',
        'PRES. DERIV. OF LIQUID DENS.': 'dDdP',
        'PRES. DERIV. OF GAS DENS.': 'dDdP',
        'TEMP. DERIV. OF LIQUID DENS.': 'dDdT',
        'TEMP. DERIV. OF GAS DENS.': 'dDdT',
        'GAS MASS FRACTION OF GAS +LIQUID': 'vapor_fraction',
        'WATER MASS FRACTION OF LIQUID': None,  # Need special handling
        'WATER MASS FRACTION OF GAS': None,  # Need special handling
        'LIQUID VISCOSITY (N S/M2)': 'viscosity',
        'GAS VISCOSITY (N S/M2)': 'viscosity',
        'LIQUID SPECIFIC HEAT (J/KG K)': 'cp',
        'GAS SPECIFIC HEAT (J/KG K)': 'cp',
        'LIQUID ENTHALPY (J/KG)': 'enthalpy',
        'GAS ENTHALPY (J/KG)': 'enthalpy',
        'LIQUID THERMAL COND. (W/M K)': 'thermal_conductivity',
        'GAS THERMAL COND. (W/M K)': 'thermal_conductivity',
        'SURFACE TENSION GAS/LIQUID (N/M)': 'surface_tension',
        'LIQUID ENTROPY (J/KG/C)': 'entropy',
        'GAS ENTROPY (J/KG/C)': 'entropy'
    }
    
    # Initialize property arrays with zeros
    property_arrays = {}
    for prop in olga_properties:
        property_arrays[prop] = np.zeros((np_grid, nt_grid))
    
    # Convert molar units to mass units if molar_mass is provided
    convert_to_mass = molar_mass is not None
    
    # Populate property arrays
    for item in results:
        # Get indices
        p_idx = item.get('p_idx', item.get('index', 0) // nt_grid)
        t_idx = item.get('t_idx', item.get('index', 0) % nt_grid)
        
        # Get phase information (liquid, vapor, two-phase)
        phase = extract_value(item, 'phase')
        if isinstance(phase, str):
            phase = phase.lower()
        else:
            # Try to determine phase from vapor fraction
            q = extract_value(item, 'vapor_fraction')
            if q == 0:
                phase = 'liquid'
            elif q == 1:
                phase = 'vapor'
            elif isinstance(q, (int, float)) and 0 < q < 1:
                phase = 'two-phase'
            else:
                phase = 'unknown'
        
        # Process each property
        for prop in olga_properties:
            # Get the corresponding field name in the JSON
            field = property_mapping.get(prop)
            if not field:
                continue
                
            # Look for gas/liquid specific properties
            if 'LIQUID' in prop and phase == 'liquid':
                value = extract_value(item, field)
            elif 'GAS' in prop and phase == 'vapor':
                value = extract_value(item, field)
            elif 'SURFACE TENSION' in prop and phase == 'two-phase':
                value = extract_value(item, field)
            elif 'GAS MASS FRACTION' in prop:
                value = extract_value(item, 'vapor_fraction')
            else:
                value = extract_value(item, field)
            
            # Convert molar to mass units if needed
            if convert_to_mass and value is not None:
                if 'DENSITY' in prop:
                    value = value * molar_mass / 1000.0
                elif 'SPECIFIC HEAT' in prop:
                    value = value * 1000.0 / molar_mass
                elif 'ENTHALPY' in prop or 'ENTROPY' in prop:
                    value = value * 1000.0 / molar_mass
            
            # Store the value in the property array
            if value is not None:
                try:
                    property_arrays[prop][p_idx, t_idx] = float(value)
                except:
                    pass
    
    # Add properties to OLGA TAB
    for prop in olga_properties:
        olga_tab += f" {prop}                \n"
        olga_tab += format_grid_values(property_arrays[prop].flatten())
    
    # Write to output file
    with open(output_file, 'w') as f:
        f.write(olga_tab)
    
    print(f"Successfully converted {input_file} to OLGA TAB format: {output_file}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert JSON API response to OLGA TAB format")
    parser.add_argument("input", help="Input JSON file")
    parser.add_argument("--output", "-o", default="output.tab", help="Output TAB file")
    parser.add_argument("--composition", "-c", help="Composition string (e.g., 'WATER-1.0')")
    parser.add_argument("--molar-mass", "-m", type=float, help="Molar mass in g/mol for unit conversions")
    args = parser.parse_args()
    
    # Convert the file
    convert_json_to_olga_tab(args.input, args.output, args.composition, args.molar_mass)

if __name__ == "__main__":
    main()