#!/usr/bin/env python3
import requests
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from matplotlib.lines import Line2D
import json  # For debugging API responses

# API endpoint for PT-flash calculations
PT_FLASH_URL = "http://localhost:5051/pt_flash"

# Create output directory for plots if it doesn't exist
OUTPUT_DIR = "property_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define the composition to use for all calculations
COMPOSITION = [
    {"fluid": "CO2", "fraction": 0.95},
    {"fluid": "NITROGEN", "fraction": 0.05}
]

# Define the temperatures and pressure range to analyze
TEMPERATURES = [5, 10, 20, 30, 40, 50]  # in °C
PRESSURE_RANGE = {
    "min": 1,    # minimum pressure in bar
    "max": 200,  # maximum pressure in bar
    "step": 2    # pressure step in bar
}

# Define the properties to calculate and plot
# Adding alternative property names that might be used in the API
PROPERTIES = [
    {"name": "density", "alt_names": ["density"], "unit": "mol/L", "title": "Density", "y_label": "Density (mol/L)"},
    {"name": "enthalpy", "alt_names": ["enthalpy", "h"], "unit": "J/mol", "title": "Enthalpy", "y_label": "Enthalpy (J/mol)"},
    {"name": "internal_energy", "alt_names": ["internal_energy", "e"], "unit": "J/mol", "title": "Internal Energy", "y_label": "Internal Energy (J/mol)"},
    {"name": "entropy", "alt_names": ["entropy", "s"], "unit": "J/(mol·K)", "title": "Entropy", "y_label": "Entropy (J/(mol·K))"},
    {"name": "Cv", "alt_names": ["Cv", "cv", "isochoric_heat_capacity"], "unit": "J/(mol·K)", "title": "Isochoric Heat Capacity (Cv)", "y_label": "Cv (J/(mol·K))"},
    {"name": "Cp", "alt_names": ["Cp", "cp", "isobaric_heat_capacity"], "unit": "J/(mol·K)", "title": "Isobaric Heat Capacity (Cp)", "y_label": "Cp (J/(mol·K))"},
    {"name": "sound_speed", "alt_names": ["sound_speed", "w", "speed_of_sound"], "unit": "m/s", "title": "Sound Speed", "y_label": "Sound Speed (m/s)"},
    {"name": "viscosity", "alt_names": ["viscosity", "visc"], "unit": "μPa·s", "title": "Viscosity", "y_label": "Viscosity (μPa·s)"},
    {"name": "thermal_conductivity", "alt_names": ["thermal_conductivity", "tcx"], "unit": "W/(m·K)", "title": "Thermal Conductivity", "y_label": "Thermal Conductivity (W/(m·K))"},
    {"name": "compressibility_factor", "alt_names": ["compressibility_factor", "Z", "z_factor"], "unit": "dimensionless", "title": "Compressibility Factor", "y_label": "Z-factor"}
]

def get_property_value(data_point, prop_info):
    """Helper function to get property value trying all alternative names"""
    for name in prop_info["alt_names"]:
        if name in data_point:
            return data_point[name]["value"]
    return np.nan

def calculate_properties_at_temperature(temperature):
    """
    Calculate properties at a specific temperature across the pressure range
    """
    print(f"Calculating properties at {temperature}°C...")
    
    # Define pressures
    pressures = np.arange(
        PRESSURE_RANGE["min"], 
        PRESSURE_RANGE["max"] + PRESSURE_RANGE["step"], 
        PRESSURE_RANGE["step"]
    )
    
    # Create property list for API request
    property_names = []
    for prop in PROPERTIES:
        property_names.extend(prop["alt_names"])
    # Remove duplicates while preserving order
    property_names = list(dict.fromkeys(property_names))
    
    # Prepare payload for API request
    payload = {
        "composition": COMPOSITION,
        "pressure_range": {
            "from": PRESSURE_RANGE["min"],
            "to": PRESSURE_RANGE["max"]
        },
        "temperature_range": {
            "from": temperature,
            "to": temperature
        },
        "pressure_resolution": PRESSURE_RANGE["step"],
        "temperature_resolution": 1,
        "properties": property_names + ["phase", "vapor_fraction"],
        "units_system": "SI"
    }
    
    try:
        # Make API request
        response = requests.post(PT_FLASH_URL, json=payload, timeout=60)
        response.raise_for_status()
        data = response.json()
        
        if "results" not in data or not data["results"]:
            print(f"  No data returned for temperature {temperature}°C")
            return None
        
        # Debug: Check what properties are actually in the response for the first point
        if temperature == TEMPERATURES[0]:
            print("\nProperties in API response:")
            sample_point = data["results"][0]
            property_keys = [key for key in sample_point.keys() if key not in ["index", "temperature", "pressure", "phase", "vapor_fraction"]]
            print(f"  {', '.join(property_keys)}")
            
        return data["results"]
        
    except Exception as e:
        print(f"  Error calculating properties at {temperature}°C: {str(e)}")
        return None

def organize_data_by_phase(results):
    """
    Organize results into liquid and vapor phases
    """
    if not results:
        return None
        
    data = {
        "pressure": [],
        "liquid": {prop["name"]: [] for prop in PROPERTIES},
        "vapor": {prop["name"]: [] for prop in PROPERTIES},
        "phase": []
    }
    
    for point in results:
        pressure = point["pressure"]["value"]
        data["pressure"].append(pressure)
        
        # Store phase information
        phase = point["phase"]["value"]
        data["phase"].append(phase)
        
        # For each property, store values by phase
        for prop in PROPERTIES:
            prop_value = get_property_value(point, prop)
            
            # Determine whether to store as liquid or vapor based on phase
            if phase == "Liquid" or (phase == "Two-Phase" and point["vapor_fraction"]["value"] == 0):
                data["liquid"][prop["name"]].append(prop_value)
                data["vapor"][prop["name"]].append(np.nan)
            elif phase == "Vapor" or (phase == "Two-Phase" and point["vapor_fraction"]["value"] == 1):
                data["vapor"][prop["name"]].append(prop_value)
                data["liquid"][prop["name"]].append(np.nan)
            elif phase == "Two-Phase":
                # For two-phase, store the same value in both liquid and vapor
                # since we can't distinguish them from the results
                data["liquid"][prop["name"]].append(prop_value)
                data["vapor"][prop["name"]].append(prop_value)
            else:
                # For supercritical or unknown, store as vapor phase
                data["vapor"][prop["name"]].append(prop_value)
                data["liquid"][prop["name"]].append(np.nan)
    
    return data

def plot_property(property_info, data_by_temperature):
    """
    Create a plot for a specific property across all temperatures
    """
    prop_name = property_info["name"]
    plt.figure(figsize=(12, 8))
    
    # Create a colormap for different temperatures
    viridis = plt.cm.viridis
    colors = [viridis(i/len(TEMPERATURES)) for i in range(len(TEMPERATURES))]
    
    # Create custom legend for phases
    legend_elements = [
        Line2D([0], [0], color='k', linestyle='-', label='Liquid'),
        Line2D([0], [0], color='k', linestyle='--', label='Vapor'),
    ]
    
    # Check if we have any valid data for this property
    has_data = False
    
    # Plot for each temperature
    for i, temp in enumerate(TEMPERATURES):
        if temp not in data_by_temperature or data_by_temperature[temp] is None:
            continue
            
        data = data_by_temperature[temp]
        
        # Convert lists to numpy arrays
        pressures = np.array(data["pressure"])
        
        # Plot liquid phase - use numpy arrays to properly handle missing data
        liquid_values = np.array(data["liquid"][prop_name], dtype=float)
        # Create a mask for valid data points
        valid_liquid = ~np.isnan(liquid_values)
        
        if np.any(valid_liquid):
            has_data = True
            plt.plot(
                pressures[valid_liquid], 
                liquid_values[valid_liquid], 
                '-', 
                color=colors[i], 
                label=f"{temp}°C (Liquid)" if i == 0 else None,
                linewidth=2
            )
        
        # Plot vapor phase - use numpy arrays to properly handle missing data
        vapor_values = np.array(data["vapor"][prop_name], dtype=float)
        # Create a mask for valid data points
        valid_vapor = ~np.isnan(vapor_values)
        
        if np.any(valid_vapor):
            has_data = True
            plt.plot(
                pressures[valid_vapor], 
                vapor_values[valid_vapor], 
                '--', 
                color=colors[i], 
                label=f"{temp}°C (Vapor)" if i == 0 else None,
                linewidth=2
            )
    
    if not has_data:
        plt.text(0.5, 0.5, f"No data available for {property_info['title']}", 
                 horizontalalignment='center', verticalalignment='center',
                 transform=plt.gca().transAxes, fontsize=14)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        print(f"Warning: No data available for {property_info['title']}")
    else:
        # Add temperature legend
        temp_legend_elements = [
            Line2D([0], [0], color=colors[i], linestyle='-', label=f"{temp}°C")
            for i, temp in enumerate(TEMPERATURES)
            if temp in data_by_temperature and data_by_temperature[temp] is not None
        ]
        
        # Add grid and labels
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.xlabel("Pressure (bar)")
        plt.ylabel(property_info["y_label"])
        
        # Add two separate legends: one for phases, one for temperatures
        phase_legend = plt.legend(handles=legend_elements, loc='upper right')
        plt.gca().add_artist(phase_legend)
        plt.legend(handles=temp_legend_elements, loc='upper left')
    
    plt.title(f"{property_info['title']} vs Pressure for CO2/N2 (95%/5%)")
    
    # Save plot
    filename = f"{OUTPUT_DIR}/{prop_name}_vs_pressure.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to {filename}")

def main():
    print(f"Calculating properties for {len(TEMPERATURES)} temperatures and {len(PROPERTIES)} properties...")
    
    # Calculate properties for each temperature
    data_by_temperature = {}
    for temp in TEMPERATURES:
        results = calculate_properties_at_temperature(temp)
        if results:
            data_by_temperature[temp] = organize_data_by_phase(results)
            print(f"  Finished calculations for {temp}°C")
        time.sleep(0.5)  # Small delay to avoid overwhelming the API
    
    # Plot each property
    for property_info in PROPERTIES:
        print(f"Plotting {property_info['title']}...")
        plot_property(property_info, data_by_temperature)
    
    print(f"All plots have been generated in the '{OUTPUT_DIR}' directory.")

if __name__ == "__main__":
    main()