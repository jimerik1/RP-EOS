#!/usr/bin/env python3
import requests
import numpy as np
import pandas as pd
import os
import time
from tabulate import tabulate
import matplotlib.pyplot as plt
import seaborn as sns

# API endpoint for PT-flash calculations
PT_FLASH_URL = "http://localhost:5051/pt_flash"

# Create output directory for results
OUTPUT_DIR = "thermal_conductivity_data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define the composition to use for all calculations
COMPOSITION = [
    {"fluid": "CO2", "fraction": 0.95},
    {"fluid": "NITROGEN", "fraction": 0.05}
]

# Define the temperatures and pressure grid
TEMPERATURES = [-30, -10, 0, 10, 20, 30, 40, 50, 75, 100]  # in °C
PRESSURES = [1, 5, 10, 20, 40, 60, 80, 100, 120, 150, 200]  # in bar

# Define properties we want to extract
PROPERTIES = ["thermal_conductivity", "phase", "vapor_fraction", "density"]

def get_thermal_conductivity_grid():
    """
    Calculate thermal conductivity and related properties for a grid of P-T values
    """
    print(f"Calculating thermal conductivity for {len(TEMPERATURES)} temperatures and {len(PRESSURES)} pressures...")
    
    # Create result arrays
    tc_values = np.zeros((len(TEMPERATURES), len(PRESSURES)))
    phases = np.empty((len(TEMPERATURES), len(PRESSURES)), dtype=object)
    densities = np.zeros((len(TEMPERATURES), len(PRESSURES)))
    
    # Loop through temperatures
    for i, temperature in enumerate(TEMPERATURES):
        print(f"Processing temperature {temperature}°C...")
        
        # Prepare payload for API request
        payload = {
            "composition": COMPOSITION,
            "temperature_range": {
                "from": temperature,
                "to": temperature
            },
            "pressure_range": {
                "from": min(PRESSURES),
                "to": max(PRESSURES)
            },
            "temperature_resolution": 1,
            "pressure_resolution": 1,  # We'll filter the exact pressure points later
            "properties": PROPERTIES,
            "units_system": "SI"
        }
        
        try:
            # Make API request
            response = requests.post(PT_FLASH_URL, json=payload, timeout=60)
            response.raise_for_status()
            data = response.json()
            
            if "results" not in data or not data["results"]:
                print(f"  No data returned for temperature {temperature}°C")
                continue
                
            results = data["results"]
            
            # Extract the values for each pressure
            for j, pressure in enumerate(PRESSURES):
                # Find the closest point to our target pressure
                closest_point = min(results, 
                                   key=lambda x: abs(x["pressure"]["value"] - pressure))
                
                # Check if the pressure is close enough (within 1% or 1 bar, whichever is smaller)
                tolerance = min(1.0, pressure * 0.01)
                if abs(closest_point["pressure"]["value"] - pressure) > tolerance:
                    print(f"  Warning: No exact match for P={pressure} bar at T={temperature}°C")
                    # Skip this point or use None/NaN
                    tc_values[i, j] = np.nan
                    phases[i, j] = "N/A"
                    densities[i, j] = np.nan
                    continue
                
                # Extract thermal conductivity
                if "thermal_conductivity" in closest_point:
                    tc_values[i, j] = closest_point["thermal_conductivity"]["value"]
                else:
                    # Try alternative names
                    if "tcx" in closest_point:
                        tc_values[i, j] = closest_point["tcx"]["value"]
                    else:
                        tc_values[i, j] = np.nan
                
                # Extract phase
                phases[i, j] = closest_point["phase"]["value"]
                
                # Extract density
                if "density" in closest_point:
                    densities[i, j] = closest_point["density"]["value"]
                else:
                    densities[i, j] = np.nan
                
        except Exception as e:
            print(f"  Error for temperature {temperature}°C: {str(e)}")
            # Fill with NaN for this temperature row
            tc_values[i, :] = np.nan
            phases[i, :] = "Error"
            densities[i, :] = np.nan
        
        # Small delay to avoid overwhelming the API
        time.sleep(0.5)
    
    return tc_values, phases, densities

def create_tables_and_plots(tc_values, phases, densities):
    """
    Create tables and visualization plots from the calculated data
    """
    # 1. Create pandas DataFrame for the thermal conductivity values
    tc_df = pd.DataFrame(tc_values, 
                        index=TEMPERATURES, 
                        columns=PRESSURES)
    tc_df.index.name = "Temperature (°C)"
    tc_df.columns.name = "Pressure (bar)"
    
    # 2. Create pandas DataFrame for the phases
    phases_df = pd.DataFrame(phases, 
                            index=TEMPERATURES, 
                            columns=PRESSURES)
    phases_df.index.name = "Temperature (°C)"
    phases_df.columns.name = "Pressure (bar)"
    
    # 3. Create pandas DataFrame for the densities
    density_df = pd.DataFrame(densities, 
                             index=TEMPERATURES, 
                             columns=PRESSURES)
    density_df.index.name = "Temperature (°C)"
    density_df.columns.name = "Pressure (bar)"
    
    # Save the DataFrames to CSV files
    tc_df.to_csv(f"{OUTPUT_DIR}/thermal_conductivity_table.csv")
    phases_df.to_csv(f"{OUTPUT_DIR}/phases_table.csv")
    density_df.to_csv(f"{OUTPUT_DIR}/density_table.csv")
    
    # Create pretty tables for the report
    print("\nThermal Conductivity [W/(m·K)] Table:")
    tc_table = tc_df.copy()
    # Format the values for better readability
    tc_table = tc_table.applymap(lambda x: f"{x:.5f}" if not np.isnan(x) else "N/A")
    print(tabulate(tc_table, headers=tc_table.columns, tablefmt="grid"))
    
    # Save pretty table to text file
    with open(f"{OUTPUT_DIR}/thermal_conductivity_table.txt", "w") as f:
        f.write("Thermal Conductivity [W/(m·K)] for CO2/N2 (95%/5%)\n\n")
        f.write(tabulate(tc_table, headers=tc_table.columns, tablefmt="grid"))
    
    print("\nPhase Table:")
    print(tabulate(phases_df, headers=phases_df.columns, tablefmt="grid"))
    
    # Create heatmap visualizations
    create_heatmaps(tc_df, phases_df, density_df)
    
    return tc_df, phases_df, density_df

def create_heatmaps(tc_df, phases_df, density_df):
    """
    Create heatmap visualizations from the data
    """
    # Set the style
    sns.set(style="white")
    
    # 1. Thermal conductivity heatmap
    plt.figure(figsize=(12, 8))
    ax = sns.heatmap(tc_df, annot=True, fmt=".5f", cmap="viridis", 
                    linewidths=.5, cbar_kws={'label': 'Thermal Conductivity [W/(m·K)]'})
    plt.title(f"Thermal Conductivity [W/(m·K)] for CO2/N2 (95%/5%)")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/thermal_conductivity_heatmap.png", dpi=300)
    plt.close()
    
    # 2. Phase visualization (categorical heatmap)
    plt.figure(figsize=(12, 8))
    # Convert phases to numeric values for coloring
    phase_map = {"Liquid": 0, "Vapor": 1, "Two-Phase": 2, "Supercritical": 3, "N/A": 4, "Error": 5}
    phase_numeric = phases_df.applymap(lambda x: phase_map.get(x, 4))
    
    # Create custom colormap for phases
    colors = ["blue", "red", "green", "purple", "gray", "black"]
    cmap = plt.cm.colors.ListedColormap(colors[:len(set(phase_map.values()))])
    
    ax = sns.heatmap(phase_numeric, cmap=cmap, cbar=False, linewidths=.5)
    
    # Add text annotations with the actual phase names
    for i in range(len(phases_df.index)):
        for j in range(len(phases_df.columns)):
            ax.text(j + 0.5, i + 0.5, phases_df.iloc[i, j],
                    ha="center", va="center", color="white")
    
    # Add a custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=colors[i], label=phase)
                      for phase, i in phase_map.items() if phase != "N/A" and phase != "Error"]
    plt.legend(handles=legend_elements, loc='upper right')
    
    plt.title(f"Phase Diagram for CO2/N2 (95%/5%)")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/phase_diagram.png", dpi=300)
    plt.close()
    
    # 3. Density heatmap
    plt.figure(figsize=(12, 8))
    ax = sns.heatmap(density_df, annot=True, fmt=".2f", cmap="rocket_r",
                    linewidths=.5, cbar_kws={'label': 'Density [mol/L]'})
    plt.title(f"Density [mol/L] for CO2/N2 (95%/5%)")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/density_heatmap.png", dpi=300)
    plt.close()
    
    print(f"\nVisualization plots saved to {OUTPUT_DIR}/")

def main():
    print("Generating thermal conductivity tables and visualizations...")
    
    # Get thermal conductivity and phase data
    tc_values, phases, densities = get_thermal_conductivity_grid()
    
    # Create tables and plots
    tc_df, phases_df, density_df = create_tables_and_plots(tc_values, phases, densities)
    
    print(f"\nAll results have been saved to the '{OUTPUT_DIR}' directory.")
    print(f"- thermal_conductivity_table.csv: Raw thermal conductivity values")
    print(f"- thermal_conductivity_table.txt: Formatted table")
    print(f"- phases_table.csv: Phase at each P-T point")
    print(f"- density_table.csv: Density at each P-T point")
    print(f"- thermal_conductivity_heatmap.png: Visualization of thermal conductivity")
    print(f"- phase_diagram.png: Visualization of phases")
    print(f"- density_heatmap.png: Visualization of density")

if __name__ == "__main__":
    main()