import requests
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter
import json
import time
import sys

# API endpoint
api_url = "http://localhost:5051/calculate"  # Modify this if your API is hosted elsewhere

# Configure the request payload for CO2
payload = {
    "composition": [
    {"fluid": "CO2", "fraction": 0.9000},
    {"fluid": "NITROGEN", "fraction": 0.1}

    ],
    "pressure_range": {
        "from": 10,   # Start at 10 bar
        "to": 300     # Up to 300 bar
    },
    "temperature_range": {
        "from": -60,  # Start at -60°C
        "to": 150     # Up to 150°C
    },
    "pressure_resolution": 2,   # Increased resolution for smoother lines
    "temperature_resolution": 5, # Temperature step in °C
    "properties": [
        "enthalpy",
        "vapor_fraction",
        "phase",
        "critical_temperature",
        "critical_pressure",
        "critical_density"  # Added critical density to compare all critical properties
    ],
    "units_system": "SI"  # Use SI units
}

# Function to send request to API with retry
def get_fluid_data(max_retries=3, retry_delay=2):
    for attempt in range(max_retries):
        try:
            print(f"API request attempt {attempt+1}/{max_retries}...")
            response = requests.post(api_url, json=payload, timeout=30)
            response.raise_for_status()  # Raise exception for HTTP errors
            return response.json()["results"]
        except requests.exceptions.RequestException as e:
            print(f"API request error: {e}")
            if attempt < max_retries - 1:
                print(f"Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
            else:
                print("Max retries reached. Using sample data instead.")
                return generate_sample_data()

# Function to generate sample data for testing when API is unavailable
def generate_sample_data():
    print("Generating sample CO2 data for demonstration...")
    
    # CO2 critical properties
    crit_temp_c = 31.1  # Critical temperature in °C
    crit_press_bar = 73.9  # Critical pressure in bar
    
    # Generate a grid of temperature and pressure points
    temps = np.arange(-60, 151, 5)  # °C
    pressures = np.logspace(1, np.log10(300), 30)  # 10 to 300 bar
    
    # Molecular weight of CO2 for unit conversion
    mw_co2 = 44.01  # g/mol
    
    results = []
    idx = 0
    
    for temp in temps:
        for press in pressures:
            # Determine phase based on critical point
            if temp > crit_temp_c and press > crit_press_bar:
                phase = "Supercritical"
                vf = 999  # Conventional value for supercritical state
            elif temp > crit_temp_c:
                phase = "Vapor"  # Supercritical temperature but subcritical pressure
                vf = 1
            else:
                # Simplified phase determination based on CO2 properties
                # This is a very rough approximation
                sat_press = 10 * np.exp(0.05 * (temp + 60))  # Approximated saturation curve
                
                if press < sat_press:
                    phase = "Vapor"
                    vf = 1
                elif press > sat_press * 1.1:
                    phase = "Liquid"
                    vf = 0
                else:
                    phase = "Two-Phase"
                    # Calculate vapor fraction based on where we are in the two-phase region
                    # This is just an approximation for visualization
                    vf = max(0, min(1, 1.1 - press/sat_press))
            
            # Generate sample enthalpy value (very approximate)
            # For CO2, enthalpy increases with temperature and decreases with pressure
            base_enthalpy = 150 + 2 * (temp + 60)  # Base component from temperature
            pressure_effect = -0.1 * np.log(press/10)  # Pressure component
            
            # Add phase-specific adjustments
            if phase == "Vapor":
                enthalpy_adjustment = 50
            elif phase == "Two-Phase":
                enthalpy_adjustment = 20 * vf
            else:
                enthalpy_adjustment = 0
                
            enthalpy = base_enthalpy + pressure_effect + enthalpy_adjustment
            
            # Create a sample data point
            results.append({
                "index": idx,
                "temperature": {"value": temp, "unit": "°C"},
                "pressure": {"value": press, "unit": "bar"},
                "enthalpy": {"value": enthalpy * mw_co2 / 1000, "unit": "J/mol"},  # Convert from kJ/kg to J/mol
                "vapor_fraction": {"value": vf, "unit": "dimensionless"},
                "phase": {"value": phase, "unit": None},
                "critical_temperature": {"value": crit_temp_c + 273.15, "unit": "K"},
                "critical_pressure": {"value": crit_press_bar * 100, "unit": "kPa"}  # Convert bar to kPa
            })
            idx += 1
    
    return results

# Get data from API or generate sample data
print("Fetching data from REFPROP API...")
data = get_fluid_data()

if not data:
    print("Failed to retrieve or generate data")
    sys.exit(1)

# Extract data for plotting
print("Processing data for plotting...")
pressures = []
enthalpies = []
vapor_fractions = []
phases = []

# Convert J/mol to kJ/kg for CO2 (molecular weight = 44.01 g/mol)
mw_co2 = 44.01  # g/mol
conversion_factor = 1000 / mw_co2  # (1000 g/kg) / (g/mol) = 1000/mw_co2

for point in data:
    # Extract values from the returned data structure
    pressure = point["pressure"]["value"]
    enthalpy_jmol = point["enthalpy"]["value"]
    
    # Skip any points with invalid enthalpy values (if any)
    if enthalpy_jmol is None or np.isnan(enthalpy_jmol):
        continue
    
    # Convert enthalpy from J/mol to kJ/kg
    enthalpy_kjkg = enthalpy_jmol * conversion_factor / 1000  # J/mol * (mol/g) * (1000g/kg) * (1kJ/1000J)
        
    vapor_fraction = point["vapor_fraction"]["value"]
    phase = point["phase"]["value"]
    
    pressures.append(pressure)
    enthalpies.append(enthalpy_kjkg)
    vapor_fractions.append(vapor_fraction)
    phases.append(phase)

# Convert to numpy arrays for easier manipulation
pressures = np.array(pressures)
enthalpies = np.array(enthalpies)
vapor_fractions = np.array(vapor_fractions)
phases = np.array(phases)

# Check if we have data
if len(pressures) == 0:
    print("No valid data points available for plotting")
    sys.exit(1)

print(f"Plotting {len(pressures)} data points...")

# Set up the figure with appropriate size and DPI
plt.figure(figsize=(12, 8), dpi=100)

# Create different datasets for different phases
liquid_mask = (vapor_fractions == 0)
vapor_mask = (vapor_fractions == 1)
twophase_mask = (vapor_fractions > 0) & (vapor_fractions < 1) 
supercritical_mask = np.array([p == "Supercritical" for p in phases])

# Create main plot with logarithmic y-axis
plt.semilogy()

# Plot isotherms (constant temperature lines)
# Group by temperature and sort by pressure
unique_temps = set()
for point in data:
    unique_temps.add(point["temperature"]["value"])

print(f"Plotting {len(unique_temps)} isotherms...")
for temp in sorted(unique_temps):
    temp_points = [(p["pressure"]["value"], 
                   p["enthalpy"]["value"] * conversion_factor / 1000) 
                   for p in data if p["temperature"]["value"] == temp]
    if temp_points:
        temp_points.sort(key=lambda x: x[0])  # Sort by pressure
        temp_pressures, temp_enthalpies = zip(*temp_points)
        
        # Use thicker lines for reference temperatures
        if temp % 20 == 0:  # Every 20°C
            plt.plot(temp_enthalpies, temp_pressures, 'k-', linewidth=0.8, alpha=0.7)
            # Add temperature label
            midpoint = len(temp_pressures) // 2
            if midpoint > 0:
                plt.text(temp_enthalpies[midpoint], temp_pressures[midpoint], 
                         f"{temp}°C", fontsize=8, ha='center', va='bottom')
        else:
            plt.plot(temp_enthalpies, temp_pressures, 'k-', linewidth=0.4, alpha=0.5)

# Find and highlight the saturation curve (vapor fraction = 0 or 1 boundary)
# This requires identifying the phase envelope
print("Highlighting saturation curve...")
sat_liquid_points = []
sat_vapor_points = []

# Improve the detection of saturation points for a cleaner boundary
for i in range(len(vapor_fractions)):
    # Find points at the liquid boundary (where vapor fraction transitions from 0)
    if vapor_fractions[i] == 0:
        # Check nearby points for transition to two-phase
        nearby_points = range(max(0, i-10), min(len(vapor_fractions), i+10))
        if any(0 < vapor_fractions[j] < 1 for j in nearby_points):
            sat_liquid_points.append((enthalpies[i], pressures[i]))
    
    # Find points at the vapor boundary (where vapor fraction transitions to 1)
    elif vapor_fractions[i] == 1:
        # Check nearby points for transition from two-phase
        nearby_points = range(max(0, i-10), min(len(vapor_fractions), i+10))
        if any(0 < vapor_fractions[j] < 1 for j in nearby_points):
            sat_vapor_points.append((enthalpies[i], pressures[i]))
            
    # Alternatively, can directly use two-phase points adjacent to single-phase
    elif 0 < vapor_fractions[i] < 1:
        # Find edges of two-phase region where it meets liquid or vapor
        nearby_points = range(max(0, i-5), min(len(vapor_fractions), i+5))
        if any(vapor_fractions[j] == 0 for j in nearby_points):
            # This is near the liquid boundary
            sat_liquid_points.append((enthalpies[i], pressures[i]))
        if any(vapor_fractions[j] == 1 for j in nearby_points):
            # This is near the vapor boundary
            sat_vapor_points.append((enthalpies[i], pressures[i]))

# Sort saturation points by pressure for proper curve plotting
sat_liquid_points.sort(key=lambda x: x[1])
sat_vapor_points.sort(key=lambda x: x[1])

# Extract sorted coordinates
if sat_liquid_points:
    sat_liquid_h, sat_liquid_p = zip(*sat_liquid_points)
    # Add a thicker, more visible line for the saturated liquid boundary
    plt.plot(sat_liquid_h, sat_liquid_p, 'r-', linewidth=3.5)
    # Add a white "glow" effect to make it stand out
    plt.plot(sat_liquid_h, sat_liquid_p, 'w-', linewidth=5, alpha=0.3)

if sat_vapor_points:
    sat_vapor_h, sat_vapor_p = zip(*sat_vapor_points)
    # Add a thicker, more visible line for the saturated vapor boundary
    plt.plot(sat_vapor_h, sat_vapor_p, 'r-', linewidth=3.5)
    # Add a white "glow" effect to make it stand out
    plt.plot(sat_vapor_h, sat_vapor_p, 'w-', linewidth=5, alpha=0.3)

# Set labels and title
plt.xlabel('Enthalpy (kJ/kg)')
plt.ylabel('Pressure (bar)')
plt.title('Pressure-Enthalpy Diagram for CO2')

# Set reasonable axis limits based on data
enthalpy_min = max(min(enthalpies), 100)  # Don't go below 100 kJ/kg
enthalpy_max = min(max(enthalpies), 600)  # Don't go above 600 kJ/kg
plt.xlim(enthalpy_min, enthalpy_max)
plt.ylim(10, 300)  # Pressure from 10 to 300 bar

# Add grid
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

# Add annotations for phases
plt.text(enthalpy_min + (enthalpy_max - enthalpy_min) * 0.2, 50, 'LIQUID', fontsize=12)
plt.text(enthalpy_min + (enthalpy_max - enthalpy_min) * 0.7, 20, 'VAPOR', fontsize=12)
plt.text(enthalpy_min + (enthalpy_max - enthalpy_min) * 0.7, 200, 'SUPERCRITICAL', fontsize=12)
plt.text(enthalpy_min + (enthalpy_max - enthalpy_min) * 0.5, 30, 'TWO-PHASE', fontsize=12)

# Extract and verify the critical properties directly from REFPROP
print("Analyzing critical properties from REFPROP...")
critical_properties = {}

# Extract the critical values - should be the same for all data points
for point in data:
    if "critical_temperature" in point and "critical_pressure" in point and "critical_density" in point:
        crit_T = point["critical_temperature"]["value"]
        crit_p = point["critical_pressure"]["value"]
        crit_d = point["critical_density"]["value"]
        
        # Convert units to standard values
        crit_T_C = crit_T - 273.15  # K to °C
        crit_p_bar = crit_p / 100    # kPa to bar
        
        # Store for comparison
        critical_properties = {
            "temperature": crit_T_C,
            "pressure": crit_p_bar,
            "density": crit_d
        }
        break

# Compare with literature values for verification
literature_values = {
    "temperature": 31.1,  # °C
    "pressure": 73.9,     # bar
    "density": 10.6       # mol/L (approximate)
}

if critical_properties:
    print("\nCritical properties comparison between REFPROP and literature values:")
    print(f"Property    | REFPROP      | Literature   | Difference (%)")
    print(f"------------|--------------|--------------|---------------")
    
    for prop, lit_value in literature_values.items():
        if prop in critical_properties:
            refprop_value = critical_properties[prop]
            diff_percent = abs(refprop_value - lit_value) / lit_value * 100
            print(f"{prop.capitalize():12}| {refprop_value:.4f} | {lit_value:.4f} | {diff_percent:.4f}%")
    
    # Find and mark the critical point on the diagram
    crit_T_C = critical_properties["temperature"]
    crit_p_bar = critical_properties["pressure"]
    
    # Need to find the enthalpy at or near the critical point
    critical_points = []
    for p in data:
        if (abs(p["temperature"]["value"] - crit_T_C) < 5 and 
            abs(p["pressure"]["value"] - crit_p_bar) < 5):
            # Convert enthalpy to kJ/kg
            h_kjkg = p["enthalpy"]["value"] * conversion_factor / 1000
            critical_points.append((h_kjkg, p["pressure"]["value"]))
    
    if critical_points:
        # Use the average if multiple points near the critical state
        crit_h = sum(h for h, _ in critical_points) / len(critical_points)
        crit_p_bar_avg = sum(p for _, p in critical_points) / len(critical_points)
        
        # Plot the critical point
        plt.plot(crit_h, crit_p_bar_avg, 'ro', markersize=8)
        
        # Add a label with critical properties
        plt.text(crit_h, crit_p_bar_avg * 1.1, 
                f'Critical Point (REFPROP)\n{crit_T_C:.2f}°C, {crit_p_bar:.2f} bar', 
                fontsize=10, ha='center', va='bottom', 
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='red', boxstyle='round,pad=0.5'))
else:
    print("Warning: Critical properties not found in the data.")
    # Plot a fallback critical point using literature values
    critical_points = []
    for p in data:
        if (abs(p["temperature"]["value"] - literature_values["temperature"]) < 5 and 
            abs(p["pressure"]["value"] - literature_values["pressure"]) < 5):
            h_kjkg = p["enthalpy"]["value"] * conversion_factor / 1000
            critical_points.append((h_kjkg, p["pressure"]["value"]))
    
    if critical_points:
        crit_h = sum(h for h, _ in critical_points) / len(critical_points)
        crit_p_bar_avg = sum(p for _, p in critical_points) / len(critical_points)
        
        plt.plot(crit_h, crit_p_bar_avg, 'ro', markersize=8)
        plt.text(crit_h, crit_p_bar_avg * 1.1, 
                f'Critical Point (Literature)\n{literature_values["temperature"]:.2f}°C, {literature_values["pressure"]:.2f} bar', 
                fontsize=10, ha='center', va='bottom', 
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='red', boxstyle='round,pad=0.5'))

# Add a custom legend with the phase envelope highlighted
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='red', lw=3.5)]
plt.legend(custom_lines, ['Phase Envelope'], loc='lower right')

# Show the plot
plt.tight_layout()
plt.savefig('co2_ph_diagram_kjkg.png', dpi=300)
plt.show()

print("Phase diagram generation complete.")