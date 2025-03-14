import requests
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import time
import sys
import json

# API endpoint
api_url = "http://localhost:5051/ph_flash"

# Configure the request payload for CO2
payload = {
    "composition": [
        {"fluid": "CO2", "fraction": 1.0}
    ],
    "pressure_range": {
        "from": 10,      # Start at 10 bar
        "to": 200        # Up to 200 bar
    },
    "enthalpy_range": {
        "from": 5000,    # Lower range for liquid region
        "to": 40000      # Upper range for vapor region
    },
    "pressure_resolution": 5,      # Pressure step in bar
    "enthalpy_resolution": 2000,   # Enthalpy step in J/mol
    "properties": [
        "temperature",
        "density", 
        "vapor_fraction",
        "phase",
        "entropy",
        "sound_speed",
        "critical_temperature",
        "critical_pressure"
    ],
    "units_system": "SI"
}

# Function to send request to API with retry
def get_fluid_data(max_retries=3, retry_delay=2):
    for attempt in range(max_retries):
        try:
            print(f"API request attempt {attempt+1}/{max_retries}...")
            response = requests.post(api_url, json=payload, timeout=60)
            response.raise_for_status()  # Raise exception for HTTP errors
            return response.json()["results"]
        except requests.exceptions.RequestException as e:
            print(f"API request error: {e}")
            if attempt < max_retries - 1:
                print(f"Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
            else:
                print("Max retries reached. Exiting.")
                sys.exit(1)

# Get data from API
print("Fetching data from REFPROP API using PH flash...")
data = get_fluid_data()

if not data:
    print("Failed to retrieve data")
    sys.exit(1)

print(f"Processing {len(data)} data points for plotting...")

# Extract data for plotting
pressures = []
enthalpies = []
temperatures = []
vapor_fractions = []
phases = []

for point in data:
    pressure = point["pressure"]["value"]
    enthalpy = point["enthalpy"]["value"]
    temperature = point["temperature"]["value"]
    vapor_fraction = point["vapor_fraction"]["value"]
    phase = point["phase"]["value"]
    
    pressures.append(pressure)
    enthalpies.append(enthalpy)
    temperatures.append(temperature)
    vapor_fractions.append(vapor_fraction)
    phases.append(phase)

# Convert to numpy arrays for easier manipulation
pressures = np.array(pressures)
enthalpies = np.array(enthalpies)
temperatures = np.array(temperatures)
vapor_fractions = np.array(vapor_fractions)
phases = np.array(phases)

# Create figure with a specific axes - important for colorbar
fig, ax = plt.subplots(figsize=(14, 10), dpi=100)

# Create different datasets for different phases
liquid_mask = (vapor_fractions == 0)
vapor_mask = (vapor_fractions == 1)
twophase_mask = (vapor_fractions > 0) & (vapor_fractions < 1)

# Create main plot with logarithmic y-axis for pressure
ax.set_yscale('log')

# Plot the data points categorized by phase
if np.any(liquid_mask):
    sc_liquid = ax.scatter(enthalpies[liquid_mask], pressures[liquid_mask], 
                           c=temperatures[liquid_mask], 
                           cmap='Blues', marker='o', alpha=0.6, label='Liquid')

if np.any(vapor_mask):
    sc_vapor = ax.scatter(enthalpies[vapor_mask], pressures[vapor_mask], 
                          c=temperatures[vapor_mask], 
                          cmap='Reds', marker='s', alpha=0.6, label='Vapor')

if np.any(twophase_mask):
    sc_twophase = ax.scatter(enthalpies[twophase_mask], pressures[twophase_mask], 
                             c=temperatures[twophase_mask], 
                             cmap='Greens', marker='*', alpha=0.8, label='Two-Phase')

# Identify temperature contours (isotherms)
unique_temps = sorted(set([round(t / 10) * 10 for t in temperatures]))  # Round to nearest 10°C

for temp in unique_temps:
    # Find points with this temperature
    temp_mask = np.abs(temperatures - temp) < 5  # Within 5°C of the target temperature
    
    if np.sum(temp_mask) > 2:  # Only plot if enough points
        points = np.column_stack((enthalpies[temp_mask], pressures[temp_mask]))
        # Sort points by enthalpy to get a smoother curve
        points = points[points[:, 0].argsort()]
        
        if temp % 20 == 0:  # Every 20°C, highlight the isotherm
            ax.plot(points[:, 0], points[:, 1], 'k-', linewidth=1.5, alpha=0.8)
            # Label the isotherm at the midpoint
            mid_idx = len(points) // 2
            if mid_idx > 0:
                ax.text(points[mid_idx, 0], points[mid_idx, 1], f"{temp}°C", 
                         fontsize=9, color='black', ha='center', va='bottom')
        else:
            ax.plot(points[:, 0], points[:, 1], 'k--', linewidth=0.7, alpha=0.5)

# Identify saturation curve (boundary of two-phase region)
def find_phase_boundary_points():
    boundary_points = []
    
    # Identify all two-phase points
    tp_points = np.column_stack((enthalpies[twophase_mask], pressures[twophase_mask]))
    
    # Find points that are adjacent to single-phase regions
    for h in np.unique(enthalpies):
        # Get all points at this enthalpy
        h_mask = np.abs(enthalpies - h) < 1e-3
        p_values = pressures[h_mask]
        vf_values = vapor_fractions[h_mask]
        
        if len(p_values) > 1:
            # Check if there are both two-phase and single-phase points
            has_tp = np.any((vf_values > 0) & (vf_values < 1))
            has_liquid = np.any(vf_values == 0)
            has_vapor = np.any(vf_values == 1)
            
            if has_tp and (has_liquid or has_vapor):
                # Find transition points
                if has_liquid:
                    liq_p = np.max(p_values[vf_values == 0])
                    tp_p = np.min(p_values[(vf_values > 0) & (vf_values < 1)])
                    # Take average as boundary
                    boundary_points.append((h, (liq_p + tp_p) / 2))
                
                if has_vapor:
                    vap_p = np.min(p_values[vf_values == 1])
                    tp_p = np.max(p_values[(vf_values > 0) & (vf_values < 1)])
                    # Take average as boundary
                    boundary_points.append((h, (vap_p + tp_p) / 2))
    
    return np.array(boundary_points) if boundary_points else np.array([])

# Find and plot phase boundary
boundary_points = find_phase_boundary_points()
if len(boundary_points) > 0:
    # Sort by enthalpy
    boundary_points = boundary_points[boundary_points[:, 0].argsort()]
    ax.plot(boundary_points[:, 0], boundary_points[:, 1], 'r-', linewidth=3, alpha=0.9, label='Saturation Curve')

# Extract critical point info - check if available
critical_temp = None
critical_pressure = None

# Extract from the first data point that has the critical properties
for point in data:
    if "critical_temperature" in point and "critical_pressure" in point:
        if point["critical_temperature"]["value"] is not None and point["critical_pressure"]["value"] is not None:
            critical_temp = point["critical_temperature"]["value"] - 273.15  # K to °C
            critical_pressure = point["critical_pressure"]["value"] / 100  # kPa to bar
            break

if critical_temp is None or critical_pressure is None:
    print("Warning: Critical point data not available")

# Find a point close to the critical point to mark on the diagram (if critical data available)
if critical_temp is not None and critical_pressure is not None:
    closest_points = []
    for i, (t, p) in enumerate(zip(temperatures, pressures)):
        if abs(t - critical_temp) < 10 and abs(p - critical_pressure) < 10:
            closest_points.append((enthalpies[i], p))

    if closest_points:
        # Average the closest points
        crit_h = sum(p[0] for p in closest_points) / len(closest_points)
        crit_p = sum(p[1] for p in closest_points) / len(closest_points)
        
        # Plot critical point
        ax.plot(crit_h, crit_p, 'ko', markersize=10, label=f'Critical Point: {critical_temp:.1f}°C, {critical_pressure:.1f} bar')

# Set labels and title
ax.set_xlabel('Enthalpy (J/mol)')
ax.set_ylabel('Pressure (bar)')
ax.set_title('Pressure-Enthalpy Diagram for CO2 Generated with REFPROP')

# Add grid with minor gridlines
ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

# Create a proper mappable object for the colorbar
norm = plt.Normalize(min(temperatures), max(temperatures))
sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=norm)
sm.set_array([])  # You need to set an array to avoid warnings, empty is fine

# Add colorbar for temperature reference with explicit axes reference
cbar = fig.colorbar(sm, ax=ax)
cbar.set_label('Temperature (°C)')

# Add legend with best position
ax.legend(loc='best')

# Set reasonable axis limits
ax.set_ylim(payload["pressure_range"]["from"], payload["pressure_range"]["to"])
ax.set_xlim(payload["enthalpy_range"]["from"], payload["enthalpy_range"]["to"])

# Add annotations for phases
ax.text(payload["enthalpy_range"]["from"]*1.1, payload["pressure_range"]["to"]*0.7, 'LIQUID', 
         fontsize=14, color='blue', fontweight='bold')
ax.text(payload["enthalpy_range"]["to"]*0.8, payload["pressure_range"]["from"]*1.5, 'VAPOR', 
         fontsize=14, color='red', fontweight='bold')
ax.text((payload["enthalpy_range"]["from"] + payload["enthalpy_range"]["to"])/2, 
         payload["pressure_range"]["from"]*3, 'TWO-PHASE', 
         fontsize=14, color='green', fontweight='bold')

plt.tight_layout()
plt.savefig('co2_ph_diagram.png', dpi=300)
plt.show()

print("Pressure-enthalpy diagram generated successfully.")