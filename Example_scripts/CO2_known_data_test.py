import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Directly define the reference data with correct values
# This bypasses the parsing issues with the original table format
reference_data = [
    {'Temperature': -20, 'Density': 2.1371, 'Viscosity': 1.27E-05, 'SpecificHeat': 0.8080, 'ThermalConductivity': 0.013178, 'PrandtlNumber': 0.780099},
    {'Temperature': -15, 'Density': 2.0946, 'Viscosity': 1.30E-05, 'SpecificHeat': 0.8126, 'ThermalConductivity': 0.013544, 'PrandtlNumber': 0.778229},
    {'Temperature': -10, 'Density': 2.0538, 'Viscosity': 1.32E-05, 'SpecificHeat': 0.8173, 'ThermalConductivity': 0.013916, 'PrandtlNumber': 0.776295},
    {'Temperature': -5, 'Density': 2.0145, 'Viscosity': 1.35E-05, 'SpecificHeat': 0.8220, 'ThermalConductivity': 0.014293, 'PrandtlNumber': 0.77439},
    {'Temperature': 0, 'Density': 1.9768, 'Viscosity': 1.37E-05, 'SpecificHeat': 0.8268, 'ThermalConductivity': 0.014674, 'PrandtlNumber': 0.772587},
    {'Temperature': 5, 'Density': 1.9405, 'Viscosity': 1.40E-05, 'SpecificHeat': 0.8316, 'ThermalConductivity': 0.01506, 'PrandtlNumber': 0.770675},
    {'Temperature': 10, 'Density': 1.9055, 'Viscosity': 1.42E-05, 'SpecificHeat': 0.8364, 'ThermalConductivity': 0.01545, 'PrandtlNumber': 0.768821},
    {'Temperature': 15, 'Density': 1.8718, 'Viscosity': 1.44E-05, 'SpecificHeat': 0.8412, 'ThermalConductivity': 0.015844, 'PrandtlNumber': 0.767022},
    {'Temperature': 20, 'Density': 1.8393, 'Viscosity': 1.47E-05, 'SpecificHeat': 0.8460, 'ThermalConductivity': 0.016242, 'PrandtlNumber': 0.765163},
    {'Temperature': 25, 'Density': 1.8080, 'Viscosity': 1.49E-05, 'SpecificHeat': 0.8508, 'ThermalConductivity': 0.016643, 'PrandtlNumber': 0.763378},
    {'Temperature': 30, 'Density': 1.7777, 'Viscosity': 1.52E-05, 'SpecificHeat': 0.8556, 'ThermalConductivity': 0.017046, 'PrandtlNumber': 0.761664},
    {'Temperature': 35, 'Density': 1.7484, 'Viscosity': 1.54E-05, 'SpecificHeat': 0.8604, 'ThermalConductivity': 0.017452, 'PrandtlNumber': 0.759974},
    {'Temperature': 40, 'Density': 1.7201, 'Viscosity': 1.57E-05, 'SpecificHeat': 0.8651, 'ThermalConductivity': 0.017861, 'PrandtlNumber': 0.758278},
    {'Temperature': 45, 'Density': 1.6926, 'Viscosity': 1.59E-05, 'SpecificHeat': 0.8698, 'ThermalConductivity': 0.018272, 'PrandtlNumber': 0.756691},
    {'Temperature': 50, 'Density': 1.6661, 'Viscosity': 1.61E-05, 'SpecificHeat': 0.8745, 'ThermalConductivity': 0.018685, 'PrandtlNumber': 0.755133},
    {'Temperature': 55, 'Density': 1.6404, 'Viscosity': 1.64E-05, 'SpecificHeat': 0.8791, 'ThermalConductivity': 0.019099, 'PrandtlNumber': 0.753649},
    {'Temperature': 60, 'Density': 1.6155, 'Viscosity': 1.66E-05, 'SpecificHeat': 0.8838, 'ThermalConductivity': 0.019515, 'PrandtlNumber': 0.752192},
    {'Temperature': 65, 'Density': 1.5913, 'Viscosity': 1.68E-05, 'SpecificHeat': 0.8883, 'ThermalConductivity': 0.019932, 'PrandtlNumber': 0.750791},
    {'Temperature': 70, 'Density': 1.5679, 'Viscosity': 1.71E-05, 'SpecificHeat': 0.8929, 'ThermalConductivity': 0.02035, 'PrandtlNumber': 0.749482},
    {'Temperature': 75, 'Density': 1.5452, 'Viscosity': 1.73E-05, 'SpecificHeat': 0.8974, 'ThermalConductivity': 0.020769, 'PrandtlNumber': 0.748182},
    {'Temperature': 80, 'Density': 1.5231, 'Viscosity': 1.75E-05, 'SpecificHeat': 0.9018, 'ThermalConductivity': 0.02119, 'PrandtlNumber': 0.746922},
    {'Temperature': 85, 'Density': 1.5016, 'Viscosity': 1.78E-05, 'SpecificHeat': 0.9063, 'ThermalConductivity': 0.02161, 'PrandtlNumber': 0.745774},
    {'Temperature': 90, 'Density': 1.4807, 'Viscosity': 1.80E-05, 'SpecificHeat': 0.9107, 'ThermalConductivity': 0.022032, 'PrandtlNumber': 0.744623},
    {'Temperature': 95, 'Density': 1.4604, 'Viscosity': 1.82E-05, 'SpecificHeat': 0.9150, 'ThermalConductivity': 0.022453, 'PrandtlNumber': 0.743557},
    {'Temperature': 100, 'Density': 1.4407, 'Viscosity': 1.85E-05, 'SpecificHeat': 0.9193, 'ThermalConductivity': 0.022875, 'PrandtlNumber': 0.742521},
    {'Temperature': 105, 'Density': 1.4215, 'Viscosity': 1.87E-05, 'SpecificHeat': 0.9236, 'ThermalConductivity': 0.023298, 'PrandtlNumber': 0.741497},
    {'Temperature': 110, 'Density': 1.4028, 'Viscosity': 1.89E-05, 'SpecificHeat': 0.9278, 'ThermalConductivity': 0.02372, 'PrandtlNumber': 0.74059},
    {'Temperature': 115, 'Density': 1.3846, 'Viscosity': 1.92E-05, 'SpecificHeat': 0.9320, 'ThermalConductivity': 0.024142, 'PrandtlNumber': 0.739686},
    {'Temperature': 120, 'Density': 1.3669, 'Viscosity': 1.94E-05, 'SpecificHeat': 0.9361, 'ThermalConductivity': 0.024565, 'PrandtlNumber': 0.738829}
]

# Create dataframe
ref_df = pd.DataFrame(reference_data)
print(f"Successfully loaded {len(ref_df)} reference data points")

# API endpoint
API_URL = "http://127.0.0.1:5051/pt_flash"

# Standard atmospheric pressure in bar (1 atm ≈ 1.01325 bar)
PRESSURE_BAR = 1.01325

# Function to convert mol/L to kg/m³ for CO2
def convert_molL_to_kgm3(density_molL):
    """Convert density from mol/L to kg/m³ for CO2"""
    molar_mass_CO2 = 44.01  # g/mol
    return density_molL * molar_mass_CO2  # mol/L * g/mol = g/L = kg/m³

# Create a dataframe to store the API results
api_results = []

# Make API calls for each temperature point
for temp in ref_df['Temperature']:
    # Replace this payload structure
    payload = {
        "composition": [
            {"fluid": "CO2", "fraction": 1.0}
        ],
        "variables": {
            "pressure": {
                "range": {"from": PRESSURE_BAR, "to": PRESSURE_BAR},
                "resolution": 1
            },
            "temperature": {
                "range": {"from": temp, "to": temp},
                "resolution": 1
            }
        },
        "calculation": {
            "properties": [
                "density", 
                "viscosity", 
                "cp", 
                "thermal_conductivity", 
                "prandtl_number",
                "phase"
            ],
            "units_system": "SI"
        }
    }
    
    try:
        print(f"Making API call for temperature {temp}°C...")
        response = requests.post(API_URL, json=payload)
        if response.status_code == 200:
            data = response.json()['results'][0]
            
            # Extract properties and convert units if needed
            density_molL = data['density']['value']
            density_kgm3 = convert_molL_to_kgm3(density_molL)
            
            viscosity = data['viscosity']['value'] * 1e-6  # Convert μPa·s to Pa·s
            specific_heat = data['cp']['value'] / 44.01  # Convert J/(mol·K) to kJ/(kg·K) for CO2
            thermal_conductivity = data['thermal_conductivity']['value']
            prandtl_number = data['prandtl_number']['value']
            phase = data['phase']['value']
            
            api_results.append({
                'Temperature': temp,
                'Density': density_kgm3,
                'Viscosity': viscosity,
                'SpecificHeat': specific_heat,
                'ThermalConductivity': thermal_conductivity,
                'PrandtlNumber': prandtl_number,
                'Phase': phase
            })
            print(f"Success: T={temp}°C, Phase={phase}, Density={density_kgm3:.4f} kg/m³")
        else:
            print(f"API error for temperature {temp}°C: {response.text}")
    except Exception as e:
        print(f"Exception for temperature {temp}°C: {str(e)}")

# Convert to DataFrame
api_df = pd.DataFrame(api_results)

if len(api_df) == 0:
    print("No valid API results were obtained. Check if the API is running correctly.")
    exit(1)

# Calculate percentage differences
comparison_df = pd.merge(ref_df, api_df, on='Temperature', suffixes=('_ref', '_api'))
for prop in ['Density', 'Viscosity', 'SpecificHeat', 'ThermalConductivity', 'PrandtlNumber']:
    comparison_df[f'{prop}_diff_pct'] = (comparison_df[f'{prop}_api'] - comparison_df[f'{prop}_ref']) / comparison_df[f'{prop}_ref'] * 100

# Print comparison results
print("\nComparison Results (% difference):\n")
summary_columns = ['Temperature']
for prop in ['Density', 'Viscosity', 'SpecificHeat', 'ThermalConductivity', 'PrandtlNumber']:
    if f'{prop}_diff_pct' in comparison_df.columns:
        summary_columns.append(f'{prop}_diff_pct')

if 'Phase_api' in comparison_df.columns:
    summary_columns.append('Phase_api')

summary_df = comparison_df[summary_columns].round(2)
print(summary_df)

# Calculate average absolute difference for each property
print("\nAverage Absolute Percentage Differences:")
for prop in ['Density', 'Viscosity', 'SpecificHeat', 'ThermalConductivity', 'PrandtlNumber']:
    if f'{prop}_diff_pct' in comparison_df.columns:
        avg_diff = comparison_df[f'{prop}_diff_pct'].abs().mean()
        print(f"{prop}: {avg_diff:.2f}%")

# Plot the comparison
fig, axs = plt.subplots(3, 2, figsize=(15, 15))
fig.suptitle('CO2 Properties: Reference vs API Values (1 atm)', fontsize=16)

properties = [
    ('Density', 'kg/m³'), 
    ('Viscosity', 'Pa·s'), 
    ('SpecificHeat', 'kJ/(kg·K)'),
    ('ThermalConductivity', 'W/(m·K)'), 
    ('PrandtlNumber', '')
]

for i, (prop, unit) in enumerate(properties):
    row, col = i // 2, i % 2
    ax = axs[row, col]
    
    if all(f'{prop}_{suffix}' in comparison_df.columns for suffix in ['ref', 'api']):
        ax.plot(comparison_df['Temperature'], comparison_df[f'{prop}_ref'], 'b-', label='Reference')
        ax.plot(comparison_df['Temperature'], comparison_df[f'{prop}_api'], 'r--', label='API')
        
        ax.set_xlabel('Temperature (°C)')
        ax.set_ylabel(f'{prop} ({unit})')
        ax.set_title(f'{prop} vs Temperature')
        ax.grid(True)
        ax.legend()
    else:
        ax.text(0.5, 0.5, f"Data for {prop} not available", 
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes)
        ax.set_title(f"{prop} - Missing Data")

# Plot percentage differences
ax = axs[2, 1]
diff_plotted = False
for prop in ['Density', 'Viscosity', 'SpecificHeat', 'ThermalConductivity', 'PrandtlNumber']:
    if f'{prop}_diff_pct' in comparison_df.columns:
        ax.plot(comparison_df['Temperature'], comparison_df[f'{prop}_diff_pct'], label=prop)
        diff_plotted = True

if diff_plotted:
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Difference (%)')
    ax.set_title('Percentage Difference vs Temperature')
    ax.grid(True)
    ax.legend()
else:
    ax.text(0.5, 0.5, "No difference data available", 
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
    ax.set_title("Differences - Missing Data")

plt.tight_layout(rect=[0, 0, 1, 0.97])  # Adjust layout to make room for the title
plt.savefig('co2_property_comparison.png')
plt.close()

# Export results to CSV
comparison_df.to_csv('co2_property_comparison.csv', index=False)

print("\nResults saved to 'co2_property_comparison.csv'")
print("Plots saved to 'co2_property_comparison.png'")