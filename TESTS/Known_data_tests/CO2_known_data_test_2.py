import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Create dataframe from reference data - Saturated CO2 Properties (Table D1)
# Temperature in °C, Pressure in MPa, specific volumes in m³/kg, energies in kJ/kg, entropy in kJ/(kg·K)
reference_data = [
    # T, P, vf, vg, uf, ug, hf, hg, sf, sg
    [-50, 0.68234, 0.000866, 0.055789, 92.35, 394.61, 92.94, 432.68, 0.5794, 2.1018],
    [-48, 0.73949, 0.000872, 0.051618, 96.26, 395.12, 96.90, 433.29, 0.5968, 2.0909],
    [-46, 0.80015, 0.000878, 0.047819, 100.18, 395.60, 100.88, 433.86, 0.6142, 2.0801],
    [-44, 0.86445, 0.000883, 0.044352, 104.11, 396.05, 104.87, 434.39, 0.6314, 2.0694],
    [-42, 0.93252, 0.000889, 0.041184, 108.05, 396.47, 108.88, 434.88, 0.6486, 2.0589],
    [-40, 1.00450, 0.000896, 0.038284, 112.00, 396.87, 112.90, 435.32, 0.6656, 2.0485],
    [-38, 1.08051, 0.000902, 0.035624, 115.97, 397.23, 116.95, 435.72, 0.6826, 2.0382],
    [-36, 1.16071, 0.000909, 0.033181, 119.96, 397.56, 121.01, 436.07, 0.6995, 2.0281],
    [-34, 1.24522, 0.000915, 0.030935, 123.96, 397.85, 125.10, 436.37, 0.7163, 2.0180],
    [-32, 1.33419, 0.000922, 0.028865, 127.97, 398.11, 129.20, 436.62, 0.7331, 2.0079],
    [-30, 1.42776, 0.000930, 0.026956, 132.01, 398.33, 133.34, 436.82, 0.7498, 1.9980],
    [-28, 1.52607, 0.000937, 0.025192, 136.07, 398.52, 137.50, 436.96, 0.7665, 1.9880],
    [-26, 1.62926, 0.000945, 0.023560, 140.15, 398.65, 141.69, 437.04, 0.7831, 1.9781],
    [-24, 1.73749, 0.000953, 0.022048, 144.25, 398.75, 145.91, 437.06, 0.7997, 1.9683],
    [-22, 1.85089, 0.000961, 0.020645, 148.38, 398.80, 150.16, 437.01, 0.8163, 1.9584],
    [-20, 1.96963, 0.000969, 0.019343, 152.54, 398.79, 154.45, 436.89, 0.8328, 1.9485],
    [-18, 2.09384, 0.000978, 0.018131, 156.73, 398.74, 158.78, 436.70, 0.8494, 1.9387],
    [-16, 2.22370, 0.000987, 0.017002, 160.95, 398.63, 163.14, 436.44, 0.8659, 1.9287],
    [-14, 2.35935, 0.000997, 0.015950, 165.20, 398.46, 167.55, 436.09, 0.8825, 1.9187],
    [-12, 2.50095, 0.001007, 0.014967, 169.49, 398.23, 172.01, 435.66, 0.8991, 1.9087],
    [-10, 2.64868, 0.001017, 0.014048, 173.83, 397.93, 176.52, 435.14, 0.9157, 1.8985],
    [-8, 2.80269, 0.001028, 0.013188, 178.20, 397.55, 181.09, 434.51, 0.9324, 1.8882],
    [-6, 2.96316, 0.001040, 0.012381, 182.63, 397.10, 185.71, 433.79, 0.9491, 1.8778],
    [-4, 3.13027, 0.001052, 0.011624, 187.11, 396.57, 190.40, 432.95, 0.9660, 1.8672],
    [-2, 3.30420, 0.001065, 0.010912, 191.64, 395.93, 195.16, 431.99, 0.9829, 1.8564],
    [0, 3.48514, 0.001078, 0.010241, 196.24, 395.20, 200.00, 430.89, 1.0000, 1.8453],
    [2, 3.67329, 0.001093, 0.009609, 200.91, 394.36, 204.93, 429.65, 1.0172, 1.8340],
    [4, 3.86884, 0.001108, 0.009011, 205.66, 393.39, 209.95, 428.25, 1.0346, 1.8223],
    [6, 4.07202, 0.001124, 0.008445, 210.50, 392.28, 215.08, 426.67, 1.0523, 1.8103],
    [8, 4.28306, 0.001142, 0.007909, 215.44, 391.02, 220.34, 424.89, 1.0702, 1.7977],
    [10, 4.50218, 0.001161, 0.007399, 220.50, 389.57, 225.73, 422.88, 1.0884, 1.7847],
    [12, 4.72966, 0.001182, 0.006912, 225.69, 387.92, 231.29, 420.62, 1.1070, 1.7710],
    [14, 4.96577, 0.001205, 0.006447, 231.05, 386.03, 237.03, 418.05, 1.1261, 1.7565],
    [16, 5.21080, 0.001231, 0.006000, 236.59, 383.86, 243.01, 415.12, 1.1458, 1.7411],
    [18, 5.46511, 0.001260, 0.005569, 242.38, 381.33, 249.26, 411.77, 1.1663, 1.7244],
    [20, 5.72905, 0.001293, 0.005149, 248.46, 378.37, 255.87, 407.87, 1.1877, 1.7062],
    [22, 6.00308, 0.001332, 0.004737, 254.93, 374.83, 262.93, 403.27, 1.2105, 1.6860],
    [24, 6.28773, 0.001379, 0.004327, 261.94, 370.50, 270.61, 397.70, 1.2352, 1.6629],
    [26, 6.58368, 0.001440, 0.003908, 269.78, 364.98, 279.26, 390.71, 1.2628, 1.6353],
    [28, 6.89182, 0.001526, 0.003459, 279.11, 357.36, 289.62, 381.20, 1.2958, 1.5999],
    [30, 7.21369, 0.001685, 0.002898, 292.40, 344.23, 304.55, 365.13, 1.3435, 1.5433],
    [30.9782, 7.37730, 0.002139, 0.002139, 316.47, 316.47, 332.25, 332.25, 1.4336, 1.4336]
]

# Create DataFrame from reference data
columns = ['Temperature', 'Pressure', 'vf', 'vg', 'uf', 'ug', 'hf', 'hg', 'sf', 'sg']
ref_df = pd.DataFrame(reference_data, columns=columns)

# Convert pressures from MPa to bar
ref_df['Pressure_bar'] = ref_df['Pressure'] * 10

# API endpoint
API_URL = "http://127.0.0.1:5051"
PT_FLASH_ENDPOINT = f"{API_URL}/pt_flash"  # Main endpoint for property calculations

# Constants for CO2
MOLAR_MASS_CO2 = 44.01  # g/mol
MOLAR_MASS_CO2_KG = MOLAR_MASS_CO2 / 1000  # kg/mol

# For debugging - make one test call to see the exact format of returned properties
def print_test_call_results():
    # Use a condition that we know works well
    payload = {
        "composition": [
            {"fluid": "CO2", "fraction": 1.0}
        ],
        "variables": {
            "pressure": {
                "range": {"from": 10, "to": 10},
                "resolution": 1
            },
            "temperature": {
                "range": {"from": 0, "to": 0},
                "resolution": 1
            }
        },
        "calculation": {
            "properties": [
                "density", 
                "liquid_density",
                "vapor_density",
                "internal_energy",
                "enthalpy",
                "entropy",
                "vapor_fraction",
                "phase"
            ],
            "units_system": "SI"
        }
    }
    
    try:
        print("Making test API call to verify returned units...")
        response = requests.post(PT_FLASH_ENDPOINT, json=payload)
        if response.status_code == 200:
            data = response.json()['results'][0]
            print("\nSample API Response Structure:")
            for key, value in data.items():
                if isinstance(value, dict) and 'value' in value and 'unit' in value:
                    print(f"{key}: {value['value']} {value['unit']}")
                else:
                    print(f"{key}: {value}")
            return data
        else:
            print(f"API error for test call: {response.text}")
            return None
    except Exception as e:
        print(f"Exception for test call: {str(e)}")
        return None

# Run the test call to see what units the API returns
test_response = print_test_call_results()

# Functions to convert reference data to match API units
def convert_specific_volume_to_density(volume_m3kg):
    """Convert specific volume (m³/kg) to density (mol/L)"""
    # m³/kg -> kg/m³
    density_kgm3 = 1 / volume_m3kg
    # kg/m³ -> mol/L: divide by molar mass (kg/mol) and multiply by 0.001 (m³ to L)
    density_molL = density_kgm3 * 0.001 / MOLAR_MASS_CO2_KG
    return density_molL

def convert_energy_kJkg_to_Jmol(energy_kJkg):
    """Convert energy (kJ/kg) to energy (J/mol)"""
    # kJ/kg -> J/kg: multiply by 1000
    energy_Jkg = energy_kJkg * 1000
    # J/kg -> J/mol: multiply by molar mass (kg/mol)
    energy_Jmol = energy_Jkg * MOLAR_MASS_CO2_KG
    return energy_Jmol

def convert_entropy_kJkgK_to_JmolK(entropy_kJkgK):
    """Convert entropy (kJ/(kg·K)) to entropy (J/(mol·K))"""
    # kJ/(kg·K) -> J/(kg·K): multiply by 1000
    entropy_JkgK = entropy_kJkgK * 1000
    # J/(kg·K) -> J/(mol·K): multiply by molar mass (kg/mol)
    entropy_JmolK = entropy_JkgK * MOLAR_MASS_CO2_KG
    return entropy_JmolK

# Convert reference data to the units returned by the API
ref_df['density_f_api_units'] = ref_df['vf'].apply(convert_specific_volume_to_density)
ref_df['density_g_api_units'] = ref_df['vg'].apply(convert_specific_volume_to_density)
ref_df['internal_energy_f_api_units'] = ref_df['uf'].apply(convert_energy_kJkg_to_Jmol)
ref_df['internal_energy_g_api_units'] = ref_df['ug'].apply(convert_energy_kJkg_to_Jmol)
ref_df['enthalpy_f_api_units'] = ref_df['hf'].apply(convert_energy_kJkg_to_Jmol)
ref_df['enthalpy_g_api_units'] = ref_df['hg'].apply(convert_energy_kJkg_to_Jmol)
ref_df['entropy_f_api_units'] = ref_df['sf'].apply(convert_entropy_kJkgK_to_JmolK)
ref_df['entropy_g_api_units'] = ref_df['sg'].apply(convert_entropy_kJkgK_to_JmolK)

# Create a dataframe to store the API results
api_results = []

# For saturation properties, we need to run separate liquid and vapor calculations
# We'll do this by setting T and P at saturation, and using a phase-specific calculation
print("\nGathering data for saturation properties comparison...")

# Sample a subset of temperature points to reduce API calls
sample_temps = ref_df['Temperature'].iloc[::3].tolist()  # Take every third point

for temp_idx, temp in enumerate(sample_temps):
    # Lookup corresponding pressure from reference data
    pressure_bar = ref_df.loc[ref_df['Temperature'] == temp, 'Pressure_bar'].values[0]
    
    # LIQUID PHASE CALCULATION
    try:
        print(f"Calculating liquid properties at T={temp}°C, P={pressure_bar:.2f} bar...")
        
        # Use PT flash with a very small temperature offset to ensure liquid phase
        # For CO2, lowering temperature slightly should favor liquid phase
        liquid_payload = {
            "composition": [
                {"fluid": "CO2", "fraction": 1.0}
            ],
            "variables": {
                "pressure": {
                    "range": {"from": pressure_bar, "to": pressure_bar},
                    "resolution": 1
                },
                "temperature": {
                    "range": {"from": temp-0.1, "to": temp-0.1},  # Slight offset for liquid phase
                    "resolution": 1
                }
            },
            "calculation": {
                "properties": [
                    "density", 
                    "internal_energy",
                    "enthalpy",
                    "entropy",
                    "phase"
                ],
                "units_system": "SI"
            }
        }
        
        response = requests.post(PT_FLASH_ENDPOINT, json=liquid_payload)
        if response.status_code == 200:
            data = response.json()['results'][0]
            
            phase = data['phase']['value']
            if 'Liquid' in phase:
                # Extract properties
                liquid_density = data['density']['value']
                liquid_internal_energy = data['internal_energy']['value'] if 'internal_energy' in data else None
                liquid_enthalpy = data['enthalpy']['value'] if 'enthalpy' in data else None
                liquid_entropy = data['entropy']['value'] if 'entropy' in data else None
                
                print(f"  Success: Found liquid properties at T={temp-0.1}°C")
            else:
                print(f"  Warning: Expected Liquid but got {phase} phase at T={temp-0.1}°C")
                liquid_density = None
                liquid_internal_energy = None
                liquid_enthalpy = None
                liquid_entropy = None
        else:
            print(f"  API error for liquid properties: {response.text}")
            liquid_density = None
            liquid_internal_energy = None
            liquid_enthalpy = None
            liquid_entropy = None
    except Exception as e:
        print(f"  Exception for liquid properties: {str(e)}")
        liquid_density = None
        liquid_internal_energy = None
        liquid_enthalpy = None
        liquid_entropy = None
        
    # VAPOR PHASE CALCULATION
    try:
        print(f"Calculating vapor properties at T={temp}°C, P={pressure_bar:.2f} bar...")
        
        # Use PT flash with a very small temperature offset to ensure vapor phase
        # For CO2, raising temperature slightly should favor vapor phase
        vapor_payload = {
            "composition": [
                {"fluid": "CO2", "fraction": 1.0}
            ],
            "variables": {
                "pressure": {
                    "range": {"from": pressure_bar, "to": pressure_bar},
                    "resolution": 1
                },
                "temperature": {
                    "range": {"from": temp+0.1, "to": temp+0.1},  # Slight offset for vapor phase
                    "resolution": 1
                }
            },
            "calculation": {
                "properties": [
                    "density", 
                    "internal_energy",
                    "enthalpy",
                    "entropy",
                    "phase"
                ],
                "units_system": "SI"
            }
        }
        
        response = requests.post(PT_FLASH_ENDPOINT, json=vapor_payload)
        if response.status_code == 200:
            data = response.json()['results'][0]
            
            phase = data['phase']['value']
            if 'Vapor' in phase:
                # Extract properties
                vapor_density = data['density']['value']
                vapor_internal_energy = data['internal_energy']['value'] if 'internal_energy' in data else None
                vapor_enthalpy = data['enthalpy']['value'] if 'enthalpy' in data else None
                vapor_entropy = data['entropy']['value'] if 'entropy' in data else None
                
                print(f"  Success: Found vapor properties at T={temp+0.1}°C")
            else:
                print(f"  Warning: Expected Vapor but got {phase} phase at T={temp+0.1}°C")
                vapor_density = None
                vapor_internal_energy = None
                vapor_enthalpy = None
                vapor_entropy = None
        else:
            print(f"  API error for vapor properties: {response.text}")
            vapor_density = None
            vapor_internal_energy = None
            vapor_enthalpy = None
            vapor_entropy = None
    except Exception as e:
        print(f"  Exception for vapor properties: {str(e)}")
        vapor_density = None
        vapor_internal_energy = None
        vapor_enthalpy = None
        vapor_entropy = None
    
    # Store results
    api_results.append({
        'Temperature': temp,
        'Pressure': pressure_bar / 10,  # Convert back to MPa for reference
        'Pressure_bar': pressure_bar,
        'liquid_density_api': liquid_density,
        'vapor_density_api': vapor_density,
        'liquid_internal_energy_api': liquid_internal_energy,
        'vapor_internal_energy_api': vapor_internal_energy,
        'liquid_enthalpy_api': liquid_enthalpy,
        'vapor_enthalpy_api': vapor_enthalpy,
        'liquid_entropy_api': liquid_entropy,
        'vapor_entropy_api': vapor_entropy
    })

# Convert API results to DataFrame
api_df = pd.DataFrame(api_results)

# Create comparison dataframe by merging reference and API data
comparison_df = pd.merge(ref_df, api_df, on=['Temperature', 'Pressure'], how='inner')

# Calculate percentage differences (directly comparing in API units)
if not comparison_df.empty:
    # Liquid density
    if 'liquid_density_api' in comparison_df.columns and 'density_f_api_units' in comparison_df.columns:
        comparison_df['liquid_density_diff_pct'] = (comparison_df['liquid_density_api'] - comparison_df['density_f_api_units']) / comparison_df['density_f_api_units'] * 100
    
    # Vapor density
    if 'vapor_density_api' in comparison_df.columns and 'density_g_api_units' in comparison_df.columns:
        comparison_df['vapor_density_diff_pct'] = (comparison_df['vapor_density_api'] - comparison_df['density_g_api_units']) / comparison_df['density_g_api_units'] * 100
    
    # Liquid internal energy
    if 'liquid_internal_energy_api' in comparison_df.columns and 'internal_energy_f_api_units' in comparison_df.columns:
        comparison_df['liquid_internal_energy_diff_pct'] = (comparison_df['liquid_internal_energy_api'] - comparison_df['internal_energy_f_api_units']) / comparison_df['internal_energy_f_api_units'] * 100
    
    # Vapor internal energy
    if 'vapor_internal_energy_api' in comparison_df.columns and 'internal_energy_g_api_units' in comparison_df.columns:
        comparison_df['vapor_internal_energy_diff_pct'] = (comparison_df['vapor_internal_energy_api'] - comparison_df['internal_energy_g_api_units']) / comparison_df['internal_energy_g_api_units'] * 100
    
    # Liquid enthalpy
    if 'liquid_enthalpy_api' in comparison_df.columns and 'enthalpy_f_api_units' in comparison_df.columns:
        comparison_df['liquid_enthalpy_diff_pct'] = (comparison_df['liquid_enthalpy_api'] - comparison_df['enthalpy_f_api_units']) / comparison_df['enthalpy_f_api_units'] * 100
    
    # Vapor enthalpy
    if 'vapor_enthalpy_api' in comparison_df.columns and 'enthalpy_g_api_units' in comparison_df.columns:
        comparison_df['vapor_enthalpy_diff_pct'] = (comparison_df['vapor_enthalpy_api'] - comparison_df['enthalpy_g_api_units']) / comparison_df['enthalpy_g_api_units'] * 100
    
    # Liquid entropy
    if 'liquid_entropy_api' in comparison_df.columns and 'entropy_f_api_units' in comparison_df.columns:
        comparison_df['liquid_entropy_diff_pct'] = (comparison_df['liquid_entropy_api'] - comparison_df['entropy_f_api_units']) / comparison_df['entropy_f_api_units'] * 100
    
    # Vapor entropy
    if 'vapor_entropy_api' in comparison_df.columns and 'entropy_g_api_units' in comparison_df.columns:
        comparison_df['vapor_entropy_diff_pct'] = (comparison_df['vapor_entropy_api'] - comparison_df['entropy_g_api_units']) / comparison_df['entropy_g_api_units'] * 100

    # Print comparison header
    print("\n" + "="*80)
    print("Comparison Results: Reference Data vs. API Calculations")
    print("="*80)

    # Print liquid properties comparison
    print("\nLiquid Phase Properties Comparison:")
    diff_columns = [col for col in comparison_df.columns if 'liquid' in col and 'diff_pct' in col]
    if diff_columns:
        summary_cols = ['Temperature', 'Pressure'] + diff_columns
        print(comparison_df[summary_cols].round(2))
        
        # Print average absolute differences
        print("\nAverage Absolute Percentage Differences (Liquid Phase):")
        for col in diff_columns:
            if not comparison_df[col].isnull().all():  # Skip if all values are NaN
                avg_diff = comparison_df[col].abs().mean()
                print(f"{col.replace('_diff_pct', '')}: {avg_diff:.2f}%")
    else:
        print("No liquid phase comparison data available.")

    # Print vapor properties comparison
    print("\nVapor Phase Properties Comparison:")
    diff_columns = [col for col in comparison_df.columns if 'vapor' in col and 'diff_pct' in col]
    if diff_columns:
        summary_cols = ['Temperature', 'Pressure'] + diff_columns
        print(comparison_df[summary_cols].round(2))
        
        # Print average absolute differences
        print("\nAverage Absolute Percentage Differences (Vapor Phase):")
        for col in diff_columns:
            if not comparison_df[col].isnull().all():  # Skip if all values are NaN
                avg_diff = comparison_df[col].abs().mean()
                print(f"{col.replace('_diff_pct', '')}: {avg_diff:.2f}%")
    else:
        print("No vapor phase comparison data available.")

    # Create plots
    # Plot 1: Liquid Density Comparison
    if 'liquid_density_api' in comparison_df.columns and 'density_f_api_units' in comparison_df.columns:
        plt.figure(figsize=(12, 8))
        plt.plot(comparison_df['Temperature'], comparison_df['density_f_api_units'], 'bo-', label='Reference')
        plt.plot(comparison_df['Temperature'], comparison_df['liquid_density_api'], 'ro--', label='API')
        plt.xlabel('Temperature (°C)')
        plt.ylabel('Liquid Density (mol/L)')
        plt.title('Liquid Density: Reference vs API')
        plt.grid(True)
        plt.legend()
        plt.savefig('co2_liquid_density_comparison.png')
        plt.close()
        
    # Plot 2: Vapor Density Comparison
    if 'vapor_density_api' in comparison_df.columns and 'density_g_api_units' in comparison_df.columns:
        plt.figure(figsize=(12, 8))
        plt.plot(comparison_df['Temperature'], comparison_df['density_g_api_units'], 'bo-', label='Reference')
        plt.plot(comparison_df['Temperature'], comparison_df['vapor_density_api'], 'ro--', label='API')
        plt.xlabel('Temperature (°C)')
        plt.ylabel('Vapor Density (mol/L)')
        plt.title('Vapor Density: Reference vs API')
        plt.grid(True)
        plt.legend()
        plt.savefig('co2_vapor_density_comparison.png')
        plt.close()
        
    # Plot 3: Liquid Enthalpy Comparison
    if 'liquid_enthalpy_api' in comparison_df.columns and 'enthalpy_f_api_units' in comparison_df.columns:
        plt.figure(figsize=(12, 8))
        plt.plot(comparison_df['Temperature'], comparison_df['enthalpy_f_api_units'], 'bo-', label='Reference')
        plt.plot(comparison_df['Temperature'], comparison_df['liquid_enthalpy_api'], 'ro--', label='API')
        plt.xlabel('Temperature (°C)')
        plt.ylabel('Liquid Enthalpy (J/mol)')
        plt.title('Liquid Enthalpy: Reference vs API')
        plt.grid(True)
        plt.legend()
        plt.savefig('co2_liquid_enthalpy_comparison.png')
        plt.close()
        
    # Plot 4: Vapor Enthalpy Comparison
    if 'vapor_enthalpy_api' in comparison_df.columns and 'enthalpy_g_api_units' in comparison_df.columns:
        plt.figure(figsize=(12, 8))
        plt.plot(comparison_df['Temperature'], comparison_df['enthalpy_g_api_units'], 'bo-', label='Reference')
        plt.plot(comparison_df['Temperature'], comparison_df['vapor_enthalpy_api'], 'ro--', label='API')
        plt.xlabel('Temperature (°C)')
        plt.ylabel('Vapor Enthalpy (J/mol)')
        plt.title('Vapor Enthalpy: Reference vs API')
        plt.grid(True)
        plt.legend()
        plt.savefig('co2_vapor_enthalpy_comparison.png')
        plt.close()
    
    # Export results to CSV
    comparison_df.to_csv('co2_properties_comparison.csv', index=False)
    print("\nResults saved to 'co2_properties_comparison.csv'")
    print("Individual property comparison plots saved as PNG files.")
    
else:
    print("No matching data points for comparison. Check if the API is returning data correctly.")

# Now test superheated properties with direct PT-Flash
print("\n" + "="*80)
print("Testing Superheated CO2 Properties")
print("="*80)

# Create a subset of data points from Table D2 (Superheated CO2)
superheated_data = [
    # P = 1.0 MPa (10 bar), Sat. T = -40.12 °C
    {'T': -20, 'P': 1.0, 'v': 0.04342, 'h': 455.21, 's': 2.1312},
    {'T': 0, 'P': 1.0, 'v': 0.04799, 'h': 474.04, 's': 2.2028},
    {'T': 20, 'P': 1.0, 'v': 0.05236, 'h': 492.53, 's': 2.2681},
    {'T': 40, 'P': 1.0, 'v': 0.05660, 'h': 510.96, 's': 2.3289},
    {'T': 60, 'P': 1.0, 'v': 0.06074, 'h': 529.47, 's': 2.3862},
    
    # P = 2.0 MPa (20 bar), Sat. T = -19.50 °C
    {'T': 0, 'P': 2.0, 'v': 0.02193, 'h': 460.00, 's': 2.0341},
    {'T': 20, 'P': 2.0, 'v': 0.02453, 'h': 481.32, 's': 2.1095},
    {'T': 40, 'P': 2.0, 'v': 0.02693, 'h': 501.65, 's': 2.1766},
    {'T': 60, 'P': 2.0, 'v': 0.02922, 'h': 521.54, 's': 2.2381},
    {'T': 80, 'P': 2.0, 'v': 0.03143, 'h': 541.27, 's': 2.2956}
]
# Create a DataFrame for superheated CO2 reference data
superheated_df = pd.DataFrame(superheated_data)

# Convert reference specific volume to density in mol/L for API comparison
superheated_df['density_api_units'] = superheated_df['v'].apply(convert_specific_volume_to_density)
superheated_df['enthalpy_api_units'] = superheated_df['h'].apply(convert_energy_kJkg_to_Jmol)
superheated_df['entropy_api_units'] = superheated_df['s'].apply(convert_entropy_kJkgK_to_JmolK)

# Create a dataframe to store the API results for superheated CO2
superheated_api_results = []

# Make API calls for each superheated point
for _, row in superheated_df.iterrows():
    pressure_MPa = row['P']
    pressure_bar = pressure_MPa * 10
    temp = row['T']
    
    # Create payload for PT-Flash
    payload = {
        "composition": [
            {"fluid": "CO2", "fraction": 1.0}
        ],
        "variables": {
            "pressure": {
                "range": {"from": pressure_bar, "to": pressure_bar},
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
                "enthalpy",
                "entropy",
                "phase"
            ],
            "units_system": "SI"
        }
    }
    
    try:
        print(f"Making API call for superheated CO2: T={temp}°C, P={pressure_bar:.2f} bar...")
        response = requests.post(PT_FLASH_ENDPOINT, json=payload)
        if response.status_code == 200:
            data = response.json()['results'][0]
            
            # Extract properties
            phase = data['phase']['value']
            density = data['density']['value'] if 'density' in data else None
            enthalpy = data['enthalpy']['value'] if 'enthalpy' in data else None
            entropy = data['entropy']['value'] if 'entropy' in data else None
            
            # Store results
            superheated_api_results.append({
                'T': temp,
                'P': pressure_MPa,
                'density_api': density,
                'enthalpy_api': enthalpy,
                'entropy_api': entropy,
                'phase': phase
            })
            
            print(f"  Success: T={temp}°C, P={pressure_MPa} MPa, Phase={phase}")
        else:
            print(f"  API error for superheated T={temp}°C, P={pressure_MPa} MPa: {response.text}")
    except Exception as e:
        print(f"  Exception for superheated T={temp}°C, P={pressure_MPa} MPa: {str(e)}")

# Check if we have any superheated API results
if not superheated_api_results:
    print("\nNo valid superheated CO2 API results obtained.")
else:
    # Convert API results to DataFrame
    superheated_api_df = pd.DataFrame(superheated_api_results)
    
    # Merge with reference data
    superheated_comparison = pd.merge(superheated_df, superheated_api_df, on=['T', 'P'])
    
    # Calculate percentage differences
    superheated_comparison['density_diff_pct'] = (superheated_comparison['density_api'] - superheated_comparison['density_api_units']) / superheated_comparison['density_api_units'] * 100
    superheated_comparison['enthalpy_diff_pct'] = (superheated_comparison['enthalpy_api'] - superheated_comparison['enthalpy_api_units']) / superheated_comparison['enthalpy_api_units'] * 100
    superheated_comparison['entropy_diff_pct'] = (superheated_comparison['entropy_api'] - superheated_comparison['entropy_api_units']) / superheated_comparison['entropy_api_units'] * 100
    
    # Print comparison results
    print("\nSuperheated CO2 Comparison Results:")
    print(superheated_comparison[['T', 'P', 'density_api_units', 'density_api', 'density_diff_pct', 'enthalpy_api_units', 'enthalpy_api', 'enthalpy_diff_pct', 'entropy_api_units', 'entropy_api', 'entropy_diff_pct']].round(6))
    
    # Calculate average absolute percentage differences
    print("\nSuperheated CO2 Average Absolute Percentage Differences:")
    density_avg_diff = superheated_comparison['density_diff_pct'].abs().mean()
    enthalpy_avg_diff = superheated_comparison['enthalpy_diff_pct'].abs().mean()
    entropy_avg_diff = superheated_comparison['entropy_diff_pct'].abs().mean()
    
    print(f"Density: {density_avg_diff:.2f}%")
    print(f"Enthalpy: {enthalpy_avg_diff:.2f}%")
    print(f"Entropy: {entropy_avg_diff:.2f}%")
    
    # Create plots for superheated data
    fig, axs = plt.subplots(3, 1, figsize=(12, 15))
    
    # Group by pressure for plotting
    for pressure, group in superheated_comparison.groupby('P'):
        # Plot 1: Density
        axs[0].plot(group['T'], group['density_api_units'], 'o-', label=f'Reference ({pressure} MPa)')
        axs[0].plot(group['T'], group['density_api'], 'x--', label=f'API ({pressure} MPa)')
        
        # Plot 2: Enthalpy
        axs[1].plot(group['T'], group['enthalpy_api_units'], 'o-', label=f'Reference ({pressure} MPa)')
        axs[1].plot(group['T'], group['enthalpy_api'], 'x--', label=f'API ({pressure} MPa)')
        
        # Plot 3: Entropy
        axs[2].plot(group['T'], group['entropy_api_units'], 'o-', label=f'Reference ({pressure} MPa)')
        axs[2].plot(group['T'], group['entropy_api'], 'x--', label=f'API ({pressure} MPa)')
    
    # Set labels and titles
    axs[0].set_xlabel('Temperature (°C)')
    axs[0].set_ylabel('Density (mol/L)')
    axs[0].set_title('Density vs Temperature')
    axs[0].grid(True)
    axs[0].legend()
    
    axs[1].set_xlabel('Temperature (°C)')
    axs[1].set_ylabel('Enthalpy (J/mol)')
    axs[1].set_title('Enthalpy vs Temperature')
    axs[1].grid(True)
    axs[1].legend()
    
    axs[2].set_xlabel('Temperature (°C)')
    axs[2].set_ylabel('Entropy (J/(mol·K))')
    axs[2].set_title('Entropy vs Temperature')
    axs[2].grid(True)
    axs[2].legend()
    
    # Overall title
    fig.suptitle('Superheated CO₂ Properties: Reference vs API Values', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.97])  # Adjust layout to make room for the title
    plt.savefig('co2_superheated_comparison.png')
    plt.close()

    # Export results to CSV
    superheated_comparison.to_csv('co2_superheated_comparison.csv', index=False)
    print("\nSuperheated results saved to 'co2_superheated_comparison.csv'")
    print("Superheated plots saved to 'co2_superheated_comparison.png'")

# Test phase envelope endpoint if available
print("\n" + "="*80)
print("Testing Other API Endpoints")
print("="*80)

# Test phase envelope endpoint for CO2
try:
    print("\nTesting phase_envelope_pt endpoint...")
    payload = {
        "composition": [
            {"fluid": "CO2", "fraction": 1.0}
        ],
        "variables": {
            "temperature": {
                "range": {"from": -50, "to": 31},
                "resolution": 5
            }
        },
        "calculation": {
            "curve_type": "both"
        }
    }
    
    response = requests.post(f"{API_URL}/phase_envelope_pt", json=payload)
    if response.status_code == 200:
        data = response.json()
        
        if 'bubble_curve' in data and 'dew_curve' in data:
            # Extract data for bubble and dew points
            bubble_points = pd.DataFrame(data['bubble_curve'])
            dew_points = pd.DataFrame(data['dew_curve'])
            
            print(f"  Success: Retrieved {len(bubble_points)} bubble points and {len(dew_points)} dew points")
            
            # Plot phase envelope
            plt.figure(figsize=(10, 6))
            
            # Convert pressures to MPa for plotting (if needed)
            if 'pressure' in bubble_points.columns:
                bubble_points['pressure_MPa'] = bubble_points['pressure'] / 10
                dew_points['pressure_MPa'] = dew_points['pressure'] / 10
                
                plt.plot(bubble_points['temperature'], bubble_points['pressure_MPa'], 'bo-', label='Bubble Line (API)')
                plt.plot(dew_points['temperature'], dew_points['pressure_MPa'], 'ro-', label='Dew Line (API)')
                
                # Add reference data from Table D1
                plt.plot(ref_df['Temperature'], ref_df['Pressure'], 'ko--', label='Reference Data')
                
                plt.xlabel('Temperature (°C)')
                plt.ylabel('Pressure (MPa)')
                plt.title('CO₂ Phase Envelope: API vs Reference')
                plt.grid(True)
                plt.legend()
                plt.savefig('co2_phase_envelope.png')
                plt.close()
                
                print("  Phase envelope plot saved to 'co2_phase_envelope.png'")
        else:
            print("  API response doesn't contain expected phase envelope data")
    else:
        print(f"  API error from phase_envelope_pt: {response.text}")
except Exception as e:
    print(f"  Exception when testing phase_envelope_pt: {str(e)}")

# Test other available endpoints
try:
    print("\nTesting additional endpoints (if available)...")
    
    # Example: Testing phase_envelope_ph endpoint
    payload = {
        "composition": [
            {"fluid": "CO2", "fraction": 1.0}
        ],
        "variables": {
            "pressure": {
                "range": {"from": 10, "to": 100},
                "resolution": 10
            }
        },
        "calculation": {
            "curve_type": "both"
        }
    }
    
    response = requests.post(f"{API_URL}/phase_envelope_ph", json=payload)
    if response.status_code == 200:
        print("  phase_envelope_ph endpoint test: Success")
    else:
        print(f"  phase_envelope_ph endpoint test: Failed - {response.status_code}")
        
    # Test PH-Flash endpoint
    payload = {
        "composition": [
            {"fluid": "CO2", "fraction": 1.0}
        ],
        "variables": {
            "pressure": {
                "range": {"from": 20, "to": 20},
                "resolution": 1
            },
            "enthalpy": {
                "range": {"from": 300, "to": 300},
                "resolution": 1
            }
        },
        "calculation": {
            "properties": ["temperature", "density", "entropy", "phase"],
            "units_system": "SI"
        }
    }
    
    response = requests.post(f"{API_URL}/ph_flash", json=payload)
    if response.status_code == 200:
        print("  ph_flash endpoint test: Success")
    else:
        print(f"  ph_flash endpoint test: Failed - {response.status_code}")
    
except Exception as e:
    print(f"  Exception when testing additional endpoints: {str(e)}")

print("\n" + "="*80)
print("Testing Complete")
print("="*80)
