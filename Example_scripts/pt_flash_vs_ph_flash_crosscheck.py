import requests
import pandas as pd
import numpy as np
import time
from tabulate import tabulate

# API base URL - adjust if your API is hosted elsewhere
BASE_URL = "http://localhost:5051"

def sanity_check_eos():
    """
    Comprehensive sanity check of the Span-Wagner EOS API by comparing PT-flash and PH-flash results.
    
    This function:
    1. Calculates thermodynamic properties using PT-flash at various pressure and temperature points
    2. Uses the resulting enthalpy with the same pressure in PH-flash
    3. Compares all properties between the two calculations to verify consistency
    
    Properties compared include:
    - Temperature and density (primary validation)
    - Internal energy, entropy, heat capacities (Cv, Cp)
    - Sound speed, compressibility factor
    - Transport properties (viscosity, thermal conductivity)
    
    This comprehensive validation ensures that the EOS implementation is thermodynamically 
    consistent between different flash calculation routes.
    """
    print("Starting REFPROP API Sanity Check")
    print("=" * 80)
    
    # Define test points - expanded to cover a wider range of conditions
    # Generate a grid of points from 1 to 300 bar and -50°C to 100°C
    pressures = np.linspace(1, 300, 6)  # 6 pressure points
    temperatures = np.linspace(-50, 100, 5)  # 5 temperature points
    
    test_points = []
    for pressure in pressures:
        for temperature in temperatures:
            test_points.append({
                "pressure": round(pressure, 1),
                "temperature": round(temperature, 1)
            })
    
    # Add a few specific points of interest
    additional_points = [
        {"pressure": 73.9, "temperature": 31.1},  # Near CO2 critical point
        {"pressure": 50, "temperature": -30},     # Low temperature region
        {"pressure": 120, "temperature": 40},     # High pressure region
        {"pressure": 5, "temperature": -40},      # Very low pressure, low temperature
        {"pressure": 250, "temperature": 80}      # Very high pressure, high temperature
    ]
    
    test_points.extend(additional_points)
    
    print(f"Testing {len(test_points)} pressure-temperature points")
    
    # Fluid composition: 90% CO2, 10% Methane
    composition = [
        {"fluid": "CO2", "fraction": 0.90},
        {"fluid": "METHANE", "fraction": 0.10}
    ]
    
    # Part 1: PT-flash calculations
    # Create batched PT-flash requests to optimize API calls
    batch_size = 10  # Adjust based on API performance
    pt_batches = [test_points[i:i + batch_size] for i in range(0, len(test_points), batch_size)]
    
    pt_results = []
    print("Running PT-flash calculations...")
    
    for batch_idx, batch in enumerate(pt_batches, 1):
        print(f"Processing PT-flash batch {batch_idx}/{len(pt_batches)}...")
        
        # For each point in the batch, prepare a single PT-flash call with narrow ranges
        for point in batch:
            # Prepare the PT-flash request payload
            pt_payload = {
                "composition": composition,
                "pressure_range": {
                    "from": point["pressure"],
                    "to": point["pressure"]
                },
                "temperature_range": {
                    "from": point["temperature"],
                    "to": point["temperature"]
                },
                "pressure_resolution": 1,
                "temperature_resolution": 1,
                "properties": [
                    "density",
                    "enthalpy",
                    "internal_energy",
                    "entropy",
                    "Cv",
                    "Cp",
                    "sound_speed",
                    "compressibility_factor",
                    "phase",
                    "vapor_fraction",
                    "viscosity",
                    "thermal_conductivity"
                ],
                "units_system": "SI"
            }
            
            # Send the PT-flash request
            try:
                response = requests.post(f"{BASE_URL}/pt_flash", json=pt_payload, timeout=30)
                response.raise_for_status()
                
                # Extract results
                result = response.json()["results"][0]
                
                # Extract all available properties from PT-flash
                pt_result = {
                    "pressure": point["pressure"],
                    "temperature": point["temperature"],
                    "density": result["density"]["value"],
                    "enthalpy": result["enthalpy"]["value"],
                    "phase": result["phase"]["value"],
                    "vapor_fraction": result["vapor_fraction"]["value"]
                }
                
                # Add all additional properties that were requested
                for prop in ["internal_energy", "entropy", "Cv", "Cp", "sound_speed", 
                            "compressibility_factor", "viscosity", "thermal_conductivity"]:
                    if prop in result:
                        pt_result[prop] = result[prop]["value"]
                
                pt_results.append(pt_result)
                print(f"PT-flash completed for P={point['pressure']} bar, T={point['temperature']}°C, Phase={pt_result['phase']}")
                
            except Exception as e:
                print(f"Error in PT-flash calculation for P={point['pressure']} bar, T={point['temperature']}°C: {str(e)}")
                continue
                
        # Brief pause between batches to avoid overwhelming the API
        if batch_idx < len(pt_batches):
            time.sleep(0.5)
    
    # Pause to avoid overloading the API
    time.sleep(1)
    
    # Part 2: PH-flash calculations using the enthalpy from PT-flash
    ph_results = []
    print("\nRunning PH-flash calculations with enthalpy from PT-flash...")
    
    # Create batches for PH-flash as well
    ph_batches = [pt_results[i:i + batch_size] for i in range(0, len(pt_results), batch_size)]
    
    for batch_idx, batch in enumerate(ph_batches, 1):
        print(f"Processing PH-flash batch {batch_idx}/{len(ph_batches)}...")
        
        for pt_result in batch:
            # Prepare the PH-flash request payload
            ph_payload = {
                "composition": composition,
                "pressure_range": {
                    "from": pt_result["pressure"],
                    "to": pt_result["pressure"]
                },
                "enthalpy_range": {
                    "from": pt_result["enthalpy"],
                    "to": pt_result["enthalpy"]
                },
                "pressure_resolution": 1,
                "enthalpy_resolution": 1,
                "properties": [
                    "temperature",
                    "density",
                    "internal_energy",
                    "entropy",
                    "Cv",
                    "Cp",
                    "sound_speed",
                    "compressibility_factor",
                    "phase",
                    "vapor_fraction",
                    "viscosity",
                    "thermal_conductivity"
                ],
                "units_system": "SI"
            }
            
            # Send the PH-flash request
            try:
                response = requests.post(f"{BASE_URL}/ph_flash", json=ph_payload, timeout=30)
                response.raise_for_status()
                
                # Extract results
                result = response.json()["results"][0]
                
                # Initialize with basic properties
                ph_result = {
                    "pressure": pt_result["pressure"],
                    "enthalpy": pt_result["enthalpy"],
                    "temperature": result["temperature"]["value"],
                    "density": result["density"]["value"],
                    "phase": result["phase"]["value"],
                    "vapor_fraction": result["vapor_fraction"]["value"],
                    "original_temperature": pt_result["temperature"],
                    "original_density": pt_result["density"]
                }
                
                # Add all additional properties for comparison
                for prop in ["internal_energy", "entropy", "Cv", "Cp", "sound_speed", 
                            "compressibility_factor", "viscosity", "thermal_conductivity"]:
                    if prop in result and prop in pt_result:
                        ph_result[prop] = result[prop]["value"]
                        ph_result[f"original_{prop}"] = pt_result[prop]
                
                ph_results.append(ph_result)
                print(f"PH-flash completed for P={pt_result['pressure']} bar, H={pt_result['enthalpy']:.2f} J/mol")
                
            except Exception as e:
                print(f"Error in PH-flash calculation for P={pt_result['pressure']} bar, H={pt_result['enthalpy']:.2f} J/mol: {str(e)}")
                continue
        
        # Brief pause between batches
        if batch_idx < len(ph_batches):
            time.sleep(0.5)
    
    # Part 3: Compare results and calculate differences
    print("\nComparing PT-flash and PH-flash results:")
    print("=" * 80)
    
    comparison_results = []
    
    for ph_result in ph_results:
        # Initialize comparison dictionary with basic properties
        comparison = {
            "pressure": ph_result["pressure"],
            "original_T": ph_result["original_temperature"],
            "calculated_T": ph_result["temperature"],
            "phase": ph_result["phase"],
            "vapor_fraction": ph_result["vapor_fraction"]
        }
        
        # Process each property with proper error handling
        properties_to_compare = [
            {"name": "temperature", "short": "T", "original_key": "original_temperature", "calculated_key": "temperature", "abs_scale": 1.0, "use_kelvin": True},
            {"name": "density", "short": "D", "original_key": "original_density", "calculated_key": "density", "abs_scale": 1.0, "use_kelvin": False},
            {"name": "internal_energy", "short": "U", "original_key": "original_internal_energy", "calculated_key": "internal_energy", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "entropy", "short": "S", "original_key": "original_entropy", "calculated_key": "entropy", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "Cv", "short": "Cv", "original_key": "original_Cv", "calculated_key": "Cv", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "Cp", "short": "Cp", "original_key": "original_Cp", "calculated_key": "Cp", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "sound_speed", "short": "W", "original_key": "original_sound_speed", "calculated_key": "sound_speed", "abs_scale": 0.1, "use_kelvin": False},
            {"name": "compressibility_factor", "short": "Z", "original_key": "original_compressibility_factor", "calculated_key": "compressibility_factor", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "viscosity", "short": "Visc", "original_key": "original_viscosity", "calculated_key": "viscosity", "abs_scale": 0.01, "use_kelvin": False},
            {"name": "thermal_conductivity", "short": "TC", "original_key": "original_thermal_conductivity", "calculated_key": "thermal_conductivity", "abs_scale": 0.001, "use_kelvin": False}
        ]
        
        # Calculate differences for each property
        for prop in properties_to_compare:
            if prop["original_key"] in ph_result and prop["calculated_key"] in ph_result:
                original_val = ph_result[prop["original_key"]]
                calculated_val = ph_result[prop["calculated_key"]]
                
                # Store original and calculated values
                comparison[f"original_{prop['short']}"] = original_val
                comparison[f"calculated_{prop['short']}"] = calculated_val
                
                # Calculate absolute difference (scaled for better readability)
                diff_abs = (calculated_val - original_val) * prop["abs_scale"]
                comparison[f"{prop['short']}_diff_abs"] = diff_abs
                
                # Calculate relative difference (special handling for temperature)
                if prop["use_kelvin"]:
                    # For temperature, convert to Kelvin for relative comparison
                    diff_rel = diff_abs / ((original_val + 273.15) * prop["abs_scale"]) * 100
                else:
                    # Prevent division by zero
                    if original_val != 0:
                        diff_rel = diff_abs / (original_val * prop["abs_scale"]) * 100
                    else:
                        diff_rel = np.nan
                        
                comparison[f"{prop['short']}_diff_rel"] = diff_rel
        
        comparison_results.append(comparison)
    
    # Create a pandas DataFrame for better display
    df = pd.DataFrame(comparison_results)
    
    # Identify all property columns for rounding and statistics
    property_columns = {}
    for col in df.columns:
        if col.endswith('_diff_abs') or col.endswith('_diff_rel'):
            property_columns[col] = 6
        elif col.startswith('original_') or col.startswith('calculated_'):
            property_columns[col] = 4
    
    # Round values for better readability
    df = df.round(property_columns)
    
    # Calculate summary statistics for each property
    print("\nSummary Statistics:")
    property_shorts = ["T", "D", "U", "S", "Cv", "Cp", "W", "Z", "Visc", "TC"]
    property_names = [
        "Temperature", "Density", "Internal Energy", "Entropy", 
        "Cv", "Cp", "Sound Speed", "Compressibility", 
        "Viscosity", "Thermal Conductivity"
    ]
    property_units = [
        "°C", "mol/L", "J/mol", "J/(mol·K)", 
        "J/(mol·K)", "J/(mol·K)", "m/s", "-", 
        "μPa·s", "W/(m·K)"
    ]
    
    # Table for summary statistics
    summary_table = []
    
    for short, name, unit in zip(property_shorts, property_names, property_units):
        abs_col = f"{short}_diff_abs"
        rel_col = f"{short}_diff_rel"
        
        if abs_col in df.columns and rel_col in df.columns:
            # Calculate statistics
            abs_mean = df[abs_col].abs().mean()
            abs_max = df[abs_col].abs().max()
            rel_mean = df[rel_col].abs().mean()
            rel_max = df[rel_col].abs().max()
            
            # Add to summary table
            summary_table.append([
                name, unit, f"{abs_mean:.6f}", f"{abs_max:.6f}", 
                f"{rel_mean:.4f}%", f"{rel_max:.4f}%"
            ])
    
    # Display the summary statistics table
    summary_headers = ["Property", "Unit", "Mean Abs Diff", "Max Abs Diff", "Mean Rel Diff", "Max Rel Diff"]
    print(tabulate(summary_table, headers=summary_headers, tablefmt="grid"))
    
    # Display the temperature and density comparison as the main table (most important properties)
    print("\nDetailed Comparison of Temperature and Density:")
    headers = [
        "P (bar)", "Original T (°C)", "Calc T (°C)", "T diff", "T diff %", 
        "Original D (mol/L)", "Calc D (mol/L)", "D diff", "D diff %", "Phase"
    ]
    
    # Use only columns that definitely exist
    display_columns = ['pressure']
    for col in ['original_T', 'calculated_T', 'T_diff_abs', 'T_diff_rel',
                'original_D', 'calculated_D', 'D_diff_abs', 'D_diff_rel', 'phase']:
        if col in df.columns:
            display_columns.append(col)
    
    table_data = df[display_columns].values.tolist()
    print(tabulate(table_data, headers=headers[:len(display_columns)], tablefmt="grid"))
    
    # Check for NaN values in specific columns (T and D)
    for col in ['T_diff_rel', 'D_diff_rel']:
        if col in df.columns:
            df_stats = df[col].abs().describe()
            if col == 'T_diff_rel':
                t_diff_rel_mean = df_stats['mean']
            elif col == 'D_diff_rel':
                d_diff_rel_mean = df_stats['mean']
    
    # Add options to export detailed results
    export_csv = False  # Set to True if you want to export results to CSV
    if export_csv:
        csv_filename = "eos_sanity_check_results.csv"
        df.to_csv(csv_filename, index=False)
        print(f"\nDetailed results exported to {csv_filename}")
    
    # Display detailed property comparison for a specific test point (e.g., the first one)
    if len(df) > 0:
        print("\nDetailed comparison for first test point:")
        first_point = df.iloc[0]
        detailed_comparison = []
        
        for short, name, unit in zip(property_shorts, property_names, property_units):
            orig_col = f"original_{short}"
            calc_col = f"calculated_{short}"
            abs_col = f"{short}_diff_abs"
            rel_col = f"{short}_diff_rel"
            
            if orig_col in first_point and calc_col in first_point:
                detailed_comparison.append([
                    name,
                    unit,
                    f"{first_point[orig_col]:.6g}",
                    f"{first_point[calc_col]:.6g}",
                    f"{first_point.get(abs_col, 'N/A'):.6g}" if abs_col in first_point else "N/A",
                    f"{first_point.get(rel_col, 'N/A'):.4g}%" if rel_col in first_point else "N/A"
                ])
        
        detail_headers = ["Property", "Unit", "PT Value", "PH Value", "Abs Diff", "Rel Diff"]
        print(tabulate(detailed_comparison, headers=detail_headers, tablefmt="grid"))
    
    # Calculate overall consistency scores
    property_scores = {}
    for short in property_shorts:
        rel_col = f"{short}_diff_rel"
        if rel_col in df.columns:
            # Calculate percentage of points with relative difference < 0.5%
            good_points = (df[rel_col].abs() < 0.5).sum()
            acceptable_points = ((df[rel_col].abs() >= 0.5) & (df[rel_col].abs() < 2.0)).sum()
            
            if len(df) > 0:
                good_percent = good_points / len(df) * 100
                acceptable_percent = acceptable_points / len(df) * 100
                property_scores[short] = (good_percent, acceptable_percent)

    # Evaluate overall sanity check result
    if 'T' in property_scores and 'D' in property_scores:
        t_good, t_acceptable = property_scores['T']
        d_good, d_acceptable = property_scores['D']
        
        if t_good > 80 and d_good > 80:
            print("\nSANITY CHECK PASSED: PT-flash and PH-flash results are highly consistent.")
        elif (t_good + t_acceptable) > 80 and (d_good + d_acceptable) > 80:
            print("\nSANITY CHECK PARTIALLY PASSED: Some discrepancies exist but are within reasonable limits.")
        else:
            print("\nSANITY CHECK FAILED: Significant discrepancies between PT-flash and PH-flash results.")
    else:
        # Fallback to simpler criteria
        try:
            if t_diff_rel_mean < 0.5 and d_diff_rel_mean < 0.5:
                print("\nSANITY CHECK PASSED: PT-flash and PH-flash results are consistent.")
            elif t_diff_rel_mean < 2.0 and d_diff_rel_mean < 2.0:
                print("\nSANITY CHECK PARTIALLY PASSED: Some discrepancies exist but are within reasonable limits.")
            else:
                print("\nSANITY CHECK FAILED: Significant discrepancies between PT-flash and PH-flash results.")
        except:
            print("\nCould not evaluate overall sanity check result due to missing or invalid comparison data.")
    
    # Check for any specific property issues
    print("\nProperty-specific observations:")
    for short, name in zip(property_shorts, property_names):
        rel_col = f"{short}_diff_rel"
        if rel_col in df.columns and not df[rel_col].isnull().all():
            mean_diff = df[rel_col].abs().mean()
            
            if mean_diff > 5.0:
                print(f"- Warning: {name} shows large discrepancies (mean {mean_diff:.2f}%)")
            elif mean_diff > 2.0:
                print(f"- Note: {name} shows moderate discrepancies (mean {mean_diff:.2f}%)")
    
    print("\nPossible reasons for discrepancies:")
    print("1. Numerical precision in the solvers")
    print("2. Iteration tolerance settings in REFPROP")
    print("3. Different solution paths in PT vs PH flash calculations")
    print("4. Properties near phase boundaries or critical region")
    print("5. Implementation differences between PT-flash and PH-flash algorithms")
    print("6. Thermodynamic consistency of the underlying equation of state")
    
    return df

if __name__ == "__main__":
    try:
        results = sanity_check_eos()
        print("\nSanity check completed successfully.")
    except Exception as e:
        print(f"\nError during sanity check: {str(e)}")