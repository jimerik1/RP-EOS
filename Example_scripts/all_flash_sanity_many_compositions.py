import requests
import pandas as pd
import numpy as np
import time
from tabulate import tabulate

# API base URL - adjust if your API is hosted elsewhere
BASE_URL = "http://localhost:5051"

def sanity_check_eos_with_composition(composition, description):
    """
    Comprehensive sanity check of the Span-Wagner EOS API by comparing PT-flash, PH-flash, and TS-flash results
    for a specific composition.
    
    Args:
        composition: List of dictionaries with fluid name and fraction
        description: Description of the test for reporting
    
    Returns:
        Dictionary with test results and statistics
    """
    print(f"\n\n{'=' * 80}")
    print(f"Starting REFPROP API Sanity Check for: {description}")
    print(f"{'=' * 80}")
    
    # Display composition
    comp_table = [[c["fluid"], f"{c['fraction']*100:.1f}%"] for c in composition]
    print("\nTesting composition:")
    print(tabulate(comp_table, headers=["Component", "Fraction"], tablefmt="grid"))
    
    # Define test points - expanded to cover a wider range of conditions
    # Generate a grid of points from 1 to 300 bar and -50°C to 100°C
    pressures = np.linspace(1, 300, 10)  # 10 pressure points
    temperatures = np.linspace(-50, 100, 10)  # 10 temperature points
    
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
                    "entropy": result["entropy"]["value"],
                    "phase": result["phase"]["value"],
                    "vapor_fraction": result["vapor_fraction"]["value"]
                }
                
                # Add all additional properties that were requested
                for prop in ["internal_energy", "Cv", "Cp", "sound_speed", 
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
                    "entropy": result["entropy"]["value"],
                    "phase": result["phase"]["value"],
                    "vapor_fraction": result["vapor_fraction"]["value"],
                    "original_temperature": pt_result["temperature"],
                    "original_density": pt_result["density"],
                    "original_entropy": pt_result["entropy"]
                }
                
                # Add all additional properties for comparison
                for prop in ["internal_energy", "Cv", "Cp", "sound_speed", 
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
    
    # Pause to avoid overloading the API
    time.sleep(1)
    
    # Part 3: TS-flash calculations using the entropy from PT-flash
    ts_results = []
    print("\nRunning TS-flash calculations with entropy from PT-flash...")
    
    # Create batches for TS-flash as well
    ts_batches = [pt_results[i:i + batch_size] for i in range(0, len(pt_results), batch_size)]
    
    for batch_idx, batch in enumerate(ts_batches, 1):
        print(f"Processing TS-flash batch {batch_idx}/{len(ts_batches)}...")
        
        for pt_result in batch:
            # Prepare the TS-flash request payload
            ts_payload = {
                "composition": composition,
                "temperature_range": {
                    "from": pt_result["temperature"],
                    "to": pt_result["temperature"]
                },
                "entropy_range": {
                    "from": pt_result["entropy"],
                    "to": pt_result["entropy"]
                },
                "temperature_resolution": 1,
                "entropy_resolution": 1,
                "properties": [
                    "pressure",
                    "density",
                    "internal_energy",
                    "enthalpy",
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
            
            # Send the TS-flash request
            try:
                response = requests.post(f"{BASE_URL}/ts_flash", json=ts_payload, timeout=30)
                response.raise_for_status()
                
                # Extract results
                result = response.json()["results"][0]
                
                # Initialize with basic properties
                ts_result = {
                    "temperature": pt_result["temperature"],
                    "entropy": pt_result["entropy"],
                    "pressure": result["pressure"]["value"],
                    "density": result["density"]["value"],
                    "enthalpy": result["enthalpy"]["value"],
                    "phase": result["phase"]["value"],
                    "vapor_fraction": result["vapor_fraction"]["value"],
                    "original_pressure": pt_result["pressure"],
                    "original_density": pt_result["density"],
                    "original_enthalpy": pt_result["enthalpy"]
                }
                
                # Add all additional properties for comparison
                for prop in ["internal_energy", "Cv", "Cp", "sound_speed", 
                            "compressibility_factor", "viscosity", "thermal_conductivity"]:
                    if prop in result and prop in pt_result:
                        ts_result[prop] = result[prop]["value"]
                        ts_result[f"original_{prop}"] = pt_result[prop]
                
                ts_results.append(ts_result)
                print(f"TS-flash completed for T={pt_result['temperature']}°C, S={pt_result['entropy']:.2f} J/(mol·K)")
                
            except Exception as e:
                print(f"Error in TS-flash calculation for T={pt_result['temperature']}°C, S={pt_result['entropy']:.2f} J/(mol·K): {str(e)}")
                continue
        
        # Brief pause between batches
        if batch_idx < len(ts_batches):
            time.sleep(0.5)
    
    # Part 4: Compare results and calculate differences
    print("\nComparing PT-flash, PH-flash, and TS-flash results:")
    print("=" * 80)
    
    # 4.1: First comparison - PT vs PH
    print("\nComparison 1: PT-flash vs PH-flash")
    pt_ph_comparison_results = []
    
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
            {"name": "entropy", "short": "S", "original_key": "original_entropy", "calculated_key": "entropy", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "internal_energy", "short": "U", "original_key": "original_internal_energy", "calculated_key": "internal_energy", "abs_scale": 0.001, "use_kelvin": False},
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
        
        pt_ph_comparison_results.append(comparison)
    
    # 4.2: Second comparison - PT vs TS
    print("\nComparison 2: PT-flash vs TS-flash")
    pt_ts_comparison_results = []
    
    for ts_result in ts_results:
        # Initialize comparison dictionary with basic properties
        comparison = {
            "temperature": ts_result["temperature"],
            "original_P": ts_result["original_pressure"],
            "calculated_P": ts_result["pressure"],
            "phase": ts_result["phase"],
            "vapor_fraction": ts_result["vapor_fraction"]
        }
        
        # Process each property with proper error handling
        properties_to_compare = [
            {"name": "pressure", "short": "P", "original_key": "original_pressure", "calculated_key": "pressure", "abs_scale": 1.0, "use_kelvin": False},
            {"name": "density", "short": "D", "original_key": "original_density", "calculated_key": "density", "abs_scale": 1.0, "use_kelvin": False},
            {"name": "enthalpy", "short": "H", "original_key": "original_enthalpy", "calculated_key": "enthalpy", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "internal_energy", "short": "U", "original_key": "original_internal_energy", "calculated_key": "internal_energy", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "Cv", "short": "Cv", "original_key": "original_Cv", "calculated_key": "Cv", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "Cp", "short": "Cp", "original_key": "original_Cp", "calculated_key": "Cp", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "sound_speed", "short": "W", "original_key": "original_sound_speed", "calculated_key": "sound_speed", "abs_scale": 0.1, "use_kelvin": False},
            {"name": "compressibility_factor", "short": "Z", "original_key": "original_compressibility_factor", "calculated_key": "compressibility_factor", "abs_scale": 0.001, "use_kelvin": False},
            {"name": "viscosity", "short": "Visc", "original_key": "original_viscosity", "calculated_key": "viscosity", "abs_scale": 0.01, "use_kelvin": False},
            {"name": "thermal_conductivity", "short": "TC", "original_key": "original_thermal_conductivity", "calculated_key": "thermal_conductivity", "abs_scale": 0.001, "use_kelvin": False}
        ]
        
        # Calculate differences for each property
        for prop in properties_to_compare:
            if prop["original_key"] in ts_result and prop["calculated_key"] in ts_result:
                original_val = ts_result[prop["original_key"]]
                calculated_val = ts_result[prop["calculated_key"]]
                
                # Store original and calculated values
                comparison[f"original_{prop['short']}"] = original_val
                comparison[f"calculated_{prop['short']}"] = calculated_val
                
                # Calculate absolute difference (scaled for better readability)
                diff_abs = (calculated_val - original_val) * prop["abs_scale"]
                comparison[f"{prop['short']}_diff_abs"] = diff_abs
                
                # Calculate relative difference
                if original_val != 0:
                    diff_rel = diff_abs / (original_val * prop["abs_scale"]) * 100
                else:
                    diff_rel = np.nan
                        
                comparison[f"{prop['short']}_diff_rel"] = diff_rel
        
        pt_ts_comparison_results.append(comparison)
    
    # 4.3: Create DataFrames for better display
    df_pt_ph = pd.DataFrame(pt_ph_comparison_results)
    df_pt_ts = pd.DataFrame(pt_ts_comparison_results)
    
    # Identify all property columns for rounding and statistics
    property_columns = {}
    for col in df_pt_ph.columns:
        if col.endswith('_diff_abs') or col.endswith('_diff_rel'):
            property_columns[col] = 6
        elif col.startswith('original_') or col.startswith('calculated_'):
            property_columns[col] = 4
    
    # Round values for better readability
    df_pt_ph = df_pt_ph.round(property_columns)
    df_pt_ts = df_pt_ts.round(property_columns)
    
    # 5: Calculate summary statistics for each comparison
    print("\nSummary Statistics for PT-flash vs PH-flash:")
    property_shorts_pt_ph = ["T", "D", "S", "U", "Cv", "Cp", "W", "Z", "Visc", "TC"]
    property_names_pt_ph = [
        "Temperature", "Density", "Entropy", "Internal Energy", 
        "Cv", "Cp", "Sound Speed", "Compressibility", 
        "Viscosity", "Thermal Conductivity"
    ]
    property_units_pt_ph = [
        "°C", "mol/L", "J/(mol·K)", "J/mol", 
        "J/(mol·K)", "J/(mol·K)", "m/s", "-", 
        "μPa·s", "W/(m·K)"
    ]
    
    # Table for PT vs PH summary statistics
    summary_table_pt_ph = []
    pt_ph_stats = {}
    
    for short, name, unit in zip(property_shorts_pt_ph, property_names_pt_ph, property_units_pt_ph):
        abs_col = f"{short}_diff_abs"
        rel_col = f"{short}_diff_rel"
        
        if abs_col in df_pt_ph.columns and rel_col in df_pt_ph.columns:
            # Calculate statistics
            abs_mean = df_pt_ph[abs_col].abs().mean()
            abs_max = df_pt_ph[abs_col].abs().max()
            rel_mean = df_pt_ph[rel_col].abs().mean()
            rel_max = df_pt_ph[rel_col].abs().max()
            
            # Store stats for return value
            pt_ph_stats[short] = {
                'abs_mean': abs_mean,
                'abs_max': abs_max,
                'rel_mean': rel_mean,
                'rel_max': rel_max
            }
            
            # Add to summary table
            summary_table_pt_ph.append([
                name, unit, f"{abs_mean:.6f}", f"{abs_max:.6f}", 
                f"{rel_mean:.4f}%", f"{rel_max:.4f}%"
            ])
    
    # Display the PT vs PH summary statistics table
    summary_headers = ["Property", "Unit", "Mean Abs Diff", "Max Abs Diff", "Mean Rel Diff", "Max Rel Diff"]
    print(tabulate(summary_table_pt_ph, headers=summary_headers, tablefmt="grid"))
    
    print("\nSummary Statistics for PT-flash vs TS-flash:")
    property_shorts_pt_ts = ["P", "D", "H", "U", "Cv", "Cp", "W", "Z", "Visc", "TC"]
    property_names_pt_ts = [
        "Pressure", "Density", "Enthalpy", "Internal Energy", 
        "Cv", "Cp", "Sound Speed", "Compressibility", 
        "Viscosity", "Thermal Conductivity"
    ]
    property_units_pt_ts = [
        "bar", "mol/L", "J/mol", "J/mol", 
        "J/(mol·K)", "J/(mol·K)", "m/s", "-", 
        "μPa·s", "W/(m·K)"
    ]
    
    # Table for PT vs TS summary statistics
    summary_table_pt_ts = []
    pt_ts_stats = {}
    
    for short, name, unit in zip(property_shorts_pt_ts, property_names_pt_ts, property_units_pt_ts):
        abs_col = f"{short}_diff_abs"
        rel_col = f"{short}_diff_rel"
        
        if abs_col in df_pt_ts.columns and rel_col in df_pt_ts.columns:
            # Calculate statistics
            abs_mean = df_pt_ts[abs_col].abs().mean()
            abs_max = df_pt_ts[abs_col].abs().max()
            rel_mean = df_pt_ts[rel_col].abs().mean()
            rel_max = df_pt_ts[rel_col].abs().max()
            
            # Store stats for return value
            pt_ts_stats[short] = {
                'abs_mean': abs_mean,
                'abs_max': abs_max,
                'rel_mean': rel_mean,
                'rel_max': rel_max
            }
            
            # Add to summary table
            summary_table_pt_ts.append([
                name, unit, f"{abs_mean:.6f}", f"{abs_max:.6f}", 
                f"{rel_mean:.4f}%", f"{rel_max:.4f}%"
            ])
    
    # Display the PT vs TS summary statistics table
    print(tabulate(summary_table_pt_ts, headers=summary_headers, tablefmt="grid"))
    
    # 6: Display key comparison details
    # Display the temperature and pressure comparison as main tables (most important properties)
    print("\nDetailed PT-flash vs PH-flash Comparison (Temperature and Density):")
    headers_pt_ph = [
        "P (bar)", "Original T (°C)", "Calc T (°C)", "T diff", "T diff %", 
        "Original D (mol/L)", "Calc D (mol/L)", "D diff", "D diff %", "Phase"
    ]
    
    # Use only columns that definitely exist
    display_columns_pt_ph = ['pressure']
    for col in ['original_T', 'calculated_T', 'T_diff_abs', 'T_diff_rel',
               'original_D', 'calculated_D', 'D_diff_abs', 'D_diff_rel', 'phase']:
        if col in df_pt_ph.columns:
            display_columns_pt_ph.append(col)
    
    table_data_pt_ph = df_pt_ph[display_columns_pt_ph].head(10).values.tolist()  # Show first 10 rows
    print(tabulate(table_data_pt_ph, headers=headers_pt_ph[:len(display_columns_pt_ph)], tablefmt="grid"))
    
    print("\nDetailed PT-flash vs TS-flash Comparison (Pressure and Density):")
    headers_pt_ts = [
        "T (°C)", "Original P (bar)", "Calc P (bar)", "P diff", "P diff %", 
        "Original D (mol/L)", "Calc D (mol/L)", "D diff", "D diff %", "Phase"
    ]
    
    # Use only columns that definitely exist
    display_columns_pt_ts = ['temperature']
    for col in ['original_P', 'calculated_P', 'P_diff_abs', 'P_diff_rel',
               'original_D', 'calculated_D', 'D_diff_abs', 'D_diff_rel', 'phase']:
        if col in df_pt_ts.columns:
            display_columns_pt_ts.append(col)
    
    table_data_pt_ts = df_pt_ts[display_columns_pt_ts].head(10).values.tolist()  # Show first 10 rows
    print(tabulate(table_data_pt_ts, headers=headers_pt_ts[:len(display_columns_pt_ts)], tablefmt="grid"))
    
    # 7: Calculate overall consistency scores
    property_scores_pt_ph = {}
    for short in property_shorts_pt_ph:
        rel_col = f"{short}_diff_rel"
        if rel_col in df_pt_ph.columns:
            # Calculate percentage of points with relative difference < 0.5%
            good_points = (df_pt_ph[rel_col].abs() < 0.5).sum()
            acceptable_points = ((df_pt_ph[rel_col].abs() >= 0.5) & (df_pt_ph[rel_col].abs() < 2.0)).sum()
            
            if len(df_pt_ph) > 0:
                good_percent = good_points / len(df_pt_ph) * 100
                acceptable_percent = acceptable_points / len(df_pt_ph) * 100
                property_scores_pt_ph[short] = (good_percent, acceptable_percent)
    
    property_scores_pt_ts = {}
    for short in property_shorts_pt_ts:
        rel_col = f"{short}_diff_rel"
        if rel_col in df_pt_ts.columns:
            # Calculate percentage of points with relative difference < 0.5%
            good_points = (df_pt_ts[rel_col].abs() < 0.5).sum()
            acceptable_points = ((df_pt_ts[rel_col].abs() >= 0.5) & (df_pt_ts[rel_col].abs() < 2.0)).sum()
            
            if len(df_pt_ts) > 0:
                good_percent = good_points / len(df_pt_ts) * 100
                acceptable_percent = acceptable_points / len(df_pt_ts) * 100
                property_scores_pt_ts[short] = (good_percent, acceptable_percent)
    
    # 8: Evaluate overall sanity check result
    print("\nOverall Sanity Check Results:")
    
    # PT vs PH evaluation
    if 'T' in property_scores_pt_ph and 'D' in property_scores_pt_ph:
        t_good, t_acceptable = property_scores_pt_ph['T']
        d_good, d_acceptable = property_scores_pt_ph['D']
        
        if t_good > 80 and d_good > 80:
            pt_ph_result = "PASSED: PT-flash and PH-flash results are highly consistent."
        elif (t_good + t_acceptable) > 80 and (d_good + d_acceptable) > 80:
            pt_ph_result = "PARTIALLY PASSED: Some discrepancies exist but are within reasonable limits."
        else:
            pt_ph_result = "FAILED: Significant discrepancies between PT-flash and PH-flash results."
    else:
        # Fallback to simpler criteria if property scores can't be calculated
        try:
            t_diff_rel_mean = df_pt_ph['T_diff_rel'].abs().mean()
            d_diff_rel_mean = df_pt_ph['D_diff_rel'].abs().mean()
            
            if t_diff_rel_mean < 0.5 and d_diff_rel_mean < 0.5:
                pt_ph_result = "PASSED: PT-flash and PH-flash results are consistent."
            elif t_diff_rel_mean < 2.0 and d_diff_rel_mean < 2.0:
                pt_ph_result = "PARTIALLY PASSED: Some discrepancies exist but are within reasonable limits."
            else:
                pt_ph_result = "FAILED: Significant discrepancies between PT-flash and PH-flash results."
        except:
            pt_ph_result = "UNABLE TO EVALUATE: Missing or invalid comparison data."
    
    # PT vs TS evaluation
    if 'P' in property_scores_pt_ts and 'D' in property_scores_pt_ts:
        p_good, p_acceptable = property_scores_pt_ts['P']
        d_good, d_acceptable = property_scores_pt_ts['D']
        
        if p_good > 80 and d_good > 80:
            pt_ts_result = "PASSED: PT-flash and TS-flash results are highly consistent."
        elif (p_good + p_acceptable) > 80 and (d_good + d_acceptable) > 80:
            pt_ts_result = "PARTIALLY PASSED: Some discrepancies exist but are within reasonable limits."
        else:
            pt_ts_result = "FAILED: Significant discrepancies between PT-flash and TS-flash results."
    else:
        # Fallback to simpler criteria if property scores can't be calculated
        try:
            p_diff_rel_mean = df_pt_ts['P_diff_rel'].abs().mean()
            d_diff_rel_mean = df_pt_ts['D_diff_rel'].abs().mean()
            
            if p_diff_rel_mean < 0.5 and d_diff_rel_mean < 0.5:
                pt_ts_result = "PASSED: PT-flash and TS-flash results are consistent."
            elif p_diff_rel_mean < 2.0 and d_diff_rel_mean < 2.0:
                pt_ts_result = "PARTIALLY PASSED: Some discrepancies exist but are within reasonable limits."
            else:
                pt_ts_result = "FAILED: Significant discrepancies between PT-flash and TS-flash results."
        except:
            pt_ts_result = "UNABLE TO EVALUATE: Missing or invalid comparison data."
    
    # Show overall results
    print("\nPT-flash vs PH-flash: " + pt_ph_result)
    print("PT-flash vs TS-flash: " + pt_ts_result)
    
    # Overall cross-consistency result
    if pt_ph_result.startswith("PASSED") and pt_ts_result.startswith("PASSED"):
        overall_result = "PASSED"
        overall_message = "All three calculation methods (PT, PH, TS) are consistently producing the same thermodynamic properties."
    elif pt_ph_result.startswith("PARTIALLY") and pt_ts_result.startswith("PARTIALLY"):
        overall_result = "PARTIALLY PASSED"
        overall_message = "The three calculation methods show moderate consistency with acceptable discrepancies."
    elif "FAILED" in pt_ph_result or "FAILED" in pt_ts_result:
        overall_result = "FAILED"
        overall_message = "Significant discrepancies exist between the different calculation methods."
    else:
        overall_result = "INCONCLUSIVE"
        overall_message = "Could not properly evaluate consistency across all three calculation methods."
    
    print(f"\nOVERALL SANITY CHECK: {overall_result}")
    print(overall_message)
    
    # 9: Check for any specific property issues
    print("\nProperty-specific observations:")
    property_warnings = []
    
    # PT vs PH property issues
    for short, name in zip(property_shorts_pt_ph, property_names_pt_ph):
        rel_col = f"{short}_diff_rel"
        if rel_col in df_pt_ph.columns and not df_pt_ph[rel_col].isnull().all():
            mean_diff = df_pt_ph[rel_col].abs().mean()
            
            if mean_diff > 5.0:
                issue = f"- Warning: PT vs PH {name} shows large discrepancies (mean {mean_diff:.2f}%)"
                print(issue)
                property_warnings.append(issue)
            elif mean_diff > 2.0:
                issue = f"- Note: PT vs PH {name} shows moderate discrepancies (mean {mean_diff:.2f}%)"
                print(issue)
                property_warnings.append(issue)
    
    # PT vs TS property issues
    for short, name in zip(property_shorts_pt_ts, property_names_pt_ts):
        rel_col = f"{short}_diff_rel"
        if rel_col in df_pt_ts.columns and not df_pt_ts[rel_col].isnull().all():
            mean_diff = df_pt_ts[rel_col].abs().mean()
            
            if mean_diff > 5.0:
                issue = f"- Warning: PT vs TS {name} shows large discrepancies (mean {mean_diff:.2f}%)"
                print(issue)
                property_warnings.append(issue)
            elif mean_diff > 2.0:
                issue = f"- Note: PT vs TS {name} shows moderate discrepancies (mean {mean_diff:.2f}%)"
                print(issue)
                property_warnings.append(issue)
    
    reasons = [
        "1. Numerical precision in the solvers",
        "2. Iteration tolerance settings in REFPROP",
        "3. Different solution paths in PT vs PH vs TS flash calculations",
        "4. Properties near phase boundaries or critical region",
        "5. Implementation differences between flash algorithms",
        "6. Thermodynamic consistency of the underlying equation of state"
    ]
    
    print("\nPossible reasons for discrepancies:")
    for reason in reasons:
        print(reason)
    
    # Return results
    results = {
        "description": description,
        "composition": composition,
        "pt_ph_result": pt_ph_result,
        "pt_ts_result": pt_ts_result,
        "overall_result": overall_result,
        "overall_message": overall_message,
        "property_warnings": property_warnings,
        "pt_ph_stats": pt_ph_stats,
        "pt_ts_stats": pt_ts_stats,
        "raw_data": {
            "pt_ph_comparison": df_pt_ph,
            "pt_ts_comparison": df_pt_ts
        }
    }
    
    return results

def run_composition_tests():
    """
    Run multiple sanity checks with different compositions
    """
    # Define compositions to test
    compositions = [
        {
            "description": "90% CO2, 10% Methane",
            "components": [
                {"fluid": "CO2", "fraction": 0.90},
                {"fluid": "METHANE", "fraction": 0.10}
            ]
        },
        {
            "description": "80% CO2, 10% Methane, 10% Ethane",
            "components": [
                {"fluid": "CO2", "fraction": 0.80},
                {"fluid": "METHANE", "fraction": 0.10},
                {"fluid": "ETHANE", "fraction": 0.10}
            ]
        },
        {
            "description": "70% CO2, 10% Methane, 10% Ethane, 10% Propane",
            "components": [
                {"fluid": "CO2", "fraction": 0.70},
                {"fluid": "METHANE", "fraction": 0.10},
                {"fluid": "ETHANE", "fraction": 0.10},
                {"fluid": "PROPANE", "fraction": 0.10}
            ]
        },
        {
            "description": "60% CO2, 10% Methane, 10% Ethane, 10% Propane, 10% Nitrogen",
            "components": [
                {"fluid": "CO2", "fraction": 0.60},
                {"fluid": "METHANE", "fraction": 0.10},
                {"fluid": "ETHANE", "fraction": 0.10},
                {"fluid": "PROPANE", "fraction": 0.10},
                {"fluid": "NITROGEN", "fraction": 0.10}
            ]
        },
        {
            "description": "50% CO2, 10% Methane, 10% Ethane, 10% Propane, 10% Nitrogen, 10% Oxygen",
            "components": [
                {"fluid": "CO2", "fraction": 0.50},
                {"fluid": "METHANE", "fraction": 0.10},
                {"fluid": "ETHANE", "fraction": 0.10},
                {"fluid": "PROPANE", "fraction": 0.10},
                {"fluid": "NITROGEN", "fraction": 0.10},
                {"fluid": "OXYGEN", "fraction": 0.10}
            ]
        },
        # New test cases with 7-10 components
        {
            "description": "45% CO2, 10% Methane, 10% Ethane, 10% Propane, 10% Nitrogen, 10% Oxygen, 5% Argon",
            "components": [
                {"fluid": "CO2", "fraction": 0.45},
                {"fluid": "METHANE", "fraction": 0.10},
                {"fluid": "ETHANE", "fraction": 0.10},
                {"fluid": "PROPANE", "fraction": 0.10},
                {"fluid": "NITROGEN", "fraction": 0.10},
                {"fluid": "OXYGEN", "fraction": 0.10},
                {"fluid": "ARGON", "fraction": 0.05}
            ]
        },
        {
            "description": "40% CO2, 10% Methane, 10% Ethane, 10% Propane, 10% Nitrogen, 5% Oxygen, 5% Argon, 10% Butane",
            "components": [
                {"fluid": "CO2", "fraction": 0.40},
                {"fluid": "METHANE", "fraction": 0.10},
                {"fluid": "ETHANE", "fraction": 0.10},
                {"fluid": "PROPANE", "fraction": 0.10},
                {"fluid": "NITROGEN", "fraction": 0.10},
                {"fluid": "OXYGEN", "fraction": 0.05},
                {"fluid": "ARGON", "fraction": 0.05},
                {"fluid": "BUTANE", "fraction": 0.10}
            ]
        },
        {
            "description": "35% CO2, 10% Methane, 10% Ethane, 10% Propane, 5% Nitrogen, 5% Oxygen, 5% Argon, 10% Butane, 10% Hydrogen",
            "components": [
                {"fluid": "CO2", "fraction": 0.35},
                {"fluid": "METHANE", "fraction": 0.10},
                {"fluid": "ETHANE", "fraction": 0.10},
                {"fluid": "PROPANE", "fraction": 0.10},
                {"fluid": "NITROGEN", "fraction": 0.05},
                {"fluid": "OXYGEN", "fraction": 0.05},
                {"fluid": "ARGON", "fraction": 0.05},
                {"fluid": "BUTANE", "fraction": 0.10},
                {"fluid": "HYDROGEN", "fraction": 0.10}
            ]
        },
        {
            "description": "30% CO2, 10% Methane, 10% Ethane, 10% Propane, 5% Nitrogen, 5% Oxygen, 5% Argon, 5% Butane, 10% Hydrogen, 10% Helium",
            "components": [
                {"fluid": "CO2", "fraction": 0.30},
                {"fluid": "METHANE", "fraction": 0.10},
                {"fluid": "ETHANE", "fraction": 0.10},
                {"fluid": "PROPANE", "fraction": 0.10},
                {"fluid": "NITROGEN", "fraction": 0.05},
                {"fluid": "OXYGEN", "fraction": 0.05},
                {"fluid": "ARGON", "fraction": 0.05},
                {"fluid": "BUTANE", "fraction": 0.05},
                {"fluid": "HYDROGEN", "fraction": 0.10},
                {"fluid": "HELIUM", "fraction": 0.10}
            ]
        }
    ]
    
    # Run tests for each composition
    results = []
    for comp in compositions:
        print(f"\n\n{'#' * 100}")
        print(f"# Testing composition: {comp['description']}")
        print(f"{'#' * 100}")
        
        try:
            test_result = sanity_check_eos_with_composition(comp["components"], comp["description"])
            results.append(test_result)
        except Exception as e:
            print(f"Error running test for {comp['description']}: {str(e)}")
            # Add a failed result
            results.append({
                "description": comp["description"],
                "composition": comp["components"],
                "pt_ph_result": "ERROR: Test failed to complete",
                "pt_ts_result": "ERROR: Test failed to complete",
                "overall_result": "ERROR",
                "overall_message": f"Error running test: {str(e)}",
                "property_warnings": [],
                "pt_ph_stats": {},
                "pt_ts_stats": {},
                "raw_data": None
            })
    
    # Generate overall summary
    print(f"\n\n{'*' * 100}")
    print("* SUMMARY OF ALL COMPOSITION TESTS")
    print(f"{'*' * 100}")
    
    summary_table = []
    for result in results:
        co2_pct = next((c["fraction"] * 100 for c in result["composition"] if c["fluid"] == "CO2"), 0)
        num_components = len(result["composition"])
        
        # Extract key metrics for primary properties
        pt_ph_temp_diff = result["pt_ph_stats"].get("T", {}).get("rel_mean", float('nan'))
        pt_ph_dens_diff = result["pt_ph_stats"].get("D", {}).get("rel_mean", float('nan'))
        pt_ts_pres_diff = result["pt_ts_stats"].get("P", {}).get("rel_mean", float('nan'))
        pt_ts_dens_diff = result["pt_ts_stats"].get("D", {}).get("rel_mean", float('nan'))
        
        # Determine if there are any severe discrepancies
        severe_issues = []
        for stat_key, stat_dict in result["pt_ph_stats"].items():
            if stat_dict.get("rel_mean", 0) > 5.0:
                severe_issues.append(f"PT-PH {stat_key}")
        
        for stat_key, stat_dict in result["pt_ts_stats"].items():
            if stat_dict.get("rel_mean", 0) > 5.0:
                severe_issues.append(f"PT-TS {stat_key}")
        
        severe_str = ", ".join(severe_issues[:3])
        if len(severe_issues) > 3:
            severe_str += f" and {len(severe_issues) - 3} more"
        
        summary_table.append([
            f"{co2_pct:.0f}% CO2",
            num_components,
            result["description"],
            result["pt_ph_result"].split(":")[0],
            result["pt_ts_result"].split(":")[0],
            result["overall_result"],
            f"{pt_ph_temp_diff:.4f}%" if not np.isnan(pt_ph_temp_diff) else "N/A",
            f"{pt_ph_dens_diff:.4f}%" if not np.isnan(pt_ph_dens_diff) else "N/A",
            f"{pt_ts_pres_diff:.4f}%" if not np.isnan(pt_ts_pres_diff) else "N/A",
            f"{pt_ts_dens_diff:.4f}%" if not np.isnan(pt_ts_dens_diff) else "N/A",
            severe_str if severe_issues else "None"
        ])
    
    summary_headers = [
        "Composition", 
        "Components", 
        "Description", 
        "PT-PH Result", 
        "PT-TS Result", 
        "Overall Result",
        "T Diff (PT-PH)",
        "D Diff (PT-PH)",
        "P Diff (PT-TS)",
        "D Diff (PT-TS)",
        "Severe Issues"
    ]
    
    print("\nComposition Test Summary:")
    print(tabulate(summary_table, headers=summary_headers, tablefmt="grid"))
    
    # Draw conclusions about the best composition ranges
    print("\nConclusions:")
    
    # Count results by type
    passed_count = sum(1 for r in results if r["overall_result"] == "PASSED")
    partial_count = sum(1 for r in results if r["overall_result"] == "PARTIALLY PASSED")
    failed_count = sum(1 for r in results if r["overall_result"] == "FAILED")
    error_count = sum(1 for r in results if r["overall_result"] == "ERROR")
    
    print(f"- {passed_count} compositions PASSED the sanity check")
    print(f"- {partial_count} compositions PARTIALLY PASSED the sanity check")
    print(f"- {failed_count} compositions FAILED the sanity check")
    print(f"- {error_count} compositions had ERRORS during testing")
    
    # Analyze by CO2 percentage
    passed_co2 = [next((c["fraction"] * 100 for c in r["composition"] if c["fluid"] == "CO2"), 0) 
                  for r in results if r["overall_result"] == "PASSED"]
    
    if passed_co2:
        print(f"\nThe API performs best with compositions containing:")
        print(f"- CO2 percentages: {', '.join(f'{x:.0f}%' for x in passed_co2)}")
        
        # Component count analysis
        passed_comp_counts = [len(r["composition"]) for r in results if r["overall_result"] == "PASSED"]
        if passed_comp_counts:
            min_comps = min(passed_comp_counts)
            max_comps = max(passed_comp_counts)
            comp_range = f"{min_comps}" if min_comps == max_comps else f"{min_comps}-{max_comps}"
            print(f"- Component count: {comp_range}")
    
    # Recommendations
    print("\nRecommendations:")
    if passed_count > 0:
        print("- For most accurate results, use compositions that passed the sanity check")
    elif partial_count > 0:
        print("- Use compositions that at least partially pass the sanity check")
        print("- Exercise caution with calculations that involve transport properties, which tend to show greater inconsistencies")
    else:
        print("- The API appears to have significant consistency issues with all tested compositions")
        print("- Consider simplifying to binary mixtures or improving the underlying modeling approach")
    
    print("\nFor future work:")
    print("- Implement validation checks in the API to warn users about potential inconsistencies")
    print("- Improve the numerical methods for complex mixtures, especially for transport properties")
    print("- Consider adding composition complexity limitations in the API documentation")
    
    return results

if __name__ == "__main__":
    try:
        results = run_composition_tests()
        print("\nAll composition tests completed.")
    except Exception as e:
        print(f"\nError during composition testing: {str(e)}")