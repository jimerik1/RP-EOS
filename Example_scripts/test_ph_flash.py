import requests
import json
import matplotlib.pyplot as plt
import numpy as np

# Base URL where the API is exposed
BASE_URL = "http://localhost:5051"

def test_ph_flash():
    """Test the /ph_flash endpoint with CO2."""
    url = f"{BASE_URL}/ph_flash"
    payload = {
        "composition": [
            {"fluid": "CO2", "fraction": 1.0}
        ],
        "pressure_range": {
            "from": 10,    # in bar
            "to": 100
        },
        "enthalpy_range": {
            "from": 10000,   # in J/mol
            "to": 25000
        },
        "pressure_resolution": 5,   # Step size in bar
        "enthalpy_resolution": 1000, # Step size in J/mol
        "properties": [
            "temperature",
            "density",
            "vapor_fraction",
            "phase",
            "entropy",
            "internal_energy"
        ],
        "units_system": "SI"
    }
    
    print("Sending request to PH flash endpoint...")
    response = requests.post(url, json=payload)
    print("Status code:", response.status_code)
    
    if response.status_code == 200:
        data = response.json()
        if "results" in data:
            results = data["results"]
            print(f"Received {len(results)} data points")
            print("First result:", json.dumps(results[0], indent=2))
            
            # Extract data for plotting
            temperatures = [point["temperature"]["value"] for point in results]
            pressures = [point["pressure"]["value"] for point in results]
            enthalpies = [point["enthalpy"]["value"] for point in results]
            vapor_fractions = [point["vapor_fraction"]["value"] for point in results]
            
            # Plot a simple P-h diagram with temperature as color
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Define colors for different phases
            colors = []
            for vf in vapor_fractions:
                if vf == 0:  # Liquid
                    colors.append('blue')
                elif vf == 1:  # Vapor
                    colors.append('red')
                elif 0 < vf < 1:  # Two-phase
                    colors.append('green')
                else:  # Supercritical
                    colors.append('purple')
            
            # Create a scatter plot
            sc = ax.scatter(enthalpies, pressures, c=temperatures, cmap='coolwarm')
            
            ax.set_xlabel('Enthalpy (J/mol)')
            ax.set_ylabel('Pressure (bar)')
            ax.set_title('P-h Diagram with Temperature Contours for CO2')
            
            # Create colorbar from the scatter plot - use the figure and axes explicitly
            cbar = fig.colorbar(sc, ax=ax)
            cbar.set_label('Temperature (°C)')
            
            ax.grid(True)
            
            # Save the plot
            plt.savefig('ph_flash_test.png')
            plt.show()
            
            print("Test completed successfully!")
        else:
            print("No results in response:", data)
    else:
        print("Error response:", response.text)

if __name__ == '__main__':
    test_ph_flash()