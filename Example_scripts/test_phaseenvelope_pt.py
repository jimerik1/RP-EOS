#!/usr/bin/env python3
import requests
import matplotlib.pyplot as plt
import numpy as np
import time

def query_phase_envelope_pt(n2_fraction, T_min=-50, T_max=50, dT=5):
    """
    Sends a POST request to your local REFPROP-based endpoint that computes 
    the PT phase envelope (bubble + dew) for a CO2/N2 mixture.

    :param n2_fraction: Molar fraction of N2 in the mixture (0 to 1).
    :param T_min: Minimum temperature in Celsius (for the calculation).
    :param T_max: Maximum temperature in Celsius (for the calculation).
    :param dT: Temperature step in Celsius for the envelope calculation.
    :return: (bubble_curve, dew_curve) data from the JSON response, or None if error.
    """
    url = "http://127.0.0.1:5051/phase_envelope_pt"

    # Ensure we pass the correct composition (rest = CO2)
    co2_fraction = 1.0 - n2_fraction
    # Ensure proper rounding to avoid floating point issues
    co2_fraction = round(co2_fraction, 6)
    n2_fraction = round(n2_fraction, 6)
    
    payload = {
        "composition": [
            {"fluid": "CO2", "fraction": co2_fraction},
            {"fluid": "N2",  "fraction": n2_fraction}
        ],
        "temperature_range": {
            "from": T_min,
            "to":   T_max
        },
        "temperature_resolution": dT,
        "desired_curve": "both"
    }

    try:
        print(f"Querying phase envelope for xN2={n2_fraction:.3f}...")
        response = requests.post(url, json=payload, timeout=60)
        response.raise_for_status()  # will raise an HTTPError if the status is 4xx,5xx
        data = response.json()
        
        bubble_curve = data.get("bubble_curve", [])
        dew_curve = data.get("dew_curve", [])
        
        print(f"  Got {len(bubble_curve)} bubble points and {len(dew_curve)} dew points")
        return bubble_curve, dew_curve
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] Request failed for n2_fraction={n2_fraction}: {e}")
        if hasattr(e, 'response') and e.response is not None:
            try:
                error_info = e.response.json()
                print(f"Server error details: {error_info}")
            except:
                print(f"Status code: {e.response.status_code}, Content: {e.response.text[:200]}...")
        return None


def main():
    # Hardcode a single composition - e.g., 5% N2, 95% CO2
    n2_fraction = 0.1
    
    # Prepare a figure
    fig, ax = plt.subplots(figsize=(10,7))
    
    print(f"Calculating phase envelope for CO2/N2 mixture with xN2={n2_fraction:.3f}...")
    result = query_phase_envelope_pt(n2_fraction, T_min=-50, T_max=50, dT=1)
    if not result:
        print("Failed to calculate phase envelope.")
        return

    bubble_curve, dew_curve = result

    # Convert the results into arrays we can plot
    Tbub = [pt["temperature"] for pt in bubble_curve]
    Pbub = [pt["pressure"] for pt in bubble_curve]

    Tdew = [pt["temperature"] for pt in dew_curve]
    Pdew = [pt["pressure"] for pt in dew_curve]

    # Plot bubble and dew lines with different styles
    ax.plot(Tbub, Pbub, 'o-', color='blue', markersize=5, label="Bubble Line")
    ax.plot(Tdew, Pdew, 's--', color='red', markersize=5, label="Dew Line")

    # Set axis labels and title
    ax.set_xlabel("Temperature [°C]", fontsize=12)
    ax.set_ylabel("Pressure [bar]", fontsize=12)
    ax.set_title(f"CO2/N2 Phase Envelope (xN2={n2_fraction:.3f})", fontsize=14)
    
    # Set y-axis to logarithmic scale
    ax.set_yscale('log')
    
    # Add grid and legend
    ax.grid(True, which='both', linestyle='--', alpha=0.6)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'co2_n2_phase_envelope_{int(n2_fraction*100)}pct.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    main()