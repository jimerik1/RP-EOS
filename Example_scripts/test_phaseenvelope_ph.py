#!/usr/bin/env python3
import requests
import matplotlib.pyplot as plt
import numpy as np
import time

def query_phase_envelope_ph(n2_fraction, p_min=5, p_max=150, dp=5):
    """
    Sends a POST request to get the PH phase envelope (bubble + dew) for a CO2/N2 mixture.

    :param n2_fraction: Molar fraction of N2 in the mixture (0 to 1)
    :param p_min: Minimum pressure in bar
    :param p_max: Maximum pressure in bar
    :param dp: Pressure step in bar
    :return: (bubble_curve, dew_curve) data from the JSON response, or None if error
    """
    url = "http://127.0.0.1:5051/phase_envelope_ph"

    # Ensure proper composition fractions (avoid floating point issues)
    co2_fraction = round(1.0 - n2_fraction, 6)
    n2_fraction = round(n2_fraction, 6)
    
    payload = {
        "composition": [
            {"fluid": "CO2", "fraction": co2_fraction},
            {"fluid": "N2",  "fraction": n2_fraction}
        ],
        "pressure_range": {
            "from": p_min,
            "to": p_max
        },
        "pressure_resolution": dp,
        "desired_curve": "both"
    }

    try:
        print(f"Querying PH phase envelope for xN2={n2_fraction:.3f}...")
        response = requests.post(url, json=payload, timeout=60)
        response.raise_for_status()
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
    # Using a fixed N2 fraction of 5%
    n2_fraction = 0.05  # 5% nitrogen, 95% CO2
    
    # Pressure range for the envelope
    p_min = 5  # bar
    p_max = 150  # bar
    dp = 0.2  # bar

    # Set up the figure
    fig, ax = plt.subplots(figsize=(10, 7))
    
    print(f"Creating PH phase envelope for CO2 with {n2_fraction*100:.1f}% N2...")
    
    result = query_phase_envelope_ph(n2_fraction, p_min, p_max, dp)
    if not result:
        print("Failed to get phase envelope data.")
        return

    bubble_curve, dew_curve = result
    
    # Extract pressure and enthalpy values
    h_bubble = [pt["enthalpy"] for pt in bubble_curve]
    p_bubble = [pt["pressure"] for pt in bubble_curve]
    
    h_dew = [pt["enthalpy"] for pt in dew_curve]
    p_dew = [pt["pressure"] for pt in dew_curve]
    
    # Plot the envelope
    ax.plot(h_bubble, p_bubble, 'o-', color='blue', markersize=5, 
            label=f"Bubble Line (Liquid)", linewidth=2)
    ax.plot(h_dew, p_dew, 's--', color='red', markersize=5,
            label=f"Dew Line (Vapor)", linewidth=2)
    
    # Add temperature labels
    if len(bubble_curve) > 0 and len(dew_curve) > 0:
        # Add temperature annotations at regular intervals
        for curve, points, marker, color, offset in [
            (bubble_curve, 5, 'left', 'blue', -8), 
            (dew_curve, 5, 'right', 'red', 8)
        ]:
            if len(curve) >= points:
                step = max(1, len(curve) // points)
                for j in range(0, len(curve), step):
                    pt = curve[j]
                    temp = pt["temperature"]
                    h = pt["enthalpy"]
                    p = pt["pressure"]
                    ax.annotate(f"{temp:.1f}°C", 
                               xy=(h, p),
                               xytext=(offset, 0), 
                               textcoords="offset points",
                               fontsize=9, color=color, 
                               ha=marker, va='center',
                               bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec=color))

    # Configure the plot
    ax.set_xlabel("Enthalpy [J/mol]", fontsize=12)
    ax.set_ylabel("Pressure [bar]", fontsize=12)
    ax.set_title(f"CO2 with {n2_fraction*100:.1f}% N2 Phase Envelope (P-H Diagram)", fontsize=14)
    
    # Set y-axis to logarithmic scale
    ax.set_yscale('log')
    
    # Add grid
    ax.grid(True, which='both', linestyle='--', alpha=0.6)
    
    # Pad the enthalpy range
    h_min = min(min(h_bubble or [0]), min(h_dew or [0]))
    h_max = max(max(h_bubble or [0]), max(h_dew or [0]))
    h_range = h_max - h_min
    ax.set_xlim([h_min - 0.05 * h_range, h_max + 0.05 * h_range])
    
    # Add legend with clear labels
    ax.legend(loc='best')
    
    # Add annotation explaining the diagram
    ax.text(0.02, 0.02, 
            "Two-phase region exists between bubble and dew lines.\n"
            "Left of bubble line: subcooled liquid\n"
            "Right of dew line: superheated vapor", 
            transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.8))
    
    plt.tight_layout()
    # Save the figure
    output_filename = f'co2_n2_{int(n2_fraction*100)}pct_ph_phase_envelope.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Phase envelope saved to {output_filename}")
    plt.show()


if __name__ == "__main__":
    main()