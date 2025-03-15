import requests
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Patch
import time
import sys

# -----------------------------
# 1) API endpoints & payloads
# -----------------------------
phase_boundaries_url = "http://localhost:5051/phase_boundaries"
critical_point_url    = "http://localhost:5051/critical_point"

phase_boundaries_payload = {
    "composition": [
        {"fluid": "CO2", "fraction": 1.0},
    ],
    "variables": {
        "temperature": {
            "range": {"from": -140, "to": 60},
            "resolution": 1
        }
    },
    "calculation": {
        "boundary_types": ["melting", "sublimation", "vaporization"],
        "units_system": "SI"  # or "SI" or "ENG", depending on your API
    }
}

critical_point_payload = {
    "composition": [
        {"fluid": "CO2", "fraction": 1.0}
    ],
    "units_system": "SI"
}

# -----------------------------
# 2) Helper: GET data w/ retries
# -----------------------------
def get_fluid_data(url, payload, max_retries=3, retry_delay=2):
    for attempt in range(max_retries):
        try:
            print(f"Request {attempt+1}/{max_retries} to {url}...")
            response = requests.post(url, json=payload, timeout=120)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")
            if attempt < max_retries - 1:
                print(f"Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
            else:
                print(f"Max retries reached for {url}.")
                return None

# -----------------------------
# 3) Fetch Phase Boundaries & Critical Point
# -----------------------------
print("Fetching CO2 phase boundaries...")
phase_data = get_fluid_data(phase_boundaries_url, phase_boundaries_payload)
print("Fetching CO2 critical point data...")
cp_data = get_fluid_data(critical_point_url, critical_point_payload)

if not phase_data:
    print("Failed to get phase boundary data. Exiting.")
    sys.exit(1)

melting_curve = phase_data.get('melting_curve', [])
sublimation_curve = phase_data.get('sublimation_curve', [])
vaporization_curve = phase_data.get('vaporization_curve', [])
triple_point_data = phase_data.get('triple_point', None)

# Default triple point if not returned
Ttp_C = -56.6
Ptp_atm = 5.11

if triple_point_data:
    # temperature already in °C
    Ttp_C = triple_point_data["temperature"]["value"]
    # convert bar -> atm
    Ptp_atm = triple_point_data["pressure"]["value"] * 0.986923

critical_point = None
if cp_data and 'critical_point' in cp_data:
    critical_point = cp_data['critical_point']

# Default critical point if not returned
Tcrit_C = 31.0
Pcrit_atm = 72.8

if critical_point:
    Tcrit_C = critical_point['critical_temperature']['value'] - 273.15  # K -> °C
    Pcrit_atm = critical_point['critical_pressure']['value'] * 0.986923 # bar -> atm

# -----------------------------
# 4) Convert each boundary to (T, P) in lists, sorted by T
# -----------------------------
def extract_sorted_TP(curve):
    if not curve:
        return [], []
    T = [pt["temperature"]["value"] for pt in curve]
    P = [pt["pressure"]["value"] for pt in curve]  # in bar
    # sort by T
    sorted_pts = sorted(zip(T, P), key=lambda x: x[0])
    T_sorted, P_sorted_bar = zip(*sorted_pts)
    # convert bar->atm
    P_sorted_atm = [pb * 0.986923 for pb in P_sorted_bar]
    return list(T_sorted), list(P_sorted_atm)

T_sub, P_sub = extract_sorted_TP(sublimation_curve)
T_melt, P_melt = extract_sorted_TP(melting_curve)
T_vap, P_vap = extract_sorted_TP(vaporization_curve)

# -----------------------------
# 5) Set up the Plot
# -----------------------------
plt.figure(figsize=(12, 9), dpi=100)
ax = plt.subplot(111)
ax.set_yscale('log')

temp_min, temp_max = -140, 60
press_min, press_max = 1e-3, 1e4

plt.xlim(temp_min, temp_max)
plt.ylim(press_min, press_max)

plt.xlabel('Temperature (°C)', fontsize=14, fontweight='bold')
plt.ylabel('Pressure (atm)', fontsize=14, fontweight='bold')
plt.title('Carbon Dioxide (CO₂) Phase Diagram', fontsize=18, fontweight='bold', pad=20)

plt.grid(True, which='both', linestyle='--', alpha=0.3)
formatter = ScalarFormatter()
formatter.set_scientific(False)
ax.yaxis.set_minor_formatter(formatter)

# Colors
color_solid       = '#e63946'
color_liquid      = '#4895ef'
color_gas         = '#999999'
color_supercrit   = '#90e0ef'

# -----------------------------
# 6) Build & fill Polygons
# -----------------------------

# -------------------------------------
# FILL THE SOLID REGION CORRECTLY
# -------------------------------------
# We want all T < Ttp.  The "lower" boundary is the sublimation curve,
# the "upper" boundary is the melting curve.  That means:
#
#  (A) gather all sub-curve points with T <= Ttp
#  (B) gather all melt-curve points with T <= Ttp
#  (C) build a polygon by going up the sub-curve in ascending T,
#      then back along the melt-curve in descending T,
#      ensuring the triple point is included in both.
#
# If your sub-curve or melt-curve does not start as far left as temp_min,
# you can optionally "clamp" the polygon to (temp_min, some pressure).
# But typically you'd at least start from sub_x[0].

solid_poly_x = []
solid_poly_y = []

# 1) Sublimation boundary up to Ttp
subX = []
subY = []
for tx, py in zip(T_sub, P_sub):
    if tx <= Ttp_C:
        subX.append(tx)
        subY.append(py)
# Ensure we end exactly at Ttp
if subX and subX[-1] < Ttp_C:
    subX.append(Ttp_C)
    subY.append(Ptp_atm)

# 2) Melting boundary up to Ttp
meltX = []
meltY = []
for tx, py in zip(T_melt, P_melt):
    if tx <= Ttp_C:
        meltX.append(tx)
        meltY.append(py)
# Ensure we end exactly at Ttp
if meltX and meltX[-1] < Ttp_C:
    meltX.append(Ttp_C)
    meltY.append(Ptp_atm)

# 3) Build the polygon
#    ascending T along the sublimation curve, 
#    then descending T along the melting curve.
if subX and meltX:
    # We assume subX is already in ascending T order,
    # and meltX is in ascending T order too,
    # so we reverse meltX to come back in descending T
    solid_poly_x = subX + list(reversed(meltX))
    solid_poly_y = subY + list(reversed(meltY))

    plt.fill(solid_poly_x, solid_poly_y, color='#e63946', alpha=0.8, zorder=1)
    
# B) GAS region
# from (temp_min, press_min) up sublimation to triple point,
# then vaporization up to critical, then across to (temp_max, Pcrit_atm),
# then down to (temp_max, press_min), back to (temp_min, press_min).
gas_poly_x = [temp_min,]
gas_poly_y = [press_min,]

sub_x2, sub_y2 = [], []
for tx, py in zip(T_sub, P_sub):
    if tx <= Ttp_C:
        sub_x2.append(tx)
        sub_y2.append(py)

if sub_x2:
    gas_poly_x.extend(sub_x2)
    gas_poly_y.extend(sub_y2)
    if sub_x2[-1] < Ttp_C:
        gas_poly_x.append(Ttp_C)
        gas_poly_y.append(Ptp_atm)

# vaporization from triple point to critical
vap_x = []
vap_y = []
for tx, py in zip(T_vap, P_vap):
    if Ttp_C <= tx <= Tcrit_C:
        vap_x.append(tx)
        vap_y.append(py)
if vap_x and vap_x[0] > Ttp_C:
    vap_x.insert(0, Ttp_C)
    vap_y.insert(0, Ptp_atm)
if vap_x and vap_x[-1] < Tcrit_C:
    vap_x.append(Tcrit_C)
    vap_y.append(Pcrit_atm)
gas_poly_x.extend(vap_x)
gas_poly_y.extend(vap_y)

# across to (temp_max, Pcrit_atm)
gas_poly_x.append(temp_max)
gas_poly_y.append(Pcrit_atm)
# down to (temp_max, press_min)
gas_poly_x.append(temp_max)
gas_poly_y.append(press_min)
# close
gas_poly_x.append(temp_min)
gas_poly_y.append(press_min)

plt.fill(gas_poly_x, gas_poly_y, color=color_gas, alpha=0.8, zorder=0)

# C) LIQUID region
# bounded by melting curve + vaporization curve from Ttp to Tcrit.
liquid_poly_x = []
liquid_poly_y = []

# melting portion from Ttp -> Tcrit
melt_x2 = []
melt_y2 = []
for tx, py in zip(T_melt, P_melt):
    if Ttp_C <= tx <= Tcrit_C:
        melt_x2.append(tx)
        melt_y2.append(py)
if melt_x2 and melt_x2[0] > Ttp_C:
    melt_x2.insert(0, Ttp_C)
    melt_y2.insert(0, Ptp_atm)
if melt_x2 and melt_x2[-1] < Tcrit_C:
    melt_x2.append(Tcrit_C)
    melt_y2.append(Pcrit_atm)

liquid_poly_x.extend(melt_x2)
liquid_poly_y.extend(melt_y2)

# vaporization portion in reverse from Tcrit -> Ttp
vap_x2 = []
vap_y2 = []
for tx, py in zip(T_vap, P_vap):
    if Ttp_C <= tx <= Tcrit_C:
        vap_x2.append(tx)
        vap_y2.append(py)
if vap_x2 and vap_x2[-1] < Tcrit_C:
    vap_x2.append(Tcrit_C)
    vap_y2.append(Pcrit_atm)
# Reverse so we go from Tcrit down to Ttp
liquid_poly_x.extend(reversed(vap_x2))
liquid_poly_y.extend(reversed(vap_y2))

plt.fill(liquid_poly_x, liquid_poly_y, color=color_liquid, alpha=0.8, zorder=2)

# D) SUPERCRITICAL region
sc_x = [Tcrit_C, temp_max, temp_max, Tcrit_C]
sc_y = [Pcrit_atm, Pcrit_atm, press_max, press_max]
plt.fill(sc_x, sc_y, color=color_supercrit, alpha=0.8, zorder=3)

# -----------------------------
# 7) Plot boundary lines
# -----------------------------
plt.plot(T_sub, P_sub, 'k-', linewidth=2)   # Sublimation
plt.plot(T_melt, P_melt, 'k-', linewidth=2) # Melting
plt.plot(T_vap, P_vap, 'k-', linewidth=2)   # Vaporization

# Mark triple & critical
plt.plot(Ttp_C, Ptp_atm, 'wo', markeredgecolor='k', markersize=8)
plt.plot(Tcrit_C, Pcrit_atm, 'ko', markersize=8)

# Annotation examples
plt.text(-120, 1e-2, "Solid phase", ha='center', va='center',
         color='white', fontsize=12, fontweight='bold')
plt.text(-25, 10, "Liquid\nphase", ha='center', va='center',
         color='white', fontsize=12, fontweight='bold')
plt.text(0, 1e-2, "Gas\nphase", ha='center', va='center',
         color='white', fontsize=12, fontweight='bold')
plt.text(45, 100, "Supercritical\nphase", ha='center', va='center',
         color='black', fontsize=12, fontweight='bold')

# Add dashed lines from critical point
plt.plot([Tcrit_C, Tcrit_C], [Pcrit_atm, press_max], 'k--', alpha=0.7)
plt.plot([Tcrit_C, temp_max], [Pcrit_atm, Pcrit_atm], 'k--', alpha=0.7)

# Legend
legend_patches = [
    Patch(facecolor=color_solid,     label='Solid'),
    Patch(facecolor=color_liquid,    label='Liquid'),
    Patch(facecolor=color_gas,       label='Gas'),
    Patch(facecolor=color_supercrit, label='Supercritical'),
]
plt.legend(handles=legend_patches, loc='lower right', fontsize=12)

plt.tight_layout()
plt.savefig('co2_phase_diagram_filled.png', dpi=300, bbox_inches='tight')
plt.show()