import numpy as np
from typing import Any

def get_phase(q):
    """Determine the phase of the fluid based on quality value."""
    if q == 0:
        return 'Liquid'
    elif q == 1:
        return 'Vapor'
    elif q == -998:
        return 'Liquid'  # Subcooled liquid
    elif q == 998:
        return 'Vapor'   # Superheated vapor
    elif q == 999:
        return 'Supercritical'
    elif q == -999:
        return 'Solid'   # Added support for solid phase
    elif q == -997:
        return 'Solid-Liquid'  # Solid-liquid equilibrium (melting)
    elif q == 997:
        return 'Solid-Vapor'   # Solid-vapor equilibrium (sublimation)
    elif q == 996:
        return 'Triple-Point'  # Triple point (solid-liquid-vapor)
    elif 0 < q < 1:
        return 'Two-Phase'
    elif q < 0:
        return 'Liquid'  # Subcooled liquid
    else:
        return 'Vapor'   # Superheated vapor

def validate_composition(composition):
    """Validate composition data"""
    if not composition:
        return False
    total = sum(comp.get('fraction', 0) for comp in composition)
    return abs(total - 1.0) < 1e-6

def trim(s: bytes) -> str:
    """Trim NULL characters and decode."""
    return s.replace(b'\x00', b'').strip().decode("utf-8")

def convert_for_json(obj):
    """Convert numpy types to JSON serializable Python types."""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f'Object of type {type(obj)} is not JSON serializable')