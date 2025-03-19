"""
Configuration for OLGA TAB format property mappings and behavior.
This file defines how REFPROP properties map to OLGA TAB properties,
including unit conversions and fallback strategies.
"""

# Define the minimum set of properties required for OLGA TAB format
OLGA_REQUIRED_PROPERTIES = [
    "temperature", "pressure", "density", "liquid_density", "vapor_density",
    "viscosity", "thermal_conductivity", "enthalpy", "entropy", 
    "cp", "vapor_fraction", "surface_tension", "dDdP", "dDdT", 
    "phase", "x", "y"
]

# Define standard OLGA TAB properties and their mappings to REFPROP properties
OLGA_PROPERTY_MAPPINGS = [
    {
        'name': 'GAS DENSITY (KG/M3)',
        'key': 'vapor_density',
        'fallbacks': ['density'],
        'condition': lambda phase: phase in ['vapor', 'two-phase', 'supercritical'],
        'converter': lambda x, wmm: x * wmm if x is not None else 0.0,  # mol/L * g/mol = g/L = kg/m³
        'description': 'Density of the vapor phase'
    },
    {
        'name': 'LIQUID DENSITY (KG/M3)',
        'key': 'liquid_density',
        'fallbacks': ['density'],
        'condition': lambda phase: phase in ['liquid', 'two-phase'],
        'converter': lambda x, wmm: x * wmm if x is not None else 0.0,  # mol/L * g/mol = g/L = kg/m³
        'description': 'Density of the liquid phase'
    },
    {
        'name': 'WATER DENSITY (KG/M3)',
        'key': 'water_density',
        'fallbacks': [],
        'condition': lambda phase: True,  # Always include, zero if not applicable
        'converter': lambda x, wmm: x if x is not None else 0.0,
        'description': 'Density of water component'
    },
    {
        'name': 'DRHOG/DP (S2/M2)',
        'key': 'dDdP_vapor',
        'fallbacks': ['dDdP'],
        'condition': lambda phase: phase in ['vapor', 'two-phase'],
        'converter': lambda x, wmm: x * wmm if x is not None else 0.0,  # (mol/L)/kPa * g/mol = (kg/m³)/kPa
        'description': 'Derivative of gas density with respect to pressure'
    },
    {
        'name': 'DRHOL/DP (S2/M2)',
        'key': 'dDdP_liquid',
        'fallbacks': ['dDdP'],
        'condition': lambda phase: phase in ['liquid', 'two-phase'],
        'converter': lambda x, wmm: x * wmm if x is not None else 0.0,  # (mol/L)/kPa * g/mol = (kg/m³)/kPa
        'description': 'Derivative of liquid density with respect to pressure'
    },
    {
        'name': 'DRHOG/DT (KG/M3/K)',
        'key': 'dDdT_vapor',
        'fallbacks': ['dDdT'],
        'condition': lambda phase: phase in ['vapor', 'two-phase'],
        'converter': lambda x, wmm: x * wmm if x is not None else 0.0,  # (mol/L)/K * g/mol = (kg/m³)/K
        'description': 'Derivative of gas density with respect to temperature'
    },
    {
        'name': 'DRHOL/DT (KG/M3/K)',
        'key': 'dDdT_liquid',
        'fallbacks': ['dDdT'],
        'condition': lambda phase: phase in ['liquid', 'two-phase'],
        'converter': lambda x, wmm: x * wmm if x is not None else 0.0,  # (mol/L)/K * g/mol = (kg/m³)/K
        'description': 'Derivative of liquid density with respect to temperature'
    },
    {
        'name': 'GAS MASS FRACTION (-)',
        'key': 'vapor_fraction',
        'fallbacks': [],
        'condition': lambda phase: True,  # Always include
        'converter': lambda x, wmm: x if x is not None and 0 <= x <= 1 else (1.0 if phase == 'vapor' else 0.0),
        'description': 'Gas mass fraction'
    },
    {
        'name': 'GAS VISCOSITY (NS/M2)',
        'key': 'vapor_viscosity',
        'fallbacks': ['viscosity'],
        'condition': lambda phase: phase in ['vapor', 'two-phase'],
        'converter': lambda x, wmm: x * 1e-3 if x is not None else 0.0,  # μPa·s * 1e-3 = mPa·s = N·s/m²
        'description': 'Gas dynamic viscosity'
    },
    {
        'name': 'LIQUID VISCOSITY (NS/M2)',
        'key': 'liquid_viscosity',
        'fallbacks': ['viscosity'],
        'condition': lambda phase: phase in ['liquid', 'two-phase'],
        'converter': lambda x, wmm: x * 1e-3 if x is not None else 0.0,  # μPa·s * 1e-3 = mPa·s = N·s/m²
        'description': 'Liquid dynamic viscosity'
    },
    {
        'name': 'GAS HEAT CAPACITY (J/KG/K)',
        'key': 'vapor_cp',
        'fallbacks': ['cp'],
        'condition': lambda phase: phase in ['vapor', 'two-phase'],
        'converter': lambda x, wmm: x * 1000.0 / wmm if x is not None else 0.0,  # J/(mol·K) * 1000 / (g/mol) = J/(kg·K)
        'description': 'Gas specific heat capacity at constant pressure'
    },
    {
        'name': 'LIQUID HEAT CAPACITY (J/KG/K)',
        'key': 'liquid_cp',
        'fallbacks': ['cp'],
        'condition': lambda phase: phase in ['liquid', 'two-phase'],
        'converter': lambda x, wmm: x * 1000.0 / wmm if x is not None else 0.0,  # J/(mol·K) * 1000 / (g/mol) = J/(kg·K)
        'description': 'Liquid specific heat capacity at constant pressure'
    },
    {
        'name': 'GAS ENTHALPY (J/KG)',
        'key': 'vapor_enthalpy',
        'fallbacks': ['enthalpy'],
        'condition': lambda phase: phase in ['vapor', 'two-phase'],
        'converter': lambda x, wmm: x * 1000.0 / wmm if x is not None else 0.0,  # J/mol * 1000 / (g/mol) = J/kg
        'description': 'Gas specific enthalpy'
    },
    {
        'name': 'LIQUID ENTHALPY (J/KG)',
        'key': 'liquid_enthalpy',
        'fallbacks': ['enthalpy'],
        'condition': lambda phase: phase in ['liquid', 'two-phase'],
        'converter': lambda x, wmm: x * 1000.0 / wmm if x is not None else 0.0,  # J/mol * 1000 / (g/mol) = J/kg
        'description': 'Liquid specific enthalpy'
    },
    {
        'name': 'GAS THERMAL CONDUCTIVITY (W/M/K)',
        'key': 'vapor_thermal_conductivity',
        'fallbacks': ['thermal_conductivity'],
        'condition': lambda phase: phase in ['vapor', 'two-phase'],
        'converter': lambda x, wmm: x if x is not None else 0.0,
        'description': 'Gas thermal conductivity'
    },
    {
        'name': 'LIQUID THERMAL CONDUCTIVITY (W/M/K)',
        'key': 'liquid_thermal_conductivity',
        'fallbacks': ['thermal_conductivity'],
        'condition': lambda phase: phase in ['liquid', 'two-phase'],
        'converter': lambda x, wmm: x if x is not None else 0.0,
        'description': 'Liquid thermal conductivity'
    },
    {
        'name': 'SURFACE TENSION (N/M)',
        'key': 'surface_tension',
        'fallbacks': [],
        'condition': lambda phase: phase == 'two-phase',
        'converter': lambda x, wmm: x if x is not None else 0.0,
        'description': 'Surface tension between liquid and vapor phases'
    },
    {
        'name': 'GAS ENTROPY (J/KG/C)',
        'key': 'vapor_entropy',
        'fallbacks': ['entropy'],
        'condition': lambda phase: phase in ['vapor', 'two-phase'],
        'converter': lambda x, wmm: x * 1000.0 / wmm if x is not None else 0.0,  # J/(mol·K) * 1000 / (g/mol) = J/(kg·K)
        'description': 'Gas specific entropy'
    },
    {
        'name': 'LIQUID ENTROPY (J/KG/C)',
        'key': 'liquid_entropy',
        'fallbacks': ['entropy'],
        'condition': lambda phase: phase in ['liquid', 'two-phase'],
        'converter': lambda x, wmm: x * 1000.0 / wmm if x is not None else 0.0,  # J/(mol·K) * 1000 / (g/mol) = J/(kg·K)
        'description': 'Liquid specific entropy'
    }
]

# Define phase mappings between REFPROP phase strings and numeric values for OLGA
PHASE_MAPPING = {
    'liquid': 0.0,
    'vapor': 1.0,
    'two-phase': 0.5,  # This will typically use the actual vapor fraction
    'supercritical': 1.0,  # Treat as vapor in OLGA
    'unknown': 0.5  # Default to two-phase if unknown
}

# Special value indicators
SPECIAL_VALUES = {
    'not_applicable': 9980.0,  # Value to use when property is not applicable
    'error_value': -2426408E+08  # Value to use when calculation fails
}

# Default options for OLGA TAB formatting
DEFAULT_OLGA_OPTIONS = {
    'values_per_line': 5,  # Number of values per line in the TAB file
    'indent_spaces': 5,    # Number of spaces for indentation
    'debug_level': 1,      # Debug level: 0=none, 1=warnings, 2=info
    'use_fallbacks': True  # Whether to use fallback calculations for missing values
}