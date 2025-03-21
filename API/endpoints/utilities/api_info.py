from flask import Blueprint, jsonify, request
import sys
import traceback
from pathlib import Path
import os

# Create a blueprint for the api_info endpoint
api_info_bp = Blueprint('api_info', __name__)

# Dictionary of supported flash calculation types
FLASH_CALCULATIONS = [
    {
        "id": "pt_flash",
        "name": "Pressure-Temperature Flash",
        "description": "Calculate thermodynamic properties based on pressure and temperature.",
        "endpoint": "/pt_flash",
        "primary_variables": ["pressure", "temperature"],
        "example_request": {
            "composition": [
                {"fluid": "METHANE", "fraction": 0.8},
                {"fluid": "ETHANE", "fraction": 0.2}
            ],
            "variables": {
                "pressure": {
                    "range": {"from": 10, "to": 50},
                    "resolution": 10
                },
                "temperature": {
                    "range": {"from": -20, "to": 30},
                    "resolution": 5
                }
            },
            "calculation": {
                "properties": [
                    "density", "enthalpy", "entropy", "sound_speed", "viscosity", "phase"
                ],
                "units_system": "SI",
                "response_format": "json",
                "grid_type": "adaptive"
            }
        }
    },
    {
        "id": "ph_flash",
        "name": "Pressure-Enthalpy Flash",
        "description": "Calculate thermodynamic properties based on pressure and enthalpy.",
        "endpoint": "/ph_flash",
        "primary_variables": ["pressure", "enthalpy"],
        "example_request": {
            "composition": [
                {"fluid": "METHANE", "fraction": 0.85},
                {"fluid": "ETHANE", "fraction": 0.1},
                {"fluid": "PROPANE", "fraction": 0.05}
            ],
            "variables": {
                "pressure": {
                    "range": {"from": 5, "to": 50},
                    "resolution": 5
                },
                "enthalpy": {
                    "range": {"from": 300, "to": 800},
                    "resolution": 50
                }
            },
            "calculation": {
                "properties": [
                    "temperature", "density", "vapor_fraction", "entropy", "cp", "cv", "sound_speed", "phase"
                ],
                "units_system": "SI",
                "response_format": "json",
                "grid_type": "adaptive"
            }
        }
    },
    {
        "id": "ts_flash",
        "name": "Temperature-Entropy Flash",
        "description": "Calculate thermodynamic properties based on temperature and entropy.",
        "endpoint": "/ts_flash",
        "primary_variables": ["temperature", "entropy"],
        "example_request": {
            "composition": [
                {"fluid": "METHANE", "fraction": 0.85},
                {"fluid": "ETHANE", "fraction": 0.1},
                {"fluid": "PROPANE", "fraction": 0.05}
            ],
            "variables": {
                "temperature": {
                    "range": {"from": -50, "to": 50},
                    "resolution": 10
                },
                "entropy": {
                    "range": {"from": 160, "to": 240},
                    "resolution": 20
                }
            },
            "calculation": {
                "properties": [
                    "pressure", "density", "vapor_fraction", "enthalpy", "cp", "cv", "sound_speed", "phase"
                ],
                "units_system": "SI",
                "response_format": "json"
            }
        }
    },
    {
        "id": "vt_flash",
        "name": "Volume-Temperature Flash",
        "description": "Calculate thermodynamic properties based on specific volume and temperature.",
        "endpoint": "/vt_flash",
        "primary_variables": ["specific_volume", "temperature"],
        "example_request": {
            "composition": [
                {"fluid": "METHANE", "fraction": 1}
            ],
            "variables": {
                "temperature": {
                    "range": {"from": -20, "to": 30},
                    "resolution": 5
                },
                "specific_volume": {
                    "range": {"from": 0.02, "to": 0.1},
                    "resolution": 0.01
                }
            },
            "calculation": {
                "properties": [
                    "temperature", "pressure", "density", "enthalpy", "entropy", "phase"
                ],
                "units_system": "SI"
            }
        }
    },
    {
        "id": "uv_flash",
        "name": "Energy-Volume Flash",
        "description": "Calculate thermodynamic properties based on internal energy and specific volume.",
        "endpoint": "/uv_flash",
        "primary_variables": ["internal_energy", "specific_volume"],
        "example_request": {
            "composition": [
                {"fluid": "METHANE", "fraction": 1}
            ],
            "variables": {
                "internal_energy": {
                    "range": {"from": -2000, "to": 3000},
                    "resolution": 500
                },
                "specific_volume": {
                    "range": {"from": 0.02, "to": 0.1},
                    "resolution": 0.01
                }
            },
            "calculation": {
                "properties": [
                    "temperature", "pressure", "density", "enthalpy", "entropy", "phase"
                ],
                "units_system": "SI"
            }
        }
    }
]

# Other endpoints
OTHER_ENDPOINTS = [
    {
        "id": "phase_envelope_pt",
        "name": "Phase Envelope PT",
        "description": "Calculate phase envelope (bubble and dew curves) in the P-T plane.",
        "endpoint": "/phase_envelope_pt",
        "example_request": {
            "composition": [
                {"fluid": "CO2", "fraction": 0.5},
                {"fluid": "METHANE", "fraction": 0.5}
            ],
            "variables": {
                "temperature": {
                    "range": {"from": -80, "to": 50},
                    "resolution": 2
                }
            },
            "calculation": {
                "curve_type": "both"
            }
        }
    },
    {
        "id": "phase_envelope_ph",
        "name": "Phase Envelope PH",
        "description": "Calculate phase envelope (bubble and dew curves) in the P-H plane.",
        "endpoint": "/phase_envelope_ph",
        "example_request": {
            "composition": [
                {"fluid": "NITROGEN", "fraction": 0.2},
                {"fluid": "OXYGEN", "fraction": 0.8}
            ],
            "variables": {
                "pressure": {
                    "range": {"from": 10, "to": 200},
                    "resolution": 5
                }
            },
            "calculation": {
                "curve_type": "both"
            }
        }
    },
    {
        "id": "critical_point",
        "name": "Critical Point",
        "description": "Calculate critical point parameters for a mixture.",
        "endpoint": "/critical_point",
        "example_request": {
            "composition": [
                {"fluid": "CO2", "fraction": 0.7},
                {"fluid": "NITROGEN", "fraction": 0.3}
            ],
            "units_system": "SI"
        }
    },
    {
        "id": "models_info",
        "name": "Models Information",
        "description": "Get information about thermodynamic models used for a composition.",
        "endpoint": "/models_info",
        "example_request": {
            "composition": [
                {"fluid": "CO2", "fraction": 0.5},
                {"fluid": "N2", "fraction": 0.2},
                {"fluid": "METHANE", "fraction": 0.3}
            ]
        }
    }
]

# Utility endpoints
UTILITY_ENDPOINTS = [
    {
        "id": "healthz",
        "name": "Health Check",
        "description": "Check the health status of the API.",
        "endpoint": "/healthz",
        "parameters": [
            {
                "name": "verbose",
                "type": "boolean",
                "description": "Include detailed system information",
                "required": False,
                "default": False
            }
        ]
    },
    {
        "id": "available_fluids",
        "name": "Available Fluids",
        "description": "Get a list of all available fluids.",
        "endpoint": "/available_fluids",
        "parameters": [
            {
                "name": "short_only",
                "type": "boolean",
                "description": "Return only basic fluid information",
                "required": False,
                "default": False
            },
            {
                "name": "search",
                "type": "string",
                "description": "Filter fluids by name, formula, or CAS number",
                "required": False
            }
        ]
    },
    {
        "id": "available_fluids_by_id",
        "name": "Fluid Information",
        "description": "Get detailed information about a specific fluid.",
        "endpoint": "/available_fluids/{fluid_id}",
        "parameters": [
            {
                "name": "fluid_id",
                "type": "string",
                "description": "ID of the fluid to retrieve",
                "required": True
            }
        ]
    },
    {
        "id": "available_properties",
        "name": "Available Properties",
        "description": "Get a list of all available properties.",
        "endpoint": "/available_properties",
        "parameters": [
            {
                "name": "flash_type",
                "type": "string",
                "description": "Filter properties by compatible flash type",
                "required": False
            },
            {
                "name": "input_only",
                "type": "boolean",
                "description": "Return only input properties",
                "required": False,
                "default": False
            },
            {
                "name": "output_only",
                "type": "boolean",
                "description": "Return only output properties",
                "required": False,
                "default": False
            }
        ]
    }
]

# Configuration options
CONFIG_OPTIONS = {
    "unit_systems": [
        {
            "id": "SI",
            "name": "SI Units",
            "description": "International System of Units",
            "examples": {
                "temperature": "°C",
                "pressure": "bar",
                "density": "mol/L",
                "enthalpy": "J/mol"
            }
        },
        {
            "id": "CGS",
            "name": "CGS Units",
            "description": "Centimeter-Gram-Second system",
            "examples": {
                "temperature": "°C",
                "pressure": "dyn/cm²",
                "density": "g/cm³",
                "enthalpy": "erg/g"
            }
        }
    ],
    "grid_types": [
        {
            "id": "equidistant",
            "name": "Equidistant Grid",
            "description": "Regular grid with equal spacing between points (default)"
        },
        {
            "id": "adaptive",
            "name": "Adaptive Grid",
            "description": "Higher resolution near phase boundaries where fluid properties change rapidly",
            "parameters": [
                {
                    "name": "enhancement_factor",
                    "type": "float",
                    "description": "Controls how much to increase resolution near boundaries",
                    "default": 5.0
                },
                {
                    "name": "boundary_zone_width",
                    "type": "float",
                    "description": "Width of zone around boundary where resolution is enhanced",
                    "default": "10% of total range (null)"
                }
            ]
        },
        {
            "id": "logarithmic",
            "name": "Logarithmic Grid",
            "description": "More points at lower values, fewer at higher values. Useful for ranges spanning multiple orders of magnitude."
        },
        {
            "id": "exponential",
            "name": "Exponential Grid",
            "description": "More points at higher values, fewer at lower values. Opposite of logarithmic grid."
        }
    ],
    "response_formats": [
        {
            "id": "json",
            "name": "JSON",
            "description": "Standard JSON format with property values and units"
        },
        {
            "id": "olga_tab",
            "name": "OLGA TAB",
            "description": "Specialized format used by multiphase flow simulators"
        }
    ]
}

@api_info_bp.route('/api_info', methods=['GET'])
def api_info():
    """
    Get comprehensive information about the API, its endpoints, and capabilities.
    
    Query parameters:
    - section: Filter to show only a specific section ('flash', 'other', 'utility', 'config')
    """
    try:
        # Parse query parameters
        section = request.args.get('section')
        
        # Build response based on section filter
        response = {
            "api_name": "Span-Wagner EOS API",
            "version": "1.0.0",
            "description": "A RESTful API for thermodynamic property calculations using the Span-Wagner Equation of State"
        }
        
        if not section or section == 'flash':
            response["flash_calculations"] = FLASH_CALCULATIONS
            
        if not section or section == 'other':
            response["other_endpoints"] = OTHER_ENDPOINTS
            
        if not section or section == 'utility':
            response["utility_endpoints"] = UTILITY_ENDPOINTS
            
        if not section or section == 'config':
            response["configuration_options"] = CONFIG_OPTIONS
        
        # Add documentation links
        response["documentation"] = {
            "readme": "/api_info/readme",
            "postman_collection": "/api_info/postman_collection"
        }
        
        return jsonify(response)
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@api_info_bp.route('/api_info/readme', methods=['GET'])
def api_info_readme():
    """Get the README.md file content"""
    try:
        # Base directory is the parent of the API folder
        base_dir = Path(__file__).resolve().parent.parent.parent.parent
        readme_path = base_dir / 'README.md'
        
        if not readme_path.exists():
            return jsonify({'error': 'README.md file not found'}), 404
            
        with open(readme_path, 'r', encoding='utf-8') as f:
            readme_content = f.read()
            
        return jsonify({
            "filename": "README.md",
            "content": readme_content
        })
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@api_info_bp.route('/api_info/postman_collection', methods=['GET'])
def api_info_postman_collection():
    """Get the Postman collection file"""
    try:
        # Base directory is the parent of the API folder
        base_dir = Path(__file__).resolve().parent.parent.parent.parent
        postman_paths = [
            base_dir / 'Span & Wagner EOS.postman_collection.json',
            base_dir / 'Span-Wagner-EOS.postman_collection.json',
            base_dir / 'postman_collection.json'
        ]
        
        for path in postman_paths:
            if path.exists():
                with open(path, 'r', encoding='utf-8') as f:
                    postman_content = f.read()
                    
                return jsonify({
                    "filename": path.name,
                    "content": postman_content
                })
                
        return jsonify({'error': 'Postman collection file not found'}), 404
        
    except Exception as e:
        print("Error processing request:", file=sys.stderr)
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500