from flask import Blueprint

# Import blueprints from endpoint modules
from API.endpoints.utilities.healthz import healthz_bp
from API.endpoints.utilities.available_fluids import available_fluids_bp
from API.endpoints.utilities.available_properties import available_properties_bp
from API.endpoints.utilities.api_info import api_info_bp

# List of all utility blueprints
utility_blueprints = [
    healthz_bp,
    available_fluids_bp,
    available_properties_bp,
    api_info_bp
]