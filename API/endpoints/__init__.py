from flask import Blueprint

# Create blueprints - only define them here
pt_flash_bp = Blueprint('pt_flash', __name__)
ph_flash_bp = Blueprint('ph_flash', __name__)

# List of all blueprints to register with the app
blueprints = [
    pt_flash_bp,
    ph_flash_bp,
    # Add new blueprints here as you create them
]

# Import views AFTER defining blueprints to avoid circular imports
# These imports will apply the route decorators to the blueprints
from API.endpoints import pt_flash
from API.endpoints import ph_flash