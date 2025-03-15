from flask import Blueprint

pt_flash_bp = Blueprint('pt_flash', __name__)
ph_flash_bp = Blueprint('ph_flash', __name__)
ts_flash_bp = Blueprint('ts_flash', __name__)
phase_envelope_pt_bp = Blueprint('phase_envelope_pt', __name__)  
phase_envelope_ph_bp = Blueprint('phase_envelope_ph', __name__) 

blueprints = [
    pt_flash_bp,
    ph_flash_bp,
    ts_flash_bp,
    phase_envelope_pt_bp,
    phase_envelope_ph_bp,
    # Add more here ...
]

# Import the endpoints AFTER defining the blueprint
from API.endpoints import pt_flash
from API.endpoints import ph_flash
from API.endpoints import ts_flash
from API.endpoints import phase_envelope_pt  
from API.endpoints import phase_envelope_ph  