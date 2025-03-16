from flask import Blueprint

pt_flash_bp = Blueprint('pt_flash', __name__)
ph_flash_bp = Blueprint('ph_flash', __name__)
ts_flash_bp = Blueprint('ts_flash', __name__)
phase_envelope_pt_bp = Blueprint('phase_envelope_pt', __name__)  
phase_envelope_ph_bp = Blueprint('phase_envelope_ph', __name__)
critical_point_bp = Blueprint('critical_point', __name__)
extended_pt_flash_bp = Blueprint('extended_pt_flash', __name__)
phase_boundaries_bp = Blueprint('phase_boundaries', __name__)
vt_flash_bp = Blueprint('vt_flash', __name__)  
uv_flash_bp = Blueprint('uv_flash', __name__)  

blueprints = [
    pt_flash_bp,
    ph_flash_bp,
    ts_flash_bp,
    phase_envelope_pt_bp,
    phase_envelope_ph_bp,
    critical_point_bp,
    extended_pt_flash_bp,
    phase_boundaries_bp,
    vt_flash_bp,    
    uv_flash_bp,    
    # Add more here ...
]

# Import the endpoints AFTER defining the blueprint
from API.endpoints import pt_flash
from API.endpoints import ph_flash
from API.endpoints import ts_flash
from API.endpoints import phase_envelope_pt  
from API.endpoints import phase_envelope_ph
from API.endpoints import critical_point
from API.endpoints import extended_pt_flash
from API.endpoints import phase_boundaries
from API.endpoints import vt_flash    
from API.endpoints import uv_flash    