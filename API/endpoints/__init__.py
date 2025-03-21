from flask import Blueprint

# Create blueprints for all endpoints
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
models_info_bp = Blueprint('models_info', __name__)  

# New dedicated OLGA TAB format endpoints
pt_flash_olga_bp = Blueprint('pt_flash_olga', __name__)
ph_flash_olga_bp = Blueprint('ph_flash_olga', __name__)
ts_flash_olga_bp = Blueprint('ts_flash_olga', __name__)
vt_flash_olga_bp = Blueprint('vt_flash_olga', __name__)
uv_flash_olga_bp = Blueprint('uv_flash_olga', __name__)

# List of all blueprints to register
blueprints = [
    # Standard JSON endpoints
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
    models_info_bp,
    
    # OLGA TAB format endpoints
    pt_flash_olga_bp,
    ph_flash_olga_bp,
    ts_flash_olga_bp,
    vt_flash_olga_bp,
    uv_flash_olga_bp,
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
from API.endpoints import models_info

# Import the new OLGA TAB format endpoints
from API.endpoints import pt_flash_olga
from API.endpoints import ph_flash_olga
from API.endpoints import ts_flash_olga
from API.endpoints import vt_flash_olga
from API.endpoints import uv_flash_olga