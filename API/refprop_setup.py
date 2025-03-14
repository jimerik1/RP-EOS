import os
import sys
from pathlib import Path
import traceback

# Base directory is the parent of the API folder
BASE_DIR = Path(__file__).resolve().parent.parent
REFPROP_ROOT = BASE_DIR

# Initialize REFPROP
def initialize_refprop():
    try:
        from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
        os.environ['RPPREFIX'] = str(REFPROP_ROOT)
        rp = REFPROPFunctionLibrary(str(REFPROP_ROOT))
        rp.SETPATHdll(str(REFPROP_ROOT))
        print("REFPROP initialized successfully.")
        return rp
    except Exception as e:
        print("Error during REFPROP initialization:", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

# Global REFPROP instance
RP = initialize_refprop()