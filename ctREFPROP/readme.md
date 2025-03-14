# ctREFPROP

This directory contains Python bindings for REFPROP (Reference Fluid Thermodynamic and Transport Properties Database).

## Files

- **ctREFPROP.py**: Python wrapper for the REFPROP dynamic library, providing ctypes-based bindings to access REFPROP functionality from Python

## Purpose

The ctREFPROP module provides a bridge between the Python API and the REFPROP FORTRAN library, allowing thermodynamic property calculations using the NIST REFPROP database within Python applications.

## Usage

```python
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

# Initialize REFPROP
RP = REFPROPFunctionLibrary('/path/to/REFPROP/fortran/directory')

# Example usage for fluid property calculation
result = RP.REFPROP('WATER', 'TP', 'D', 300, 101.325)
```

Refer to the Example_scripts directory for more detailed usage examples.