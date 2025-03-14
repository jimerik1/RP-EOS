# FORTRAN

This directory contains FORTRAN libraries and source code for REFPROP integration.

## Purpose

The FORTRAN directory stores the compiled REFPROP FORTRAN libraries and potentially any custom FORTRAN extensions or modifications. These libraries provide the core thermodynamic calculation capabilities used by the API.

## Contents

- REFPROP dynamic libraries (.dll, .so, or .dylib files)
- Custom FORTRAN modules (if applicable)
- Configuration files for FORTRAN libraries

## Integration

The ctREFPROP Python module interacts with these FORTRAN libraries using ctypes bindings, allowing the API to leverage the high-performance thermodynamic calculations implemented in FORTRAN while providing a modern API interface.

## Installation Note

When deploying the application:
1. Ensure the appropriate REFPROP FORTRAN libraries for your platform (Windows, Linux, macOS) are placed in this directory
2. Verify file permissions allow execution of the libraries
3. The Docker configuration handles the appropriate setup of these libraries in the container environment