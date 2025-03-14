# FLUIDS

This directory contains fluid property files used by REFPROP for thermodynamic calculations.

## Purpose

The FLUIDS directory stores the necessary data files for pure fluids that REFPROP requires to perform calculations. These files contain equation of state parameters, transport property correlations, and other thermodynamic data.

## Usage

When adding new fluids:
1. Place the fluid property file (usually with .FLD extension) in this directory
2. The fluid will be automatically recognized by REFPROP when referenced by name

## File Format

Fluid property files follow the NIST REFPROP format for fluid definitions, including:
- Critical parameters
- Equation of state coefficients
- Reference state definitions
- Transport property correlations

## Default Fluids

REFPROP includes many common fluids by default, including:
- Water
- Common refrigerants
- Natural gas components
- Industrial fluids

Custom fluids can be added as needed for specific applications.