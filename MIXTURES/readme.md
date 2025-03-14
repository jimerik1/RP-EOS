# MIXTURES

This directory contains mixture definition files used by REFPROP for multi-component thermodynamic calculations.

## Purpose

The MIXTURES directory stores predefined mixture compositions and interaction parameters that REFPROP uses when performing calculations on fluid mixtures. These files define how components interact with each other in mixtures.

## Usage

When working with mixtures:
1. Predefined mixtures can be referenced by name in the API
2. Custom mixture files (usually with .MIX extension) can be added to this directory
3. Binary interaction parameters are defined in these files for accurate mixture property prediction

## File Format

Mixture files follow the NIST REFPROP format for mixture definitions, including:
- Component names and proportions
- Binary interaction parameters
- Mixing rules
- Optional: specialized models for specific mixture types

## Common Mixtures

REFPROP includes definitions for common mixtures like:
- Natural gas compositions
- Refrigerant blends
- Air components
- Industrial gas mixtures

Custom mixtures can be defined for specific applications when needed.