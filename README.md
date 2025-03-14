# Create README.md file (if needed)
cat > README.md << 'EOF'
# Span-Wagner EOS API

A RESTful API for thermodynamic property calculations using the Span-Wagner Equation of State, implemented with REFPROP and Flask.

## Overview

This project provides a containerized web API for calculating thermodynamic properties of fluids and mixtures. It leverages the REFPROP Fortran library for high-accuracy calculations and exposes the functionality through a modern REST API.

## Features

- RESTful API for thermodynamic calculations
- Support for various flash calculations:
  - PT-Flash (Pressure-Temperature)
- Support for pure fluids and mixtures
- Property calculations including:
  - Density, enthalpy, entropy
  - Phase determination
  - Transport properties
  - Critical properties
- Docker containerization for easy deployment
- Modular Flask architecture for extensibility

## API Endpoints

### `/calculate` (PT-Flash)
Calculate thermodynamic properties based on pressure and temperature.

## Installation

### Prerequisites
- Docker and Docker Compose

### Setup
1. Clone this repository
2. Run `docker-compose up -d`
3. Access the API at `http://localhost:5051`

## Usage Examples

See the included test scripts for examples of how to interact with the API and visualize the results.
EOF