# Span-Wagner EOS API

A RESTful API for thermodynamic property calculations using the Span-Wagner Equation of State, implemented with REFPROP and Flask.

## Overview

This project provides a containerized web API for calculating thermodynamic properties of fluids and mixtures. It leverages the REFPROP Fortran library for high-accuracy calculations and exposes the functionality through a modern REST API.

## Features

- RESTful API for thermodynamic calculations
- Support for various flash calculations:
  - PT-Flash (Pressure-Temperature)
  - PH-Flash (Pressure-Enthalpy)
- Phase envelope calculations:
  - P-T (Pressure-Temperature) plane
  - P-H (Pressure-Enthalpy) plane
- Support for pure fluids and mixtures
- Property calculations including:
  - Density, enthalpy, entropy
  - Phase determination
  - Transport properties
  - Critical properties
- Multiple unit systems (SI and CGS)
- Docker containerization for easy deployment
- Modular Flask architecture for extensibility

## Installation

### Prerequisites
- Docker and Docker Compose

### Setup
1. Clone this repository
2. Run `docker-compose up -d`
3. Access the API at `http://localhost:5051`

## API Endpoints

The API provides the following endpoints for thermodynamic calculations:

### 1. PT-Flash: `/pt_flash`

Calculate thermodynamic properties based on pressure and temperature.

#### Request Format

```json
{
  "composition": [
    {"fluid": "FLUID_NAME1", "fraction": X1},
    {"fluid": "FLUID_NAME2", "fraction": X2},
    ...
  ],
  "pressure_range": {
    "from": P_MIN,
    "to": P_MAX
  },
  "temperature_range": {
    "from": T_MIN,
    "to": T_MAX
  },
  "pressure_resolution": P_STEP,
  "temperature_resolution": T_STEP,
  "properties": [
    "property1",
    "property2",
    ...
  ],
  "units_system": "SI" or "CGS"
}
```

#### Example

```json
{
  "composition": [
    {"fluid": "CO2", "fraction": 0.7},
    {"fluid": "N2", "fraction": 0.3}
  ],
  "pressure_range": {
    "from": 10,
    "to": 50
  },
  "temperature_range": {
    "from": -20,
    "to": 30
  },
  "pressure_resolution": 10,
  "temperature_resolution": 5,
  "properties": [
    "density",
    "enthalpy",
    "entropy",
    "sound_speed",
    "viscosity",
    "phase"
  ],
  "units_system": "SI"
}
```

### 2. PH-Flash: `/ph_flash`

Calculate thermodynamic properties based on pressure and enthalpy.

#### Request Format

```json
{
  "composition": [
    {"fluid": "FLUID_NAME1", "fraction": X1},
    {"fluid": "FLUID_NAME2", "fraction": X2},
    ...
  ],
  "pressure_range": {
    "from": P_MIN,
    "to": P_MAX
  },
  "enthalpy_range": {
    "from": H_MIN,
    "to": H_MAX
  },
  "pressure_resolution": P_STEP,
  "enthalpy_resolution": H_STEP,
  "properties": [
    "property1",
    "property2",
    ...
  ],
  "units_system": "SI" or "CGS"
}
```

#### Example

```json
{
  "composition": [
    {"fluid": "METHANE", "fraction": 0.85},
    {"fluid": "ETHANE", "fraction": 0.1},
    {"fluid": "PROPANE", "fraction": 0.05}
  ],
  "pressure_range": {
    "from": 5,
    "to": 50
  },
  "enthalpy_range": {
    "from": 300,
    "to": 800
  },
  "pressure_resolution": 5,
  "enthalpy_resolution": 50,
  "properties": [
    "temperature",
    "density",
    "vapor_fraction",
    "entropy",
    "cp",
    "cv",
    "sound_speed",
    "phase"
  ],
  "units_system": "SI"
}
```

### 3. Phase Envelope PT: `/phase_envelope_pt`

Calculate phase envelope (bubble and dew curves) in the pressure-temperature plane.

#### Request Format

```json
{
  "composition": [
    {"fluid": "FLUID_NAME1", "fraction": X1},
    {"fluid": "FLUID_NAME2", "fraction": X2},
    ...
  ],
  "temperature_range": {
    "from": T_MIN,
    "to": T_MAX
  },
  "temperature_resolution": T_STEP,
  "desired_curve": "bubble", "dew", or "both"
}
```

#### Example

```json
{
  "composition": [
    {"fluid": "CO2", "fraction": 0.5},
    {"fluid": "METHANE", "fraction": 0.5}
  ],
  "temperature_range": {
    "from": -80,
    "to": 50
  },
  "temperature_resolution": 2,
  "desired_curve": "both"
}
```

### 4. Phase Envelope PH: `/phase_envelope_ph`

Calculate phase envelope (bubble and dew curves) in the pressure-enthalpy plane.

#### Request Format

```json
{
  "composition": [
    {"fluid": "FLUID_NAME1", "fraction": X1},
    {"fluid": "FLUID_NAME2", "fraction": X2},
    ...
  ],
  "pressure_range": {
    "from": P_MIN,
    "to": P_MAX
  },
  "pressure_resolution": P_STEP,
  "desired_curve": "bubble", "dew", or "both"
}
```

#### Example

```json
{
  "composition": [
    {"fluid": "NITROGEN", "fraction": 0.2},
    {"fluid": "OXYGEN", "fraction": 0.8}
  ],
  "pressure_range": {
    "from": 10,
    "to": 200
  },
  "pressure_resolution": 5,
  "desired_curve": "both"
}
```

## Available Properties

The following properties can be requested in the API calls:

| Property                    | Description                                      | SI Unit                | CGS Unit                 |
|-----------------------------|--------------------------------------------------|------------------------|--------------------------|
| `temperature`               | Temperature                                      | °C                     | °C                       |
| `pressure`                  | Pressure                                         | bar                    | dyn/cm²                  |
| `density`                   | Bulk density                                     | mol/L                  | g/cm³                    |
| `liquid_density`            | Liquid phase density                             | mol/L                  | g/cm³                    |
| `vapor_density`             | Vapor phase density                              | mol/L                  | g/cm³                    |
| `critical_density`          | Critical density                                 | mol/L                  | g/cm³                    |
| `critical_pressure`         | Critical pressure                                | bar                    | dyn/cm²                  |
| `critical_temperature`      | Critical temperature                             | K                      | K                        |
| `internal_energy`           | Specific internal energy                         | J/mol                  | erg/g                    |
| `enthalpy`                  | Specific enthalpy                                | J/mol                  | erg/g                    |
| `entropy`                   | Specific entropy                                 | J/(mol·K)              | erg/(g·K)                |
| `cv`                        | Isochoric specific heat capacity                 | J/(mol·K)              | erg/(g·K)                |
| `cp`                        | Isobaric specific heat capacity                  | J/(mol·K)              | erg/(g·K)                |
| `sound_speed`               | Speed of sound                                   | m/s                    | cm/s                     |
| `viscosity`                 | Dynamic viscosity                                | μPa·s                  | poise                    |
| `thermal_conductivity`      | Thermal conductivity                             | W/(m·K)                | cal/(s·cm·K)             |
| `surface_tension`           | Surface tension (two-phase only)                 | N/m                    | dyn/cm                   |
| `vapor_fraction`            | Vapor quality/fraction                           | dimensionless          | dimensionless            |
| `compressibility_factor`    | Compressibility factor (Z)                       | dimensionless          | dimensionless            |
| `prandtl_number`            | Prandtl number                                   | dimensionless          | dimensionless            |
| `isothermal_compressibility`| Isothermal compressibility                       | 1/kPa                  | cm²/dyn                  |
| `volume_expansivity`        | Volume expansivity                               | 1/K                    | 1/K                      |
| `dp_dt_saturation`          | dP/dT along saturation curve                     | kPa/K                  | dyn/(cm²·K)              |
| `joule_thomson_coefficient` | Joule-Thomson coefficient                        | K/bar                  | K·cm²/dyn                |
| `kinematic_viscosity`       | Kinematic viscosity                              | cm²/s                  | cm²/s                    |
| `thermal_diffusivity`       | Thermal diffusivity                              | cm²/s                  | cm²/s                    |
| `phase`                     | Phase state (Liquid, Vapor, etc.)                | text                   | text                     |

## Available Fluids

The API supports all fluids available in the REFPROP database. The most commonly used fluids include:

### Pure Compounds
- Water (`WATER`)
- Carbon dioxide (`CO2`)
- Nitrogen (`NITROGEN`)
- Oxygen (`OXYGEN`)
- Hydrogen (`HYDROGEN`)
- Methane (`METHANE`)
- Ethane (`ETHANE`)
- Propane (`PROPANE`)
- Butane (`BUTANE`)
- Isobutane (`ISOBUTANE`)
- Pentane (`PENTANE`)
- Isopentane (`ISOPENTANE`)
- Hexane (`HEXANE`)
- Heptane (`HEPTANE`)
- Octane (`OCTANE`)
- Ammonia (`AMMONIA`)
- Argon (`ARGON`)
- R134a (`R134A`)
- And many more refrigerants, hydrocarbons, and industrial fluids

### Predefined Mixtures
- Air (`AIR`)
- Natural gas compositions (custom)

## Unit Systems

The API supports two unit systems:

### SI Units (International System of Units)
- Temperature: °C (Celsius) for temperature, K (Kelvin) for critical temperature
- Pressure: bar
- Density: mol/L
- Energy: J/mol
- Entropy: J/(mol·K)
- Viscosity: μPa·s
- Thermal conductivity: W/(m·K)
- Surface tension: N/m

### CGS Units (Centimeter-Gram-Second system)
- Temperature: °C (Celsius) for temperature, K (Kelvin) for critical temperature
- Pressure: dyn/cm²
- Density: g/cm³
- Energy: erg/g
- Entropy: erg/(g·K)
- Viscosity: poise
- Thermal conductivity: cal/(s·cm·K)
- Surface tension: dyn/cm

## Response Format

The API returns a JSON response with the following structure:

```json
{
  "results": [
    {
      "index": 0,
      "temperature": {"value": X, "unit": "°C"},
      "pressure": {"value": Y, "unit": "bar"},
      "property1": {"value": Z1, "unit": "unit1"},
      "property2": {"value": Z2, "unit": "unit2"},
      ...
    },
    ...
  ]
}
```

## Error Handling

The API returns standard HTTP status codes:
- 200: Success
- 400: Bad request (e.g., invalid input)
- 500: Internal server error

Error responses include a detailed error message:

```json
{
  "error": "Error message describing the issue"
}
```

## Health Check

The API provides a health check endpoint:
- URL: `/health`
- Method: GET
- Response: `{"status": "ok", "message": "REFPROP API is running"}`

## Usage Examples

See the included test scripts for examples of how to interact with the API and visualize the results.

## Notes on Fluid Composition

When specifying fluid compositions:
- The fractions must sum to 1.0
- For pure fluids, use a single component with fraction 1.0
- Fluid names are case-insensitive but typically uppercase in REFPROP

## Development and Contributions

Contributions to the API are welcome. Please follow the standard GitHub workflow:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

This project uses REFPROP, which requires a license from NIST. Ensure you have the appropriate licensing to use REFPROP for your application.