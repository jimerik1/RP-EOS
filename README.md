# Span-Wagner EOS API

A RESTful API for thermodynamic property calculations using the Span-Wagner Equation of State, implemented with REFPROP and Flask.

## Overview

This project provides a containerized web API for calculating thermodynamic properties of fluids and mixtures. It leverages the REFPROP Fortran library for high-accuracy calculations and exposes the functionality through a modern REST API.

## Features

- RESTful API for thermodynamic calculations
- Support for various flash calculations:
  - PT-Flash (Pressure-Temperature)
  - PH-Flash (Pressure-Enthalpy)
  - TS-Flash (Temperature-Entropy)
- Phase envelope calculations:
  - P-T (Pressure-Temperature) plane
  - P-H (Pressure-Enthalpy) plane
- Support for pure fluids and mixtures
- Property calculations including:
  - Density, enthalpy, entropy
  - Phase determination
  - Transport properties
  - Critical properties
- **Smart grid generation with adaptive resolution near phase boundaries**
- Multiple unit systems (SI and CGS)
- Multiple response formats (JSON and OLGA TAB)
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
  "variables": {
    "pressure": {
      "range": {
        "from": P_MIN,
        "to": P_MAX
      },
      "resolution": P_STEP
    },
    "temperature": {
      "range": {
        "from": T_MIN,
        "to": T_MAX
      },
      "resolution": T_STEP
    }
  },
  "calculation": {
    "properties": [
      "property1",
      "property2",
      ...
    ],
    "units_system": "SI" or "CGS",
    "response_format": "json" or "olga_tab",
    "grid_type": "equidistant" or "adaptive" or "logarithmic" or "exponential",
    "enhancement_factor": 5.0,
    "boundary_zone_width": null
  }
}
```

#### Example

```json
{
  "composition": [
    {"fluid": "CO2", "fraction": 0.7},
    {"fluid": "N2", "fraction": 0.3}
  ],
  "variables": {
    "pressure": {
      "range": {
        "from": 10,
        "to": 50
      },
      "resolution": 10
    },
    "temperature": {
      "range": {
        "from": -20,
        "to": 30
      },
      "resolution": 5
    }
  },
  "calculation": {
    "properties": [
      "density",
      "enthalpy",
      "entropy",
      "sound_speed",
      "viscosity",
      "phase"
    ],
    "units_system": "SI",
    "response_format": "json",
    "grid_type": "adaptive"
  }
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
  "variables": {
    "pressure": {
      "range": {
        "from": P_MIN,
        "to": P_MAX
      },
      "resolution": P_STEP
    },
    "enthalpy": {
      "range": {
        "from": H_MIN,
        "to": H_MAX
      },
      "resolution": H_STEP
    }
  },
  "calculation": {
    "properties": [
      "property1",
      "property2",
      ...
    ],
    "units_system": "SI" or "CGS",
    "response_format": "json" or "olga_tab",
    "grid_type": "equidistant" or "adaptive" or "logarithmic" or "exponential",
    "enhancement_factor": 5.0,
    "boundary_zone_width": null
  }
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
  "variables": {
    "pressure": {
      "range": {
        "from": 5,
        "to": 50
      },
      "resolution": 5
    },
    "enthalpy": {
      "range": {
        "from": 300,
        "to": 800
      },
      "resolution": 50
    }
  },
  "calculation": {
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
    "units_system": "SI",
    "response_format": "json",
    "grid_type": "adaptive"
  }
}
```

### 3. TS-Flash: `/ts_flash`

Calculate thermodynamic properties based on temperature and entropy.

#### Request Format

```json
{
  "composition": [
    {"fluid": "FLUID_NAME1", "fraction": X1},
    {"fluid": "FLUID_NAME2", "fraction": X2},
    ...
  ],
  "variables": {
    "temperature": {
      "range": {
        "from": T_MIN,
        "to": T_MAX
      },
      "resolution": T_STEP
    },
    "entropy": {
      "range": {
        "from": S_MIN,
        "to": S_MAX
      },
      "resolution": S_STEP
    }
  },
  "calculation": {
    "properties": [
      "property1",
      "property2",
      ...
    ],
    "units_system": "SI" or "CGS",
    "response_format": "json" or "olga_tab",
    "grid_type": "equidistant" or "adaptive" or "logarithmic" or "exponential",
    "enhancement_factor": 5.0,
    "boundary_zone_width": null
  }
}
```

### 4. Phase Envelope PT: `/phase_envelope_pt`

Calculate phase envelope (bubble and dew curves) in the pressure-temperature plane.

#### Request Format

```json
{
  "composition": [
    {"fluid": "FLUID_NAME1", "fraction": X1},
    {"fluid": "FLUID_NAME2", "fraction": X2},
    ...
  ],
  "variables": {
    "temperature": {
      "range": {
        "from": T_MIN,
        "to": T_MAX
      },
      "resolution": T_STEP
    }
  },
  "calculation": {
    "curve_type": "bubble", "dew", or "both"
  }
}
```

#### Example

```json
{
  "composition": [
    {"fluid": "CO2", "fraction": 0.5},
    {"fluid": "METHANE", "fraction": 0.5}
  ],
  "variables": {
    "temperature": {
      "range": {
        "from": -80,
        "to": 50
      },
      "resolution": 2
    }
  },
  "calculation": {
    "curve_type": "both"
  }
}
```

### 5. Phase Envelope PH: `/phase_envelope_ph`

Calculate phase envelope (bubble and dew curves) in the pressure-enthalpy plane.

#### Request Format

```json
{
  "composition": [
    {"fluid": "FLUID_NAME1", "fraction": X1},
    {"fluid": "FLUID_NAME2", "fraction": X2},
    ...
  ],
  "variables": {
    "pressure": {
      "range": {
        "from": P_MIN,
        "to": P_MAX
      },
      "resolution": P_STEP
    }
  },
  "calculation": {
    "curve_type": "bubble", "dew", or "both"
  }
}
```

#### Example

```json
{
  "composition": [
    {"fluid": "NITROGEN", "fraction": 0.2},
    {"fluid": "OXYGEN", "fraction": 0.8}
  ],
  "variables": {
    "pressure": {
      "range": {
        "from": 10,
        "to": 200
      },
      "resolution": 5
    }
  },
  "calculation": {
    "curve_type": "both"
  }
}
```

## Grid Generation Types

The API supports several types of grid generation for calculating properties:

### 1. Equidistant (Default)
Regular uniform grid with equal spacing between points. This is the default mode and matches the original behavior.

### 2. Adaptive
Higher resolution near phase boundaries where fluid properties change rapidly. The API automatically detects phase transitions and increases the grid density in these regions.

Parameters:
- `enhancement_factor`: Controls how much to increase resolution near boundaries (default: 5.0)
- `boundary_zone_width`: Width of zone around boundary where resolution is enhanced (if null, calculated as 10% of the total range)

### 3. Logarithmic
More points at lower values, fewer at higher values. Useful for pressure ranges spanning multiple orders of magnitude.

### 4. Exponential
More points at higher values, fewer at lower values. Provides the opposite distribution to logarithmic grids.

## Response Formats

The API supports multiple response formats:

### 1. JSON (Default)

JSON responses include calculated properties with their units:

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
  ],
  "grid_info": {
    "type": "adaptive",
    "pressure_points": 12,
    "temperature_points": 15,
    "total_points": 180
  }
}
```

### 2. OLGA TAB

OLGA TAB is a specialized format used by multiphase flow simulators like OLGA. It represents fluid properties as a structured table across a grid.

To request OLGA TAB format, add `"response_format": "olga_tab"` to the calculation section of your request.

Example:

```json
{
  "calculation": {
    "properties": [...],
    "units_system": "SI",
    "response_format": "olga_tab",
    "grid_type": "adaptive"
  }
}
```

The response will be a plain text file in OLGA TAB format with the following structure:

```
'Span-Wagner EOS FLUID_NAME-FRACTION'
   NP  NT    .276731E-08
    P1    P2    P3    P4    P5
    ...
    T1    T2    T3    T4    T5
    ...
 LIQUID DENSITY (KG/M3)                
    D11    D12    D13    D14    D15
    ...
 GAS DENSITY (KG/M3)                
    D21    D22    D23    D24    D25
    ...
```

OLGA TAB format includes standard property sets typically needed for multiphase flow simulations:

- Density (liquid and gas)
- Viscosity (liquid and gas)
- Specific heat (liquid and gas)
- Thermal conductivity (liquid and gas)
- Enthalpy (liquid and gas)
- Entropy (liquid and gas)
- Surface tension
- Derivative properties (density derivatives with respect to pressure and temperature)
- Phase fractions (gas mass fraction)
- Component distributions (water mass fractions in liquid and gas phases)

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
| `dDdP`                      | Pressure derivative of density                   | (mol/L)/kPa            | (g/cm³)/(dyn/cm²)        |
| `dDdT`                      | Temperature derivative of density                | (mol/L)/K              | (g/cm³)/K                |
| `x`                         | Liquid phase composition                         | mole fraction          | mole fraction            |
| `y`                         | Vapor phase composition                          | mole fraction          | mole fraction            |

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

See the included examples in the EXAMPLES directory:

- `generate_olga_tab.py`: Example of generating an OLGA TAB file directly from the API
- `json_to_olga_tab.py`: Utility to convert an existing JSON response to OLGA TAB format

## Notes on Fluid Composition

When specifying fluid compositions:
- The fractions must sum to 1.0
- For pure fluids, use a single component with fraction 1.0
- Fluid names are case-insensitive but typically uppercase in REFPROP

## Benefits of Adaptive Grid Generation

The adaptive grid feature offers several key advantages:

1. **Improved Accuracy**: Higher resolution near phase boundaries where properties change rapidly
2. **Computational Efficiency**: Fewer points needed in regions where properties change slowly
3. **Better Interpolation**: More accurate interpolation of properties between calculated points
4. **Optimal Resource Usage**: Concentrates computational effort where it's most needed

This is particularly valuable for:
- Phase transitions (liquid-vapor interfaces)
- Critical regions
- Retrograde condensation areas
- Anywhere fluid properties change non-linearly

## Development and Contributions

Contributions to the API are welcome. Please follow the standard GitHub workflow:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

This project uses REFPROP, which requires a license from NIST. Ensure you have the appropriate licensing to use REFPROP for your application.