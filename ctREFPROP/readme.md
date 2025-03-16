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

# REFPROP Python Wrapper Documentation

This document provides detailed documentation for the functions in the Python wrapper for NIST REFPROP (Reference Fluid Thermodynamic and Transport Properties Database). These comments can be added to the top of each function in the wrapper code to improve its usability and make it easier to understand what each function does.

## ABFL1dll
```python
"""
Calculate temperature, pressure, and density from a and b property pair.

Parameters
----------
a : float
    First property value
b : float
    Second property value
z : list
    Composition (mole fractions)
kph : int
    Phase flag: 1 = liquid, 2 = vapor
ab : str
    Property pair code (two characters, e.g., 'TP', 'DH', etc.)
Dmin : float
    Minimum density bound [mol/L]
Dmax : float
    Maximum density bound [mol/L]

Returns
-------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
D : float
    Density [mol/L]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## ABFL2dll
```python
"""
Calculate properties for a 2-phase state including saturation information.

Parameters
----------
a : float
    First property value
b : float
    Second property value
z : list
    Composition (mole fractions)
kq : int
    Flag for calculation type
ksat : int
    Flag for saturation state
ab : str
    Property pair code (two characters)

Returns
-------
Tbub : float
    Bubble-point temperature [K]
Tdew : float
    Dew-point temperature [K]
Pbub : float
    Bubble-point pressure [kPa]
Pdew : float
    Dew-point pressure [kPa]
Dlbub : float
    Bubble-point liquid density [mol/L]
Dvdew : float
    Dew-point vapor density [mol/L]
ybub : list
    Bubble-point vapor composition [mole fractions]
xdew : list
    Dew-point liquid composition [mole fractions]
T : float
    Temperature [K]
P : float
    Pressure [kPa]
Dl : float
    Liquid density [mol/L]
Dv : float
    Vapor density [mol/L]
x : list
    Liquid composition [mole fractions]
y : list
    Vapor composition [mole fractions]
q : float
    Vapor quality [mole basis]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## ABFLASHdll
```python
"""
Comprehensive flash calculation given a property pair.

Parameters
----------
ab : str
    Property pair code (two characters)
a : float
    First property value
b : float
    Second property value
z : list
    Composition (mole fractions)
iFlag : int
    Flag for special options

Returns
-------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
D : float
    Overall density [mol/L]
Dl : float
    Liquid density [mol/L]
Dv : float
    Vapor density [mol/L]
x : list
    Liquid composition [mole fractions]
y : list
    Vapor composition [mole fractions]
q : float
    Vapor quality [mole basis]
e : float
    Internal energy [J/mol]
h : float
    Enthalpy [J/mol]
s : float
    Entropy [J/(mol·K)]
Cv : float
    Isochoric heat capacity [J/(mol·K)]
Cp : float
    Isobaric heat capacity [J/(mol·K)]
w : float
    Speed of sound [m/s]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## ABFLSHdll
```python
"""
Flash calculation given a property pair.

Parameters
----------
ab : str
    Property pair code (two characters)
a : float
    First property value
b : float
    Second property value
z : list
    Composition (mole fractions)
iFlag : int
    Flag for special options

Returns
-------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
D : float
    Overall density [mol/L]
Dl : float
    Liquid density [mol/L]
Dv : float
    Vapor density [mol/L]
x : list
    Liquid composition [mole fractions]
y : list
    Vapor composition [mole fractions]
q : float
    Vapor quality [mole basis]
e : float
    Internal energy [J/mol]
h : float
    Enthalpy [J/mol]
s : float
    Entropy [J/(mol·K)]
Cv : float
    Isochoric heat capacity [J/(mol·K)]
Cp : float
    Isobaric heat capacity [J/(mol·K)]
w : float
    Speed of sound [m/s]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## AGdll
```python
"""
Calculate Helmholtz and Gibbs free energies.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
a : float
    Helmholtz free energy [J/mol]
g : float
    Gibbs free energy [J/mol]
"""
```

## ALLPROPS0dll
```python
"""
Calculates all thermodynamic and transport properties at once.

Parameters
----------
iIn : int
    Input code for calculation type
iOut : list
    Array of output codes
iFlag : int
    Flag for special handling
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
Output : list
    Array of calculated properties
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## ALLPROPS1dll
```python
"""
Calculate a single property given a string identifier.

Parameters
----------
hOut : str
    Property identifier string
iUnits : int
    Unit system flag
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
c : float
    Calculated property value
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## ALLPROPS20dll
```python
"""
Calculate multiple properties (up to 20) given by a string identifier.

Parameters
----------
hOut : str
    Property identifier string
iUnits : int
    Unit system flag
iMass : int
    Flag for mass or molar basis
iFlag : int
    Flag for special handling
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
Output : list
    Array of calculated properties (up to 20 values)
hUnitsArray : str
    String of unit descriptions
iUCodeArray : list
    Array of unit codes
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## ALLPROPSdll
```python
"""
Calculate all available properties given by a string identifier.

Parameters
----------
hOut : str
    Property identifier string
iUnits : int
    Unit system flag
iMass : int
    Flag for mass or molar basis
iFlag : int
    Flag for special handling
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
Output : list
    Array of calculated properties (up to 200 values)
hUnitsArray : str
    String of unit descriptions
iUCodeArray : list
    Array of unit codes
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## B12dll
```python
"""
Calculate second virial coefficient.

Parameters
----------
T : float
    Temperature [K]
z : list
    Composition (mole fractions)

Returns
-------
B : float
    Second virial coefficient [L/mol]
"""
```

## BLCRVdll
```python
"""
Calculate temperature on the boiling line at a specified density.

Parameters
----------
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)
T : float
    Initial estimate of temperature [K]

Returns
-------
T : float
    Temperature [K]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## CCRITdll
```python
"""
Calculate critical flow parameters.

Parameters
----------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
v : float
    Velocity [m/s]
z : list
    Composition (mole fractions)

Returns
-------
Cs : float
    Speed of sound at critical point [m/s]
Ts : float
    Temperature at critical point [K]
Ds : float
    Density at critical point [mol/L]
Ps : float
    Pressure at critical point [kPa]
ws : float
    Speed of sound at throat [m/s]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## CHEMPOTdll
```python
"""
Calculate chemical potentials for each component.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
u : list
    Chemical potentials [J/mol]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## CP0dll
```python
"""
Calculate ideal gas isobaric heat capacity.

Parameters
----------
T : float
    Temperature [K]
z : list
    Composition (mole fractions)

Returns
-------
Cp : float
    Ideal gas isobaric heat capacity [J/(mol·K)]
"""
```

## CRITPdll
```python
"""
Calculate critical parameters for a mixture.

Parameters
----------
z : list
    Composition (mole fractions)

Returns
-------
Tc : float
    Critical temperature [K]
Pc : float
    Critical pressure [kPa]
Dc : float
    Critical density [mol/L]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## CRTPNTdll
```python
"""
Calculate critical point parameters with initial guesses.

Parameters
----------
z : list
    Composition (mole fractions)
Tc : float
    Initial guess for critical temperature [K]
Pc : float
    Initial guess for critical pressure [kPa]

Returns
-------
Tc : float
    Critical temperature [K]
Pc : float
    Critical pressure [kPa]
Dc : float
    Critical density [mol/L]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## CSATKdll
```python
"""
Calculate saturation properties for a pure fluid.

Parameters
----------
icomp : int
    Component number (1-based)
T : float
    Temperature [K]
kph : int
    Phase flag: 1 = liquid, 2 = vapor

Returns
-------
P : float
    Saturation pressure [kPa]
D : float
    Saturation density [mol/L]
Csat : float
    dp/dT along saturation line [kPa/K]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## CSTARdll
```python
"""
Calculate critical flow parameters at throat conditions.

Parameters
----------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
v : float
    Velocity [m/s]
z : list
    Composition (mole fractions)

Returns
-------
Cs : float
    Flow parameter at throat [m/s]
Ts : float
    Temperature at throat [K]
Ds : float
    Density at throat [mol/L]
Ps : float
    Pressure at throat [kPa]
ws : float
    Speed of sound at throat [m/s]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## CV2PKdll
```python
"""
Calculate 2-phase isochoric heat capacity for a pure fluid.

Parameters
----------
icomp : int
    Component number (1-based)
T : float
    Temperature [K]
D : float
    Density [mol/L]

Returns
-------
Cv2p : float
    2-phase isochoric heat capacity [J/(mol·K)]
Csat : float
    dp/dT along saturation line [kPa/K]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## CVCPKdll
```python
"""
Calculate isochoric and isobaric heat capacities for a pure fluid.

Parameters
----------
icomp : int
    Component number (1-based)
T : float
    Temperature [K]
D : float
    Density [mol/L]

Returns
-------
Cv : float
    Isochoric heat capacity [J/(mol·K)]
Cp : float
    Isobaric heat capacity [J/(mol·K)]
"""
```

## CVCPdll
```python
"""
Calculate isochoric and isobaric heat capacities.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
Cv : float
    Isochoric heat capacity [J/(mol·K)]
Cp : float
    Isobaric heat capacity [J/(mol·K)]
"""
```

## DBDTdll
```python
"""
Calculate dB/dT for a mixture, where B is the second virial coefficient.

Parameters
----------
T : float
    Temperature [K]
z : list
    Composition (mole fractions)

Returns
-------
dBT : float
    Temperature derivative of second virial coefficient [L/(mol·K)]
"""
```

## DBFL1dll
```python
"""
Calculate temperature and pressure for a specified density and property.

Parameters
----------
D : float
    Density [mol/L]
b : float
    Property value
z : list
    Composition (mole fractions)
hab : str
    Property identifier code (e.g., 'E', 'H', 'S')

Returns
-------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DBFL2dll
```python
"""
Calculate two-phase state with specified density, property and quality.

Parameters
----------
D : float
    Density [mol/L]
b : float
    Property value
z : list
    Composition (mole fractions)
kq : int
    Phase identifier
ab : str
    Property identifier code

Returns
-------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
Dl : float
    Liquid density [mol/L]
Dv : float
    Vapor density [mol/L]
x : list
    Liquid composition [mole fractions]
y : list
    Vapor composition [mole fractions]
q : float
    Vapor quality [mole basis]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DDDPdll
```python
"""
Calculate dD/dP at constant temperature.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
dDdP : float
    dD/dP at constant T [(mol/L)/kPa]
"""
```

## DDDTdll
```python
"""
Calculate dD/dT at constant pressure.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
dDdT : float
    dD/dT at constant P [(mol/L)/K]
"""
```

## DEFL1dll
```python
"""
Calculate temperature from density and energy.

Parameters
----------
D : float
    Density [mol/L]
e : float
    Internal energy [J/mol]
z : list
    Composition (mole fractions)

Returns
-------
T : float
    Temperature [K]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DEFLSHdll
```python
"""
Flash calculation given density and internal energy.

Parameters
----------
D : float
    Density [mol/L]
e : float
    Internal energy [J/mol]
z : list
    Composition (mole fractions)

Returns
-------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
Dl : float
    Liquid density [mol/L]
Dv : float
    Vapor density [mol/L]
x : list
    Liquid composition [mole fractions]
y : list
    Vapor composition [mole fractions]
q : float
    Vapor quality [mole basis]
h : float
    Enthalpy [J/mol]
s : float
    Entropy [J/(mol·K)]
Cv : float
    Isochoric heat capacity [J/(mol·K)]
Cp : float
    Isobaric heat capacity [J/(mol·K)]
w : float
    Speed of sound [m/s]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DERVPVTdll
```python
"""
Calculate thermodynamic derivatives from temperature, density, and composition.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
dPdD : float
    dP/dD at constant T [(kPa)/(mol/L)]
dPdT : float
    dP/dT at constant D [kPa/K]
d2PdD2 : float
    d²P/dD² at constant T [(kPa)/(mol/L)²]
d2PdT2 : float
    d²P/dT² at constant D [kPa/K²]
d2PdTD : float
    d²P/dTdD [(kPa)/((mol/L)·K)]
dDdP : float
    dD/dP at constant T [(mol/L)/kPa]
dDdT : float
    dD/dT at constant P [(mol/L)/K]
d2DdP2 : float
    d²D/dP² at constant T [(mol/L)/(kPa)²]
d2DdT2 : float
    d²D/dT² at constant P [(mol/L)/K²]
d2DdPT : float
    d²D/dPdT [(mol/L)/(kPa·K)]
dTdP : float
    dT/dP at constant D [K/kPa]
dTdD : float
    dT/dD at constant P [K/(mol/L)]
d2TdP2 : float
    d²T/dP² at constant D [K/(kPa)²]
d2TdD2 : float
    d²T/dD² at constant P [K/(mol/L)²]
d2TdPD : float
    d²T/dPdD [K/(kPa·(mol/L))]
"""
```

## DHD1dll
```python
"""
Calculate various enthalpy derivatives.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
dhdt_d : float
    dh/dT at constant D [J/(mol·K)]
dhdt_p : float
    dh/dT at constant P [J/(mol·K)]
dhdd_t : float
    dh/dD at constant T [(J/mol)/(mol/L)]
dhdd_p : float
    dh/dD at constant P [(J/mol)/(mol/L)]
dhdp_t : float
    dh/dP at constant T [(J/mol)/kPa]
dhdp_d : float
    dh/dP at constant D [(J/mol)/kPa]
"""
```

## DHFL1dll
```python
"""
Calculate temperature from density and enthalpy.

Parameters
----------
D : float
    Density [mol/L]
h : float
    Enthalpy [J/mol]
z : list
    Composition (mole fractions)

Returns
-------
T : float
    Temperature [K]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DHFLSHdll
```python
"""
Flash calculation given density and enthalpy.

Parameters
----------
D : float
    Density [mol/L]
h : float
    Enthalpy [J/mol]
z : list
    Composition (mole fractions)

Returns
-------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
Dl : float
    Liquid density [mol/L]
Dv : float
    Vapor density [mol/L]
x : list
    Liquid composition [mole fractions]
y : list
    Vapor composition [mole fractions]
q : float
    Vapor quality [mole basis]
e : float
    Internal energy [J/mol]
s : float
    Entropy [J/(mol·K)]
Cv : float
    Isochoric heat capacity [J/(mol·K)]
Cp : float
    Isobaric heat capacity [J/(mol·K)]
w : float
    Speed of sound [m/s]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DIELECdll
```python
"""
Calculate dielectric constant.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
de : float
    Dielectric constant [-]
"""
```

## DLSATKdll
```python
"""
Calculate saturated liquid density for a pure fluid.

Parameters
----------
icomp : int
    Component number (1-based)
T : float
    Temperature [K]

Returns
-------
D : float
    Saturated liquid density [mol/L]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DPDD2dll
```python
"""
Calculate d²P/dD² at constant temperature.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
d2PdD2 : float
    d²P/dD² at constant T [(kPa)/(mol/L)²]
"""
```

## DPDDdll
```python
"""
Calculate dP/dD at constant temperature.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
dPdD : float
    dP/dD at constant T [(kPa)/(mol/L)]
"""
```

## DPDTdll
```python
"""
Calculate dP/dT at constant density.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
dPdT : float
    dP/dT at constant D [kPa/K]
"""
```

## DPTSATKdll
```python
"""
Calculate saturation properties and dP/dT for a pure fluid.

Parameters
----------
icomp : int
    Component number (1-based)
T : float
    Temperature [K]
kph : int
    Phase flag: 1 = liquid, 2 = vapor

Returns
-------
P : float
    Saturation pressure [kPa]
D : float
    Saturation density [mol/L]
Csat : float
    Saturation heat capacity [J/(mol·K)]
dPdT : float
    dP/dT along saturation line [kPa/K]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DQFL2dll
```python
"""
Calculate two-phase state given density, quality, and composition.

Parameters
----------
D : float
    Density [mol/L]
q : float
    Vapor quality [mole basis]
z : list
    Composition (mole fractions)
kq : int
    Flag for quality definition

Returns
-------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
Dl : float
    Liquid density [mol/L]
Dv : float
    Vapor density [mol/L]
x : list
    Liquid composition [mole fractions]
y : list
    Vapor composition [mole fractions]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DSD1dll
```python
"""
Calculate various entropy derivatives.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
dsdt_d : float
    ds/dT at constant D [J/(mol·K²)]
dsdt_p : float
    ds/dT at constant P [J/(mol·K²)]
dsdd_t : float
    ds/dD at constant T [(J/(mol·K))/(mol/L)]
dsdd_p : float
    ds/dD at constant P [(J/(mol·K))/(mol/L)]
dsdp_t : float
    ds/dP at constant T [(J/(mol·K))/kPa]
dsdp_d : float
    ds/dP at constant D [(J/(mol·K))/kPa]
"""
```

## DSFL1dll
```python
"""
Calculate temperature from density and entropy.

Parameters
----------
D : float
    Density [mol/L]
s : float
    Entropy [J/(mol·K)]
z : list
    Composition (mole fractions)

Returns
-------
T : float
    Temperature [K]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DSFLSHdll
```python
"""
Flash calculation given density and entropy.

Parameters
----------
D : float
    Density [mol/L]
s : float
    Entropy [J/(mol·K)]
z : list
    Composition (mole fractions)

Returns
-------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
Dl : float
    Liquid density [mol/L]
Dv : float
    Vapor density [mol/L]
x : list
    Liquid composition [mole fractions]
y : list
    Vapor composition [mole fractions]
q : float
    Vapor quality [mole basis]
e : float
    Internal energy [J/mol]
h : float
    Enthalpy [J/mol]
Cv : float
    Isochoric heat capacity [J/(mol·K)]
Cp : float
    Isobaric heat capacity [J/(mol·K)]
w : float
    Speed of sound [m/s]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## DVSATKdll
```python
"""
Calculate saturated vapor density for a pure fluid.

Parameters
----------
icomp : int
    Component number (1-based)
T : float
    Temperature [K]

Returns
-------
D : float
    Saturated vapor density [mol/L]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## ENTHALdll
```python
"""
Calculate enthalpy.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
h : float
    Enthalpy [J/mol]
"""
```

## ENTROdll
```python
"""
Calculate entropy.

Parameters
----------
T : float
    Temperature [K]
D : float
    Density [mol/L]
z : list
    Composition (mole fractions)

Returns
-------
s : float
    Entropy [J/(mol·K)]
"""
```

## ERRMSGdll
```python
"""
Retrieve error message for a specified error code.

Parameters
----------
ierr : int
    Error code

Returns
-------
herr : str
    Error message
"""
```

## ESFLSHdll
```python
"""
Flash calculation given internal energy and entropy.

Parameters
----------
e : float
    Internal energy [J/mol]
s : float
    Entropy [J/(mol·K)]
z : list
    Composition (mole fractions)

Returns
-------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
D : float
    Overall density [mol/L]
Dl : float
    Liquid density [mol/L]
Dv : float
    Vapor density [mol/L]
x : list
    Liquid composition [mole fractions]
y : list
    Vapor composition [mole fractions]
q : float
    Vapor quality [mole basis]
h : float
    Enthalpy [J/mol]
Cv : float
    Isochoric heat capacity [J/(mol·K)]
Cp : float
    Isobaric heat capacity [J/(mol·K)]
w : float
    Speed of sound [m/s]
ierr : int
    Error code (0 = no error)
herr : str
    Error message
"""
```

## EXCESSdll
```python
"""
Calculate excess properties.

Parameters
----------
T : float
    Temperature [K]
P : float
    Pressure [kPa]
z : list
    Composition (mole fractions)
kph : int
    Phase flag: 1 = liquid, 2 = vapor
D : float
    Density [mol/L] (initial estimate)

Returns
-------
D : float
    Density [mol/L]
vE : float
    Excess volume [L/mol]
eE : float
    Excess internal energy [J/mol]
hE : float
    Excess enthalpy [J/mol]
sE : float
    Excess entropy [J/(mol·K)]