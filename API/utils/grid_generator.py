"""
Utility functions for generating calculation grids with different distribution strategies.
"""

import numpy as np
from typing import List, Dict, Tuple, Union, Optional, Any

def generate_grid(
    range_min: float, 
    range_max: float, 
    resolution: float, 
    grid_type: str = "equidistant",
    boundaries: Optional[List[float]] = None, 
    enhancement_factor: float = 5.0,
    boundary_zone_width: Optional[float] = None,
    **kwargs
) -> np.ndarray:
    """
    Generate a calculation grid with the specified distribution strategy.
    
    Args:
        range_min: Minimum value of the range
        range_max: Maximum value of the range  
        resolution: Base resolution (spacing between points)
        grid_type: Type of grid to generate:
                  - "equidistant": Regular spacing (default)
                  - "adaptive": Higher resolution near phase boundaries
                  - "logarithmic": Logarithmic spacing (more points at lower values)
                  - "exponential": Exponential spacing (more points at higher values)
        boundaries: List of values where phase boundaries occur (for adaptive grid)
        enhancement_factor: How much to increase resolution near boundaries
        boundary_zone_width: Width of zone around boundary where resolution is enhanced
                            (If None, it's calculated as 10% of the total range)
        
    Returns:
        numpy array of grid points
    """
    if grid_type.lower() == "equidistant":
        return np.arange(range_min, range_max + resolution, resolution)
    
    elif grid_type.lower() == "adaptive":
        if not boundaries:
            # If no boundaries provided, fall back to equidistant
            return np.arange(range_min, range_max + resolution, resolution)
        
        return generate_adaptive_grid(
            range_min, range_max, resolution, boundaries, 
            enhancement_factor, boundary_zone_width
        )
    
    elif grid_type.lower() == "logarithmic":
        # Ensure range doesn't include zero or negative values for logarithmic grid
        if range_min <= 0:
            min_val = 0.001 * range_max if range_max > 0 else 0.001
        else:
            min_val = range_min
        
        # Generate logarithmically spaced points
        log_min = np.log10(min_val)
        log_max = np.log10(range_max)
        num_points = int(np.ceil((log_max - log_min) * range_max / resolution))
        
        # Ensure at least 2 points
        num_points = max(num_points, 2)
        
        return np.logspace(log_min, log_max, num_points)
    
    elif grid_type.lower() == "exponential":
        # Generate exponentially spaced points (more points at higher values)
        # Calculate number of points to get approximately the requested resolution
        range_span = range_max - range_min
        num_points = int(np.ceil(range_span / resolution))
        
        # Ensure at least 2 points
        num_points = max(num_points, 2)
        
        # Exponential spacing factor - higher values create more bias towards end
        exponent = kwargs.get('exponent', 2)
        
        # Generate normalized points between 0 and 1 with exponential distribution
        normalized = np.power(np.linspace(0, 1, num_points), exponent)
        
        # Scale to the actual range
        return range_min + normalized * range_span
    
    else:
        # Default to equidistant grid for unknown grid_type
        return np.arange(range_min, range_max + resolution, resolution)

def generate_adaptive_grid(
    range_min: float, 
    range_max: float, 
    base_resolution: float, 
    boundaries: List[float], 
    enhancement_factor: float = 5.0, 
    boundary_zone_width: Optional[float] = None
) -> np.ndarray:
    """
    Generate an adaptive grid with higher resolution near phase boundaries.
    
    Args:
        range_min: Minimum value of the range
        range_max: Maximum value of the range
        base_resolution: Base spacing between points
        boundaries: List of values where phase boundaries occur
        enhancement_factor: How much to increase resolution near boundaries
        boundary_zone_width: Width of zone around boundary where resolution is enhanced
                            (If None, it's calculated as 10% of the total range)
        
    Returns:
        numpy array of grid points
    """
    # Calculate boundary zone width if not specified
    if boundary_zone_width is None:
        boundary_zone_width = 0.1 * (range_max - range_min)
    
    # Start with basic grid
    basic_grid = np.arange(range_min, range_max + base_resolution, base_resolution)
    
    # For each phase boundary, add more points in its vicinity
    enhanced_points = []
    for boundary in boundaries:
        if range_min <= boundary <= range_max:
            zone_min = max(range_min, boundary - boundary_zone_width)
            zone_max = min(range_max, boundary + boundary_zone_width)
            
            # Create higher resolution grid in this zone
            fine_resolution = base_resolution / enhancement_factor
            zone_grid = np.arange(zone_min, zone_max + fine_resolution, fine_resolution)
            
            enhanced_points.extend(zone_grid)
    
    # Combine all points
    all_points = np.concatenate([basic_grid, enhanced_points])
    
    # Sort and remove duplicates (within a small tolerance)
    tolerance = base_resolution / (enhancement_factor * 10)
    
    # Sort the points
    sorted_points = np.sort(all_points)
    
    # Remove near-duplicates
    result = [sorted_points[0]]
    for i in range(1, len(sorted_points)):
        if sorted_points[i] - result[-1] >= tolerance:
            result.append(sorted_points[i])
    
    return np.array(result)

def get_phase_boundaries_pt(rp, z: List[float], t_range: Dict[str, float], p_range: Dict[str, float]) -> Tuple[List[float], List[float]]:
    """
    Determine phase boundaries in P-T space for a given composition.
    
    Args:
        rp: REFPROP instance
        z: Composition array
        t_range: Temperature range dictionary {'from': min, 'to': max} in °C
        p_range: Pressure range dictionary {'from': min, 'to': max} in bar
        
    Returns:
        t_boundaries: List of temperatures (°C) at phase boundaries
        p_boundaries: List of pressures (bar) at phase boundaries
    """
    t_boundaries = []
    p_boundaries = []
    
    try:
        # Get range values
        t_min, t_max = t_range['from'], t_range['to']
        p_min, p_max = p_range['from'], p_range['to']
        
        # Get critical point
        Tc, Pc, Dc, ierr, herr = rp.CRITPdll(z)
        
        if ierr == 0:
            Tc_C = Tc - 273.15  # Convert to Celsius
            Pc_bar = Pc / 100   # Convert to bar
            
            # Check if critical point is within range
            if t_min <= Tc_C <= t_max and p_min <= Pc_bar <= p_max:
                t_boundaries.append(Tc_C)
                p_boundaries.append(Pc_bar)
        
        # Create temperature points for calculating saturation curve
        num_t_points = min(50, int((t_max - t_min) / 2) + 1)  # Reasonable number of points
        t_sat_points = np.linspace(t_min, min(Tc_C if ierr == 0 else t_max, t_max), num_t_points)
        
        # Calculate saturation points at each temperature
        for T_C in t_sat_points:
            T_K = T_C + 273.15  # Convert to Kelvin
            
            try:
                # Get saturation pressure at this temperature (bubble point)
                result = rp.SATTdll(T_K, z, 1)
                
                if result.ierr == 0:
                    P_bar = result.P / 100  # Convert kPa to bar
                    
                    # Check if within pressure range
                    if p_min <= P_bar <= p_max:
                        t_boundaries.append(T_C)
                        p_boundaries.append(P_bar)
            except Exception:
                # Skip failed calculations
                pass
        
    except Exception as e:
        print(f"Error determining PT phase boundaries: {str(e)}")
    
    return t_boundaries, p_boundaries

def get_phase_boundaries_ph(rp, z: List[float], p_range: Dict[str, float], h_range: Dict[str, float]) -> Tuple[List[float], List[float]]:
    """
    Determine phase boundaries in P-H space for a given composition.
    
    Args:
        rp: REFPROP instance
        z: Composition array
        p_range: Pressure range dictionary {'from': min, 'to': max} in bar
        h_range: Enthalpy range dictionary {'from': min, 'to': max} in J/mol
        
    Returns:
        p_boundaries: List of pressures (bar) at phase boundaries
        h_boundaries: List of enthalpies (J/mol) at phase boundaries
    """
    p_boundaries = []
    h_boundaries = []
    
    try:
        # Get range values
        p_min, p_max = p_range['from'], p_range['to']
        h_min, h_max = h_range['from'], h_range['to']
        
        # Get critical point
        Tc, Pc, Dc, ierr, herr = rp.CRITPdll(z)
        
        if ierr == 0:
            Pc_bar = Pc / 100  # Convert to bar
            
            # Check if critical pressure is within range
            if p_min <= Pc_bar <= p_max:
                # Calculate enthalpy at critical point
                try:
                    # Use TPFLSHdll to get properties at critical point
                    result = rp.TPFLSHdll(Tc, Pc, z)
                    if result.ierr == 0:
                        h_crit = result.h
                        
                        # Check if critical enthalpy is within range
                        if h_min <= h_crit <= h_max:
                            p_boundaries.append(Pc_bar)
                            h_boundaries.append(h_crit)
                except Exception:
                    pass
        
        # Create pressure points for calculating the phase envelope
        num_p_points = min(50, int((p_max - p_min) / 5) + 1)
        p_points = np.linspace(p_min, p_max, num_p_points)
        
        # For each pressure, find bubble and dew point enthalpies
        for P_bar in p_points:
            P_kpa = P_bar * 100  # Convert to kPa for REFPROP
            
            try:
                # Calculate saturation properties at this pressure
                result = rp.SATPdll(P_kpa, z, 1)  # Bubble point
                
                if result.ierr == 0:
                    T_sat = result.T
                    Dl = result.Dl
                    
                    # Calculate enthalpy at bubble point
                    props = rp.THERMdll(T_sat, Dl, z)
                    h_bubble = props.h
                    
                    # Check if within enthalpy range
                    if h_min <= h_bubble <= h_max:
                        p_boundaries.append(P_bar)
                        h_boundaries.append(h_bubble)
            except Exception:
                pass
                
            try:
                # Calculate dew point properties
                result = rp.SATPdll(P_kpa, z, 2)  # Dew point
                
                if result.ierr == 0:
                    T_sat = result.T
                    Dv = result.Dv
                    
                    # Calculate enthalpy at dew point
                    props = rp.THERMdll(T_sat, Dv, z)
                    h_dew = props.h
                    
                    # Check if within enthalpy range
                    if h_min <= h_dew <= h_max:
                        p_boundaries.append(P_bar)
                        h_boundaries.append(h_dew)
            except Exception:
                pass
        
    except Exception as e:
        print(f"Error determining PH phase boundaries: {str(e)}")
    
    return p_boundaries, h_boundaries

def get_phase_boundaries_ts(rp, z: List[float], t_range: Dict[str, float], s_range: Dict[str, float]) -> Tuple[List[float], List[float]]:
    """
    Determine phase boundaries in T-S space for a given composition.
    
    Args:
        rp: REFPROP instance
        z: Composition array
        t_range: Temperature range dictionary {'from': min, 'to': max} in °C
        s_range: Entropy range dictionary {'from': min, 'to': max} in J/(mol·K)
        
    Returns:
        t_boundaries: List of temperatures (°C) at phase boundaries
        s_boundaries: List of entropies (J/(mol·K)) at phase boundaries
    """
    t_boundaries = []
    s_boundaries = []
    
    try:
        # Get range values
        t_min, t_max = t_range['from'], t_range['to']
        s_min, s_max = s_range['from'], s_range['to']
        
        # Get critical point
        Tc, Pc, Dc, ierr, herr = rp.CRITPdll(z)
        
        if ierr == 0:
            Tc_C = Tc - 273.15  # Convert to Celsius
            
            # Check if critical temperature is within range
            if t_min <= Tc_C <= t_max:
                # Calculate entropy at critical point
                try:
                    # Use TPFLSHdll to get properties at critical point
                    result = rp.TPFLSHdll(Tc, Pc, z)
                    if result.ierr == 0:
                        s_crit = result.s
                        
                        # Check if critical entropy is within range
                        if s_min <= s_crit <= s_max:
                            t_boundaries.append(Tc_C)
                            s_boundaries.append(s_crit)
                except Exception:
                    pass
        
        # Create temperature points for calculating the phase envelope
        num_t_points = min(50, int((t_max - t_min) / 2) + 1)
        t_points = np.linspace(t_min, min(Tc_C if ierr == 0 else t_max, t_max), num_t_points)
        
        # For each temperature, find bubble and dew point entropies
        for T_C in t_points:
            T_K = T_C + 273.15  # Convert to Kelvin
            
            try:
                # Calculate saturation properties at this temperature
                result = rp.SATTdll(T_K, z, 1)  # Bubble point
                
                if result.ierr == 0:
                    # Get properties for liquid at saturation
                    Dl = result.Dl
                    props = rp.ENTROdll(T_K, Dl, z)
                    s_bubble = props  # ENTROdll returns entropy directly
                    
                    # Check if within entropy range
                    if s_min <= s_bubble <= s_max:
                        t_boundaries.append(T_C)
                        s_boundaries.append(s_bubble)
            except Exception:
                pass
                
            try:
                # Calculate dew point properties
                result = rp.SATTdll(T_K, z, 2)  # Dew point
                
                if result.ierr == 0:
                    # Get properties for vapor at saturation
                    Dv = result.Dv
                    props = rp.ENTROdll(T_K, Dv, z)
                    s_dew = props  # ENTROdll returns entropy directly
                    
                    # Check if within entropy range
                    if s_min <= s_dew <= s_max:
                        t_boundaries.append(T_C)
                        s_boundaries.append(s_dew)
            except Exception:
                pass
        
    except Exception as e:
        print(f"Error determining TS phase boundaries: {str(e)}")
    
    return t_boundaries, s_boundaries