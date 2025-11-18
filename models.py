"""
TIMESAT-like Phenology Models
=============================
Function fitting models for vegetation time series analysis:
- Asymmetric Gaussian
- Double Logistic
- Robust optimization using least_squares
"""

import numpy as np
from scipy.optimize import least_squares


# ============================================================================
# FUNCTION FITTING MODELS
# ============================================================================

def asymmetric_gaussian(t, a, b, c, d1, d2):
    """
    Asymmetric Gaussian function (TIMESAT method)
    
    Parameters:
    -----------
    t : array-like
        Time indices
    a : float
        Amplitude
    b : float
        Position of peak
    c : float
        Base level
    d1 : float
        Width parameter (left side)
    d2 : float
        Width parameter (right side)
    
    Returns:
    --------
    array-like
        Fitted values
    """
    result = np.zeros_like(t, dtype=float)
    left_mask = t <= b
    right_mask = t > b
    
    result[left_mask] = c + a * np.exp(-((t[left_mask] - b) ** 2) / (2 * d1 ** 2))
    result[right_mask] = c + a * np.exp(-((t[right_mask] - b) ** 2) / (2 * d2 ** 2))
    
    return result


def double_logistic(t, a, b, c, m1, m2, n1, n2):
    """
    Double Logistic function (TIMESAT method)
    
    Parameters:
    -----------
    t : array-like
        Time indices
    a : float
        Amplitude
    b : float
        Base level
    c : float
        Position of inflection
    m1, m2 : float
        Rate parameters for increase/decrease
    n1, n2 : float
        Shape parameters
    
    Returns:
    --------
    array-like
        Fitted values
    """
    with np.errstate(over='ignore', invalid='ignore'):
        left = 1 / (1 + np.exp(m1 * (t - c + n1)))
        right = 1 / (1 + np.exp(-m2 * (t - c - n2)))
        fitted = b + a * left * right
    
    return fitted


def fit_asymmetric_gaussian(t_indices, y_values, max_iterations=50):
    """
    Fit Asymmetric Gaussian to data using robust least_squares optimization
    
    Parameters:
    -----------
    t_indices : array-like
        Time indices
    y_values : array-like
        Observed values
    max_iterations : int
        Maximum iterations for optimization
    
    Returns:
    --------
    tuple
        (fitted_values, parameters, success)
    """
    try:
        # Initial parameter estimation
        y_min, y_max = np.nanmin(y_values), np.nanmax(y_values)
        peak_idx = np.nanargmax(y_values)
        
        # Initial guess
        p0 = [
            y_max - y_min,              # amplitude
            t_indices[peak_idx],        # peak position
            y_min,                      # base level
            len(t_indices) / 6,         # left width
            len(t_indices) / 6          # right width
        ]
        
        # Bounds
        bounds = (
            [0, 0, 0, 1, 1],
            [2 * (y_max - y_min), len(t_indices), y_max, 
             len(t_indices), len(t_indices)]
        )
        
        # Residual function for optimization
        def residual(p):
            return asymmetric_gaussian(t_indices, *p) - y_values
        
        # Robust fitting using least_squares
        res = least_squares(
            residual, 
            p0, 
            bounds=bounds, 
            max_nfev=max_iterations * 100
        )
        
        if res.success:
            fitted = asymmetric_gaussian(t_indices, *res.x)
            return fitted, res.x, True
        else:
            raise ValueError("Optimization failed")
    
    except Exception as e:
        print(f"      ⚠ Asymmetric Gaussian fit failed: {e}")
        return y_values, None, False


def fit_double_logistic(t_indices, y_values, max_iterations=50):
    """
    Fit Double Logistic to data using robust least_squares optimization
    
    Parameters:
    -----------
    t_indices : array-like
        Time indices
    y_values : array-like
        Observed values
    max_iterations : int
        Maximum iterations for optimization
    
    Returns:
    --------
    tuple
        (fitted_values, parameters, success)
    """
    try:
        y_min, y_max = np.nanmin(y_values), np.nanmax(y_values)
        peak_idx = np.nanargmax(y_values)
        
        # Initial guess
        p0 = [
            y_max - y_min,              # amplitude
            y_min,                      # base level
            t_indices[peak_idx],        # inflection point
            0.1, 0.1,                   # rate parameters
            len(t_indices) / 6,         # shape parameters
            len(t_indices) / 6
        ]
        
        # Bounds
        bounds = (
            [0, 0, 0, 0.01, 0.01, 1, 1],
            [2 * (y_max - y_min), y_max, len(t_indices), 
             1, 1, len(t_indices), len(t_indices)]
        )
        
        # Residual function
        def residual(p):
            return double_logistic(t_indices, *p) - y_values
        
        # Robust fitting
        res = least_squares(
            residual, 
            p0, 
            bounds=bounds, 
            max_nfev=max_iterations * 100
        )
        
        if res.success:
            fitted = double_logistic(t_indices, *res.x)
            return fitted, res.x, True
        else:
            raise ValueError("Optimization failed")
    
    except Exception as e:
        print(f"      ⚠ Double Logistic fit failed: {e}")
        return y_values, None, False


def fit_function_to_timeseries(t_indices, y_values, method='asymmetric_gaussian'):
    """
    Wrapper function for selecting fitting method
    
    Parameters:
    -----------
    t_indices : array-like
        Time indices
    y_values : array-like
        Observed values
    method : str
        'asymmetric_gaussian' or 'double_logistic'
    
    Returns:
    --------
    tuple
        (fitted_values, parameters, success)
    """
    if method == 'asymmetric_gaussian':
        return fit_asymmetric_gaussian(t_indices, y_values)
    elif method == 'double_logistic':
        return fit_double_logistic(t_indices, y_values)
    else:
        return y_values, None, False
