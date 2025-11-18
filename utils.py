"""
TIMESAT Utility Functions
=========================
Helper functions for:
- Data loading and preprocessing
- Vegetation indices calculation
- Time series smoothing
- Region extraction
- Plotting and visualization
"""

import os
import re
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
from scipy.signal import savgol_filter, find_peaks
from scipy import ndimage
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


# ============================================================================
# FILE HANDLING UTILITIES
# ============================================================================

def default_time_parser_from_filename(fname):
    """
    Extract date from filename using common patterns
    
    Parameters:
    -----------
    fname : str
        Filename to parse
    
    Returns:
    --------
    datetime
        Extracted date
    """
    patterns = [
        r'(\d{4})(\d{2})(\d{2})',      # YYYYMMDD
        r'(\d{4})-(\d{2})-(\d{2})',    # YYYY-MM-DD
        r'(\d{4})_(\d{2})_(\d{2})'     # YYYY_MM_DD
    ]
    
    for p in patterns:
        m = re.search(p, fname)
        if m:
            y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
            return datetime(y, mo, d)
    
    # Fallback to file modification time if no date found
    try:
        mtime = os.path.getmtime(fname)
        return datetime.fromtimestamp(mtime)
    except:
        raise ValueError(f"No date found in filename: {fname}")


def read_band_from_multiband_tif(fp, band_index=1):
    """
    Read specific band from multi-band GeoTIFF
    
    Parameters:
    -----------
    fp : str
        File path
    band_index : int
        Band index (1-based)
    
    Returns:
    --------
    xarray.DataArray
        Band data
    """
    da = rxr.open_rasterio(fp, masked=True)
    
    if 'band' in da.dims:
        try:
            return da.isel(band=band_index - 1)
        except Exception:
            return da.isel(band=0)
    
    return da


def stack_multiband_files(file_list, band_index=1):
    """
    Stack multiple multi-band files along time dimension
    
    Parameters:
    -----------
    file_list : list
        List of file paths
    band_index : int
        Band index to extract
    
    Returns:
    --------
    xarray.DataArray
        Stacked data with time dimension
    """
    das = []
    
    for fp in sorted(file_list):
        try:
            # Read band
            da = read_band_from_multiband_tif(fp, band_index)
            
            # Extract time from filename
            t = pd.to_datetime(default_time_parser_from_filename(os.path.basename(fp)))
            
            # Ensure standard dimension names
            if 'y' not in da.dims or 'x' not in da.dims:
                dims = list(da.dims)
                if len(dims) >= 2:
                    new = dims[:]
                    new[-2] = 'y'
                    new[-1] = 'x'
                    da = da.rename(dict(zip(dims, new)))
            
            das.append(da.expand_dims(time=[t]))
            
        except Exception as e:
            print(f"  ⚠ Could not process {fp}: {e}")
    
    if len(das) == 0:
        raise ValueError("No readable files")
    
    # Combine along time dimension
    combined = xr.concat(das, dim="time").sortby('time')
    
    # Remove redundant band dimension if exists
    if 'band' in combined.dims and combined.sizes['band'] == 1:
        combined = combined.isel(band=0)
    
    return combined


# ============================================================================
# VEGETATION INDICES CALCULATION
# ============================================================================

def compute_ndvi(nir_da, red_da):
    """
    Compute Normalized Difference Vegetation Index
    
    NDVI = (NIR - RED) / (NIR + RED)
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        ndvi = (nir_da - red_da) / (nir_da + red_da)
    return ndvi.clip(-1, 1)


def compute_savi(nir_da, red_da, L=0.5):
    """
    Compute Soil Adjusted Vegetation Index
    
    SAVI = ((NIR - RED) / (NIR + RED + L)) * (1 + L)
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        savi = ((nir_da - red_da) / (nir_da + red_da + L)) * (1 + L)
    return savi.clip(-1, 1)


def compute_msavi2(nir_da, red_da):
    """
    Compute Modified Soil Adjusted Vegetation Index 2
    
    MSAVI2 = 0.5 * (2*NIR + 1 - sqrt((2*NIR + 1)^2 - 8*(NIR - RED)))
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        msavi2 = 0.5 * (2 * nir_da + 1 - np.sqrt((2 * nir_da + 1) ** 2 - 8 * (nir_da - red_da)))
    return msavi2.clip(-1, 1)


def compute_gndvi(nir_da, green_da):
    """
    Compute Green Normalized Difference Vegetation Index
    
    GNDVI = (NIR - GREEN) / (NIR + GREEN)
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        gndvi = (nir_da - green_da) / (nir_da + green_da)
    return gndvi.clip(-1, 1)


def compute_rvi(nir_da, red_da):
    """
    Compute Ratio Vegetation Index
    
    RVI = NIR / RED
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        rvi = nir_da / red_da
    return rvi.clip(0, np.inf)


# ============================================================================
# TIME SERIES SMOOTHING
# ============================================================================

def remove_outliers(arr, z_thresh=2.5):
    """
    Remove outliers using Modified Z-score (MAD-based)
    
    Parameters:
    -----------
    arr : array-like
        Input array
    z_thresh : float
        Z-score threshold
    
    Returns:
    --------
    array-like
        Array with outliers replaced by NaN
    """
    a = np.array(arr, dtype=float)
    
    if np.all(np.isnan(a)):
        return a
    
    # Compute median and MAD
    med = np.nanmedian(a)
    mad = np.nanmedian(np.abs(a - med))
    
    if mad == 0 or np.isnan(mad):
        return a
    
    # Modified Z-score
    z = 0.6745 * (a - med) / mad
    a[np.abs(z) > z_thresh] = np.nan
    
    return a


def ensure_odd_positive_window(window_length, n):
    """
    Ensure window length is odd and valid
    
    Parameters:
    -----------
    window_length : int
        Desired window length
    n : int
        Array length
    
    Returns:
    --------
    int or None
        Valid window length
    """
    if n < 3:
        return None
    
    w = int(window_length)
    
    # Ensure window fits in array
    if w >= n:
        w = n - 1
    
    # Ensure window is odd
    if w % 2 == 0:
        w -= 1
    
    # Ensure minimum window size
    if w < 3:
        w = 3 if n >= 3 else None
    
    return w


def double_savgol_smooth(ts, window_length=7, polyorder=2):
    """
    Apply double Savitzky-Golay smoothing with outlier removal
    
    Parameters:
    -----------
    ts : array-like
        Time series to smooth
    window_length : int
        Window length (must be odd)
    polyorder : int
        Polynomial order
    
    Returns:
    --------
    array-like
        Smoothed time series
    """
    s = pd.Series(ts)
    
    # Check for all NaN
    if s.isna().all():
        return np.full_like(ts, np.nan, dtype="float64")
    
    # Interpolate missing values
    s_interp = s.interpolate(limit=3, limit_direction="both").ffill().bfill()
    arr = s_interp.values.astype("float64")
    
    # Remove outliers
    arr = remove_outliers(arr)
    
    # Validate window
    w = ensure_odd_positive_window(window_length, len(arr))
    if w is None or w < 3:
        return arr
    
    try:
        # First smoothing pass
        smooth1 = savgol_filter(arr, w, polyorder, mode='interp')
        
        # Second smoothing pass
        smooth2 = savgol_filter(smooth1, w, polyorder, mode='interp')
        
        return smooth2
    
    except Exception:
        return arr


# ============================================================================
# REGION EXTRACTION
# ============================================================================

def extract_regions_from_label(label_array, label_value, min_pixels=5):
    """
    Extract connected regions for specific label value
    
    Parameters:
    -----------
    label_array : array-like
        Label array
    label_value : int
        Target label value
    min_pixels : int
        Minimum pixels for valid region
    
    Returns:
    --------
    dict
        Dictionary of regions {region_id: {'mask': mask, 'n_pixels': count}}
    """
    # Create binary mask for target label
    binary_mask = (label_array == label_value)
    
    if not np.any(binary_mask):
        return {}
    
    # Label connected components
    labeled_regions, num_regions = ndimage.label(binary_mask)
    
    regions = {}
    
    for region_id in range(1, num_regions + 1):
        region_mask = (labeled_regions == region_id)
        n_pixels = int(np.sum(region_mask))
        
        # Filter by minimum size
        if n_pixels >= min_pixels:
            regions[region_id] = {
                'mask': region_mask,
                'n_pixels': n_pixels
            }
    
    return regions


def load_label_for_tile(label_folder, tile_id, match_crs_da=None):
    """
    Load label file for specific tile
    
    Parameters:
    -----------
    label_folder : str
        Path to label folder
    tile_id : str
        Tile identifier
    match_crs_da : xarray.DataArray, optional
        DataArray to match CRS and grid
    
    Returns:
    --------
    xarray.DataArray
        Label data
    """
    # Find matching label files
    candidates = []
    for f in os.listdir(label_folder):
        if not f.lower().endswith(".tif"):
            continue
        if f.startswith(f"{tile_id}_") or f"_{tile_id}_" in f or f"_{tile_id}." in f:
            candidates.append(os.path.join(label_folder, f))
    
    if len(candidates) == 0:
        raise FileNotFoundError(f"No label file for tile '{tile_id}'")
    
    if len(candidates) > 1:
        print(f"  Warning: Multiple labels found for {tile_id}, using first: {candidates[0]}")
    
    label_fp = sorted(candidates)[0]
    print(f"  Loading label: {os.path.basename(label_fp)}")
    
    # Load label
    label_da = rxr.open_rasterio(label_fp, masked=True)
    
    # Remove band dimension if exists
    if 'band' in label_da.dims:
        label_da = label_da.isel(band=0)
    
    # Reproject to match reference if provided
    if match_crs_da is not None:
        try:
            label_da = label_da.rio.reproject_match(match_crs_da)
        except Exception as e:
            print(f"  ⚠ Reproject failed: {e}")
    
    # Standardize dimension names
    if 'y' not in label_da.dims or 'x' not in label_da.dims:
        dims = list(label_da.dims)
        if len(dims) >= 2:
            new = dims[:]
            new[-2] = 'y'
            new[-1] = 'x'
            label_da = label_da.rename(dict(zip(dims, new)))
    
    # Convert to integer array
    label_arr = label_da.values
    if np.ma.is_masked(label_arr):
        lab = np.array(label_arr.filled(0), dtype=np.int32)
    else:
        lab = np.array(label_arr, dtype=np.float64)
        lab[np.isnan(lab)] = 0
        lab[np.isinf(lab)] = 0
        lab = lab.astype(np.int32)
    
    # Print unique values
    unique_vals = np.unique(lab)
    print(f"  Label unique values: {unique_vals[:10]}{'...' if len(unique_vals) > 10 else ''}")
    
    # Filter anomalous large values
    if np.any(np.abs(lab) > 10000):
        print(f"  WARNING: Filtering large values...")
        lab[np.abs(lab) > 10000] = 0
    
    return xr.DataArray(
        lab,
        coords={"y": label_da.coords["y"], "x": label_da.coords["x"]},
        dims=("y", "x")
    )


# ============================================================================
# MULTIPLE SEASONS DETECTION
# ============================================================================

def detect_multiple_seasons(ts_smooth, times, min_peak_prominence=0.1,
                            min_peak_distance=15, min_amplitude=0.05):
    """
    Detect multiple growing seasons in time series (corrected to avoid IndexError)
    
    Parameters:
    -----------
    ts_smooth : array-like
        Smoothed time series
    times : array-like
        Time points
    min_peak_prominence : float
        Minimum peak prominence
    min_peak_distance : int
        Minimum distance between peaks
    min_amplitude : float
        Minimum amplitude for valid season
    
    Returns:
    --------
    list
        List of detected seasons with metadata
    """
    if np.all(np.isnan(ts_smooth)):
        return []
    
    valid_idx = ~np.isnan(ts_smooth)
    if np.sum(valid_idx) < min_peak_distance:
        return []
    
    # Find peaks
    peaks, properties = find_peaks(
        ts_smooth,
        prominence=min_peak_prominence,
        distance=min_peak_distance
    )
    
    if len(peaks) == 0:
        return []
    
    print(f"      Found {len(peaks)} potential peak(s) at indices: {peaks}")
    
    seasons = []
    base_value = np.nanmin(ts_smooth)
    
    for i, peak_idx in enumerate(peaks):
        peak_val = ts_smooth[peak_idx]
        amplitude = peak_val - base_value
        
        # Check minimum amplitude
        if amplitude < min_amplitude:
            print(f"      Skipping peak {i + 1}: amplitude {amplitude:.3f} < {min_amplitude}")
            continue
        
        # Get previous season end (safe)
        prev_end = seasons[-1]['end_idx'] if len(seasons) > 0 else -1
        
        # Find start: scan backwards from peak
        start_idx = None
        j_start = max(peak_idx - 1, prev_end + 1)
        
        for j in range(j_start, max(prev_end, 0) - 1, -1):
            # Check boundaries
            if j - 1 < 0 or j + 1 >= len(ts_smooth):
                start_idx = max(j, 0)
                break
            
            # Local minimum heuristic
            if (not np.isnan(ts_smooth[j - 1]) and 
                not np.isnan(ts_smooth[j + 1]) and 
                not np.isnan(ts_smooth[j])):
                if (ts_smooth[j] <= ts_smooth[j - 1] and 
                    ts_smooth[j] <= ts_smooth[j + 1]):
                    start_idx = j
                    break
        
        if start_idx is None:
            start_idx = prev_end + 1 if (prev_end + 1) < peak_idx else max(0, peak_idx - 1)
        
        # Find end: scan forward from peak
        end_idx = None
        j_end_limit = len(ts_smooth) - 1
        
        for j in range(peak_idx + 1, j_end_limit):
            # Check boundaries
            if j - 1 < 0 or j + 1 >= len(ts_smooth):
                end_idx = min(j_end_limit, j)
                break
            
            # Local minimum heuristic
            if (not np.isnan(ts_smooth[j - 1]) and 
                not np.isnan(ts_smooth[j + 1]) and 
                not np.isnan(ts_smooth[j])):
                if (ts_smooth[j] <= ts_smooth[j - 1] and 
                    ts_smooth[j] <= ts_smooth[j + 1]):
                    end_idx = j
                    break
        
        if end_idx is None:
            end_idx = j_end_limit
        
        # Validate season: no overlap, minimum width
        if (end_idx > start_idx and 
            (end_idx - start_idx) >= min_peak_distance and 
            start_idx > prev_end):
            
            season = {
                'season_number': len(seasons) + 1,
                'start_idx': int(start_idx),
                'end_idx': int(end_idx),
                'peak_idx': int(peak_idx),
                'peak_value': float(peak_val),
                'amplitude': float(amplitude),
                'base_value': float(base_value)
            }
            seasons.append(season)
            print(f"      ✓ Season {len(seasons)}: idx {start_idx}-{end_idx}, peak={peak_val:.3f}")
        else:
            print(f"      Skipping peak at idx {peak_idx}: invalid range or overlap")
    
    return seasons


# ============================================================================
# PLOTTING AND VISUALIZATION
# ============================================================================

def plot_multi_season_phenology(df_mean, times, phenology_list, smooth_series,
                                label_id, region_id, tile_id, output_path):
    """
    Create advanced multi-season phenology plot
    
    Parameters:
    -----------
    df_mean : pd.DataFrame
        Mean time series data
    times : array-like
        Time points
    phenology_list : list
        List of phenology dictionaries per season
    smooth_series : array-like
        Smoothed NDVI values
    label_id : int
        Label identifier
    region_id : int
        Region identifier
    tile_id : str
        Tile identifier
    output_path : str
        Output file path
    """
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(3, 1, height_ratios=[4, 2, 1], hspace=0.35)
    
    # Main time series plot
    ax1 = fig.add_subplot(gs[0])
    times_pd = pd.to_datetime(times)
    
    # Plot raw and smoothed NDVI
    ax1.plot(times_pd, df_mean['NDVI'], 'o', alpha=0.4, markersize=5,
             color='gray', label='Raw NDVI')
    ax1.plot(times_pd, smooth_series, '-', linewidth=3,
             color='darkgreen', label='Smoothed NDVI')
    
    # Color palette for seasons
    season_colors = ['red', 'blue', 'orange', 'purple', 'brown']
    
    # Plot phenology markers for each season
    for i, ph in enumerate(phenology_list):
        season_num = ph.get('season_number', i + 1)
        color = season_colors[i % len(season_colors)]
        
        def doy_to_date(doy):
            """Convert DOY to date"""
            if doy is None or np.isnan(doy):
                return None
            doy_int = int(round(doy))
            doy_series = times_pd.dayofyear
            idx = int((np.abs(doy_series - doy_int)).argmin())
            return times_pd[idx]
        
        # Start of Season
        sos_date = doy_to_date(ph.get('SOS_DOY'))
        if sos_date:
            ax1.axvline(sos_date, color=color, linestyle='--',
                        linewidth=2, alpha=0.6, label=f'S{season_num} SOS')
        
        # Peak
        peak_date = doy_to_date(ph.get('Peak_DOY'))
        if peak_date:
            ax1.axvline(peak_date, color=color, linestyle='-',
                        linewidth=2.5, alpha=0.8)
            ax1.plot(peak_date, ph.get('Peak_value', np.nan),
                     marker='*', markersize=18, color=color,
                     markeredgecolor='black', markeredgewidth=1.5,
                     label=f'S{season_num} Peak')
        
        # End of Season
        eos_date = doy_to_date(ph.get('EOS_DOY'))
        if eos_date:
            ax1.axvline(eos_date, color=color, linestyle=':',
                        linewidth=2, alpha=0.6, label=f'S{season_num} EOS')
    
    ax1.set_ylabel('NDVI', fontsize=13, fontweight='bold')
    ax1.set_title(f'Multi-Season Phenology: Label {label_id}, Region {region_id} (Tile {tile_id})',
                  fontsize=15, fontweight='bold')
    ax1.legend(loc='upper left', framealpha=0.9, 
               ncol=min(3, len(phenology_list) + 1), fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    # Season comparison table
    ax2 = fig.add_subplot(gs[1])
    ax2.axis('off')
    
    if len(phenology_list) > 0:
        table_text = "Season-by-Season Metrics:\n\n"
        table_text += f"{'Season':<8} {'SOS':>8} {'Peak':>8} {'EOS':>8} {'Amp':>8} {'LOS':>8} {'Method':<12}\n"
        table_text += "─" * 75 + "\n"
        
        for ph in phenology_list:
            sn = ph.get('season_number', '?')
            sos = ph.get('SOS_DOY', np.nan)
            peak = ph.get('Peak_DOY', np.nan)
            eos = ph.get('EOS_DOY', np.nan)
            amp = ph.get('Amplitude', np.nan)
            los = ph.get('LOS', np.nan)
            method = ph.get('fit_method', 'none')
            
            def fmt(v):
                return f"{v:.1f}" if not np.isnan(v) else "nan"
            
            table_text += f"{sn:<8} {fmt(sos):>8} {fmt(peak):>8} {fmt(eos):>8} {fmt(amp):>8} {fmt(los):>8} {method:<12}\n"
        
        ax2.text(0.05, 0.95, table_text, transform=ax2.transAxes,
                 fontsize=9, verticalalignment='top', fontfamily='monospace',
                 bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    
    # Detailed metrics for first season
    ax3 = fig.add_subplot(gs[2])
    ax3.axis('off')
    
    if len(phenology_list) > 0:
        ph = phenology_list[0]
        
        def sfmt(val, fmt="{:.1f}"):
            try:
                return fmt.format(val) if not np.isnan(val) else "nan"
            except:
                return "nan"
        
        detail_text = f"""
        Detailed Metrics (Season 1):
        • Rate increase: {sfmt(ph.get('Rate_increase'), '{:.4f}')}    • Large integral: {sfmt(ph.get('Large_integral'), '{:.2f}')}
        • Rate decrease: {sfmt(ph.get('Rate_decrease'), '{:.4f}')}    • Small integral: {sfmt(ph.get('Small_integral'), '{:.2f}')}
        • Base value:    {sfmt(ph.get('Base_value'), '{:.3f}')}    • Harvest DOY:   {sfmt(ph.get('Harvest_DOY'))}
        """
        ax3.text(0.05, 0.5, detail_text, transform=ax3.transAxes,
                 fontsize=9, verticalalignment='center', fontfamily='monospace',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    ✓ Saved multi-season plot")