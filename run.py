"""
TIMESAT Main Processing Pipeline
=================================
Main execution script for:
- Multi-tile processing
- Per-region phenology extraction
- Parallel processing support
- Results aggregation and reporting
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from multiprocessing import Pool, cpu_count
import warnings

# Import custom modules
from models import fit_function_to_timeseries
from config import *
from utils import (
    stack_multiband_files, compute_ndvi, compute_savi, compute_msavi2,
    compute_gndvi, compute_rvi, double_savgol_smooth, load_label_for_tile,
    extract_regions_from_label, detect_multiple_seasons, plot_multi_season_phenology
)

warnings.filterwarnings('ignore')


# ============================================================================
# PHENOLOGY EXTRACTION CORE FUNCTIONS
# ============================================================================

def extract_phenology_single_season(ts_smooth, times, season_range=None,
                                    method='seasonal_amplitude',
                                    amplitude_threshold=0.2,
                                    harvest_ratio=0.95,
                                    use_function_fitting=True,
                                    fitting_method='asymmetric_gaussian'):
    """
    Extract phenology metrics for a single season
    
    Parameters:
    -----------
    ts_smooth : array-like
        Smoothed time series
    times : array-like
        Datetime objects
    season_range : tuple, optional
        (start_idx, end_idx) for season boundaries
    method : str
        Threshold method
    amplitude_threshold : float
        Threshold for SOS/EOS detection
    harvest_ratio : float
        Ratio for harvest detection
    use_function_fitting : bool
        Whether to use function fitting
    fitting_method : str
        Fitting method name
    
    Returns:
    --------
    dict
        Phenology metrics dictionary
    """
    # Restrict to season range if provided
    if season_range is not None:
        start_idx, end_idx = season_range
        ts_season = ts_smooth[start_idx:end_idx + 1]
        times_season = times[start_idx:end_idx + 1]
        idx_offset = start_idx
    else:
        ts_season = ts_smooth
        times_season = times
        idx_offset = 0
    
    # Check for valid data
    if np.all(np.isnan(ts_season)) or len(ts_season) < 3:
        return _empty_phenology_dict()
    
    # Optional function fitting
    if use_function_fitting:
        t_indices = np.arange(len(ts_season))
        fitted_values, fit_params, fit_success = fit_function_to_timeseries(
            t_indices, ts_season, method=fitting_method
        )
        if fit_success:
            print(f"      ✓ Function fitting successful ({fitting_method})")
            ts_season = fitted_values
        else:
            print(f"      ⚠ Function fitting failed, using smoothed data")
    
    # Basic statistics
    minv, maxv = np.nanmin(ts_season), np.nanmax(ts_season)
    amp = maxv - minv
    base = minv
    
    # Handle zero amplitude
    if amp == 0 or np.isnan(amp):
        peak_idx = int(np.nanargmax(ts_season))
        result = _simple_phenology_dict(ts_season, times_season, peak_idx)
        result['fit_method'] = fitting_method if use_function_fitting else 'none'
        return result
    
    # Calculate threshold
    if method == 'seasonal_amplitude':
        threshold_val = base + amplitude_threshold * amp
    else:
        threshold_val = base + amplitude_threshold * amp
    
    # Find peak
    peak_idx = int(np.nanargmax(ts_season))
    peak_val = ts_season[peak_idx]
    
    # Find Start of Season (SOS) with interpolation
    sos_idx = _find_crossing(ts_season, threshold_val, 0, peak_idx, ascending=True)
    
    # Find End of Season (EOS)
    eos_idx = _find_crossing(ts_season, threshold_val, peak_idx, 
                             len(ts_season) - 1, ascending=False)
    
    # Calculate Length of Season (LOS)
    los = (eos_idx - sos_idx) if (sos_idx is not None and eos_idx is not None) else None
    
    # Calculate rates
    rate_increase = None
    rate_decrease = None
    if sos_idx is not None:
        denom = (peak_idx - sos_idx)
        rate_increase = (peak_val - threshold_val) / denom if denom > 0 else np.nan
    if eos_idx is not None:
        denom = (eos_idx - peak_idx)
        rate_decrease = (peak_val - threshold_val) / denom if denom > 0 else np.nan
    
    # Calculate integrals
    large_integral = None
    small_integral = None
    if sos_idx is not None and eos_idx is not None:
        start_i = int(np.floor(sos_idx))
        end_i = int(np.ceil(eos_idx))
        if end_i > start_i:
            large_integral = np.trapz(ts_season[start_i:end_i + 1] - base)
            small_integral = np.trapz(
                np.maximum(ts_season[start_i:end_i + 1] - threshold_val, 0)
            )
    
    # Find harvest point
    harvest_idx = None
    threshold_harvest = peak_val * harvest_ratio
    
    # Search after peak
    for i in range(peak_idx, len(ts_season)):
        if np.isnan(ts_season[i]):
            continue
        if ts_season[i] <= threshold_harvest:
            harvest_idx = float(i)
            break
    
    # Fallback: search before peak
    if harvest_idx is None:
        for i in range(peak_idx - 1, -1, -1):
            if np.isnan(ts_season[i]):
                continue
            if ts_season[i] <= threshold_harvest:
                harvest_idx = float(i)
                break
    
    # Mid-season points (50% and 80% of amplitude)
    val_50 = base + 0.5 * amp
    val_80 = base + 0.8 * amp
    
    mid_50_start = _find_crossing(ts_season, val_50, 0, peak_idx, ascending=True)
    mid_50_end = _find_crossing(ts_season, val_50, peak_idx, 
                                len(ts_season) - 1, ascending=False)
    mid_80_start = _find_crossing(ts_season, val_80, 0, peak_idx, ascending=True)
    mid_80_end = _find_crossing(ts_season, val_80, peak_idx, 
                                len(ts_season) - 1, ascending=False)
    
    # Convert indices to DOY
    def idx_to_doy(idx):
        """Convert fractional index to Day of Year"""
        if idx is None or np.isnan(idx):
            return np.nan
        actual_idx = idx + idx_offset
        i_floor = int(np.floor(actual_idx))
        i_ceil = int(np.ceil(actual_idx))
        if i_floor >= len(times):
            i_floor = len(times) - 1
        if i_ceil >= len(times):
            i_ceil = len(times) - 1
        frac = actual_idx - i_floor
        date_floor = pd.to_datetime(times[i_floor])
        date_ceil = pd.to_datetime(times[i_ceil])
        interpolated_date = date_floor + (date_ceil - date_floor) * frac
        return float(interpolated_date.dayofyear)
    
    # Build result dictionary
    result = {
        "SOS_DOY": idx_to_doy(sos_idx),
        "EOS_DOY": idx_to_doy(eos_idx),
        "Peak_DOY": idx_to_doy(peak_idx),
        "Peak_value": float(peak_val),
        "Base_value": float(base),
        "Amplitude": float(amp),
        "LOS": float(los) if los is not None else np.nan,
        "Rate_increase": float(rate_increase) if rate_increase is not None else np.nan,
        "Rate_decrease": float(rate_decrease) if rate_decrease is not None else np.nan,
        "Large_integral": float(large_integral) if large_integral is not None else np.nan,
        "Small_integral": float(small_integral) if small_integral is not None else np.nan,
        "Mid_50_start_DOY": idx_to_doy(mid_50_start),
        "Mid_50_end_DOY": idx_to_doy(mid_50_end),
        "Mid_80_start_DOY": idx_to_doy(mid_80_start),
        "Mid_80_end_DOY": idx_to_doy(mid_80_end),
        "Harvest_DOY": idx_to_doy(harvest_idx),
        "fit_method": fitting_method if use_function_fitting else 'none'
    }
    
    # Correct LOS calculation based on actual DOY (handle year wrap)
    sos_doy = result.get("SOS_DOY", np.nan)
    eos_doy = result.get("EOS_DOY", np.nan)
    
    if not np.isnan(sos_doy) and not np.isnan(eos_doy):
        los_corrected = eos_doy - sos_doy
        if los_corrected < 0:
            los_corrected += 365  # Assuming non-leap year
        result["LOS"] = los_corrected
    
    return result


def _find_crossing(arr, threshold, start, end, ascending=True):
    """
    Find first threshold crossing with linear interpolation
    
    Parameters:
    -----------
    arr : array-like
        Values
    threshold : float
        Threshold value
    start : int
        Start index
    end : int
        End index
    ascending : bool
        Search direction
    
    Returns:
    --------
    float or None
        Interpolated crossing index
    """
    if ascending:
        i_range = range(start + 1, end + 1)
    else:
        i_range = range(end, start, -1)
    
    for i in i_range:
        prev_i = i - 1
        
        # Check boundaries
        if prev_i < 0 or i >= len(arr):
            continue
        
        # Skip NaN values
        if np.isnan(arr[prev_i]) or np.isnan(arr[i]):
            continue
        
        # Check for crossing
        if (ascending and arr[prev_i] < threshold <= arr[i]) or \
           (not ascending and arr[i] < threshold <= arr[prev_i]):
            denom = arr[i] - arr[prev_i] if ascending else arr[prev_i] - arr[i]
            if denom <= 0:
                continue  # Skip if not in expected direction
            frac = (threshold - arr[prev_i]) / denom if ascending else \
                   (arr[prev_i] - threshold) / denom
            return prev_i + frac
    
    return None


def _empty_phenology_dict():
    """Return dictionary with NaN values for all phenology metrics"""
    return {k: np.nan for k in [
        "SOS_DOY", "EOS_DOY", "Peak_DOY", "Peak_value", "Base_value", "Amplitude",
        "LOS", "Rate_increase", "Rate_decrease", "Large_integral", "Small_integral",
        "Mid_50_start_DOY", "Mid_50_end_DOY", "Mid_80_start_DOY", "Mid_80_end_DOY",
        "Harvest_DOY", "fit_method"
    ]}


def _simple_phenology_dict(ts, times, peak_idx):
    """Return simple phenology dictionary with only peak information"""
    result = _empty_phenology_dict()
    result["Peak_DOY"] = float(pd.to_datetime(times[peak_idx]).dayofyear)
    result["Peak_value"] = float(ts[peak_idx])
    result["fit_method"] = 'none'
    return result


def extract_timesat_phenology_multi_season(ts_smooth, times,
                                           detect_multiple=True,
                                           use_function_fitting=True,
                                           fitting_method='asymmetric_gaussian',
                                           amplitude_threshold=0.2,
                                           harvest_ratio=0.95,
                                           min_peak_prominence=0.1,
                                           min_peak_distance=15):
    """
    Main function for extracting phenology with multiple season detection
    
    Parameters:
    -----------
    ts_smooth : array-like
        Smoothed time series
    times : array-like
        Datetime objects
    detect_multiple : bool
        Whether to detect multiple seasons
    use_function_fitting : bool
        Whether to use function fitting
    fitting_method : str
        Fitting method name
    amplitude_threshold : float
        Threshold for SOS/EOS
    harvest_ratio : float
        Harvest detection ratio
    min_peak_prominence : float
        Minimum peak prominence
    min_peak_distance : int
        Minimum distance between peaks
    
    Returns:
    --------
    list
        List of phenology dictionaries per season
    """
    results = []
    
    # Single season mode
    if not detect_multiple:
        ph = extract_phenology_single_season(
            ts_smooth, times, season_range=None,
            amplitude_threshold=amplitude_threshold,
            harvest_ratio=harvest_ratio,
            use_function_fitting=use_function_fitting,
            fitting_method=fitting_method
        )
        ph['season_number'] = 1
        results.append(ph)
        return results
    
    # Detect multiple seasons
    seasons = detect_multiple_seasons(
        ts_smooth, times,
        min_peak_prominence=min_peak_prominence,
        min_peak_distance=min_peak_distance
    )
    
    # Fallback to single season if no seasons detected
    if len(seasons) == 0:
        print("      No distinct seasons detected, treating as single season")
        ph = extract_phenology_single_season(
            ts_smooth, times, season_range=None,
            amplitude_threshold=amplitude_threshold,
            harvest_ratio=harvest_ratio,
            use_function_fitting=use_function_fitting,
            fitting_method=fitting_method
        )
        ph['season_number'] = 1
        results.append(ph)
        return results
    
    # Process each detected season
    for season in seasons:
        season_range = (season['start_idx'], season['end_idx'])
        ph = extract_phenology_single_season(
            ts_smooth, times, season_range=season_range,
            amplitude_threshold=amplitude_threshold,
            harvest_ratio=harvest_ratio,
            use_function_fitting=use_function_fitting,
            fitting_method=fitting_method
        )
        ph['season_number'] = season['season_number']
        results.append(ph)
    
    return results


# ============================================================================
# REGION PROCESSING (PARALLEL SUPPORT)
# ============================================================================

def process_region(args):
    """
    Process single region (wrapper for parallel execution)
    
    Parameters:
    -----------
    args : tuple
        Packed arguments for region processing
    
    Returns:
    --------
    list
        Phenology records for the region
    """
    try:
        (tile_id, label_val, region_id, region_info, indices, times, times_pd,
         smoothing_method, window_length, polyorder, detect_multiple_seasons,
         use_function_fitting, fitting_method, amplitude_threshold, harvest_ratio,
         output_folder, save_timeseries_csv, save_plots) = args
        
        region_mask = region_info['mask']
        n_pixels = region_info['n_pixels']
        
        print(f"\n  Region {region_id}: {n_pixels} pixels")
        
        # Extract time series for all indices
        df_ts = pd.DataFrame(index=times_pd)
        for name, da in indices.items():
            mask_da = xr.DataArray(
                region_mask,
                coords={"y": da.coords["y"], "x": da.coords["x"]},
                dims=("y", "x")
            )
            masked = da.where(mask_da)
            series = masked.mean(dim=("y", "x"), skipna=True).values
            df_ts[name] = series
        
        # Save raw time series
        if save_timeseries_csv:
            out_csv = os.path.join(
                output_folder,
                f"timeseries_raw_tile{tile_id}_label{label_val}_region{region_id}.csv"
            )
            df_ts.to_csv(out_csv, index_label="time")
            print(f"    ✓ Raw timeseries saved")
        
        # Smooth NDVI
        ndvi_series = df_ts['NDVI'].values
        if smoothing_method == 'double_savgol':
            ndvi_smooth = double_savgol_smooth(ndvi_series, window_length, polyorder)
        else:
            ndvi_smooth = double_savgol_smooth(ndvi_series, window_length, polyorder)
        
        # Extract phenology
        print(f"    Extracting phenology...")
        phenology_list = extract_timesat_phenology_multi_season(
            ndvi_smooth, times,
            detect_multiple=detect_multiple_seasons,
            use_function_fitting=use_function_fitting,
            fitting_method=fitting_method,
            amplitude_threshold=amplitude_threshold,
            harvest_ratio=harvest_ratio,
            min_peak_prominence=MIN_PEAK_PROMINENCE,
            min_peak_distance=MIN_PEAK_DISTANCE
        )
        
        print(f"    ✓ Found {len(phenology_list)} season(s)")
        
        # Save smoothed time series
        if save_timeseries_csv:
            df_smooth = df_ts.copy()
            df_smooth['NDVI_smoothed'] = ndvi_smooth
            out_smooth = os.path.join(
                output_folder,
                f"timeseries_smoothed_tile{tile_id}_label{label_val}_region{region_id}.csv"
            )
            df_smooth.to_csv(out_smooth, index_label="time")
            print(f"    ✓ Smoothed timeseries saved")
        
        # Generate plot
        if save_plots:
            try:
                out_png = os.path.join(
                    output_folder,
                    f"phenology_tile{tile_id}_label{label_val}_region{region_id}.png"
                )
                plot_multi_season_phenology(
                    df_ts, times, phenology_list,
                    ndvi_smooth, label_val, region_id,
                    tile_id, out_png
                )
            except Exception as e:
                print(f"    ✗ Plot error: {e}")
        
        # Build result rows
        rows = []
        for ph in phenology_list:
            row = {
                "tile": tile_id,
                "label": label_val,
                "region": region_id,
                "n_pixels": n_pixels
            }
            row.update(ph)
            rows.append(row)
        
        return rows
    
    except Exception as e:
        # Log and return empty list (prevents pool crash)
        print(f"    ✗ Error in process_region: {e}")
        return []


# ============================================================================
# TILE PROCESSING
# ============================================================================

def process_tile_per_region_advanced(tile_id, input_files_for_tile, label_folder, output_folder,
                                     red_band_index=3, nir_band_index=7,
                                     blue_band_index=1, green_band_index=2,
                                     smoothing_method='double_savgol',
                                     window_length=9, polyorder=2,
                                     detect_multiple_seasons=True,
                                     use_function_fitting=True,
                                     fitting_method='asymmetric_gaussian',
                                     min_peak_prominence=0.1,
                                     min_peak_distance=15,
                                     amplitude_threshold=0.2,
                                     harvest_ratio=0.95,
                                     min_region_pixels=5,
                                     save_timeseries_csv=True,
                                     save_plots=True):
    """
    Advanced per-region processing for single tile
    
    Parameters:
    -----------
    tile_id : str
        Tile identifier
    input_files_for_tile : list
        List of input files for this tile
    label_folder : str
        Path to label folder
    output_folder : str
        Path to output folder
    red_band_index : int
        Red band index
    nir_band_index : int
        NIR band index
    blue_band_index : int
        Blue band index
    green_band_index : int
        Green band index
    smoothing_method : str
        Smoothing method name
    window_length : int
        Window length for smoothing
    polyorder : int
        Polynomial order for smoothing
    detect_multiple_seasons : bool
        Enable multiple season detection
    use_function_fitting : bool
        Enable function fitting
    fitting_method : str
        Fitting method name
    min_peak_prominence : float
        Minimum peak prominence
    min_peak_distance : int
        Minimum peak distance
    amplitude_threshold : float
        Amplitude threshold for phenology
    harvest_ratio : float
        Harvest detection ratio
    min_region_pixels : int
        Minimum pixels per region
    save_timeseries_csv : bool
        Save time series CSV files
    save_plots : bool
        Save phenology plots
    """
    if len(input_files_for_tile) == 0:
        print(f"No images for tile {tile_id}")
        return
    
    print(f"\n{'=' * 70}")
    print(f"Processing tile {tile_id} with {len(input_files_for_tile)} images")
    print(f"  Multi-season: {detect_multiple_seasons}")
    print(f"  Function fitting: {use_function_fitting} ({fitting_method})")
    print(f"{'=' * 70}")
    
    # Stack bands
    print("Stacking bands...")
    try:
        red_da = stack_multiband_files(input_files_for_tile, red_band_index)
        nir_da = stack_multiband_files(input_files_for_tile, nir_band_index)
        blue_da = stack_multiband_files(input_files_for_tile, blue_band_index)
        green_da = stack_multiband_files(input_files_for_tile, green_band_index)
    except Exception as e:
        print(f"Error stacking bands: {e}")
        return
    
    times = red_da['time'].values
    times_pd = pd.to_datetime(times)
    
    # Load label
    print("Loading label...")
    try:
        label_da = load_label_for_tile(label_folder, tile_id, red_da.isel(time=0))
    except Exception as e:
        print(f"Error loading label: {e}")
        return
    
    # Compute vegetation indices
    print("Computing indices...")
    indices = {
        'NDVI': compute_ndvi(nir_da, red_da),
        'SAVI': compute_savi(nir_da, red_da),
        'MSAVI2': compute_msavi2(nir_da, red_da),
        'GNDVI': compute_gndvi(nir_da, green_da),
        'RVI': compute_rvi(nir_da, red_da)
    }
    
    # Extract unique labels
    label_array = label_da.values
    if np.ma.is_masked(label_array):
        label_array = label_array.filled(0)
    
    unique_labels = [
        int(x) for x in np.unique(label_array)
        if not np.isnan(x) and int(x) != 0 and abs(int(x)) < 10000
    ]
    
    if len(unique_labels) == 0:
        print("No valid labels found!")
        return
    
    print(f"Found {len(unique_labels)} labels: {unique_labels}")
    
    phenology_rows = []
    total_regions = 0
    
    # Process each label
    for label_val in unique_labels:
        print(f"\n{'─' * 60}")
        print(f"Label {label_val}")
        print(f"{'─' * 60}")
        
        # Extract regions for this label
        regions = extract_regions_from_label(label_array, label_val, min_region_pixels)
        
        if len(regions) == 0:
            continue
        
        print(f"  Found {len(regions)} region(s)")
        
        # Prepare arguments for parallel processing
        pool_args = []
        for region_id, region_info in regions.items():
            pool_args.append((
                tile_id, label_val, region_id, region_info, indices, times, times_pd,
                smoothing_method, window_length, polyorder,
                detect_multiple_seasons, use_function_fitting, fitting_method,
                amplitude_threshold, harvest_ratio,
                output_folder, save_timeseries_csv, save_plots
            ))
        
        # Process regions in parallel
        with Pool(processes=min(cpu_count(), len(regions))) as pool:
            results = pool.map(process_region, pool_args)
        
        # Aggregate results
        for rows in results:
            phenology_rows.extend(rows)
            total_regions += 1
    
    # Save results
    if len(phenology_rows) > 0:
        df_ph = pd.DataFrame(phenology_rows)
        df_ph = df_ph.sort_values(['label', 'region', 'season_number'])
        
        out_csv = os.path.join(output_folder, f"phenology_per_region_tile{tile_id}.csv")
        df_ph.to_csv(out_csv, index=False)
        
        print(f"\n{'=' * 70}")
        print(f"✓ Saved {len(phenology_rows)} season records from {total_regions} regions")
        print(f"  File: {out_csv}")
        print(f"{'=' * 70}")
        
        print("\nSummary:")
        print(df_ph[['label', 'region', 'season_number', 'SOS_DOY',
                     'Peak_DOY', 'EOS_DOY', 'Amplitude']].to_string(index=False))
    
    print(f"\n{'=' * 70}")
    print(f"Completed tile {tile_id}")
    print(f"{'=' * 70}\n")


# ============================================================================
# MAIN EXECUTION FUNCTION
# ============================================================================

def main():
    """
    Main execution function for TIMESAT processing pipeline
    """
    # Validate configuration
    print("\n" + "═" * 70)
    print("TIMESAT-LIKE PHENOLOGY EXTRACTION PIPELINE")
    print("═" * 70)
    print(get_config_summary())
    validate_config()
    
    # Create output folder
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    # Find input files
    input_files = [
        os.path.join(INPUT_FOLDER, f)
        for f in os.listdir(INPUT_FOLDER)
        if f.lower().endswith(".tif") and f.startswith("S2_image")
    ]
    
    if len(input_files) == 0:
        raise FileNotFoundError(f"No files found in {INPUT_FOLDER}")
    
    print(f"\nFound {len(input_files)} input files")
    
    # Extract tile IDs
    tile_ids = set()
    for f in input_files:
        basename = os.path.basename(f)
        parts = basename.split('_')
        if len(parts) >= 3:
            tile_ids.add(parts[2])
    
    tile_ids = sorted(list(tile_ids))
    print(f"Found {len(tile_ids)} tiles: {tile_ids}\n")
    
    # Process each tile
    all_results = []
    
    for tile_id in tile_ids:
        # Get files for this tile
        files_for_tile = [
            f for f in input_files
            if f"_image_{tile_id}_" in os.path.basename(f)
        ]
        files_for_tile.sort()
        
        print(f"\n{'=' * 70}")
        print(f"Tile {tile_id}: {len(files_for_tile)} files")
        print(f"{'=' * 70}")
        
        # Process tile
        process_tile_per_region_advanced(
            tile_id,
            files_for_tile,
            LABEL_FOLDER,
            OUTPUT_FOLDER,
            red_band_index=RED_BAND_INDEX,
            nir_band_index=NIR_BAND_INDEX,
            blue_band_index=BLUE_BAND_INDEX,
            green_band_index=GREEN_BAND_INDEX,
            smoothing_method=SMOOTHING_METHOD,
            window_length=WINDOW_LENGTH,
            polyorder=POLYORDER,
            detect_multiple_seasons=DETECT_MULTIPLE_SEASONS,
            use_function_fitting=USE_FUNCTION_FITTING,
            fitting_method=FITTING_METHOD,
            min_peak_prominence=MIN_PEAK_PROMINENCE,
            min_peak_distance=MIN_PEAK_DISTANCE,
            amplitude_threshold=AMPLITUDE_THRESHOLD,
            harvest_ratio=HARVEST_RATIO,
            min_region_pixels=MIN_REGION_PIXELS,
            save_timeseries_csv=SAVE_TIMESERIES_CSV,
            save_plots=SAVE_PLOTS
        )
        
        # Load tile results
        result_file = os.path.join(OUTPUT_FOLDER, f"phenology_per_region_tile{tile_id}.csv")
        if os.path.exists(result_file):
            df = pd.read_csv(result_file)
            all_results.append(df)
    
    # Combine all results
    if len(all_results) > 0:
        df_all = pd.concat(all_results, ignore_index=True)
        out_all = os.path.join(OUTPUT_FOLDER, "phenology_ALL_TILES_COMPLETE.csv")
        df_all.to_csv(out_all, index=False)
        
        print("\n" + "═" * 70)
        print("FINAL SUMMARY")
        print("═" * 70)
        print(f"Total records: {len(df_all)}")
        print(f"Total tiles: {df_all['tile'].nunique()}")
        print(f"Total labels: {df_all['label'].nunique()}")
        print(f"Total regions: {df_all.groupby(['tile', 'label', 'region']).ngroups}")
        print(f"Total seasons: {len(df_all)}")
        
        # Season breakdown
        print("\nSeasons breakdown:")
        season_counts = df_all['season_number'].value_counts().sort_index()
        for season_num, count in season_counts.items():
            print(f"  Season {season_num}: {count} records")
        
        # Fitting method statistics
        if 'fit_method' in df_all.columns:
            print("\nFunction fitting methods:")
            fit_counts = df_all['fit_method'].value_counts()
            for method, count in fit_counts.items():
                print(f"  {method}: {count} ({100 * count / len(df_all):.1f}%)")
        
        print(f"\n✓ Combined results: {out_all}")
        
        # Overall statistics
        print("\nOverall Statistics:")
        print("─" * 70)
        for col in ['SOS_DOY', 'Peak_DOY', 'EOS_DOY', 'Amplitude', 'LOS']:
            if col in df_all.columns:
                try:
                    print(f"{col:15s}: {df_all[col].mean():7.2f} ± {df_all[col].std():6.2f}  "
                          f"(range: {df_all[col].min():.2f} - {df_all[col].max():.2f})")
                except:
                    pass
        
        print("\n" + "═" * 70)
        print("ALL PROCESSING COMPLETED!")
        print("═" * 70)
        print(f"\nOutput folder: {OUTPUT_FOLDER}")
        print("\nGenerated files:")
        print("  • phenology_per_region_tile*.csv - per tile results")
        print("  • phenology_ALL_TILES_COMPLETE.csv - combined results")
        print("  • timeseries_*.csv - time series data")
        print("  • phenology_tile*_label*_region*.png - plots")
    else:
        print("\n⚠ No results generated")


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    main()