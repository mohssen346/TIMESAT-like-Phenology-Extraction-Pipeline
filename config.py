"""
TIMESAT Configuration Module
============================
Configuration parameters for phenology extraction pipeline

Complete configuration file with all settings for:
- File paths
- Band indices
- Smoothing parameters
- Season detection
- Function fitting
- Phenology extraction
- Output options
- Preprocessing settings
"""


# ============================================================================
# FILE PATHS CONFIGURATION
# ============================================================================

# Input folder containing multi-band satellite images (e.g., Sentinel-2)
# Expected format: S2_image_{tile_id}_{date}.tif
INPUT_FOLDER = r"D:\Project\Cotton\Timesat\Data\S2_individual_images\2018"

# Folder containing rasterized label files (one per tile)
# Expected format: {tile_id}_label.tif
LABEL_FOLDER = r"D:\Project\Cotton\Timesat\Data\Labels"

# Output folder for phenology results, time series, and plots
OUTPUT_FOLDER = r"D:\Project\Cotton\Timesat\Data\Output_Complete"


# ============================================================================
# PREPROCESSING PATHS (for preprocess.py)
# ============================================================================

# Folder containing raw satellite images (before duplicate removal)
PREPROCESS_IMAGES_FOLDER = r"D:\Project\Cotton\Timesat\Data\S2_individual_images\2018"

# Path to polygons GeoJSON file (labeled areas - ground truth)
PREPROCESS_POLYGONS_PATH = r"D:\Project\Cotton\Timesat\Data\Final_Golestan.geojson"

# Path to patches/tiles GeoJSON file (defines image tiles/patches)
PREPROCESS_PATCHES_PATH = r"D:\Project\Cotton\Timesat\Data\golestan_grid_2016_Edit.geojson"

# Output directory for rasterized labels (will be used as LABEL_FOLDER)
PREPROCESS_LABELS_OUTPUT = r"D:\Project\Cotton\Timesat\Data\Labels"


# ============================================================================
# BAND INDICES CONFIGURATION
# ============================================================================

# Band indices for multi-band images (1-based indexing)
# Default values are for Sentinel-2 L2A products
RED_BAND_INDEX = 3      # Red band (B4 in S2)
NIR_BAND_INDEX = 7      # Near-infrared band (B8 in S2)
BLUE_BAND_INDEX = 1     # Blue band (B2 in S2)
GREEN_BAND_INDEX = 2    # Green band (B3 in S2)


# ============================================================================
# SMOOTHING CONFIGURATION
# ============================================================================

# Smoothing method for time series
# Options: 'double_savgol' (double Savitzky-Golay filter - recommended)
SMOOTHING_METHOD = 'double_savgol'

# Savitzky-Golay filter parameters
WINDOW_LENGTH = 7       # Window length for smoothing (must be odd, >= 3)
                        # Larger values = more smoothing
                        # Recommended: 5-11 for 16-day intervals

POLYORDER = 2           # Polynomial order for fitting (must be < WINDOW_LENGTH)
                        # Recommended: 2-3 for vegetation time series

# Outlier removal using Modified Z-score (MAD-based)
Z_THRESHOLD = 2.5       # Z-score threshold for outlier detection
                        # Lower = more aggressive outlier removal
                        # Recommended: 2.5-3.5 for NDVI data


# ============================================================================
# MULTIPLE SEASONS DETECTION CONFIGURATION
# ============================================================================

# Enable detection of multiple growing seasons per year
# Set to False for single-season crops (e.g., cotton, wheat)
# Set to True for multi-cropping systems (e.g., rice-wheat rotation)
DETECT_MULTIPLE_SEASONS = True

# Peak detection parameters (for scipy.signal.find_peaks)
MIN_PEAK_PROMINENCE = 0.1   # Minimum prominence for peak detection
                            # Higher = more selective (fewer seasons detected)
                            # Recommended: 0.05-0.15 for NDVI

MIN_PEAK_DISTANCE = 15      # Minimum distance between peaks (in time steps)
                            # For 16-day intervals: 15 steps ≈ 240 days
                            # Prevents detection of noise as separate seasons

MIN_AMPLITUDE = 0.05        # Minimum amplitude for valid season
                            # Filters out very small variations
                            # Recommended: 0.05-0.10 for NDVI


# ============================================================================
# FUNCTION FITTING CONFIGURATION
# ============================================================================

# Enable function fitting to smooth phenological curves
# Recommended: True for robust phenology metrics
# Set to False to use only Savitzky-Golay smoothed data
USE_FUNCTION_FITTING = True

# Fitting method for phenological curves
# Options:
#   - 'asymmetric_gaussian': Asymmetric Gaussian (TIMESAT default)
#     Good for most crops, handles asymmetric green-up/senescence
#   - 'double_logistic': Double logistic function
#     Good for crops with sharp transitions
FITTING_METHOD = 'asymmetric_gaussian'

# Maximum iterations for optimization
MAX_FITTING_ITERATIONS = 50


# ============================================================================
# PHENOLOGY EXTRACTION CONFIGURATION
# ============================================================================

# Threshold method for determining Start/End of Season (SOS/EOS)
# Options: 'seasonal_amplitude' (recommended)
THRESHOLD_METHOD = 'seasonal_amplitude'

# Amplitude threshold for SOS/EOS detection (0-1)
# SOS/EOS defined as points where NDVI = base + threshold * amplitude
# Lower = earlier SOS, later EOS (more conservative)
# Higher = later SOS, earlier EOS (more strict)
# Recommended: 0.2-0.3 for most crops
AMPLITUDE_THRESHOLD = 0.2

# Harvest detection ratio (0-1)
# Harvest point = when NDVI drops to (peak_value * harvest_ratio)
# Typical: 0.90-0.95 (near peak value)
HARVEST_RATIO = 0.95


# ============================================================================
# REGION PROCESSING CONFIGURATION
# ============================================================================

# Minimum number of pixels for a region to be considered valid
# Smaller regions are filtered out to avoid noise
# Recommended: 5-10 pixels (depends on image resolution)
MIN_REGION_PIXELS = 5


# ============================================================================
# OUTPUT CONFIGURATION
# ============================================================================

# Save raw and smoothed time series as CSV files
# Set to False to save disk space (only phenology metrics will be saved)
SAVE_TIMESERIES_CSV = True

# Save phenology plots as PNG files
# Plots show raw data, smoothed curve, and phenology markers
# Set to False to save processing time
SAVE_PLOTS = True

# Plot configuration
PLOT_DPI = 300              # Resolution for saved plots (dots per inch)
PLOT_FIGSIZE = (16, 10)     # Figure size in inches (width, height)


# ============================================================================
# PARALLEL PROCESSING CONFIGURATION
# ============================================================================

# Enable parallel processing for regions within each tile
# Significantly speeds up processing for tiles with many regions
# Recommended: True (unless debugging)
USE_PARALLEL_PROCESSING = True

# Number of parallel processes
# None = auto-detect CPU count (recommended)
# Integer = specific number of processes (useful for limiting CPU usage)
N_PROCESSES = None


# ============================================================================
# VALIDATION PARAMETERS
# ============================================================================

# Maximum valid label value (filter out anomalous labels)
# Labels with absolute value > MAX_LABEL_VALUE are set to 0
# Prevents processing of corrupted or invalid label values
MAX_LABEL_VALUE = 10000

# Valid NDVI range (values outside this range are clipped)
NDVI_MIN = -1.0     # Theoretical minimum NDVI
NDVI_MAX = 1.0      # Theoretical maximum NDVI


# ============================================================================
# ADVANCED SETTINGS (rarely need modification)
# ============================================================================

# Interpolation limit for gap-filling in time series
# Maximum number of consecutive NaN values to interpolate
INTERPOLATION_LIMIT = 3

# Edge handling for Savitzky-Golay filter
# Options: 'interp', 'nearest', 'mirror', 'wrap'
SAVGOL_MODE = 'interp'


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_config_summary():
    """
    Return a formatted summary of current configuration

    Returns:
    --------
    str
        Formatted configuration summary string

    Example:
    --------
    >>> print(get_config_summary())
    """
    summary = f"""
    ╔════════════════════════════════════════════════════════════════════════╗
    ║                  TIMESAT CONFIGURATION SUMMARY                         ║
    ╠════════════════════════════════════════════════════════════════════════╣
    ║ INPUT PATHS                                                            ║
    ╟────────────────────────────────────────────────────────────────────────╢
    ║ Input Folder:      {INPUT_FOLDER:<50s} ║
    ║ Label Folder:      {LABEL_FOLDER:<50s} ║
    ║ Output Folder:     {OUTPUT_FOLDER:<50s} ║
    ╠════════════════════════════════════════════════════════════════════════╣
    ║ SMOOTHING PARAMETERS                                                   ║
    ╟────────────────────────────────────────────────────────────────────────╢
    ║ Method:            {SMOOTHING_METHOD:<50s} ║
    ║ Window Length:     {WINDOW_LENGTH:<50d} ║
    ║ Polynomial Order:  {POLYORDER:<50d} ║
    ║ Z-Threshold:       {Z_THRESHOLD:<50.2f} ║
    ╠════════════════════════════════════════════════════════════════════════╣
    ║ PHENOLOGY SETTINGS                                                     ║
    ╟────────────────────────────────────────────────────────────────────────╢
    ║ Multiple Seasons:  {str(DETECT_MULTIPLE_SEASONS):<50s} ║
    ║ Function Fitting:  {str(USE_FUNCTION_FITTING):<50s} ║
    ║ Fitting Method:    {FITTING_METHOD:<50s} ║
    ║ Amplitude Thresh:  {AMPLITUDE_THRESHOLD:<50.2f} ║
    ║ Harvest Ratio:     {HARVEST_RATIO:<50.2f} ║
    ╠════════════════════════════════════════════════════════════════════════╣
    ║ PROCESSING OPTIONS                                                     ║
    ╟────────────────────────────────────────────────────────────────────────╢
    ║ Min Peak Distance: {MIN_PEAK_DISTANCE:<50d} ║
    ║ Min Region Pixels: {MIN_REGION_PIXELS:<50d} ║
    ║ Save CSV:          {str(SAVE_TIMESERIES_CSV):<50s} ║
    ║ Save Plots:        {str(SAVE_PLOTS):<50s} ║
    ║ Parallel Process:  {str(USE_PARALLEL_PROCESSING):<50s} ║
    ╚════════════════════════════════════════════════════════════════════════╝
    """
    return summary


def validate_config():
    """
    Validate configuration parameters

    Checks all configuration values for validity and raises
    ValueError if any parameters are invalid.

    Validations:
    - Window length must be odd and >= 3
    - Polynomial order must be < window length
    - Thresholds must be in valid ranges (0-1)
    - Fitting method must be valid
    - Band indices must be positive

    Raises:
    -------
    ValueError
        If any configuration parameter is invalid

    Example:
    --------
    >>> validate_config()  # Prints "✓ Configuration validation passed"
    """
    # ─────────────────────────────────────────────────────────────
    # Validate smoothing parameters
    # ─────────────────────────────────────────────────────────────
    if WINDOW_LENGTH < 3 or WINDOW_LENGTH % 2 == 0:
        raise ValueError("WINDOW_LENGTH must be odd and >= 3")

    if POLYORDER >= WINDOW_LENGTH:
        raise ValueError("POLYORDER must be < WINDOW_LENGTH")

    if Z_THRESHOLD <= 0:
        raise ValueError("Z_THRESHOLD must be positive")

    # ─────────────────────────────────────────────────────────────
    # Validate phenology thresholds
    # ─────────────────────────────────────────────────────────────
    if not (0 <= AMPLITUDE_THRESHOLD <= 1):
        raise ValueError("AMPLITUDE_THRESHOLD must be between 0 and 1")

    if not (0 <= HARVEST_RATIO <= 1):
        raise ValueError("HARVEST_RATIO must be between 0 and 1")

    # ─────────────────────────────────────────────────────────────
    # Validate function fitting method
    # ─────────────────────────────────────────────────────────────
    valid_methods = ['asymmetric_gaussian', 'double_logistic']
    if FITTING_METHOD not in valid_methods:
        raise ValueError(f"FITTING_METHOD must be one of {valid_methods}")

    # ─────────────────────────────────────────────────────────────
    # Validate band indices
    # ─────────────────────────────────────────────────────────────
    band_indices = [RED_BAND_INDEX, NIR_BAND_INDEX, BLUE_BAND_INDEX, GREEN_BAND_INDEX]
    if any(idx < 1 for idx in band_indices):
        raise ValueError("All band indices must be >= 1 (1-based indexing)")

    # ─────────────────────────────────────────────────────────────
    # Validate season detection parameters
    # ─────────────────────────────────────────────────────────────
    if MIN_PEAK_PROMINENCE <= 0:
        raise ValueError("MIN_PEAK_PROMINENCE must be positive")

    if MIN_PEAK_DISTANCE < 1:
        raise ValueError("MIN_PEAK_DISTANCE must be >= 1")

    if MIN_AMPLITUDE <= 0:
        raise ValueError("MIN_AMPLITUDE must be positive")

    # ─────────────────────────────────────────────────────────────
    # Validate region parameters
    # ─────────────────────────────────────────────────────────────
    if MIN_REGION_PIXELS < 1:
        raise ValueError("MIN_REGION_PIXELS must be >= 1")

    # ─────────────────────────────────────────────────────────────
    # Validate NDVI range
    # ─────────────────────────────────────────────────────────────
    if NDVI_MIN >= NDVI_MAX:
        raise ValueError("NDVI_MIN must be < NDVI_MAX")

    if not (-1 <= NDVI_MIN <= 1) or not (-1 <= NDVI_MAX <= 1):
        raise ValueError("NDVI range must be within [-1, 1]")

    # All validations passed
    print("✓ Configuration validation passed")


def get_preprocessing_config():
    """
    Get preprocessing configuration as dictionary

    Returns:
    --------
    dict
        Dictionary with preprocessing paths

    Example:
    --------
    >>> config = get_preprocessing_config()
    >>> print(config['images_folder'])
    """
    return {
        'images_folder': PREPROCESS_IMAGES_FOLDER,
        'polygons_path': PREPROCESS_POLYGONS_PATH,
        'patches_path': PREPROCESS_PATCHES_PATH,
        'labels_output': PREPROCESS_LABELS_OUTPUT
    }


def get_processing_config():
    """
    Get main processing configuration as dictionary

    Returns:
    --------
    dict
        Dictionary with all processing parameters

    Example:
    --------
    >>> config = get_processing_config()
    >>> print(config['window_length'])
    7
    """
    return {
        # Paths
        'input_folder': INPUT_FOLDER,
        'label_folder': LABEL_FOLDER,
        'output_folder': OUTPUT_FOLDER,

        # Bands
        'red_band': RED_BAND_INDEX,
        'nir_band': NIR_BAND_INDEX,
        'blue_band': BLUE_BAND_INDEX,
        'green_band': GREEN_BAND_INDEX,

        # Smoothing
        'smoothing_method': SMOOTHING_METHOD,
        'window_length': WINDOW_LENGTH,
        'polyorder': POLYORDER,
        'z_threshold': Z_THRESHOLD,

        # Phenology
        'detect_multiple': DETECT_MULTIPLE_SEASONS,
        'use_fitting': USE_FUNCTION_FITTING,
        'fitting_method': FITTING_METHOD,
        'amplitude_threshold': AMPLITUDE_THRESHOLD,
        'harvest_ratio': HARVEST_RATIO,

        # Season detection
        'min_peak_prominence': MIN_PEAK_PROMINENCE,
        'min_peak_distance': MIN_PEAK_DISTANCE,
        'min_amplitude': MIN_AMPLITUDE,

        # Output
        'save_csv': SAVE_TIMESERIES_CSV,
        'save_plots': SAVE_PLOTS,
        'plot_dpi': PLOT_DPI,

        # Processing
        'parallel': USE_PARALLEL_PROCESSING,
        'n_processes': N_PROCESSES,
        'min_pixels': MIN_REGION_PIXELS
    }