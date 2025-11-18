"""
TIMESAT Data Preprocessing Module
==================================
Preprocessing functions for:
1. Removing duplicate files (ignoring '_T' suffix)
2. Rasterizing vector labels to match image patches

This module handles data preparation before phenology analysis.
"""

import os
import re
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_bounds
import numpy as np
import config

# ============================================================================
# DUPLICATE FILE REMOVAL
# ============================================================================

def remove_duplicate_files_ignore_T(folder_path):
    """
    Remove duplicate files in a folder by ignoring '_T' suffix in filenames

    This function:
    - Walks through all files in the folder
    - Extracts base name by removing everything after '_T'
    - Keeps only the first occurrence of each base name
    - Removes all duplicate files

    Parameters:
    -----------
    folder_path : str
        Path to folder containing files to check

    Returns:
    --------
    None
        Prints status messages for each operation

    Example:
    --------
    Files before: S2_image_tile1_T001.tif, S2_image_tile1_T002.tif
    Files after: S2_image_tile1_T001.tif (keeps first occurrence)
    """
    # Dictionary to track first occurrence of each base filename
    seen = {}

    # List to store paths of removed files
    removed_files = []

    # Walk through all files in the folder and subdirectories
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            # Get full file path
            file_path = os.path.join(root, file)

            # Remove file extension from filename
            name_without_ext = os.path.splitext(file)[0]

            # Extract base name by removing everything after '_T'
            # Example: "S2_image_tile1_T001" -> "S2_image_tile1"
            base_name = re.split(r'_T', name_without_ext)[0]

            # Check if this base name has been seen before
            if base_name in seen:
                # This is a duplicate - remove it
                try:
                    os.remove(file_path)
                    removed_files.append(file_path)
                    print(f"Removed: {file_path}")
                except Exception as e:
                    print(f"Error removing {file_path}: {e}")
            else:
                # First occurrence - keep it
                seen[base_name] = file_path

    # Print summary
    print("\n✅ Operation completed.")
    if removed_files:
        print("Removed files:")
        for f in removed_files:
            print(" -", f)
    else:
        print("No duplicate files found.")


# ============================================================================
# LABEL RASTERIZATION
# ============================================================================

def create_patch_labels(polygons_path, patches_path, output_dir, default_width=512, default_height=512,
                        all_touched=True, nodata=0, resolution=None):
    """
    Rasterize vector labels from polygons GeoJSON for each patch in the patches GeoJSON.

    Each patch should ideally have 'width' and 'height' properties for pixel dimensions.
    If not provided, default_width and default_height are used.
    If resolution is specified, width/height will be calculated based on patch bounds.

    Skips patches that already have existing label files.
    Empty areas (no overlapping polygons) will be filled with 0.
    Output is strictly binary: 0 or 1 only.

    Returns:
    --------
    int -- Number of new label files created
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load polygons GeoJSON
    try:
        gdf_poly = gpd.read_file(polygons_path)
        if gdf_poly.empty:
            print(f"⚠️ Empty polygons GeoJSON: '{polygons_path}'. Creating all-zero labels if applicable.")
    except Exception as e:
        print(f"❌ Failed to read polygons GeoJSON: {e}")
        return 0

    # Load patches GeoJSON
    try:
        gdf_patches = gpd.read_file(patches_path)
        if gdf_patches.empty:
            print(f"⚠️ Empty patches GeoJSON: '{patches_path}'.")
            return 0
    except Exception as e:
        print(f"❌ Failed to read patches GeoJSON: {e}")
        return 0

    count_written = 0

    # Process each patch feature
    for idx, row in gdf_patches.iterrows():
        patch_geom = row.geometry
        if patch_geom is None or patch_geom.is_empty:
            print(f"⚠️ Empty geometry for patch index {idx}.")
            continue

        properties = row.to_dict()  # Ensure properties are in dictionary format

        patch_id = properties.get('id', str(idx))
        out_name = f"{patch_id}_label.tif"
        out_path = os.path.join(output_dir, out_name)

        if os.path.exists(out_path):
            print(f"⏭️ Exists, skip: {out_path}")
            continue

        crs = gdf_patches.crs

        # Get spatial bounds
        minx, miny, maxx, maxy = patch_geom.bounds

        # Determine output raster shape
        if resolution is not None:
            width = int((maxx - minx) / resolution) + 1
            height = int((maxy - miny) / resolution) + 1
        else:
            width = properties.get('width', default_width)
            height = properties.get('height', default_height)

        transform = from_bounds(minx, miny, maxx, maxy, width, height)

        # Rasterize patch geometry as mask
        mask_shapes = [(patch_geom, 1)]
        try:
            mask_arr = rasterize(
                shapes=mask_shapes,
                out_shape=(height, width),
                transform=transform,
                fill=0,
                dtype="uint8",
                all_touched=all_touched
            )
            mask_arr = (mask_arr > 0).astype("uint8")  # Ensure strictly binary
        except Exception as e:
            print(f"❌ Mask rasterize failed for patch {patch_id}: {e}")
            continue

        # If polygons GeoJSON is empty
        if 'gdf_poly' not in locals() or gdf_poly.empty:
            arr = np.zeros((height, width), dtype="uint8")
        else:
            gdf_poly_reproj = gdf_poly.to_crs(crs) if gdf_poly.crs != crs else gdf_poly

            clipped_gdf = gpd.clip(gdf_poly_reproj, patch_geom)
            shapes = ((geom, 1) for geom in clipped_gdf.geometry if geom is not None and not geom.is_empty)

            try:
                arr = rasterize(
                    shapes=shapes,
                    out_shape=(height, width),
                    transform=transform,
                    fill=0,
                    dtype="uint8",
                    all_touched=all_touched
                )
                arr = (arr > 0).astype("uint8")  # Strictly binary
            except Exception as e:
                print(f"❌ Rasterize failed for patch {patch_id}: {e}")
                continue

        # Apply mask
        arr = np.multiply(arr, mask_arr)

        arr = (arr > 0).astype("uint8")  # Final binary enforcement

        # Save profile
        profile = {
            'driver': 'GTiff',
            'height': height,
            'width': width,
            'count': 1,
            'dtype': 'uint8',
            'crs': crs,
            'transform': transform,
            'nodata': 0,
            'compress': 'lzw',
            'NBITS': 1
        }

        try:
            with rasterio.open(out_path, "w", **profile) as dst:
                dst.write(arr, 1)
            count_written += 1
            print(f"✅ Saved label: {out_path}")
        except Exception as e:
            print(f"❌ Write failed for {out_path}: {e}")

    return count_written


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def preprocess_pipeline(images_folder, polygons_path, patches_path, labels_output_dir):
    """
    Complete preprocessing pipeline:
    1. Remove duplicate images
    2. Rasterize labels

    Returns summary as dictionary
    """
    print("=" * 70)
    print("STARTING PREPROCESSING PIPELINE")
    print("=" * 70)

    print("\n" + "─" * 70)
    print("STEP 1: Removing duplicate files")
    print("─" * 70)

    if os.path.isdir(images_folder):
        remove_duplicate_files_ignore_T(images_folder)
    else:
        print(f"❌ Invalid images folder: {images_folder}")
        return {"status": "failed", "reason": "Invalid images folder"}

    print("\n" + "─" * 70)
    print("STEP 2: Rasterizing vector labels")
    print("─" * 70)

    num_created = create_patch_labels(polygons_path, patches_path, labels_output_dir)

    print("\n" + "=" * 70)
    print("PREPROCESSING COMPLETED")
    print("=" * 70)
    print(f"Total new labels created: {num_created}")

    return {
        "status": "success",
        "labels_created": num_created,
        "images_folder": images_folder,
        "labels_folder": labels_output_dir
    }

# ============================================================================
# UPDATED EXAMPLE USAGE (USE config.py PATHS)
# ============================================================================

if __name__ == "__main__":
    try:
        config.validate_config()
    except Exception as e:
        print(f"⚠️ Config validation failed: {e}")

    prep_cfg = config.get_preprocessing_config()
    images_folder = prep_cfg.get('images_folder')
    polygons_path = prep_cfg.get('polygons_path')
    patches_path = prep_cfg.get('patches_path')
    labels_output_dir = prep_cfg.get('labels_output')

    print("\nUsing preprocessing paths from config.py:")
    print(f" - images_folder    : {images_folder}")
    print(f" - polygons_path    : {polygons_path}")
    print(f" - patches_path     : {patches_path}")
    print(f" - labels_output_dir: {labels_output_dir}\n")

    result = preprocess_pipeline(
        images_folder=images_folder,
        polygons_path=polygons_path,
        patches_path=patches_path,
        labels_output_dir=labels_output_dir
    )

    print("\nPreprocessing result:")
    print(result)
