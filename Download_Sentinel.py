# =========================
# Colab setup
# =========================

!pip -q install earthengine-api

import ee, json, os, re
from google.colab import drive

drive.mount('/content/drive', force_remount=True)
ee.Authenticate()    # run once and authorize
ee.Initialize(project='ee-imageproccing11111')

# =========================
# User settings (configurable)
# =========================
START_DATE = "2018-01-01"
END_DATE   = "2019-12-31"

DRIVE_FOLDER   = "S2_individual_images"   # Output folder
IMAGE_PREFIX   = "S2_image"               # File name prefix

# Sentinel-2 bands
S2_BANDS = ['B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12']

# Path to the GeoJSON file
GEOJSON_PATH = "/content/drive/MyDrive/Cotton/golestan_grid_2016_Edit.geojson"


# =========================
# Helpers
# =========================
def sanitize(s):
    return re.sub(r'[^A-Za-z0-9_\-]+', '_', str(s))[:50]

def parse_geojson(path):
    with open(path, "r") as f:
        gj = json.load(f)

    feats = []
    t = gj.get("type", "")
    if t == "FeatureCollection":
        source = gj.get("features", [])
    elif t == "Feature":
        source = [gj]
    else:
        source = [{"type": "Feature", "geometry": gj, "properties": {}}]

    for i, feat in enumerate(source):
        geom = feat.get("geometry", None)
        if geom is None:
            continue
        gtype = geom.get("type", "")
        if gtype not in ("Polygon", "MultiPolygon"):
            continue

        props = feat.get("properties", {}) or {}
        cand = None
        for key in ["id","Id","ID","name","Name","title","TITLE","uid","UID"]:
            if key in props and props[key]:
                cand = props[key]
                break
        uid = sanitize(cand if cand is not None else f"feat_{i+1}")

        try:
            ee_geom = ee.Geometry(geom)
        except Exception:
            continue

        feats.append((uid, ee_geom))
    return feats

def s2_collection(start, end, roi):
    """Sentinel-2 collection without cloud mask"""
    return (ee.ImageCollection("COPERNICUS/S2")
            .filterDate(start, end)
            .filterBounds(roi)
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 50))
            .select(S2_BANDS))


# =========================
# Export each image individually
# =========================
def export_images_for_feature(uid, geom):
    col = s2_collection(START_DATE, END_DATE, geom)

    try:
        n = int(col.size().getInfo())
    except Exception as e:
        print(f"⚠️ {uid}: Error reading number of images: {e}")
        return

    if n == 0:
        print(f"⚠️ {uid}: No images found in the date range.")
        return

    col_list = col.toList(n)

    for i in range(n):
        try:
            img = ee.Image(col_list.get(i))
            img_id = img.get('system:index').getInfo()
            img_date = ee.Date(img.get('system:time_start')).format('yyyyMMdd').getInfo()

            fname_prefix = f"{IMAGE_PREFIX}_{uid}_{img_date}_{sanitize(img_id)}"
            desc = fname_prefix

            img_to_export = img.clip(geom)

            task = ee.batch.Export.image.toDrive(
                image=img_to_export,
                description=desc,
                folder=DRIVE_FOLDER,
                fileNamePrefix=fname_prefix,
                region=geom,
                scale=10,
                crs='EPSG:4326',
                maxPixels=1e13
            )
            task.start()
            print(f"🚀 Started export: {fname_prefix}.tif")
        except Exception as e:
            print(f"❌ Error preparing image {i+1}/{n} for {uid}: {e}")


# =========================
# Run
# =========================
features = parse_geojson(GEOJSON_PATH)
if not features:
    raise SystemExit("❌ No valid Polygon/MultiPolygon was found in the GeoJSON.")

print(f"Found {len(features)} features. Starting exports...")
for uid, geom in features:
    export_images_for_feature(uid, geom)
print("All export tasks started. Check the Tasks tab.")
