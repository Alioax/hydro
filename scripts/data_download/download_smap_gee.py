# -*- coding: utf-8 -*-
"""
Download SMAP L3 9 km basin-averaged soil moisture using Google Earth Engine.

Uses NASA/SMAP/SPL3SMP_E (L3 Radiometer Global Daily 9 km). Each day has one
AM pass (6:00 local) and one PM pass (18:00 local). For each pass we compute
basin mean soil moisture and the number of valid pixels (quality). Processed
in one-year windows to avoid timeouts.

Output CSV columns: date, sm, quality, pass
- date: YYYY-MM-DD
- sm: basin mean soil moisture (volume fraction; valid pixels only), or "nan" when no data
- quality: number of 9 km pixels with valid retrieval in the basin (0 when no data)
- pass: AM or PM

Date range: --start-date and --end-date. SMAP L3 from 2015-03-31.
Output: CSV to data/raw_data/SMAP/daily/{station_code}/
"""
import csv
import os
from datetime import datetime
from pathlib import Path
import geopandas as gpd
import ee

# Get project root (parent of scripts/)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "camels_required"
RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw_data"
SMAP_DIR = RAW_DATA_DIR / "SMAP" / "daily"

# Configuration
STATION_CODE = "K287191001"  # Default station
GEE_PROJECT_ID = os.getenv("GEE_PROJECT_ID", "surface-hydrology")
START_DATE = "2015-04-01"   # SMAP from 2015-03-31
END_DATE = "2023-12-31"     # 005 ends 2023-12-03; 006 continues
SMAP_SCALE_M = 9000         # 9 km EASE-Grid
# L3 9 km: 005 until 2023-12-03, 006 from 2023-12-04
SMAP_ASSET_005 = "NASA/SMAP/SPL3SMP_E/005"
SMAP_ASSET_006 = "NASA/SMAP/SPL3SMP_E/006"


def initialize_gee(project_id=None):
    """Initialize Google Earth Engine."""
    try:
        ee.Initialize(project=project_id or GEE_PROJECT_ID)
        print("Google Earth Engine initialized successfully")
        return True
    except Exception as e:
        print(f"Error initializing GEE: {e}")
        print("To authenticate, run: earthengine authenticate")
        return False


def get_catchment_boundary(station_code):
    """Load catchment boundary for a given station."""
    boundaries = gpd.read_file(DATA_DIR / "CAMELS_FR_geography" / "CAMELS_FR_catchment_boundaries.gpkg")
    basin = boundaries[boundaries["sta_code_h3"] == station_code].copy()
    if len(basin) == 0:
        raise ValueError(f"Station {station_code} not found")
    return basin.to_crs(epsg=4326)


def basin_to_ee_geometry(basin_gdf):
    """Convert GeoDataFrame to Earth Engine geometry."""
    geom = basin_gdf.iloc[0].geometry
    geojson_dict = geom.__geo_interface__

    geom_type = geojson_dict["type"]
    coords = geojson_dict["coordinates"]

    if geom_type == "MultiPolygon":
        from shapely.geometry import Polygon as ShapelyPolygon
        shapely_polys = [ShapelyPolygon(poly[0]) for poly in coords]
        areas = [p.area for p in shapely_polys]
        largest_idx = areas.index(max(areas))
        return ee.Geometry.Polygon(coords[largest_idx][0])
    elif geom_type == "Polygon":
        return ee.Geometry.Polygon(coords[0])
    else:
        raise ValueError(f"Unsupported geometry type: {geom_type}")


def _process_day_smap(d, collection, dummy_image, geom, scale_m):
    """
    For one day: get daily image from collection (or dummy if no data), compute AM and PM
    basin mean and pixel count. When no image exists for the day, dummy is used (fully
    masked), so reduceRegion returns null sm and 0 quality. Returns 2 features (AM, PM).
    """
    day_end = ee.Date(d).advance(1, "day")
    filtered = collection.filterDate(d, day_end)
    # Merge with dummy so we never call .mean() on empty collection (avoids no-bands error)
    combined = filtered.merge(ee.ImageCollection([dummy_image]))
    img = ee.Image(combined.mean())
    date_str = ee.Date(d).format("yyyy-MM-dd")

    # AM: use soil_moisture_am as-is (quality flag can be a bitmask; using any pixel with data)
    # quality = number of 9 km pixels with a soil moisture value in the basin
    am_sm = img.select("soil_moisture_am")
    am_clip = am_sm.clip(geom)
    am_mean = am_clip.reduceRegion(
        ee.Reducer.mean(), geometry=geom, scale=scale_m, maxPixels=1e9, bestEffort=True
    ).get("soil_moisture_am")
    am_count = am_clip.reduceRegion(
        ee.Reducer.count(), geometry=geom, scale=scale_m, maxPixels=1e9, bestEffort=True
    ).get("soil_moisture_am")
    feat_am = ee.Feature(None, {"date": date_str, "soil_moisture": am_mean, "quality": am_count, "pass": "AM"})

    # PM: same
    pm_sm = img.select("soil_moisture_pm")
    pm_clip = pm_sm.clip(geom)
    pm_mean = pm_clip.reduceRegion(
        ee.Reducer.mean(), geometry=geom, scale=scale_m, maxPixels=1e9, bestEffort=True
    ).get("soil_moisture_pm")
    pm_count = pm_clip.reduceRegion(
        ee.Reducer.count(), geometry=geom, scale=scale_m, maxPixels=1e9, bestEffort=True
    ).get("soil_moisture_pm")
    feat_pm = ee.Feature(None, {"date": date_str, "soil_moisture": pm_mean, "quality": pm_count, "pass": "PM"})

    return ee.FeatureCollection([feat_am, feat_pm])


def download_smap_basin_csv(station_code, start_date, end_date, output_dir=None):
    """
    Compute basin-averaged SMAP L3 9 km soil moisture per day and pass (AM/PM).
    Output CSV: date, sm, quality, pass. Loops by year to avoid getInfo() timeout.
    """
    if not initialize_gee():
        raise RuntimeError("Failed to initialize Google Earth Engine")

    print(f"Loading catchment boundary for station {station_code}...")
    basin = get_catchment_boundary(station_code)
    bounds = basin.total_bounds
    print(f"Catchment bounds (WGS84): West={bounds[0]:.6f}, South={bounds[1]:.6f}, East={bounds[2]:.6f}, North={bounds[3]:.6f}")
    ee_geom = basin_to_ee_geometry(basin)

    if output_dir is None:
        output_dir = SMAP_DIR / station_code
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"SMAP_basin_9km_{station_code}_{start_date}_{end_date}.csv"

    start_year = int(start_date[:4])
    end_year = int(end_date[:4])

    csv_file = open(output_file, "w", newline="", encoding="utf-8")
    writer = csv.writer(csv_file)
    writer.writerow(["date", "sm", "quality", "pass"])
    total_rows = 0

    for year in range(start_year, end_year + 1):
        y_start = f"{year}-01-01"
        y_end = f"{year}-12-31"
        if year == start_year:
            y_start = start_date
        if year == end_year:
            y_end = end_date
        d0 = datetime.strptime(y_start, "%Y-%m-%d")
        d1 = datetime.strptime(y_end, "%Y-%m-%d")
        num_days = (d1 - d0).days + 1
        if num_days <= 0:
            continue

        print(f"Processing {year} ({num_days} days)...")
        merged = ee.ImageCollection(SMAP_ASSET_005).merge(ee.ImageCollection(SMAP_ASSET_006))
        year_coll = merged.filterDate(y_start, f"{year + 1}-01-01")
        # Dummy image (same bands, fully masked) so .mean() is never called on empty collection
        dummy_image = ee.Image(merged.first()).multiply(0).updateMask(0)
        day_list = (
            ee.List.sequence(0, num_days - 1)
            .map(lambda i: ee.Date(y_start).advance(ee.Number(i), "day"))
        )
        fc = ee.FeatureCollection(day_list.map(
            lambda d: _process_day_smap(d, year_coll, dummy_image, ee_geom, SMAP_SCALE_M)
        )).flatten()

        # getInfo() on the whole FC so we get all properties; aggregate_array("soil_moisture") returns empty in GEE
        fc_info = fc.getInfo()
        features = fc_info.get("features") or []

        n_written = 0
        for feat in features:
            props = feat.get("properties") or {}
            d = props.get("date", "")
            sm = props.get("soil_moisture")
            q = props.get("quality")
            p = props.get("pass", "")
            sm_out = "nan" if (sm is None or (q is not None and q == 0)) else sm
            q_out = 0 if q is None else q
            writer.writerow([d, sm_out, q_out, p])
            n_written += 1
        total_rows += n_written
    csv_file.flush()
    csv_file.close()

    print(f"Wrote {total_rows} rows to {output_file}")
    return output_file


def main():
    """Main function."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Download SMAP L3 9 km basin soil moisture (AM/PM) using Google Earth Engine"
    )
    parser.add_argument("--station", type=str, default=STATION_CODE, help=f"Station code (default: {STATION_CODE})")
    parser.add_argument("--start-date", type=str, default=START_DATE, help=f"Start date YYYY-MM-DD (default: {START_DATE})")
    parser.add_argument("--end-date", type=str, default=END_DATE, help=f"End date YYYY-MM-DD (default: {END_DATE})")
    parser.add_argument("--project-id", type=str, default=None, help="GEE project ID (or set GEE_PROJECT_ID env var)")

    args = parser.parse_args()
    project_id = args.project_id or GEE_PROJECT_ID

    print("=" * 60)
    print("SMAP L3 9 km Daily Download - Google Earth Engine")
    print("=" * 60)
    print(f"Station: {args.station}")
    print(f"Date range: {args.start_date} to {args.end_date}")
    print(f"Project ID: {project_id}")
    print("=" * 60)

    try:
        output_file = download_smap_basin_csv(args.station, args.start_date, args.end_date)
        print("Done. Basin SMAP soil moisture (date, sm, quality, pass) saved to CSV.")
        print(f"  Output: {output_file}")
        return 0
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
