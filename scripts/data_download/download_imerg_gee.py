# -*- coding: utf-8 -*-
"""
Download daily GPM IMERG basin-averaged precipitation using Google Earth Engine.

Uses NASA/GPM_L3/IMERG_V07 (half-hourly). For each day we sum half-hourly
precipitation to get daily total per pixel, then compute area-weighted mean
over the catchment. Processed in one-year windows to avoid timeouts.

Date range: --start-date and --end-date (default 2010-01-01 to 2021-12-31).
Output: CSV to data/raw_data/GPM_IMERG/daily/{station_code}/
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
GPM_DIR = RAW_DATA_DIR / "GPM_IMERG" / "daily"

# Configuration
STATION_CODE = "K287191001"  # Default station
GEE_PROJECT_ID = os.getenv("GEE_PROJECT_ID", "surface-hydrology")
START_DATE = "2010-01-01"   # Default start (IMERG from 2000)
END_DATE = "2021-12-31"     # Default end
IMERG_ASSET = "NASA/GPM_L3/IMERG_V07"  # Half-hourly; we aggregate to daily
IMERG_SCALE_M = 11132      # ~0.1 deg in m at equator
PRECIP_BAND = "precipitation"  # V07: "precipitation" (calibrated) or "precipitationUncal"


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
        # Compute areas locally using Shapely to avoid server calls
        from shapely.geometry import Polygon as ShapelyPolygon
        shapely_polys = [ShapelyPolygon(poly[0]) for poly in coords]
        areas = [p.area for p in shapely_polys]
        largest_idx = areas.index(max(areas))
        return ee.Geometry.Polygon(coords[largest_idx][0])
    elif geom_type == "Polygon":
        return ee.Geometry.Polygon(coords[0])
    else:
        raise ValueError(f"Unsupported geometry type: {geom_type}")


def _daily_basin_feature(d, imerg_year, geom, band, scale_m):
    """
    For one day d (ee.Date): filter half-hourly to that day, sum precip (×0.5 hr → daily mm),
    then area-weighted mean over basin. Returns ee.Feature with date, precipitation_mm.
    """
    day_end = ee.Date(d).advance(1, "day")
    # Half-hourly band is mm/hr; sum over day then ×0.5 → daily mm per pixel
    daily_img = (
        imerg_year.filterDate(d, day_end)
        .select(band)
        .sum()
        .multiply(0.5)
    )
    mean_val = (
        daily_img.clip(geom)
        .reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=geom,
            scale=scale_m,
            maxPixels=1e9,
            bestEffort=True,
        )
        .get(band)
    )
    return ee.Feature(
        None,
        {
            "date": ee.Date(d).format("yyyy-MM-dd"),
            "precipitation_mm": mean_val,
        },
    )


def download_imerg_basin_csv(station_code, start_date, end_date, precip_band=None, output_dir=None):
    """
    Compute daily basin-averaged IMERG precipitation (daily sum then area-weighted mean).
    Uses IMERG V07 half-hourly; sums to daily total per pixel, then mean over basin.
    Loops by year to avoid getInfo() timeout.
    """
    if precip_band is None:
        precip_band = PRECIP_BAND
    if not initialize_gee():
        raise RuntimeError("Failed to initialize Google Earth Engine")

    print(f"Loading catchment boundary for station {station_code}...")
    basin = get_catchment_boundary(station_code)
    bounds = basin.total_bounds
    print(f"Catchment bounds (WGS84): West={bounds[0]:.6f}, South={bounds[1]:.6f}, East={bounds[2]:.6f}, North={bounds[3]:.6f}")
    ee_geom = basin_to_ee_geometry(basin)

    if output_dir is None:
        output_dir = GPM_DIR / station_code
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"IMERG_basin_precip_{station_code}_{start_date}_{end_date}.csv"

    start_year = int(start_date[:4])
    end_year = int(end_date[:4])

    csv_file = open(output_file, "w", newline="", encoding="utf-8")
    writer = csv.writer(csv_file)
    writer.writerow(["date", "precipitation_mm"])
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
        # Filter to full year for GEE so all days in range have data
        # filterDate end is exclusive; use next-year-01-01 to include full Dec 31
        imerg_year = (
            ee.ImageCollection(IMERG_ASSET)
            .filterDate(y_start, f"{year + 1}-01-01")
            .select(precip_band)
        )
        day_list = (
            ee.List.sequence(0, num_days - 1)
            .map(lambda i: ee.Date(y_start).advance(ee.Number(i), "day"))
        )
        fc = ee.FeatureCollection(
            day_list.map(
                lambda d: _daily_basin_feature(d, imerg_year, ee_geom, precip_band, IMERG_SCALE_M)
            )
        )
        dates = fc.aggregate_array("date").getInfo()
        precip_mm = fc.aggregate_array("precipitation_mm").getInfo()
        for d, p in zip(dates, precip_mm):
            writer.writerow([d, p if p is not None else ""])
        total_rows += len(dates)
    csv_file.close()

    print(f"Wrote {total_rows} rows to {output_file}")
    return output_file


def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Download daily GPM IMERG precipitation data using Google Earth Engine")
    parser.add_argument("--station", type=str, default=STATION_CODE, help=f"Station code (default: {STATION_CODE})")
    parser.add_argument("--start-date", type=str, default=START_DATE, help=f"Start date YYYY-MM-DD (default: {START_DATE})")
    parser.add_argument("--end-date", type=str, default=END_DATE, help=f"End date YYYY-MM-DD (default: {END_DATE})")
    parser.add_argument("--precip-band", type=str, default=PRECIP_BAND,
                       choices=["precipitation", "precipitationUncal"],
                       help=f"Precipitation band (default: {PRECIP_BAND}); V07 uses 'precipitation' (calibrated)")
    parser.add_argument("--project-id", type=str, default=None, help="GEE project ID (or set GEE_PROJECT_ID env var)")
    
    args = parser.parse_args()
    project_id = args.project_id or GEE_PROJECT_ID
    
    print("=" * 60)
    print("GPM IMERG Daily Download Script - Google Earth Engine")
    print("=" * 60)
    print(f"Station: {args.station}")
    print(f"Date range: {args.start_date} to {args.end_date}")
    print(f"Precipitation band: {args.precip_band}")
    print(f"Project ID: {project_id}")
    print("=" * 60)
    
    try:
        output_file = download_imerg_basin_csv(
            args.station, args.start_date, args.end_date, args.precip_band
        )
        print("Done. Basin daily precipitation (area-weighted mean) saved to CSV.")
        print(f"  Output: {output_file}")
        return 0
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
