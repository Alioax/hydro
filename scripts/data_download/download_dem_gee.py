# -*- coding: utf-8 -*-
"""
Download DEM (Digital Elevation Model) data using Google Earth Engine.

This script downloads DEM data for a specified catchment using Google Earth Engine.
Available DEM datasets:
- SRTM (30m resolution, global coverage)
- ASTER GDEM (30m resolution)
- Copernicus DEM (30m resolution, Europe)
- NASADEM (30m resolution, improved SRTM)

Output: GeoTIFF files saved to data/raw_data/DEM/
"""
import os
import json
from pathlib import Path
import geopandas as gpd
import ee

# Get project root (parent of scripts/)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "camels_required"
RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw_data"
DEM_DIR = RAW_DATA_DIR / "DEM"

# Configuration
STATION_CODE = "K287191001"  # Default station
GEE_PROJECT_ID = os.getenv("GEE_PROJECT_ID", "surface-hydrology")
DEM_DATASET = "SRTM"  # Options: "SRTM", "ASTER", "COPERNICUS", "NASADEM"
SCALE = 30  # Resolution in meters

# DEM dataset mapping
DEM_DATASETS = {
    "SRTM": "USGS/SRTMGL1_003",
    "ASTER": "NASA/ASTER_GED/AG002_003",
    "COPERNICUS": "COPERNICUS/DEM/GLO30",
    "NASADEM": "NASA/NASADEM_HGT/001"
}


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
    # Use __geo_interface__ to get GeoJSON dict directly from geometry
    geom = basin_gdf.iloc[0].geometry
    geojson_dict = geom.__geo_interface__
    
    geom_type = geojson_dict["type"]
    coords = geojson_dict["coordinates"]
    
    if geom_type == "MultiPolygon":
        # Take the largest polygon from MultiPolygon
        polygons = [ee.Geometry.Polygon(poly[0]) for poly in coords]
        areas = [p.area().getInfo() for p in polygons]
        return polygons[areas.index(max(areas))]
    elif geom_type == "Polygon":
        return ee.Geometry.Polygon(coords[0])
    else:
        raise ValueError(f"Unsupported geometry type: {geom_type}")


def download_dem(station_code, dem_dataset="SRTM", scale=30, output_dir=None):
    """Download DEM data for a catchment using Google Earth Engine."""
    if dem_dataset not in DEM_DATASETS:
        raise ValueError(f"Unknown DEM dataset: {dem_dataset}")
    
    if not initialize_gee():
        raise RuntimeError("Failed to initialize Google Earth Engine")
    
    print(f"Loading catchment boundary for station {station_code}...")
    basin = get_catchment_boundary(station_code)
    bounds = basin.total_bounds
    print(f"Catchment bounds (WGS84): West={bounds[0]:.6f}, South={bounds[1]:.6f}, East={bounds[2]:.6f}, North={bounds[3]:.6f}")
    
    ee_geom = basin_to_ee_geometry(basin)
    
    print(f"Loading {dem_dataset} DEM dataset...")
    dem = ee.Image(DEM_DATASETS[dem_dataset]).clip(ee_geom)
    
    if output_dir is None:
        output_dir = DEM_DIR / dem_dataset.lower()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"DEM_{dem_dataset}_{station_code}_{scale}m.tif"
    
    print(f"Downloading DEM to {output_file}...")
    
    task = ee.batch.Export.image.toDrive(
        image=dem,
        description=f"DEM_{dem_dataset}_{station_code}",
        folder="GEE_Exports",
        fileNamePrefix=f"DEM_{dem_dataset}_{station_code}_{scale}m",
        scale=scale,
        region=ee_geom,
        crs="EPSG:4326",
        fileFormat="GeoTIFF",
        maxPixels=1e9
    )
    
    task.start()
    print(f"Export task started. Task ID: {task.id}")
    print("Next steps:")
    print("1. Wait for export to complete in Google Earth Engine")
    print("2. Download from Google Drive (folder: GEE_Exports)")
    print(f"3. Move to: {output_file}")
    
    return output_file, task.id


def main():
    """Main function."""
    import argparse
    parser = argparse.ArgumentParser(description="Download DEM data using Google Earth Engine")
    parser.add_argument("--station", type=str, default=STATION_CODE, help=f"Station code (default: {STATION_CODE})")
    parser.add_argument("--dataset", type=str, default=DEM_DATASET, choices=list(DEM_DATASETS.keys()), help=f"DEM dataset (default: {DEM_DATASET})")
    parser.add_argument("--scale", type=int, default=SCALE, help=f"Resolution in meters (default: {SCALE})")
    parser.add_argument("--project-id", type=str, default=None, help="GEE project ID (or set GEE_PROJECT_ID env var)")
    
    args = parser.parse_args()
    project_id = args.project_id or GEE_PROJECT_ID
    
    print("=" * 60)
    print("DEM Download Script - Google Earth Engine")
    print("=" * 60)
    print(f"Station: {args.station}")
    print(f"Dataset: {args.dataset}")
    print(f"Resolution: {args.scale}m")
    print(f"Project ID: {project_id}")
    print("=" * 60)
    
    try:
        output_file, task_id = download_dem(args.station, args.dataset, args.scale)
        print(f"Download initiated successfully!")
        print(f"  Output: {output_file}")
        print(f"  Task ID: {task_id}")
        return 0
    except Exception as e:
        print(f"Error: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
