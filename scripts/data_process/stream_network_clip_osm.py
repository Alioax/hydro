# -*- coding: utf-8 -*-
"""
Build OSM stream network clipped to a CAMELS-FR basin and save as GeoPackage.

Fetches waterways from OpenStreetMap (Overpass API) for the basin bounding box,
clips them to the CAMELS-FR polygon, and saves to data/processed_data/
for use in mapping (e.g. fig_basin_outlet_map.py) and GIS.

Run from project root:  python scripts/data_process/stream_network_clip_osm.py
"""
from pathlib import Path
import urllib.request
import urllib.parse
import json

import geopandas as gpd
from shapely.geometry import LineString

# Project root (parent of scripts/)
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
GEOGRAPHY_DIR = PROJECT_ROOT / "data" / "camels_required" / "CAMELS_FR_geography"
STATION_CODE = "K287191001"
OUT_GPKG = PROJECT_ROOT / "data" / "processed_data" / "stream_network_osm_clipped.gpkg"

OVERPASS_URL = "https://overpass-api.de/api/interpreter"
BBOX_BUFFER = 0.02


def get_basin_wgs84():
    """CAMELS-FR basin as GeoDataFrame (one row) in WGS84."""
    boundaries = gpd.read_file(
        GEOGRAPHY_DIR / "CAMELS_FR_catchment_boundaries.gpkg"
    )
    row = boundaries[boundaries["sta_code_h3"] == STATION_CODE]
    if len(row) == 0:
        raise ValueError(f"Station {STATION_CODE} not found.")
    gdf = row.to_crs(epsg=4326) if row.crs and str(row.crs) != "EPSG:4326" else row
    return gdf


def fetch_osm_waterways(south, west, north, east, timeout_sec=120):
    """
    Fetch OSM ways with waterway tag in bbox via Overpass API.
    Returns list of Shapely LineString geometries.
    """
    query = f"""
    [out:json][timeout:{timeout_sec}];
    (
      way["waterway"]({south},{west},{north},{east});
    );
    out geom;
    """
    data = urllib.parse.urlencode({"data": query}).encode("utf-8")
    req = urllib.request.Request(
        OVERPASS_URL,
        data=data,
        method="POST",
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    with urllib.request.urlopen(req, timeout=timeout_sec + 10) as resp:
        result = json.loads(resp.read().decode("utf-8"))
    lines = []
    for el in result.get("elements", []):
        if el.get("type") != "way":
            continue
        geom = el.get("geometry")
        if not geom or len(geom) < 2:
            continue
        coords = [(n["lon"], n["lat"]) for n in geom]
        try:
            line = LineString(coords)
            if not line.is_valid:
                line = line.buffer(0)
            if not line.is_empty:
                lines.append(line)
        except Exception:
            continue
    return lines


def main():
    print("Loading CAMELS-FR basin...")
    basin_gdf = get_basin_wgs84()
    w, s, e, n = basin_gdf.total_bounds
    south = s - BBOX_BUFFER
    west = w - BBOX_BUFFER
    north = n + BBOX_BUFFER
    east = e + BBOX_BUFFER
    print(f"  Basin bbox (W,S,E,N): {w:.5f}, {s:.5f}, {e:.5f}, {n:.5f}")

    print("Fetching waterways from OpenStreetMap (Overpass API)...")
    lines = fetch_osm_waterways(south, west, north, east)
    if not lines:
        print("  No waterways returned. Check bbox or Overpass API.")
        return 1
    print(f"  Got {len(lines)} waterway segments")

    streams_gdf = gpd.GeoDataFrame(
        {"geometry": lines},
        crs="EPSG:4326",
    )
    print("Clipping streams to CAMELS-FR basin...")
    clipped = gpd.clip(streams_gdf, basin_gdf)
    if clipped.empty or clipped.geometry.is_empty.all():
        print("  No segments inside basin after clip.")
        return 1
    clipped = clipped[~clipped.geometry.is_empty].explode(index_parts=False).reset_index(drop=True)
    print(f"  Segments inside basin: {len(clipped)}")

    OUT_GPKG.parent.mkdir(parents=True, exist_ok=True)
    clipped.to_file(OUT_GPKG, driver="GPKG")
    print(f"  Saved: {OUT_GPKG}")
    return 0


if __name__ == "__main__":
    exit(main())
