# -*- coding: utf-8 -*-
"""
Generate the basin physiography "hero" table for the report.

Reads CAMELS-FR attributes (station, topography), optional DEM metrics from
dem_stream_processing output, and boundary geometry for perimeter and basin length.
Writes LaTeX table to reports/report 1/basin_metrics_table.tex.
"""
import csv
from pathlib import Path
import geopandas as gpd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "camels_required"
REPORT_DIR = PROJECT_ROOT / "reports" / "report 1"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed_data"
STATION_CODE = "K287191001"


def _num(x, fmt=".1f"):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "---"
    try:
        v = float(x)
        return format(v, fmt) if fmt else str(int(round(v)))
    except (ValueError, TypeError):
        return "---"


def run():
    # Station attributes
    station_row = None
    with open(DATA_DIR / "CAMELS_FR_attributes" / "static_attributes" / "CAMELS_FR_station_general_attributes.csv", "r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter=";"):
            if row["sta_code_h3"] == STATION_CODE:
                station_row = row
                break
    if not station_row:
        raise ValueError(f"Station {STATION_CODE} not found")

    # Topography attributes
    topo_row = None
    with open(DATA_DIR / "CAMELS_FR_attributes" / "static_attributes" / "CAMELS_FR_topography_general_attributes.csv", "r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter=";"):
            if row["sta_code_h3"] == STATION_CODE:
                topo_row = row
                break
    if not topo_row:
        raise ValueError(f"Topography for {STATION_CODE} not found")

    # Geometry: perimeter, centroid, basin length (max distance in basin or centroid-to-outlet * 2 approx), mean width
    boundaries = gpd.read_file(DATA_DIR / "CAMELS_FR_geography" / "CAMELS_FR_catchment_boundaries.gpkg")
    basin = boundaries[boundaries["sta_code_h3"] == STATION_CODE].copy()
    if len(basin) == 0:
        raise ValueError(f"Boundary for {STATION_CODE} not found")
    geom = basin.iloc[0].geometry
    # Reproject to metric for lengths
    basin_metric = basin.to_crs(epsg=2154)  # Lambert 93 France, meters
    geom_m = basin_metric.iloc[0].geometry
    area_km2 = station_row.get("sta_area_snap")
    try:
        area_km2 = float(area_km2)
    except (TypeError, ValueError):
        area_km2 = geom_m.area / 1e6
    perimeter_km = geom_m.length / 1000.0
    centroid = geom_m.centroid
    # Basin length: approximate as 2 * max distance from centroid to boundary (sample points)
    try:
        if geom_m.geom_type == "Polygon" and geom_m.exterior is not None:
            pts = [geom_m.exterior.interpolate(i, normalized=True) for i in [0, 0.25, 0.5, 0.75]]
        else:
            pts = [geom_m.centroid]
        basin_length_km = 2 * max(centroid.distance(p) for p in pts) / 1000.0
    except Exception:
        basin_length_km = np.sqrt(area_km2 * 1e6) * 2 / 1000.0  # rough
    if basin_length_km <= 0:
        basin_length_km = 2 * np.sqrt(area_km2 * 1e6 / np.pi) / 1000.0  # diameter of circle of same area
    mean_width_km = area_km2 / basin_length_km if basin_length_km > 0 else "---"

    # DEM-derived (optional)
    dem_metrics = {}
    npz_path = PROCESSED_DIR / "dem_metrics.npz"
    if npz_path.exists():
        with np.load(npz_path, allow_pickle=True) as z:
            dem_metrics = {k: float(v) if np.isscalar(v) or (hasattr(v, "item") and v.ndim == 0) else v for k, v in z.items()}

    # Outlet elevation: basin min cannot be below outlet
    try:
        outlet_alt_val = float(station_row.get("sta_altitude_snap") or 398)
    except (TypeError, ValueError):
        outlet_alt_val = 398.0
    # Elevation: from topography quantiles we have min/max; CAMELS top_altitude_mean
    top_alt_mean = topo_row.get("top_altitude_mean")
    try:
        top_alt_mean = float(top_alt_mean)
    except (TypeError, ValueError):
        top_alt_mean = None
    # Relief from topography quantiles (100th - 0th) or from DEM
    elev_min = dem_metrics.get("dem_elev_min")
    elev_max = dem_metrics.get("dem_elev_max")
    elev_mean = dem_metrics.get("dem_elev_mean")
    if elev_min is None or (isinstance(elev_min, float) and np.isnan(elev_min)):
        elev_min = outlet_alt_val
    else:
        elev_min = float(elev_min)
        if elev_min < outlet_alt_val:
            elev_min = outlet_alt_val  # DEM nodata/zeros can yield spurious min
    if elev_max is None or (isinstance(elev_max, float) and np.isnan(elev_max)):
        elev_max = 1630  # from quantiles 100
    else:
        elev_max = float(elev_max)
    if elev_mean is None or (isinstance(elev_mean, float) and np.isnan(elev_mean)):
        elev_mean = top_alt_mean
    relief = elev_max - elev_min

    main_channel_km = dem_metrics.get("main_channel_length_km")
    if main_channel_km is None:
        # Fallback: longest OSM stream segment in basin (main stem proxy)
        stream_gpkg = PROCESSED_DIR / "stream_network_osm_clipped.gpkg"
        if stream_gpkg.exists():
            try:
                streams = gpd.read_file(stream_gpkg)
                if not streams.empty and streams.geometry.notna().any():
                    streams_metric = streams.to_crs(epsg=2154)
                    lengths_m = streams_metric.geometry.length
                    main_channel_km = float(lengths_m.max()) / 1000.0
            except Exception:
                pass
    max_strahler = dem_metrics.get("max_strahler")
    mean_slope_deg = dem_metrics.get("mean_slope_deg")
    if mean_slope_deg is None or (isinstance(mean_slope_deg, float) and np.isnan(mean_slope_deg)):
        mean_slope_deg = topo_row.get("top_slo_mean")
    try:
        mean_slope_deg = float(mean_slope_deg)
    except (TypeError, ValueError):
        mean_slope_deg = None

    area = _num(area_km2, ".1f")
    perimeter = _num(perimeter_km, ".1f")
    basin_len = _num(basin_length_km, ".1f")
    main_ch = _num(main_channel_km, ".2f") if main_channel_km is not None else "---"
    elev_mi = _num(elev_min, ".0f")
    elev_mx = _num(elev_max, ".0f")
    elev_mu = _num(elev_mean, ".0f")
    relief_s = _num(relief, ".0f")
    slope_s = _num(mean_slope_deg, ".2f")
    compact = _num(topo_row.get("top_mor_compact_coef"), ".2f")
    circ = _num(topo_row.get("top_mor_circ_ratio"), ".2f")
    drain_dens = _num(topo_row.get("top_drainage_density"), ".4f")
    strahler_s = _num(max_strahler, ".0f") if max_strahler is not None else "---"
    outlet_alt = _num(station_row.get("sta_altitude_snap"), ".0f")

    lines = [
        "\\begin{tabular}{ll}",
        "\\toprule",
        "Metric & Value \\\\",
        "\\midrule",
        f"Area (km²) & {area} \\\\",
        f"Perimeter (km) & {perimeter} \\\\",
        f"Basin length (km) & {basin_len} \\\\",
        f"Main channel length (km) & {main_ch} \\\\",
        f"Elevation min (m\\,a.s.l.) & {elev_mi} \\\\",
        f"Elevation mean (m\\,a.s.l.) & {elev_mu} \\\\",
        f"Elevation max (m\\,a.s.l.) & {elev_mx} \\\\",
        f"Relief (m) & {relief_s} \\\\",
        f"Outlet elevation (m\\,a.s.l.) & {outlet_alt} \\\\",
        f"Mean slope (degrees) & {slope_s} \\\\",
        f"Compactness coefficient & {compact} \\\\",
        f"Circularity ratio & {circ} \\\\",
        f"Drainage density (km$^{{-1}}$) & {drain_dens} \\\\",
        f"Max Strahler order & {strahler_s} \\\\",
        "\\bottomrule",
        "\\end{tabular}",
    ]
    out_path = REPORT_DIR / "basin_metrics_table.tex"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Wrote {out_path}")
    return out_path


if __name__ == "__main__":
    run()
