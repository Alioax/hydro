# -*- coding: utf-8 -*-
"""
Hypsometric curve and altimetric curve from DEM (basin-masked).

- Hypsometric curve: normalized elevation (z-min)/(max-min) vs normalized cumulative area (0--1).
- Hypsometric integral: area under the hypsometric curve.
- Altimetric curve: cumulative area (km² or %) vs elevation (m a.s.l.).

Outputs: reports/report 1/figs/hypsometric_curve.pdf, altimetric_curve.pdf.
"""
import csv
from pathlib import Path
import numpy as np
import geopandas as gpd
import rasterio
from rasterio import features
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.crs import CRS
from rasterio import windows
import matplotlib as mpl
import matplotlib.pyplot as plt


mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
    color=["#FF5F05", "#13294B", "#009FD4", "#FCB316",
           "#006230", "#007E8E", "#5C0E41", "#7D3E13"]
)
mpl.rcParams["figure.dpi"] = 800
plt.rcParams["font.family"] = "Times New Roman"

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "camels_required"
FIG_DIR = PROJECT_ROOT / "reports" / "report 1" / "figs"
STATION_CODE = "K287191001"
DEM_CANDIDATES = [
    PROJECT_ROOT / "data" / "raw_data" / "DEM" / "srtm" / f"DEM_SRTM_{STATION_CODE}_30m.tif",
    PROJECT_ROOT / "data" / "remote sensing" / "DEM" / "srtm" / f"DEM_SRTM_{STATION_CODE}_30m.tif",
]


def find_dem():
    for p in DEM_CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError(f"DEM not found. Tried: {DEM_CANDIDATES}")


def get_dem_masked_to_basin(dem_path, basin_geom, basin_crs):
    """Read DEM, mask to basin. Returns 1D array of elevations (valid pixels only)."""
    with rasterio.open(dem_path) as src:
        dem_crs = src.crs
        data = src.read(1)
        transform = src.transform
        nodata = src.nodata if src.nodata is not None else -9999
    basin_in_dem_crs = gpd.GeoSeries([basin_geom], crs=basin_crs).to_crs(dem_crs)[0]
    mask = features.geometry_mask(
        [basin_in_dem_crs],
        out_shape=data.shape,
        transform=transform,
        invert=True,
    )
    data = np.where(mask, data, nodata)
    valid = (data != nodata) & np.isfinite(data)
    return data[valid].ravel()


def run():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    dem_path = find_dem()
    boundaries = gpd.read_file(DATA_DIR / "CAMELS_FR_geography" / "CAMELS_FR_catchment_boundaries.gpkg")
    basin = boundaries[boundaries["sta_code_h3"] == STATION_CODE]
    if len(basin) == 0:
        raise ValueError(f"Station {STATION_CODE} not found")
    basin_geom = basin.iloc[0].geometry
    basin_crs = basin.crs

    elev = get_dem_masked_to_basin(dem_path, basin_geom, basin_crs)
    if len(elev) == 0:
        raise ValueError("No valid DEM pixels in basin")
    elev = np.sort(elev)
    z_min, z_max = elev.min(), elev.max()
    n = len(elev)
    # Cumulative area (number of pixels = area in pixel units; we normalize to 0--1)
    cum_area_norm = np.arange(1, n + 1, dtype=float) / n
    # Normalized elevation
    z_norm = (elev - z_min) / (z_max - z_min) if z_max > z_min else np.zeros_like(elev)

    # Hypsometric integral = area under curve (z_norm vs cum_area_norm)
    hypsometric_integral = np.trapz(z_norm, cum_area_norm)
    # Cumulative area in km²: scale by CAMELS basin area (not pixel count)
    basin_area_km2 = 795.0
    with open(DATA_DIR / "CAMELS_FR_attributes" / "static_attributes" / "CAMELS_FR_station_general_attributes.csv", "r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter=";"):
            if row["sta_code_h3"] == STATION_CODE:
                try:
                    basin_area_km2 = float(row.get("sta_area_snap") or basin_area_km2)
                except (TypeError, ValueError):
                    pass
                break
    cum_area_km2 = cum_area_norm * basin_area_km2

    # --- Hypsometric curve ---
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 4))
    ax.spines[["top", "right"]].set_visible(False)
    ax.plot(cum_area_norm, z_norm, color="C0", linewidth=1.5)
    ax.set_xlabel("Normalized cumulative area")
    ax.set_ylabel("Normalized elevation")
    label_fontsize = plt.rcParams["axes.labelsize"]
    ax.text(0.02, 0.02, f"integral = {hypsometric_integral:.3f}", transform=ax.transAxes,
            fontsize=label_fontsize, verticalalignment="bottom")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "hypsometric_curve.pdf", bbox_inches="tight", dpi=800)
    plt.close()

    # --- Altimetric curve with hypsometric on secondary axis ---
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 4))
    ax.spines[["top", "right"]].set_visible(False)
    total_area_km2 = cum_area_km2[-1]
    ax.plot(cum_area_km2, elev, color="C0", linewidth=1.5)
    ax.set_xlim(0, total_area_km2)
    ax.set_xticks([0, 200, 400, 600, 795])
    ax.set_ylim(z_min, z_max)  # Must match right axis [0,1] for exact overlay
    ax.set_xlabel("Cumulative area (km²)")
    ax.set_ylabel("Elevation (m a.s.l.)")
    ax.grid(True, alpha=0.3)

    # Secondary y-axis (right): normalized elevation (same curve, different scale)
    ax2 = ax.twinx()
    ax2.set_ylabel("Normalized elevation")
    ax2.set_ylim(0, 1)

    # Secondary x-axis (top): normalized cumulative area
    ax_top = ax.secondary_xaxis(
        "top",
        functions=(lambda x: x / total_area_km2, lambda x: x * total_area_km2),
    )
    ax_top.set_xlim(0, 1)
    ax_top.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax_top.set_xlabel("Normalized cumulative area")

    label_fontsize = plt.rcParams["axes.labelsize"]
    ax.text(0.02, 0.02, f"integral = {hypsometric_integral:.3f}", transform=ax.transAxes,
            fontsize=label_fontsize, verticalalignment="bottom")

    plt.tight_layout()
    plt.savefig(FIG_DIR / "altimetric_curve.pdf", bbox_inches="tight", dpi=800)
    plt.close()

    # Save integral for report
    out_metrics = PROJECT_ROOT / "data" / "processed_data"
    out_metrics.mkdir(parents=True, exist_ok=True)
    np.savez(out_metrics / "hypsometry.npz", hypsometric_integral=hypsometric_integral)
    print(f"Hypsometric integral: {hypsometric_integral:.3f}")
    print(f"Figures saved to {FIG_DIR}")
    return hypsometric_integral


if __name__ == "__main__":
    run()
