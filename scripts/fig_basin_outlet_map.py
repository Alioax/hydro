# -*- coding: utf-8 -*-
"""
Generate a map figure showing the catchment boundary and outlet gauge location
with satellite imagery basemap (Google Earth-like), and a geographic elevation map
from the DEM in a right-hand subplot.

Output: PDF figure saved to reports/report 1/figs/map_basin_outlet.pdf
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib as mpl
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import contextily as ctx
import rasterio
from rasterio import features
from rasterio.warp import reproject, Resampling
from rasterio.crs import CRS
from rasterio.transform import from_bounds

mpl.rcParams["figure.dpi"] = 800
plt.rcParams["font.family"] = "Times New Roman"
mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
    color=["#FF5F05", "#13294B", "#009FD4", "#FCB316",
           "#006230", "#007E8E", "#5C0E41", "#7D3E13"]
)

# Get project root (parent of scripts/)
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "camels_required"
GEOGRAPHY_DIR = DATA_DIR / "CAMELS_FR_geography"
FIG_DIR = PROJECT_ROOT / "reports" / "report 1" / "figs"
STREAMS_GPKG = PROJECT_ROOT / "data" / "processed_data" / "stream_network_osm_clipped.gpkg"

# Station code to plot (default: K287191001 - La Dore)
STATION_CODE = "K287191001"

DEM_CANDIDATES = [
    PROJECT_ROOT / "data" / "raw_data" / "DEM" / "srtm" / f"DEM_SRTM_{STATION_CODE}_30m.tif",
    PROJECT_ROOT / "data" / "remote sensing" / "DEM" / "srtm" / f"DEM_SRTM_{STATION_CODE}_30m.tif",
]


def find_dem():
    for p in DEM_CANDIDATES:
        if p.exists():
            return p
    return None


def get_dem_raster_masked_to_basin(dem_path, basin_geom, basin_crs, padding_fraction=0.05,
                                   target_crs=None, target_bounds=None):
    """
    Read DEM windowed to basin bbox, mask to basin. Returns (data_2d, extent, crs).
    If target_crs and target_bounds are given, reproject to that CRS so the shape
    matches (e.g. Web Mercator for same aspect as left map). extent = (left, right, bottom, top).
    The read window is clipped to the raster bounds so areas outside the TIFF stay NaN (no 0 fill).
    """
    with rasterio.open(dem_path) as src:
        dem_crs = src.crs
        transform = src.transform
        nodata = src.nodata if src.nodata is not None else -9999
        raster_bounds = src.bounds  # (left, bottom, right, top)
    basin_in_dem_crs = gpd.GeoSeries([basin_geom], crs=basin_crs).to_crs(dem_crs)[0]
    minx, miny, maxx, maxy = basin_in_dem_crs.bounds
    dx, dy = maxx - minx, maxy - miny
    pad_x, pad_y = dx * padding_fraction, dy * padding_fraction
    req_left = minx - pad_x
    req_bottom = miny - pad_y
    req_right = maxx + pad_x
    req_top = maxy + pad_y
    # Clip requested bounds to actual raster extent so we never read outside the file
    # (reading outside fills with 0 and causes a false "0 elevation" strip)
    clip_left = max(req_left, raster_bounds.left)
    clip_bottom = max(req_bottom, raster_bounds.bottom)
    clip_right = min(req_right, raster_bounds.right)
    clip_top = min(req_top, raster_bounds.top)
    if clip_left >= clip_right or clip_bottom >= clip_top:
        # No overlap: return empty-like array and extent
        out_crs = CRS.from_epsg(target_crs) if target_crs else dem_crs
        if target_crs and target_bounds:
            left_t, bottom_t, right_t, top_t = target_bounds
            out_h = out_w = 256
            dst_transform = from_bounds(left_t, bottom_t, right_t, top_t, out_w, out_h)
            return np.full((out_h, out_w), np.nan), (left_t, right_t, bottom_t, top_t), out_crs
        return np.array([[np.nan]]), (0, 1, 0, 1), dem_crs
    window = rasterio.windows.from_bounds(clip_left, clip_bottom, clip_right, clip_top, transform=transform)
    with rasterio.open(dem_path) as src:
        data = src.read(1, window=window)
        win_transform = rasterio.windows.transform(window, src.transform)
    # Mask outside basin
    mask = features.geometry_mask(
        [basin_in_dem_crs],
        out_shape=data.shape,
        transform=win_transform,
        invert=True,
    )
    data = np.where(mask, data, np.nan)
    # Nodata and invalid
    data = np.where(np.isfinite(data) & (data != nodata), data, np.nan)
    # SRTM and many DEMs use 0 for nodata when metadata has no nodata; avoid false 0 m elevation
    if nodata != 0:
        data = np.where(data == 0, np.nan, data)
    out_crs = dem_crs
    left = win_transform.c
    top = win_transform.f
    right = win_transform.c + data.shape[1] * win_transform.a
    bottom = win_transform.f + data.shape[0] * win_transform.e
    extent = (left, right, bottom, top)

    if target_crs is not None and target_bounds is not None:
        # Reproject to target CRS (e.g. 3857) so shape matches the left map
        left_t, bottom_t, right_t, top_t = target_bounds
        # Use fixed output shape so full extent is covered; areas with no source stay NaN
        out_w = 400
        out_h = int(out_w * (top_t - bottom_t) / (right_t - left_t))
        out_h = max(out_h, 1)
        dst_transform = from_bounds(left_t, bottom_t, right_t, top_t, out_w, out_h)
        data_out = np.full((out_h, out_w), np.nan, dtype=np.float32)
        reproject(
            source=data.astype(np.float32),
            destination=data_out,
            src_transform=win_transform,
            src_crs=dem_crs,
            dst_transform=dst_transform,
            dst_crs=CRS.from_epsg(target_crs),
            resampling=Resampling.bilinear,
            src_nodata=np.nan,
            dst_nodata=np.nan,
        )
        data = data_out
        extent = (left_t, right_t, bottom_t, top_t)
        out_crs = CRS.from_epsg(target_crs)
    return data, extent, out_crs

# Read station label from attributes
station_label = STATION_CODE
with open(DATA_DIR / "CAMELS_FR_attributes" / "static_attributes" / "CAMELS_FR_station_general_attributes.csv", 
          'r', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter=';')
    for row in reader:
        if row['sta_code_h3'] == STATION_CODE:
            station_label = row['sta_label']
            break

# Read catchment boundaries and filter to our station
boundaries = gpd.read_file(DATA_DIR / "CAMELS_FR_geography" / "CAMELS_FR_catchment_boundaries.gpkg")
basin = boundaries[boundaries['sta_code_h3'] == STATION_CODE].copy()

# Read gauge outlets and filter to our station
outlets = gpd.read_file(DATA_DIR / "CAMELS_FR_geography" / "CAMELS_FR_gauge_outlet.gpkg")
outlet = outlets[outlets['sta_code_h3'] == STATION_CODE].copy()

if len(basin) == 0:
    raise ValueError(f"Station {STATION_CODE} not found in catchment boundaries")
if len(outlet) == 0:
    raise ValueError(f"Station {STATION_CODE} not found in gauge outlets")

# Reproject to Web Mercator (EPSG:3857) for basemap tiles
basin_web = basin.to_crs(epsg=3857)
outlet_web = outlet.to_crs(epsg=3857)

# Subplot mosaic: map (top), elevation (bottom); share x/y so axis labels don't repeat
fig, axd = plt.subplot_mosaic([['map'], ['elevation']], figsize=(3.5, 10), sharex=True, sharey=True)
ax = axd['map']
ax_elev = axd['elevation']

# Set axis limits based on basin bounds with some padding
bounds = basin_web.total_bounds
padding = (bounds[2] - bounds[0]) * 0.1  # 10% padding
ax.set_xlim(bounds[0] - padding, bounds[2] + padding)
ax.set_ylim(bounds[1] - padding, bounds[3] + padding)

# Add satellite imagery basemap (Esri World Imagery - Google Earth-like)
try:
    ctx.add_basemap(ax,
                   crs=basin_web.crs,
                   source=ctx.providers.Esri.WorldImagery,
                   zoom='auto',
                   attribution_size=2)
except Exception as e:
    print(f"Warning: Could not add Esri World Imagery basemap: {e}")
    print("Falling back to OpenStreetMap basemap...")
    try:
        ctx.add_basemap(ax,
                       crs=basin_web.crs,
                       source=ctx.providers.OpenStreetMap.Mapnik,
                       zoom='auto',
                       attribution_size=2)
    except Exception as e2:
        print(f"Warning: Could not add OpenStreetMap basemap: {e2}")
        print("Proceeding without basemap...")

# Plot catchment boundary (outline only, no fill) so satellite image is visible
basin_web.plot(ax=ax, facecolor='none', edgecolor='C0',
               linewidth=1.2, zorder=5)

# Overlay OSM stream network if available
if STREAMS_GPKG.exists():
    streams = gpd.read_file(STREAMS_GPKG).to_crs(epsg=3857)
    streams.plot(ax=ax, color='C2', linewidth=0.7, alpha=0.9, zorder=6)

# Plot outlet point (X mark, outline only, no fill)
xy = outlet_web.geometry.iloc[0]
ax.plot(xy.x, xy.y, marker='x', color='C3', markeredgewidth=2.5,
        markersize=14, fillstyle='none', linestyle='none', zorder=10)

# Axis labels: top subplot = y only; bottom will show x (vertical layout)
ax.set_ylabel('Northing (m)')
ax.tick_params(axis='both', which='major', labelsize=8)
ax.tick_params(axis='x', labelbottom=False)  # x-label and ticks on bottom subplot only
plt.setp(ax.get_xticklabels(), visible=False)

# Elevation imshow return value for colorbar (set when DEM is plotted)
im_elev = None

# Set equal aspect ratio for map
ax.set_aspect('equal')

# --- Right subplot: geographic elevation map from DEM (same CRS as left: 3857) ---
dem_path = find_dem()
map_bounds = (bounds[0] - padding, bounds[1] - padding, bounds[2] + padding, bounds[3] + padding)
if dem_path is not None:
    basin_geom = basin.iloc[0].geometry
    basin_crs = basin.crs
    try:
        dem_data, dem_extent, _ = get_dem_raster_masked_to_basin(
            dem_path, basin_geom, basin_crs,
            target_crs=3857, target_bounds=map_bounds
        )
        cmap = plt.cm.terrain
        cmap.set_bad('none', alpha=0)
        im_elev = ax_elev.imshow(dem_data, extent=dem_extent, cmap=cmap, origin='upper',
                                 aspect='equal', interpolation='bilinear')
        ax_elev.set_xlim(dem_extent[0], dem_extent[1])
        ax_elev.set_ylim(dem_extent[2], dem_extent[3])
        # Outlet marker on elevation subplot (same CRS as map)
        ax_elev.plot(xy.x, xy.y, marker='x', color='C3', markeredgewidth=2.5,
                     markersize=14, fillstyle='none', linestyle='none', zorder=10)
        ax_elev.set_xlabel('Easting (m)')
        ax_elev.set_ylabel('Northing (m)')
        ax_elev.tick_params(axis='both', which='major', labelsize=8)
    except Exception as e:
        ax_elev.text(0.5, 0.5, f'Could not load DEM: {e}', ha='center', va='center',
                    transform=ax_elev.transAxes, fontsize=8)
        ax_elev.set_xticks([])
        ax_elev.set_yticks([])
else:
    ax_elev.text(0.5, 0.5, 'DEM not found', ha='center', va='center',
                transform=ax_elev.transAxes)
    ax_elev.set_xticks([])
    ax_elev.set_yticks([])

# Remove top and right spines from both subplots
for ax_i in (ax, ax_elev):
    ax_i.spines[['top', 'right']].set_visible(False)

# Legend: outside, top center (above both subplots), one row, no frame
legend_elements = [
    Line2D([0], [0], color='C0', linewidth=1.2, label='Catchment boundary'),
    Line2D([0], [0], color='C2', linewidth=0.7, alpha=0.9, label='Stream network'),
    Line2D([0], [0], color='C3', marker='x', linestyle='None', markersize=10,
           markeredgewidth=2, label='Monitoring station')
]
if not STREAMS_GPKG.exists():
    legend_elements = [legend_elements[0], legend_elements[2]]
fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=1, frameon=False, fontsize=10)

# Tight layout first so both subplots get identical width and stay aligned
plt.tight_layout()

# Colorbar: add after layout so it doesn't resize subplots; place to the right of bottom subplot, same height
if im_elev is not None:
    pos = ax_elev.get_position()
    pad, cbar_width = 0.02, 0.02
    cax = fig.add_axes([pos.x1 + pad, pos.y0, cbar_width, pos.height])
    cbar = fig.colorbar(im_elev, cax=cax, label='Elevation (m a.s.l.)')
    cbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True, steps=[1, 2, 5, 10]))
    cbar.ax.tick_params(labelsize=8)

# Ensure output directory exists
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Save as PDF
output_file = FIG_DIR / "map_basin_outlet.pdf"
plt.savefig(output_file, format='pdf', bbox_inches='tight', dpi=800)
print(f"Map figure saved to: {output_file}")

plt.close()
