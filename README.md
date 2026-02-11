# Hydrology Research Project

Hydrology research project using the CAMELS-FR dataset and remote sensing (GPM IMERG, SMAP, DEM), supervised by Afshin Ashrafzadeh.

## Project Structure

```
.
├── data/
│   ├── camels_required/        # Minimal CAMELS-FR subset used by scripts (see data/camels_required/README.md)
│   ├── CAMELS-FR/              # Full CAMELS-FR (can be archived; scripts use camels_required)
│   └── raw_data/               # Remote sensing: GPM_IMERG, SMAP, DEM
├── scripts/
│   ├── data_download/          # GEE download scripts + requirements_gee.txt
│   ├── data_process/           # (placeholder)
│   ├── gen_basin_table.py      # LaTeX basin table
│   └── fig_basin_outlet_map.py # Basin + outlet map figure
├── reports/report 1/           # LaTeX report, figs, basin_table_rows.tex
└── .archive/                   # Optional dataset archives
```

## Data

**CAMELS-FR** — scripts use `data/camels_required/` (minimal subset; full dataset in `data/CAMELS-FR/` can be archived):

- `CAMELS_FR_attributes/` — station, topography, land cover, hydrological signatures
- `CAMELS_FR_geography/` — catchment boundaries and gauge outlets (`.gpkg`)
- `CAMELS_FR_time_series/daily/` — daily time series for station K287191001

**Remote sensing** (under `data/raw_data/`, created by scripts):

- `GPM_IMERG/daily/{station_code}/` — daily basin-averaged precipitation (CSV)
- `SMAP/daily/{station_code}/` — daily basin soil moisture AM/PM (CSV)
- `DEM/` — DEM GeoTIFFs (e.g. SRTM; GEE exports to Drive, then move here)

## Scripts

Run from project root.

### Report / figures

| Script | Purpose | Output |
|--------|---------|--------|
| `python scripts/gen_basin_table.py` | Build LaTeX table rows (basins 500–1000 km²) | `reports/report 1/basin_table_rows.tex` |
| `python scripts/fig_basin_outlet_map.py` | Map: catchment + outlet on basemap | `reports/report 1/figs/map_basin_outlet.pdf` |

### Data process (analysis and report figures)

Requires `scripts/data_process/requirements_process.txt` (geopandas, matplotlib, rasterio, numpy, pandas). DEM must be present at `data/raw_data/DEM/srtm/` or `data/remote sensing/DEM/srtm/` for station K287191001.

| Script | Purpose | Output |
|--------|---------|--------|
| `python scripts/data_process/dem_stream_processing.py` | DEM clip, flow routing, stream network, main channel | `reports/report 1/figs/` (elevation, slope, stream_network, longitudinal_profile), `data/processed_data/dem_metrics.npz` |
| `python scripts/data_process/basin_metrics_table.py` | Basin physiography “hero” table | `reports/report 1/basin_metrics_table.tex` |
| `python scripts/data_process/hypsometry.py` | Hypsometric and altimetric curves | `reports/report 1/figs/hypsometric_curve.pdf`, `altimetric_curve.pdf` |
| `python scripts/data_process/timeseries_qc_aggregate.py` | Time series QC, monthly/annual aggregates | `data/processed_data/monthly_*.csv`, `annual_*.csv`, `daily_*.csv` |
| `python scripts/data_process/water_balance_plots.py` | Time series snippet (2017–2020), climatology, FDC, annual water balance | `reports/report 1/figs/` (timeseries_flow_precip_temp_2017_2020, monthly_climatology, flow_duration_curve, annual_water_balance) |

Recommended order: run `dem_stream_processing` and `timeseries_qc_aggregate` first, then `basin_metrics_table`, `hypsometry`, and `water_balance_plots`.

### Data download (Google Earth Engine)

Requires `earthengine authenticate` and `scripts/data_download/requirements_gee.txt`.

| Script | Purpose | Output |
|--------|---------|--------|
| `python scripts/data_download/download_imerg_gee.py` | Daily basin precipitation (GPM IMERG V07) | `data/raw_data/GPM_IMERG/daily/{station}/` CSV |
| `python scripts/data_download/download_smap_gee.py` | Daily basin soil moisture AM/PM (SMAP L3 9 km) | `data/raw_data/SMAP/daily/{station}/` CSV (date, sm, quality, pass) |
| `python scripts/data_download/download_dem_gee.py` | DEM for catchment (SRTM, ASTER, etc.) | Export to GEE Drive; move to `data/raw_data/DEM/` |

**Common options:** `--station` (default `K287191001`), `--start-date`, `--end-date`, `--project-id` (or env `GEE_PROJECT_ID`).

## Report

LaTeX in `reports/report 1/`: main file `hydrology_project_haghighi.tex`, figures in `figs/`, basin metrics table `basin_metrics_table.tex` (generated), appendix from `basin_table_rows.tex`. Build with `pdflatex hydrology_project_haghighi` (run twice for references).

## Dependencies

- **Report/figures:** `geopandas`, `matplotlib`, `contextily`
- **Data process (DEM, hypsometry, time series, plots):** see `scripts/data_process/requirements_process.txt` (`geopandas`, `matplotlib`, `rasterio`, `numpy`, `pandas`)
- **GEE downloads:** see `scripts/data_download/requirements_gee.txt` (`earthengine-api`, `geopandas`, `rasterio`, `numpy`)
