# CAMELS-FR Required Data

This directory contains a minimal subset of the full CAMELS-FR dataset, copied for use by the project scripts. The full CAMELS-FR (>1 GB) can be archived.

## Contents

### Attributes
- **static_attributes/**
  - `CAMELS_FR_station_general_attributes.csv` – station metadata
  - `CAMELS_FR_topography_general_attributes.csv` – topography metrics
  - `CAMELS_FR_land_cover_attributes.csv` – land cover (used by gen_basin_table)
- **time_series_statistics/**
  - `CAMELS_FR_hydrological_signatures.csv` – hydrologic signatures (used by gen_basin_table)

### Geography
- `CAMELS_FR_catchment_boundaries.gpkg` – catchment polygons
- `CAMELS_FR_gauge_outlet.gpkg` – gauge outlet points

### Time series (station K287191001 – La Dore)
- **daily/** `CAMELS_FR_tsd_K287191001.csv` – daily discharge, precipitation, PET, temperature

## Scripts using this data

- `basin_metrics_table.py`, `fig_basin_outlet_map.py`, `hypsometry.py`, `stream_network_clip_osm.py`
- `timeseries_qc_aggregate.py`, `water_balance_plots.py` (via processed_data)
- `gen_basin_table.py`
- `download_dem_gee.py`, `download_imerg_gee.py`, `download_smap_gee.py`
