# -*- coding: utf-8 -*-
import csv
import os
from pathlib import Path

def escape_tex(s):
    if not s or s == 'NA':
        return '---'
    s = str(s)
    s = s.replace('\\', '\\textbackslash{}')
    s = s.replace('&', '\\&')
    s = s.replace('%', '\\%')
    s = s.replace('#', '\\#')
    s = s.replace('_', '\\_')
    s = s.replace('{', '\\{')
    s = s.replace('}', '\\}')
    return s

def num(s, fmt='.1f'):
    if s is None or s == '' or s == 'NA':
        return '---'
    try:
        x = float(s)
        return format(x, fmt) if fmt else str(int(round(x)))
    except (ValueError, TypeError):
        return '---'

# Get project root (parent of scripts/)
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "camels_required"
REPORT_DIR = PROJECT_ROOT / "reports" / "report 1"

base_attr = DATA_DIR / "CAMELS_FR_attributes"
base_static = base_attr / "static_attributes"
base_ts = base_attr / "time_series_statistics"

# Station general
with open(base_static / "CAMELS_FR_station_general_attributes.csv", 'r', encoding='utf-8') as f:
    r = csv.DictReader(f, delimiter=';')
    rows = list(r)

# Land cover (sta_code_h3 -> clc_2018_lvl1_2, clc_2018_lvl1_3)
land = {}
with open(base_static / "CAMELS_FR_land_cover_attributes.csv", 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter=';'):
        land[row['sta_code_h3']] = {
            'agr': row.get('clc_2018_lvl1_2'),
            'forest': row.get('clc_2018_lvl1_3'),
        }

# Hydrological signatures (sta_code_h3 -> hyd_q_mean_yr)
hyd = {}
with open(base_ts / "CAMELS_FR_hydrological_signatures.csv", 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter=';'):
        hyd[row['sta_code_h3']] = row.get('hyd_q_mean_yr')

candidates = []
for row in rows:
    try:
        a = float(row['sta_area_snap'])
        if 500 <= a <= 1000:
            start = row.get('sta_date_start', '')[:10] if row.get('sta_date_start') else '---'
            end = (row.get('sta_date_end', '') or '').strip()
            if not end or end == 'NA':
                end = 'present'
            else:
                end = end[:10]
            code = row['sta_code_h3']
            outlet_alt = row.get('sta_altitude_snap')
            lc = land.get(code, {})
            agr = lc.get('agr')
            mean_flow = hyd.get(code)
            candidates.append((a, code, row['sta_label'], start, end, outlet_alt, mean_flow, agr))
    except (ValueError, KeyError):
        pass
candidates.sort(key=lambda x: x[0])

lines = []
for i, (a, code, label, start, end, outlet_alt, mean_flow, agr) in enumerate(candidates, 1):
    label_tex = escape_tex(label)
    code_tex = code.replace('_', '\\_')
    outlet_str = num(outlet_alt, '.0f')
    flow_str = num(mean_flow, '.1f')
    agr_str = num(agr, '.1f')
    lines.append(f"{i} & {a:.1f} & \\texttt{{{code_tex}}} & {label_tex} & {start}--{end} & {outlet_str} & {flow_str} & {agr_str}\\% \\\\")

out = REPORT_DIR / "basin_table_rows.tex"
with open(out, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print(f"Wrote {len(lines)} rows to {out}")
