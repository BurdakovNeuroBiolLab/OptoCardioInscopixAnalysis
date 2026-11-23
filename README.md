# OptoCardioInscopixAnalysis

MATLAB analysis code for opto / cardio / Inscopix experiments.

## Layout

| Folder | Contents |
|---|---|
| `inscopix/` | Inscopix (IDPS) export analysis |
| `heka/` | HEKA data extraction |

## Scripts

- `inscopix/plot_cellmap_and_traces.m`: pick a folder containing one IDPS
  cell props CSV, one cell traces CSV and one cell map TIFF. Plots the cell map
  with colour-coded, labelled ROIs (left) next to the stacked traces in matching
  colours (right).
