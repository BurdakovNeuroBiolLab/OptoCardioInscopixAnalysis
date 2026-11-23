# OptoCardioInscopixAnalysis

MATLAB analysis code for opto / cardio / Inscopix experiments.

## Layout

| Folder | Contents |
|---|---|
| `inscopix/` | Inscopix (IDPS) export analysis |
| `heka/` | HEKA data extraction |

## Scripts

- `inscopix/plot_cellmap_and_traces.m`: add `inscopix/` to the MATLAB path,
  set the current folder to the recording folder and run the script. It needs
  exactly one `*celltraces.csv`, one `*celltraces-props.csv` and one cell map
  `.tif`/`.tiff` in that folder. Plots the cell map
  with colour-coded, labelled ROIs (left) next to the stacked traces in matching
  colours (right).
