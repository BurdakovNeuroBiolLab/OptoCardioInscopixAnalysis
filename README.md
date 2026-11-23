# OptoCardioInscopixAnalysis

MATLAB analysis code for opto / cardio / Inscopix experiments.

## Layout

| Path | Contents |
|---|---|
| `B_AlignHekaToCellTraces.m` | Parent script: HEKA data downsampled onto the Inscopix cell-trace frames |
| `inscopix/` | Inscopix (IDPS) export analysis |
| `heka/` | HEKA data extraction |

## Usage

Set the MATLAB current folder to the recording folder and run the script by
name. Each recording folder holds exactly one of each:

- `*celltraces.csv`: IDPS cell traces
- `*celltraces-props.csv`: IDPS cell props
- `*.tif` / `*.tiff`: cell map with white ROI outlines
- `*L.mat`: HEKA export

## Scripts

- `B_AlignHekaToCellTraces.m`: runs the two scripts below, then bins
  `extractedMatGRINcrop` into one row per cell-trace frame (`hekaDS`, saved to
  `<traces>_hekaAligned.mat`). The GRIN synch high period is taken as the full
  Inscopix recording: its rising edge is the first frame, its falling edge the
  end of the last frame. Analogue channels are averaged per frame, run edge
  events are summed (encoder ticks per frame). Frames outside the GRIN/camera
  crop are NaN. Heart rate is computed from ECG R peaks (`hrDS`, bpm). The
  `plot_cellmap_and_traces` figure gets three panels under the cell traces on
  the same time axis: running (encoder ticks/s), heart rate and laser.
- `inscopix/plot_cellmap_and_traces.m`: plots the cell map with colour-coded,
  labelled ROIs (left) next to the stacked traces in matching colours (right).
- `heka/A_Extract_HekaData_GRIN.m`: loads the HEKA sweeps, concatenates
  channels, binarises digital channels and crops to the period where both
  the GRIN and camera sync are high (`extractedMatGRINcrop`).
