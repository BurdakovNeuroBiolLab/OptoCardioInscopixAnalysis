# OptoCardioInscopixAnalysis

MATLAB analysis code for opto / cardio / Inscopix experiments.

## Files

All scripts are in the repository root. Run `A_Parent_Align_Heka_Pupil_Inscopix`;
it calls the others.

| File | Contents |
|---|---|
| `A_Parent_Align_Heka_Pupil_Inscopix.m` | Parent script: pupil and HEKA data aligned to the Inscopix cell-trace frames |
| `plot_cellmap_and_traces.m` | Inscopix (IDPS) cell map and stacked traces |
| `Extract_HekaData_GRIN.m` | HEKA data extraction |
| `pupilCircleFit.m` | DeepLabCut pupil circle fit |
| `alignPupilToHeka.m` | Pupil diameter onto the HEKA time base |

## Usage

Add the repository folder to the MATLAB path, set the current folder to the
recording folder and run `A_Parent_Align_Heka_Pupil_Inscopix`
by name. Each recording folder holds exactly one of each:

- `*celltrace*.csv` without "props" in the name: IDPS cell traces
- `*celltrace*props*.csv`: IDPS cell props (`celltrace`/`celltraces`, `-props`/`_props` all work)
- `*.tif` / `*.tiff`: cell map with white ROI outlines
- `*L.mat`: HEKA export
- `*zoompup*.csv`: DeepLabCut pupil tracking

## Scripts

- `A_Parent_Align_Heka_Pupil_Inscopix.m`: runs the pupil fit, the cell-trace plot and
  the HEKA extraction, puts the pupil diameter on the HEKA time base
  (`alignPupilToHeka`, column 8 of `extractedMatGRINcrop`), then bins
  `extractedMatGRINcrop` into one row per cell-trace frame (`hekaDS`, saved to
  `<traces>_hekaAligned.mat`). The GRIN synch high period is taken as the full
  Inscopix recording: its rising edge is the first frame, its falling edge the
  end of the last frame. Analogue channels are averaged per frame, run edge
  events are summed (encoder ticks per frame). Frames outside the GRIN/camera
  crop are NaN. Heart rate is computed from ECG R peaks (`hrDS`, bpm). The
  `plot_cellmap_and_traces` figure gets four panels under the cell traces on
  the same time axis: running (encoder ticks/s), heart rate, laser and pupil
  diameter.
- `plot_cellmap_and_traces.m`: plots the cell map with colour-coded,
  labelled ROIs (left) next to the stacked traces in matching colours (right).
- `Extract_HekaData_GRIN.m`: loads the HEKA sweeps, concatenates
  channels, binarises digital channels and crops to the period where both
  the GRIN and camera sync are high (`extractedMatGRINcrop`).
- `pupilCircleFit.m`: circle fit to the DeepLabCut pupil points per
  video frame, with blink / poor-fit / outlier removal.
- `alignPupilToHeka.m`: the Cam Synch period (first rise to last fall)
  covers the whole pupil CSV; video frames are spread evenly over it and each
  HEKA sample takes the diameter of the frame it falls in.
