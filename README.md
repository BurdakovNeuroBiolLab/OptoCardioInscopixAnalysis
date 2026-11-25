# OptoCardioInscopixAnalysis

MATLAB analysis code for opto / cardio / Inscopix experiments, plus Python
batch processing of raw Inscopix `.isxd` movies.

## Files

All scripts are in the repository root. Run `A_Parent_Align_Heka_Pupil_Inscopix`;
it calls the others. The Python files (`isx_*.py`, `isxd_numpy.py`) are
described under [Inscopix batch processing (Python)](#inscopix-batch-processing-python).

| File | Contents |
|---|---|
| `A_Parent_Align_Heka_Pupil_Inscopix.m` | Parent script: pupil and HEKA data aligned to the Inscopix cell-trace frames |
| `plot_cellmap_and_traces.m` | Inscopix (IDPS) cell map and stacked traces |
| `Extract_HekaData_GRIN.m` | HEKA data extraction |
| `pupilCircleFit.m` | DeepLabCut pupil circle fit |
| `alignPupilToHeka.m` | Pupil diameter onto the HEKA time base |
| `continuous_hrv_mouse.m` | Windowed HR / HRV from beat times (Nik and Mattia, with fixes) |
| `dfa_alpha1_local.m` | DFA alpha1 of an IBI series |
| `flag_bad_ibi.m` | Marks IBIs from missed / extra beats |
| `Nik_and_Mattia_HeartRate_code_vs_Claude.md` | Comparison of the two heart-rate methods |
| `isx_batch_gui.py` | Python GUI: select `.isxd` files, set parameters, batch process |
| `isx_batch_process.py` | Python batch pipeline (command line), IDPS 1.9.x steps |
| `isxd_numpy.py` | Reads/writes `.isxd` and re-implements the IDPS steps in NumPy/SciPy |
| `isx_viewer.py` | Python viewer for processed results: cell map, traces, cell images |

## Usage

Add the repository folder to the MATLAB path, set the current folder to the
recording folder and run `A_Parent_Align_Heka_Pupil_Inscopix`
by name. Each recording folder holds exactly one of each:

- `*celltrace.csv`: IDPS cell traces
- `*celltrace-props.csv`: IDPS cell props
- `*cellmap.tif` / `*cellmap.tiff`: cell map with white ROI outlines
- `*L.mat`: HEKA export
- `*zoompup*.csv`: DeepLabCut pupil tracking

Files starting with `._` (macOS resource files) are ignored by every script.

## Scripts

- `A_Parent_Align_Heka_Pupil_Inscopix.m`: runs the pupil fit, the cell-trace plot and
  the HEKA extraction, puts the pupil diameter on the HEKA time base
  (`alignPupilToHeka`, column 8 of `extractedMatGRINcrop`), then bins
  `extractedMatGRINcrop` into one row per cell-trace frame (`hekaDS`, saved to
  `<traces>_hekaAligned.mat`). The GRIN synch high period is taken as the full
  Inscopix recording: its rising edge is the first frame, its falling edge the
  end of the last frame. Analogue channels are averaged per frame, run edge
  events are summed (encoder ticks per frame). Frames outside the GRIN/camera
  crop are NaN. Heart rate is computed from ECG R peaks refined to
  sub-sample timing (`hrDS`, bpm), after removing missed / extra beats; HRV
  (HR, RMSSD, CVSD, CV_IBI, residual) comes from `continuous_hrv_mouse` in
  10 s windows (`hrvDS`), and DFA alpha1 over the whole crop (`alpha1Whole`).
  The `plot_cellmap_and_traces` figure gets five panels under the cell
  traces on the same time axis: running (encoder ticks/s), heart rate,
  RMSSD, laser and pupil diameter. Set `nCellsToShow` at the top to a range of cells (1 = C00, e.g.
  `1:10` or `5:20`; [] = all). The cells outside the range go to a second
  figure, "Unlikely to be cells", with their traces and outlines in colour
  and the selected cells' outlines in black.
  The main figure is saved as `<traces>_cellmap_traces.tif` (300 dpi) and
  `.svg` in the current folder (`saveFigure`).
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

## Inscopix batch processing (Python)

Processes raw `.isxd` movies with the steps of Inscopix Data Processing (IDPS)
1.9.x, for any number of files:

preprocess → spatial filter → motion correction (mean-image reference) → crop
to the valid area → DF/F → PCA-ICA or CNMFe → event detection → auto
accept/reject → export.

The exports are named so that `A_Parent_Align_Heka_Pupil_Inscopix` /
`plot_cellmap_and_traces` find them: copy them into the recording folder
next to the HEKA and pupil files.

| Output | Contents |
|---|---|
| `<name>_celltrace.csv` | Cell traces (IDPS format: names row, status row, time + one column per cell) |
| `<name>_celltrace-props.csv` | Cell props (Name, Status, ColorR/G/B, CentroidX/Y, NumComponents, Size) |
| `<name>_cellmap.tif` | Cell map with the ROI outlines in white |
| `<name>_cellimage_C00.tif` … | One image per cell |
| `<name>_events.csv` | Detected events |

`<name>` is the movie's file name (independent mode) or the series name.

Defaults are the IDPS defaults, except spatial downsampling (2).

**GUI**

```
python isx_batch_gui.py
```

In Spyder, open `isx_batch_gui.py` (or `isx_viewer.py`) and press Run (F5):
started from Spyder's console or Jupyter, the window opens as its own process,
so no "Execute in an external system terminal" setting is needed and the
console stays free. If the window cannot start, the error is printed in the
console.

Add files or folders, choose an output folder, set any parameter (every IDPS
1.9.x parameter has a field with an explanation) and press Run. Settings can be
saved/loaded as JSON; each run also saves its parameters as
`params_<date-time>.json` in the output folder. Cancel stops at once; the
interrupted step is redone on the next run.

**Rerunning with changed settings** only redoes what changed. Each output
folder keeps `<name>_steps.json`, recording which settings and earlier results
produced each step's output. A step runs again when its own settings changed,
an earlier step was redone, or its output is missing; otherwise the existing
output is reused. For example, a new number of ICs or cell diameter reruns
cell identification, event detection, accept/reject and export only; a new
spatial downsampling reruns everything. Unticked steps reuse their existing
output (the log warns if it was made with other settings). "Overwrite" redoes
every ticked step.

**Viewing results**

"View results…" in the GUI (or `python isx_viewer.py <folder>`) opens every
`*_celltrace.csv` found in the folder and its subfolders, with its props and
cell map; IDPS exports work too. Left: the cell map with all ROI outlines in
yellow. Right, bottom: the traces, 10 at a time (scroll for more), each with two
check boxes. The left box selects the cell: its trace and outline get their own
colour (several cells at once). The right box shows that cell's image above
(one at a time), or "Cell Images Not Found" if no per-cell TIFFs were saved.

**Command line**

```
python isx_batch_process.py data/*.isxd -o processed
python isx_batch_process.py seg1.isxd seg2.isxd -o out --mode series
python isx_batch_process.py data/ -o out --config params.json
python isx_batch_process.py --help
```

**Several files**

- `independent` (default): each file is its own session, in its own output folder.
- `series`: the files are segments of one recording (same field of view); they
  share the motion-correction reference, F0, the cell set and accept/reject
  decisions, and export to one CSV with real time gaps.

**Backends**

- `isx`: the Python API shipped with IDPS 1.9.x. Results identical to IDPS.
  Run the GUI from that environment, or set its "Python" field to that
  environment's interpreter.
- `numpy`: no Inscopix software needed (`pip install numpy scipy tifffile`).
  Reads the binary `.isxd` format directly and re-implements the steps, so
  results are close to, but not identical with, IDPS. No CNMFe; event sets
  are written in this code's own format (use the CSV export).
- `auto` (default): `isx` if installed, otherwise `numpy`.

Load an `.isxd` file in Python:

```python
import isxd_numpy as ix
mov = ix.IsxdMovie("rec.isxd")       # mov.shape, mov.period
frames = mov.get_frames(0, 100)      # float32 (frames, height, width)
```
