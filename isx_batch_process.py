#!/usr/bin/env python
"""
Batch-process Inscopix .isxd movies with the Inscopix Data Processing (IDPS) 1.9.x
Python API (`isx`). The steps and default parameters match the IDPS 1.9.x GUI:

    raw .isxd
      -> Preprocess (spatial/temporal downsample, defective-pixel fix)    PP
      -> Spatial band-pass filter                                          BP
      -> Mean projection (motion-correction reference)                     MEAN
      -> Motion correction                                                 MC
      -> Crop to the motion-correction crop rectangle                      CROP
      -> DF/F                                                              DFF
      -> Cell identification: PCA-ICA (on DF/F) or CNMFe (on MC movie)     PCA-ICA / CNMFE
      -> Event detection                                                   ED
      -> Auto accept/reject                                                (in place)
      -> Export, named for the MATLAB scripts in this repository:
           <name>_celltrace.csv, <name>_celltrace-props.csv, <name>_cellmap.tif,
           <name>_cellimage_C00.tif ..., <name>_events.csv

Two ways to handle several input files:
  --mode independent  (default) Each file is its own session; it is processed on its own
                      and gets its own outputs.
  --mode series       All files are one recording split into segments (same FOV), the way
                      IDPS handles a "series": each isx call gets the whole list, so motion
                      correction uses one reference and cell extraction yields one cell set
                      covering every segment.

Backends
  isx    (preferred) the Inscopix IDPS 1.9.x Python API: results identical to the GUI.
  numpy  isxd_numpy.py (next to this script): reads/writes the binary .isxd format
         directly with NumPy/SciPy; same steps, re-implemented, so close but not identical
         to IDPS. No CNMFe. Needs: pip install numpy scipy tifffile
  --backend auto (default) uses isx when it is installed, otherwise numpy.

Requirements for the isx backend
  The `isx` package that ships with IDPS 1.9.x (in <IDPS install>/Inscopix Data Processing.app/
  Contents/API/Python on macOS, or <IDPS install>\\api\\python on Windows); install it as
  described in its README, usually in a conda environment built from its environment.yml.
  The open-source `isx` 2.x on PyPI can read .isxd files but does NOT include the processing
  algorithms; with it installed, --backend auto falls back to numpy.

Examples
  python isx_batch_process.py data/*.isxd -o processed
  python isx_batch_process.py data/ -o processed --spatial-ds 4 --method cnmfe
  python isx_batch_process.py seg1.isxd seg2.isxd seg3.isxd -o out --mode series
  python isx_batch_process.py data/ -o out --steps pp bp mc dff     # stop after DF/F
"""

import argparse
import hashlib
import json
import glob
import logging
import os
import sys
import time
from dataclasses import dataclass, field, asdict
from typing import List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import isxd_numpy  # noqa: E402  (cell map writer, used with both backends)
import cell_quality  # noqa: E402  (duplicate ROI rejection on the exports)

isx = None  # set in main(): the IDPS `isx` API, or the isxd_numpy fallback
BACKEND = ""  # "isx" or "numpy"


def load_backend(name):
    """'isx' = Inscopix IDPS 1.9.x API; 'numpy' = isxd_numpy.py (reads the binary .isxd
    directly, no Inscopix software needed); 'auto' = isx if installed, else numpy."""
    if name in ("auto", "isx"):
        try:
            import isx as mod
            if hasattr(mod, "preprocess"):
                return mod, "isx"
            if name == "isx":
                sys.exit("The installed `isx` has no processing functions (it is the PyPI "
                         "reader); use the API shipped with IDPS 1.9.x, or --backend numpy.")
        except ImportError:
            if name == "isx":
                sys.exit("Could not import `isx` from IDPS 1.9.x; try --backend numpy.")
    import isxd_numpy as mod
    return mod, "numpy"

log = logging.getLogger("isx_batch")

ALL_STEPS = ["pp", "bp", "mc", "dff", "cells", "events", "accept", "export"]


# ----------------------------------------------------------------------------------------
# Parameters (defaults = IDPS 1.9.x GUI defaults unless noted)
# ----------------------------------------------------------------------------------------
@dataclass
class Params:
    # Preprocess
    temporal_downsample_factor: int = 1
    spatial_downsample_factor: int = 2          # GUI default is 1; 2 is Inscopix's recommendation for nVista
    crop_rect: Optional[Tuple[int, int, int, int]] = None   # (top, left, bottom, right)
    fix_defective_pixels: bool = True
    trim_early_frames: bool = True

    # Spatial filter
    sf_low_cutoff: float = 0.005
    sf_high_cutoff: float = 0.500
    sf_retain_mean: bool = False
    sf_subtract_global_minimum: bool = True

    # Motion correction
    mc_max_translation: int = 20
    mc_low_bandpass_cutoff: Optional[float] = None   # None: input is already spatially filtered
    mc_high_bandpass_cutoff: Optional[float] = None
    mc_use_mean_projection_reference: bool = True    # else uses reference_frame_index
    mc_reference_file: str = ""                      # custom reference image (.isxd); overrides both
    mc_reference_segment_index: int = 0
    mc_reference_frame_index: int = 0
    mc_global_registration_weight: float = 1.0
    mc_crop: bool = True        # crop MC movies to the area valid in every frame (crop rect)

    # DF/F
    f0_type: str = "mean"                        # 'mean' or 'min'

    # Cell identification
    method: str = "pca-ica"                      # 'pca-ica' or 'cnmfe'
    # PCA-ICA
    num_ics: int = 120
    num_pcs: Optional[int] = None                # None -> round(1.2 * num_ics), as in the GUI
    unmix_type: str = "spatial"
    ica_temporal_weight: float = 0.0
    max_iterations: int = 100
    convergence_threshold: float = 1e-5
    block_size: int = 1000
    auto_estimate_num_ics: bool = False
    average_cell_diameter: int = 10
    # CNMFe
    cell_diameter: int = 7
    min_corr: float = 0.8
    min_pnr: float = 10
    bg_spatial_subsampling: int = 2
    ring_size_factor: float = 1.4
    gaussian_kernel_size: int = 0
    closing_kernel_size: int = 0
    merge_threshold: float = 0.7
    processing_mode: str = "parallel_patches"
    num_threads: int = 4
    patch_size: int = 80
    patch_overlap: int = 20
    output_unit_type: str = "df_over_noise"

    # Event detection
    ed_threshold: float = 5.0                    # in median absolute deviations
    ed_tau: float = 0.2                          # seconds
    ed_event_time_ref: str = "beginning"
    ed_ignore_negative_transients: bool = True
    ed_accepted_cells_only: bool = False

    # Auto accept/reject
    accept_filters: List[Tuple[str, str, float]] = field(
        default_factory=lambda: [("# Comps", "=", 1), ("Cell Size", ">", 0), ("SNR", ">", 3)]
    )

    # Duplicate ROIs (applied to the exported CSVs, see cell_quality.py)
    reject_duplicates: bool = True
    dup_distance_frac: float = 0.25     # x average cell diameter
    dup_min_corr: float = 0.5           # traces less correlated than this = different cells

    # Export
    export_time_ref: str = "start"              # 'start', 'unix' or 'tsc'
    export_cell_images: bool = True             # one .tif per cell
    export_cell_map: bool = True                # <name>_cellmap.tif with white outlines
    export_events_sparse: bool = True           # events CSV: one row per event
    export_events_amplitude: bool = True        # events CSV: amplitude (else 1)


def params_to_json(P: "Params", path: str) -> None:
    import json
    with open(path, "w") as f:
        json.dump(asdict(P), f, indent=2)


def params_from_json(path: str) -> "Params":
    """Load parameters saved by params_to_json / the GUI. Unknown keys are ignored,
    missing keys keep their defaults."""
    import json
    with open(path) as f:
        d = json.load(f)
    P = Params()
    for k, v in d.items():
        if not hasattr(P, k):
            log.warning("Unknown parameter in %s ignored: %s", path, k)
            continue
        if k == "crop_rect" and v is not None:
            v = tuple(int(x) for x in v)
        if k == "accept_filters":
            v = [(str(m), str(o), float(x)) for m, o, x in v]
        setattr(P, k, v)
    return P


# ----------------------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------------------
def collect_inputs(paths: List[str], recursive: bool) -> List[str]:
    """Expand files, directories and glob patterns into a sorted list of .isxd files."""
    found = []
    for p in paths:
        if os.path.isdir(p):
            pattern = os.path.join(p, "**", "*.isxd") if recursive else os.path.join(p, "*.isxd")
            found.extend(glob.glob(pattern, recursive=recursive))
        elif any(c in p for c in "*?["):
            found.extend(glob.glob(p, recursive=recursive))
        elif os.path.isfile(p):
            found.append(p)
        else:
            log.warning("Not found, skipped: %s", p)
    files = sorted({os.path.abspath(f) for f in found if f.lower().endswith(".isxd")})
    # Don't re-process our own outputs if the output dir sits inside an input dir.
    suffixes = ("-PP.isxd", "-BP.isxd", "-MC.isxd", "-CROP.isxd", "-DFF.isxd", "-MEAN.isxd",
                "-PCA-ICA.isxd", "-CNMFE.isxd", "-ED.isxd")
    return [f for f in files if not f.endswith(suffixes)]


def out_paths(inputs: List[str], out_dir: str, suffix: str, ext: str = "isxd") -> List[str]:
    """<out_dir>/<input basename>-<suffix>.<ext>, like IDPS naming."""
    res = []
    for f in inputs:
        base = os.path.splitext(os.path.basename(f))[0]
        # Strip a previous step suffix so names stay short: rec-PP-BP -> rec-BP
        for s in ("-PP", "-BP", "-MC", "-CROP", "-DFF"):
            if base.endswith(s):
                base = base[: -len(s)]
                break
        res.append(os.path.join(out_dir, f"{base}-{suffix}.{ext}"))
    return res


def _marker(path: str) -> str:
    return path + ".incomplete"


def read_crop_rect(path):
    """Last line of the motion-correction crop-rect CSV holding 4 integers."""
    rect = None
    with open(path) as f:
        for line in f:
            parts = [x.strip() for x in line.replace(";", ",").split(",") if x.strip()]
            if len(parts) == 4 and all(x.lstrip("-").isdigit() for x in parts):
                rect = tuple(int(x) for x in parts)
    if rect is None:
        raise ValueError(f"Could not read a crop rectangle from {path}")
    return rect


# ----------------------------------------------------------------------------------------
# Step bookkeeping: rerun only what changed
# ----------------------------------------------------------------------------------------
# Parameters each step depends on. A step is redone when these, or the result of a step it
# builds on, change; otherwise its existing output is reused.
STEP_PARAMS = {
    "pp": ["temporal_downsample_factor", "spatial_downsample_factor", "crop_rect",
           "fix_defective_pixels", "trim_early_frames"],
    "bp": ["sf_low_cutoff", "sf_high_cutoff", "sf_retain_mean", "sf_subtract_global_minimum"],
    "mean": [],
    "mc": ["mc_max_translation", "mc_low_bandpass_cutoff", "mc_high_bandpass_cutoff",
           "mc_use_mean_projection_reference", "mc_reference_file",
           "mc_reference_segment_index", "mc_reference_frame_index",
           "mc_global_registration_weight"],
    "crop": ["mc_crop"],
    "dff": ["f0_type"],
    "pca-ica": ["num_ics", "num_pcs", "unmix_type", "ica_temporal_weight", "max_iterations",
                "convergence_threshold", "block_size", "auto_estimate_num_ics",
                "average_cell_diameter"],
    "cnmfe": ["cell_diameter", "min_corr", "min_pnr", "bg_spatial_subsampling",
              "ring_size_factor", "gaussian_kernel_size", "closing_kernel_size",
              "merge_threshold", "processing_mode", "num_threads", "patch_size",
              "patch_overlap", "output_unit_type"],
    "events": ["ed_threshold", "ed_tau", "ed_event_time_ref", "ed_ignore_negative_transients",
               "ed_accepted_cells_only"],
    "accept": ["accept_filters"],
    "export": ["export_time_ref", "export_cell_images", "export_cell_map",
               "export_events_sparse", "export_events_amplitude",
               "reject_duplicates", "dup_distance_frac", "dup_min_corr"],
}


def _file_sig(path: str):
    """Identifies an input file without reading it: name and size."""
    try:
        return [os.path.basename(path), os.path.getsize(path)]
    except OSError:
        return [os.path.basename(path), None]


class StepRecord:
    """Remembers, per output folder, which parameters and upstream results produced each
    step's outputs (<name>_steps.json), and decides what to run, reuse or skip."""

    def __init__(self, out_dir, name, P, overwrite, root):
        self.path = os.path.join(out_dir, f"{name}_steps.json")
        self.P, self.overwrite = P, overwrite
        self.keys = {"raw": self._hash(root)}
        self.ran = set()
        try:
            with open(self.path) as f:
                self.rec = json.load(f)
        except (OSError, ValueError):
            self.rec = {}

    @staticmethod
    def _hash(obj):
        return hashlib.sha1(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]

    def _params(self, step, extra=None):
        d = {k: getattr(self.P, k) for k in STEP_PARAMS.get(step, [])}
        if "accept_filters" in d:          # 1 and 1.0 are the same filter
            d["accept_filters"] = [[m, o, float(v)] for m, o, v in d["accept_filters"]]
        if step == "mc" and self.P.mc_reference_file:
            d["reference_sig"] = _file_sig(self.P.mc_reference_file)
        d.update(extra or {})
        return json.loads(json.dumps(d, default=str))   # same form as read back from file

    def _save(self):
        tmp = self.path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(self.rec, f, indent=1, default=str)
        os.replace(tmp, self.path)

    def step(self, step, label, outputs, deps, enabled, run, extra_params=None, side=()):
        """Run, reuse or skip one step. Returns True when its outputs are available
        (computed now or earlier), False when it is skipped (no outputs)."""
        params = self._params(step, extra_params)
        key = self._hash({"step": step, "params": params,
                          "deps": [self.keys[d] for d in deps]})
        for o in outputs:                    # leftovers of an interrupted run
            if os.path.exists(_marker(o)):
                log.info("  removing incomplete output from an interrupted run: %s",
                         os.path.basename(o))
                for f in (o, _marker(o)):
                    if os.path.exists(f):
                        os.remove(f)
                self.rec.pop(step, None)
        exist = all(os.path.exists(o) for o in outputs)
        rec = self.rec.get(step)
        names = ", ".join(os.path.basename(o) for o in outputs) or "in place"

        if not enabled:
            if exist and (outputs or rec):
                stale = rec is None or rec.get("key") != key
                self.keys[step] = rec["key"] if rec else key
                log.info("[%s] not selected; using existing %s%s", label, names,
                         " (made with other settings; tick the step to update it)"
                         if stale else "")
                return True
            self.keys[step] = self.keys[deps[0]] if deps else key
            return False

        reason = ("overwrite" if self.overwrite else
                  "no output yet" if not exist else
                  "no record of the settings used" if rec is None else
                  "earlier step redone" if any(d in self.ran for d in deps) else
                  "settings changed" if rec.get("key") != key else None)
        self.keys[step] = key
        if reason is None:
            log.info("[%s] unchanged, reusing %s", label, names)
            return True
        log.info("[%s] -> %s (%s)", label, names, reason)
        if rec and rec.get("params") != params and rec.get("params") is not None:
            changed = sorted(k for k in set(params) | set(rec["params"])
                             if params.get(k) != rec["params"].get(k))
            if changed:
                log.info("  changed: %s", ", ".join(changed))
        for f in list(outputs) + list(side):
            if os.path.exists(f):
                os.remove(f)
        for o in outputs:
            open(_marker(o), "w").close()
        t0 = time.time()
        run()
        for o in outputs:
            if os.path.exists(_marker(o)):
                os.remove(_marker(o))
        self.rec[step] = {"key": key, "params": params,
                          "outputs": [os.path.basename(o) for o in outputs],
                          "finished": time.strftime("%Y-%m-%d %H:%M:%S")}
        self._save()
        self.ran.add(step)
        log.info("  done in %.1f s", time.time() - t0)
        return True


# ----------------------------------------------------------------------------------------
# Pipeline for one group of files (one session, or one series)
# ----------------------------------------------------------------------------------------
def process_group(inputs: List[str], out_dir: str, name: str, P: Params,
                  steps: List[str], overwrite: bool) -> None:
    os.makedirs(out_dir, exist_ok=True)
    log.info("=== %s: %d file(s) -> %s", name, len(inputs), out_dir)
    root = {"inputs": [_file_sig(f) for f in inputs], "backend": BACKEND,
            "backend_version": getattr(isx, "READER_VERSION", getattr(isx, "__version__", ""))}
    R = StepRecord(out_dir, name, P, overwrite, root)
    current, last = inputs, "raw"

    # 1. Preprocess
    pp = out_paths(current, out_dir, "PP")
    if R.step("pp", "preprocess", pp, [last], "pp" in steps, lambda: isx.preprocess(
            input_movie_files=inputs, output_movie_files=pp,
            temporal_downsample_factor=P.temporal_downsample_factor,
            spatial_downsample_factor=P.spatial_downsample_factor,
            crop_rect=P.crop_rect, fix_defective_pixels=P.fix_defective_pixels,
            trim_early_frames=P.trim_early_frames)):
        current, last = pp, "pp"

    # 2. Spatial band-pass filter
    src = current
    bp = out_paths(src, out_dir, "BP")
    filtered = R.step("bp", "spatial filter", bp, [last], "bp" in steps, lambda: isx.spatial_filter(
            input_movie_files=src, output_movie_files=bp,
            low_cutoff=P.sf_low_cutoff, high_cutoff=P.sf_high_cutoff,
            retain_mean=P.sf_retain_mean, subtract_global_minimum=P.sf_subtract_global_minimum))
    if filtered:
        current, last = bp, "bp"
    # Without a spatially filtered movie, motion correction band-passes internally
    # (IDPS default cutoffs) unless cutoffs were given.
    mc_lo, mc_hi = P.mc_low_bandpass_cutoff, P.mc_high_bandpass_cutoff
    if not filtered and mc_lo is None and mc_hi is None:
        mc_lo, mc_hi = 0.004, 0.016

    # 3. Motion correction (reference = mean projection of the filtered movie(s)), then crop
    src = current
    mc = out_paths(src, out_dir, "MC")
    translations = out_paths(src, out_dir, "translations", "csv")
    crop_file = os.path.join(out_dir, f"{name}-crop_rect.csv")
    ref, mc_deps = P.mc_reference_file or "", [last]
    if not ref and P.mc_use_mean_projection_reference:
        ref = os.path.join(out_dir, f"{name}-MEAN.isxd")
        R.step("mean", "mean projection", [ref], [last], "mc" in steps,
               lambda: isx.project_movie(input_movie_files=src, output_image_file=ref,
                                         stat_type="mean"))
        mc_deps = [last, "mean"]
    if R.step("mc", "motion correct", mc, mc_deps, "mc" in steps, lambda: isx.motion_correct(
            input_movie_files=src, output_movie_files=mc,
            max_translation=P.mc_max_translation,
            low_bandpass_cutoff=mc_lo, high_bandpass_cutoff=mc_hi,
            reference_segment_index=P.mc_reference_segment_index,
            reference_frame_index=P.mc_reference_frame_index,
            reference_file_name=ref,
            global_registration_weight=P.mc_global_registration_weight,
            output_translation_files=translations, output_crop_rect_file=crop_file),
            side=translations + [crop_file],
            extra_params={"bandpass_used": [mc_lo, mc_hi]}):
        current, last = mc, "mc"
        if P.mc_crop and os.path.exists(crop_file):
            src = current
            cropped = out_paths(src, out_dir, "CROP")
            rect = read_crop_rect(crop_file)
            if R.step("crop", "crop", cropped, [last], "mc" in steps,
                      lambda: (log.info("  crop rect (top, left, bottom, right) = %s", rect),
                               isx.preprocess(input_movie_files=src, output_movie_files=cropped,
                                              crop_rect=rect, fix_defective_pixels=False,
                                              trim_early_frames=False))):
                current, last = cropped, "crop"
    mc_movies, mc_last = current, last

    # 4. DF/F
    src = current
    dff = out_paths(src, out_dir, "DFF")
    if R.step("dff", "dF/F", dff, [last], "dff" in steps,
              lambda: isx.dff(input_movie_files=src, output_movie_files=dff, f0_type=P.f0_type),
              extra_params={"impl": getattr(isx, "DFF_VERSION", "")}):
        current, last = dff, "dff"

    # 5. Cell identification -> one cell set per input movie (series: one per segment)
    if P.method == "pca-ica":
        src, cell_dep = current, last
        cellsets = out_paths(src, out_dir, "PCA-ICA")
        num_pcs = P.num_pcs or int(round(1.2 * P.num_ics))
        run_cells = lambda: isx.pca_ica(  # noqa: E731
            input_movie_files=src, output_cell_set_files=cellsets,
            num_pcs=num_pcs, num_ics=P.num_ics, unmix_type=P.unmix_type,
            ica_temporal_weight=P.ica_temporal_weight, max_iterations=P.max_iterations,
            convergence_threshold=P.convergence_threshold, block_size=P.block_size,
            auto_estimate_num_ics=P.auto_estimate_num_ics,
            average_cell_diameter=P.average_cell_diameter)
        label = "PCA-ICA"
    elif P.method == "cnmfe":
        # CNMFe models the background itself, so IDPS runs it on the motion-corrected
        # movie, not on DF/F.
        src, cell_dep = mc_movies, mc_last
        cellsets = out_paths(src, out_dir, "CNMFE")
        run_cells = lambda: isx.cnmfe(  # noqa: E731
            input_movie_files=src, output_cell_set_files=cellsets, output_dir=out_dir,
            cell_diameter=P.cell_diameter, min_corr=P.min_corr, min_pnr=P.min_pnr,
            bg_spatial_subsampling=P.bg_spatial_subsampling,
            ring_size_factor=P.ring_size_factor, gaussian_kernel_size=P.gaussian_kernel_size,
            closing_kernel_size=P.closing_kernel_size, merge_threshold=P.merge_threshold,
            processing_mode=P.processing_mode, num_threads=P.num_threads,
            patch_size=P.patch_size, patch_overlap=P.patch_overlap,
            output_unit_type=P.output_unit_type)
        label = "CNMFe"
    else:
        raise ValueError(f"Unknown cell identification method: {P.method}")
    impl = {"impl": getattr(isx, "PCA_ICA_VERSION", "")} if P.method == "pca-ica" else None
    if not R.step(P.method, label, cellsets, [cell_dep], "cells" in steps, run_cells,
                  extra_params=impl):
        if any(s in steps for s in ("events", "accept", "export")):
            log.warning("No cell sets available; skipping events/accept/export.")
        return

    # 6. Event detection
    events = out_paths(cellsets, out_dir, "ED")
    have_events = R.step("events", "event detection", events, [P.method], "events" in steps,
                         lambda: isx.event_detection(
                             input_cell_set_files=cellsets, output_event_set_files=events,
                             threshold=P.ed_threshold, tau=P.ed_tau,
                             event_time_ref=P.ed_event_time_ref,
                             ignore_negative_transients=P.ed_ignore_negative_transients,
                             accepted_cells_only=P.ed_accepted_cells_only))

    # 7. Auto accept/reject (changes the cell statuses inside the cell sets)
    if have_events:
        R.step("accept", "auto accept/reject", [], [P.method, "events"], "accept" in steps,
               lambda: (log.info("  filters=%s", P.accept_filters),
                        isx.auto_accept_reject(cellsets, events, filters=P.accept_filters)))
    elif "accept" in steps:
        log.warning("Auto accept/reject needs event sets; skipped.")

    # 8. Export, named so A_Parent_Align_Heka_Pupil_Inscopix.m / plot_cellmap_and_traces.m
    #    find them: <name>_celltrace.csv, <name>_celltrace-props.csv, <name>_cellmap.tif.
    #    With several cell sets (series) the CSVs cover all segments.
    traces_csv = os.path.join(out_dir, f"{name}_celltrace.csv")
    props_csv = os.path.join(out_dir, f"{name}_celltrace-props.csv")
    cellmap_tif = os.path.join(out_dir, f"{name}_cellmap.tif")
    events_csv = os.path.join(out_dir, f"{name}_events.csv")
    # per-cell images: <name>_cellimage_C00.tif ... (must not end in "cellmap.tif")
    tiff_name = os.path.join(out_dir, f"{name}_cellimage.tif")
    outputs = [traces_csv, props_csv] + ([cellmap_tif] if P.export_cell_map else []) \
        + ([events_csv] if have_events else [])
    old_images = glob.glob(os.path.join(glob.escape(out_dir), f"{glob.escape(name)}_cellimage_*.tif"))

    def export():
        isx.export_cell_set_to_csv_tiff(cellsets, traces_csv,
                                        tiff_name if P.export_cell_images else "",
                                        time_ref=P.export_time_ref, output_props_file=props_csv)
        if P.export_cell_map:
            scale = isxd_numpy.write_cell_map(cellmap_tif, read_cell_images(cellsets[0]))
            log.info("  cell map enlarged %dx", scale)
        if have_events:
            isx.export_event_set_to_csv(events, events_csv, time_ref=P.export_time_ref,
                                        sparse_output=P.export_events_sparse,
                                        write_amplitude=P.export_events_amplitude)
        if P.reject_duplicates:
            pairs = cell_quality.reject_duplicates_in_exports(
                props_csv, traces_csv, P.dup_distance_frac, P.dup_min_corr)
            log.info("  duplicate ROIs rejected: %s", ", ".join(f"{a} (dup. of {b})"
                                                          for a, b in pairs) or "none")
    R.step("export", "export", outputs,
           [P.method] + (["events", "accept"] if have_events else []),
           "export" in steps, export, side=old_images)


def read_cell_images(cellset_path: str):
    """Cell images of a cell set, through the isx API when it is the backend (reads any
    IDPS cell set), otherwise through isxd_numpy."""
    if hasattr(isx, "CellSet"):
        cs = isx.CellSet.read(cellset_path)
        try:
            return [cs.get_cell_image_data(i) for i in range(cs.num_cells)]
        finally:
            cs.close()
    cs = isxd_numpy.IsxdCellSet(cellset_path)
    return [cs.get_image(i) for i in range(cs.num_cells)]


# ----------------------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Batch-process Inscopix .isxd movies like IDPS 1.9.x.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("inputs", nargs="+", help=".isxd files, directories or glob patterns")
    p.add_argument("-o", "--out-dir", required=True, help="output directory")
    p.add_argument("--mode", choices=["independent", "series"], default="independent",
                   help="independent: each file is its own session; series: all files are "
                        "segments of one recording")
    p.add_argument("--series-name", default="series", help="name used for series outputs")
    p.add_argument("-r", "--recursive", action="store_true", help="search directories recursively")
    p.add_argument("--steps", nargs="+", choices=ALL_STEPS, default=ALL_STEPS,
                   help="steps to run")
    p.add_argument("--overwrite", action="store_true", help="redo steps whose outputs exist")
    p.add_argument("--config", help="JSON file with all parameters (as saved by the GUI); "
                                    "when given, the individual parameter options are ignored")
    p.add_argument("--backend", choices=["auto", "isx", "numpy"], default="auto",
                   help="processing engine (see load_backend)")
    p.add_argument("--continue-on-error", action="store_true",
                   help="in independent mode, keep going when a file fails")

    d = Params()
    g = p.add_argument_group("parameters")
    g.add_argument("--spatial-ds", type=int, default=d.spatial_downsample_factor)
    g.add_argument("--temporal-ds", type=int, default=d.temporal_downsample_factor)
    g.add_argument("--crop", type=int, nargs=4, metavar=("TOP", "LEFT", "BOTTOM", "RIGHT"))
    g.add_argument("--no-fix-defective-pixels", action="store_true")
    g.add_argument("--sf-low", type=float, default=d.sf_low_cutoff)
    g.add_argument("--sf-high", type=float, default=d.sf_high_cutoff)
    g.add_argument("--max-translation", type=int, default=d.mc_max_translation)
    g.add_argument("--mc-frame-reference", action="store_true",
                   help="use the first frame instead of the mean projection as MC reference")
    g.add_argument("--no-mc-crop", action="store_true",
                   help="keep the motion-corrected movie uncropped (padded edges)")
    g.add_argument("--f0-type", choices=["mean", "min"], default=d.f0_type)
    g.add_argument("--method", choices=["pca-ica", "cnmfe"], default=d.method)
    g.add_argument("--num-ics", type=int, default=d.num_ics)
    g.add_argument("--num-pcs", type=int, default=None, help="default: 1.2 x num-ics")
    g.add_argument("--auto-estimate-ics", action="store_true")
    g.add_argument("--cell-diameter", type=int, default=d.cell_diameter, help="CNMFe")
    g.add_argument("--min-corr", type=float, default=d.min_corr, help="CNMFe")
    g.add_argument("--min-pnr", type=float, default=d.min_pnr, help="CNMFe")
    g.add_argument("--num-threads", type=int, default=d.num_threads, help="CNMFe")
    g.add_argument("--ed-threshold", type=float, default=d.ed_threshold)
    g.add_argument("--ed-tau", type=float, default=d.ed_tau)
    g.add_argument("--min-snr", type=float, default=3.0, help="auto accept/reject SNR cutoff")
    g.add_argument("--time-ref", choices=["start", "unix", "tsc"], default=d.export_time_ref)
    return p.parse_args(argv)


def build_params(a) -> Params:
    if a.config:
        return params_from_json(a.config)
    P = Params()
    P.spatial_downsample_factor = a.spatial_ds
    P.temporal_downsample_factor = a.temporal_ds
    P.crop_rect = tuple(a.crop) if a.crop else None
    P.fix_defective_pixels = not a.no_fix_defective_pixels
    P.sf_low_cutoff, P.sf_high_cutoff = a.sf_low, a.sf_high
    P.mc_max_translation = a.max_translation
    P.mc_use_mean_projection_reference = not a.mc_frame_reference
    P.mc_crop = not a.no_mc_crop
    P.f0_type = a.f0_type
    P.method = a.method
    P.num_ics, P.num_pcs = a.num_ics, a.num_pcs
    P.auto_estimate_num_ics = a.auto_estimate_ics
    P.cell_diameter, P.min_corr, P.min_pnr = a.cell_diameter, a.min_corr, a.min_pnr
    P.num_threads = a.num_threads
    P.ed_threshold, P.ed_tau = a.ed_threshold, a.ed_tau
    P.accept_filters = [("# Comps", "=", 1), ("Cell Size", ">", 0), ("SNR", ">", a.min_snr)]
    P.export_time_ref = a.time_ref
    return P


def main(argv=None):
    a = parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.StreamHandler(),
                  logging.FileHandler(os.path.join(a.out_dir, "isx_batch_process.log"))])

    files = collect_inputs(a.inputs, a.recursive)
    if not files:
        sys.exit("No .isxd files found.")
    P = build_params(a)
    steps = [s for s in ALL_STEPS if s in a.steps]
    global isx, BACKEND
    isx, backend = load_backend(a.backend)
    BACKEND = backend
    log.info("backend: %s (version %s)", backend, getattr(isx, "__version__", "unknown"))
    log.info("%d input file(s), mode=%s, steps=%s", len(files), a.mode, steps)
    log.info("parameters: %s", asdict(P))

    failed = []
    if a.mode == "series":
        process_group(files, a.out_dir, a.series_name, P, steps, a.overwrite)
    else:
        for i, f in enumerate(files, 1):
            name = os.path.splitext(os.path.basename(f))[0]
            log.info("---- [%d/%d] %s", i, len(files), name)
            try:
                process_group([f], os.path.join(a.out_dir, name), name, P, steps, a.overwrite)
            except Exception:
                log.exception("FAILED: %s", f)
                failed.append(f)
                if not a.continue_on_error:
                    break

    if failed:
        log.error("%d file(s) failed: %s", len(failed), failed)
        sys.exit(1)
    log.info("All done.")


if __name__ == "__main__":
    main()
