#!/usr/bin/env python
"""
Graphical front end for isx_batch_process.py: pick .isxd files, set every processing
parameter that Inscopix Data Processing 1.9.x offers, and batch-process them.

    python isx_batch_gui.py

Needs only the Python standard library (Tkinter) plus whatever the chosen backend needs.
The processing runs in a separate Python process, so the window stays responsive and
Cancel stops the work immediately; an interrupted step is redone on the next run.

Settings can be saved to / loaded from JSON. The same JSON runs headless with:
    python isx_batch_process.py <files> -o <out> --config settings.json
"""

import json
import os
import queue
import re
import subprocess
import sys
import threading
import time
import tkinter as tk
from dataclasses import asdict, fields
from tkinter import filedialog, messagebox, ttk

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from isx_batch_process import ALL_STEPS, Params, collect_inputs  # noqa: E402

SCRIPT = os.path.join(HERE, "isx_batch_process.py")
LAST_SETTINGS = os.path.join(os.path.expanduser("~"), ".isx_batch_gui.json")

STEP_LABELS = {
    "pp": "Preprocess", "bp": "Spatial filter", "mc": "Motion correction", "dff": "DF/F",
    "cells": "Cell identification", "events": "Event detection",
    "accept": "Auto accept/reject", "export": "Export (CSV/TIFF)",
}

# --------------------------------------------------------------------------------------
# Parameter definitions: (field, label, kind, extra, help)
#   kind: int | float | bool | choice | opt_int | opt_float | file
# --------------------------------------------------------------------------------------
TABS = [
    ("Preprocess", [
        ("spatial_downsample_factor", "Spatial downsample factor", "int", None,
         "Average N x N pixel blocks. IDPS default 1; Inscopix recommends 2-4 for nVista."),
        ("temporal_downsample_factor", "Temporal downsample factor", "int", None,
         "Average N consecutive frames."),
        ("crop_rect", "Crop", "crop", None,
         "Crop before processing (pixels of the raw movie, inclusive)."),
        ("fix_defective_pixels", "Fix defective pixels", "bool", None,
         "Replace hot/dead pixels by the median of their neighbours."),
        ("trim_early_frames", "Trim early frames", "bool", None,
         "Drop frames recorded before the LED reached steady state (isx backend only)."),
    ]),
    ("Spatial filter", [
        ("sf_low_cutoff", "Low cutoff (cycles/pixel)", "float", None,
         "Removes structure larger than ~1/cutoff pixels (background)."),
        ("sf_high_cutoff", "High cutoff (cycles/pixel)", "float", None,
         "Removes structure smaller than ~1/cutoff pixels (noise). 0.5 = Nyquist."),
        ("sf_retain_mean", "Retain mean", "bool", None,
         "Add each frame's mean back after filtering."),
        ("sf_subtract_global_minimum", "Subtract global minimum", "bool", None,
         "Shift the output so its minimum is 0 (keeps DF/F well defined)."),
    ]),
    ("Motion correction", [
        ("mc_max_translation", "Maximum translation (pixels)", "int", None,
         "Largest shift searched for, in each direction."),
        ("mc_low_bandpass_cutoff", "Low band-pass cutoff", "opt_float", None,
         "Blank = no filtering (use when the spatial filter step ran). IDPS: 0.004."),
        ("mc_high_bandpass_cutoff", "High band-pass cutoff", "opt_float", None,
         "Blank = no filtering. IDPS: 0.016."),
        ("mc_use_mean_projection_reference", "Reference: mean image", "bool", None,
         "Register to the mean projection; off = a single reference frame."),
        ("mc_reference_segment_index", "Reference segment index", "int", None,
         "Series only: which file holds the reference frame (0 = first)."),
        ("mc_reference_frame_index", "Reference frame index", "int", None,
         "Frame used as reference when 'mean image' is off."),
        ("mc_reference_file", "Custom reference image", "file", None,
         "Optional .isxd image to register to; overrides the two options above."),
        ("mc_global_registration_weight", "Global registration weight", "float", None,
         "0-1, weight of whole-frame vs. local registration (isx backend only)."),
        ("mc_crop", "Crop to valid area", "bool", None,
         "Crop the corrected movie to the region covered in every frame."),
    ]),
    ("DF/F", [
        ("f0_type", "F0", "choice", ["mean", "min"],
         "Per-pixel baseline F0 used for (F - F0) / F0."),
    ]),
    ("PCA-ICA", [
        ("num_ics", "Number of ICs", "int", None, "Expected number of cells (slightly more)."),
        ("num_pcs", "Number of PCs", "opt_int", None, "Blank = 1.2 x number of ICs (IDPS default)."),
        ("auto_estimate_num_ics", "Auto-estimate number of ICs", "bool", None,
         "Estimate the number of ICs from the data."),
        ("average_cell_diameter", "Average cell diameter (pixels)", "int", None,
         "Used by auto-estimation."),
        ("unmix_type", "Unmix type", "choice", ["spatial", "temporal", "both"],
         "Which dimension ICA unmixes."),
        ("ica_temporal_weight", "ICA temporal weight", "float", None,
         "0-1, used when unmix type is 'both'."),
        ("max_iterations", "Maximum iterations", "int", None, ""),
        ("convergence_threshold", "Convergence threshold", "float", None, ""),
        ("block_size", "Block size", "int", None, "Frames per block (memory use, isx backend)."),
    ]),
    ("CNMFe", [
        ("cell_diameter", "Cell diameter (pixels)", "int", None, "Typical cell size."),
        ("min_corr", "Minimum pixel correlation", "float", None, "Seed threshold."),
        ("min_pnr", "Minimum peak-to-noise ratio", "float", None, "Seed threshold."),
        ("bg_spatial_subsampling", "Background spatial subsampling", "int", None, ""),
        ("ring_size_factor", "Ring size factor", "float", None, "Background ring radius / cell diameter."),
        ("gaussian_kernel_size", "Gaussian kernel size", "int", None, "0 = automatic."),
        ("closing_kernel_size", "Closing kernel size", "int", None, "0 = automatic."),
        ("merge_threshold", "Merge threshold", "float", None,
         "Merge cells whose temporal correlation is above this."),
        ("processing_mode", "Processing mode", "choice",
         ["all_in_memory", "sequential_patches", "parallel_patches"], ""),
        ("num_threads", "Number of threads", "int", None, "For parallel patches."),
        ("patch_size", "Patch size (pixels)", "int", None, ""),
        ("patch_overlap", "Patch overlap (pixels)", "int", None, ""),
        ("output_unit_type", "Output units", "choice", ["df", "df_over_noise"], ""),
    ]),
    ("Event detection", [
        ("ed_threshold", "Threshold (MAD)", "float", None,
         "Minimum event size in median absolute deviations of the trace."),
        ("ed_tau", "Tau (s)", "float", None, "Minimum decay time / spacing between events."),
        ("ed_event_time_ref", "Event time reference", "choice",
         ["beginning", "mid_rise", "maximum"], "Which point of the transient is the event time."),
        ("ed_ignore_negative_transients", "Ignore negative transients", "bool", None, ""),
        ("ed_accepted_cells_only", "Accepted cells only", "bool", None,
         "Only detect events of cells already accepted."),
    ]),
    ("Accept/reject", [
        ("accept_filters", "Filters", "filters", None,
         "A cell is accepted when it passes every filter, otherwise rejected."),
    ]),
    ("Export", [
        ("export_time_ref", "Time reference", "choice", ["start", "unix", "tsc"],
         "start = seconds from recording start; unix = epoch seconds; tsc = hardware counter."),
        ("export_cell_map", "Export cell map (TIFF)", "bool", None,
         "<name>_cellmap.tif: cell images with white outlines, as read by "
         "plot_cellmap_and_traces.m."),
        ("export_cell_images", "Export cell images (TIFF)", "bool", None,
         "One <name>_cellimage_Cxx.tif per cell."),
        ("export_events_sparse", "Sparse event output", "bool", None,
         "On: one row per event. Off: one row per frame, one column per cell."),
        ("export_events_amplitude", "Write event amplitude", "bool", None,
         "Off: events are written as 1."),
    ]),
]
FILTER_METRICS = ["# Comps", "Cell Size", "SNR", "Event Rate"]
FILTER_OPS = [">", ">=", "=", "<=", "<"]


# --------------------------------------------------------------------------------------
class ScrollFrame(ttk.Frame):
    """A frame with a vertical scrollbar, for tabs with many parameters."""

    def __init__(self, master):
        super().__init__(master)
        self.canvas = tk.Canvas(self, highlightthickness=0, borderwidth=0)
        sb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas, padding=8)
        self.inner.bind("<Configure>",
                        lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.canvas.bind("<Configure>",
                         lambda e: self.canvas.itemconfigure(self._win, width=e.width))
        self.canvas.configure(yscrollcommand=sb.set)
        self.canvas.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        self.inner.bind("<Enter>", lambda e: self._wheel(True))
        self.inner.bind("<Leave>", lambda e: self._wheel(False))

    def _wheel(self, on):
        if on:
            self.canvas.bind_all("<MouseWheel>", lambda e: self.canvas.yview_scroll(
                int(-e.delta / 120) if abs(e.delta) >= 120 else -int(e.delta), "units"))
            self.canvas.bind_all("<Button-4>", lambda e: self.canvas.yview_scroll(-1, "units"))
            self.canvas.bind_all("<Button-5>", lambda e: self.canvas.yview_scroll(1, "units"))
        else:
            for ev in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                self.canvas.unbind_all(ev)


class FilterEditor(ttk.Frame):
    """Editable list of auto accept/reject filters (metric, operator, value)."""

    def __init__(self, master):
        super().__init__(master)
        self.rows = []
        self.body = ttk.Frame(self)
        self.body.pack(fill="x")
        ttk.Button(self, text="Add filter", command=lambda: self.add()).pack(anchor="w", pady=4)

    def add(self, metric="SNR", op=">", value="0"):
        f = ttk.Frame(self.body)
        f.pack(fill="x", pady=1)
        mv, ov, vv = tk.StringVar(value=metric), tk.StringVar(value=op), tk.StringVar(value=str(value))
        ttk.Combobox(f, textvariable=mv, values=FILTER_METRICS, width=11, state="readonly").pack(side="left")
        ttk.Combobox(f, textvariable=ov, values=FILTER_OPS, width=4, state="readonly").pack(side="left", padx=4)
        ttk.Entry(f, textvariable=vv, width=8).pack(side="left")
        row = (f, mv, ov, vv)
        ttk.Button(f, text="Remove", command=lambda: self.remove(row)).pack(side="left", padx=4)
        self.rows.append(row)

    def remove(self, row):
        row[0].destroy()
        self.rows.remove(row)

    def get(self):
        out = []
        for _, mv, ov, vv in self.rows:
            try:
                out.append((mv.get(), ov.get(), float(vv.get())))
            except ValueError:
                raise ValueError(f"Accept/reject filter '{mv.get()} {ov.get()}': "
                                 f"'{vv.get()}' is not a number")
        return out

    def set(self, filters):
        for r in list(self.rows):
            self.remove(r)
        for m, o, v in filters:
            self.add(m, o, v)


class CropEditor(ttk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.on = tk.BooleanVar(value=False)
        ttk.Checkbutton(self, text="Crop to rectangle", variable=self.on, command=self._toggle
                        ).grid(row=0, column=0, columnspan=4, sticky="w")
        self.vars, self.entries = [], []
        for i, lab in enumerate(("top", "left", "bottom", "right")):
            ttk.Label(self, text=lab).grid(row=1 + i // 2, column=2 * (i % 2), sticky="e",
                                           padx=(8 if i % 2 else 0, 4))
            v = tk.StringVar(value="0")
            e = ttk.Entry(self, textvariable=v, width=7)
            e.grid(row=1 + i // 2, column=2 * (i % 2) + 1, sticky="w", pady=1)
            self.vars.append(v)
            self.entries.append(e)
        self._toggle()

    def _toggle(self):
        for e in self.entries:
            e.configure(state="normal" if self.on.get() else "disabled")

    def get(self):
        if not self.on.get():
            return None
        try:
            t, l, b, r = (int(v.get()) for v in self.vars)
        except ValueError:
            raise ValueError("Crop: top/left/bottom/right must be whole numbers")
        if b <= t or r <= l:
            raise ValueError("Crop: bottom must be > top and right must be > left")
        return (t, l, b, r)

    def set(self, rect):
        self.on.set(rect is not None)
        for v, x in zip(self.vars, rect or (0, 0, 0, 0)):
            v.set(str(x))
        self._toggle()


# --------------------------------------------------------------------------------------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Inscopix batch processing")
        self.geometry("1100x820")
        self.minsize(900, 640)
        self.proc = None
        self.q = queue.Queue()
        self.widgets = {}     # field -> (kind, var or editor)
        self.defaults = Params()

        self._build()
        self.load_settings(LAST_SETTINGS, quiet=True)
        self.protocol("WM_DELETE_WINDOW", self.on_close)
        self.after(100, self._poll)

    # ---------------------------------------------------------------- layout
    def _build(self):
        pw = ttk.PanedWindow(self, orient="vertical")
        pw.pack(fill="both", expand=True, padx=8, pady=8)
        top = ttk.Frame(pw)
        pw.add(top, weight=3)
        logf = ttk.LabelFrame(pw, text="Log", padding=4)
        pw.add(logf, weight=1)

        # --- left column: files, output, run options, steps
        left = ttk.Frame(top)
        left.pack(side="left", fill="both", expand=False)
        right = ttk.Frame(top)
        right.pack(side="left", fill="both", expand=True, padx=(8, 0))

        ff = ttk.LabelFrame(left, text="Input .isxd files", padding=6)
        ff.pack(fill="both", expand=True)
        lb_frame = ttk.Frame(ff)
        lb_frame.pack(fill="both", expand=True)
        self.files = tk.Listbox(lb_frame, selectmode="extended", width=40, height=10)
        sb = ttk.Scrollbar(lb_frame, orient="vertical", command=self.files.yview)
        self.files.configure(yscrollcommand=sb.set)
        self.files.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        bf = ttk.Frame(ff)
        bf.pack(fill="x", pady=(4, 0))
        ttk.Button(bf, text="Add files…", command=self.add_files).pack(side="left")
        ttk.Button(bf, text="Add folder…", command=self.add_folder).pack(side="left", padx=4)
        ttk.Button(bf, text="Remove", command=self.remove_files).pack(side="left")
        ttk.Button(bf, text="Clear", command=lambda: (self.files.delete(0, "end"), self._count())
                   ).pack(side="left", padx=4)
        ttk.Button(bf, text="▲", width=3, command=lambda: self.move(-1)).pack(side="left")
        ttk.Button(bf, text="▼", width=3, command=lambda: self.move(1)).pack(side="left")
        bf2 = ttk.Frame(ff)
        bf2.pack(fill="x")
        self.recursive = tk.BooleanVar(value=False)
        ttk.Checkbutton(bf2, text="Include subfolders", variable=self.recursive).pack(side="left")
        self.count_lbl = ttk.Label(bf2, text="0 files")
        self.count_lbl.pack(side="right")

        of = ttk.LabelFrame(left, text="Output", padding=6)
        of.pack(fill="x", pady=(8, 0))
        self.out_dir = tk.StringVar()
        ttk.Entry(of, textvariable=self.out_dir).pack(side="left", fill="x", expand=True)
        ttk.Button(of, text="Browse…", command=self.pick_out).pack(side="left", padx=(4, 0))

        rf = ttk.LabelFrame(left, text="Run options", padding=6)
        rf.pack(fill="x", pady=(8, 0))
        self.mode = tk.StringVar(value="independent")
        ttk.Radiobutton(rf, text="Each file is a separate session", value="independent",
                        variable=self.mode).grid(row=0, column=0, columnspan=3, sticky="w")
        ttk.Radiobutton(rf, text="All files are one series, name:", value="series",
                        variable=self.mode).grid(row=1, column=0, columnspan=2, sticky="w")
        self.series_name = tk.StringVar(value="series")
        ttk.Entry(rf, textvariable=self.series_name, width=14).grid(row=1, column=2, sticky="w")
        ttk.Label(rf, text="Backend").grid(row=2, column=0, sticky="w", pady=(4, 0))
        self.backend = tk.StringVar(value="auto")
        ttk.Combobox(rf, textvariable=self.backend, values=["auto", "isx", "numpy"], width=8,
                     state="readonly").grid(row=2, column=1, sticky="w", pady=(4, 0))
        ttk.Label(rf, text="Python").grid(row=3, column=0, sticky="w")
        self.python = tk.StringVar(value=sys.executable)
        ttk.Entry(rf, textvariable=self.python, width=28).grid(row=3, column=1, columnspan=2, sticky="we")
        ttk.Button(rf, text="…", width=2, command=self.pick_python).grid(row=3, column=3)
        self.overwrite = tk.BooleanVar(value=False)
        self.cont = tk.BooleanVar(value=True)
        ttk.Checkbutton(rf, text="Overwrite existing outputs", variable=self.overwrite
                        ).grid(row=4, column=0, columnspan=3, sticky="w")
        ttk.Checkbutton(rf, text="Continue with next file on error", variable=self.cont
                        ).grid(row=5, column=0, columnspan=3, sticky="w")

        # --- right: section list + parameter pages + buttons
        body = ttk.Frame(right)
        body.pack(fill="both", expand=True)
        self.nav = tk.Listbox(body, width=17, exportselection=False, activestyle="none",
                              font=("TkDefaultFont", 10))
        self.nav.pack(side="left", fill="y")
        pages = ttk.Frame(body, relief="groove", borderwidth=1)
        pages.pack(side="left", fill="both", expand=True, padx=(6, 0))
        pages.rowconfigure(0, weight=1)
        pages.columnconfigure(0, weight=1)
        self.pages = []
        self._help_labels = []
        pages.bind("<Configure>", lambda e: [h.configure(wraplength=max(200, e.width - 40))
                                             for h in self._help_labels])

        # Steps page
        sfm = ScrollFrame(pages)
        sfm.grid(row=0, column=0, sticky="nsew")
        self.pages.append(sfm)
        self.nav.insert("end", "Steps")
        g = sfm.inner
        ttk.Label(g, text="Processing steps to run", font=("TkDefaultFont", 10, "bold")
                  ).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 6))
        self.steps = {}
        for i, st in enumerate(ALL_STEPS):
            v = tk.BooleanVar(value=True)
            self.steps[st] = v
            ttk.Checkbutton(g, text=STEP_LABELS[st], variable=v).grid(row=i + 1, column=0, sticky="w")
        r = len(ALL_STEPS) + 1
        ttk.Label(g, text="Cell identification method").grid(row=r, column=0, sticky="w", pady=(10, 0))
        mv = tk.StringVar(value=self.defaults.method)
        ttk.Combobox(g, textvariable=mv, values=["pca-ica", "cnmfe"], width=10, state="readonly"
                     ).grid(row=r, column=1, sticky="w", pady=(10, 0))
        self.widgets["method"] = ("choice", mv)
        ttk.Label(g, foreground="gray40", wraplength=380, justify="left", text=(
            "Steps whose output files already exist are skipped unless 'Overwrite' is on, so "
            "you can run the first steps, check the motion-corrected movie, then continue. "
            "PCA-ICA runs on the DF/F movie; CNMFe runs on the motion-corrected movie.")
                  ).grid(row=r + 1, column=0, columnspan=2, sticky="w", pady=(12, 0))

        for tab_name, items in TABS:
            sfm = ScrollFrame(pages)
            sfm.grid(row=0, column=0, sticky="nsew")
            self.pages.append(sfm)
            self.nav.insert("end", tab_name)
            g = sfm.inner
            ttk.Label(g, text=tab_name, font=("TkDefaultFont", 10, "bold")
                      ).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 6))
            for r, (name, label, kind, extra, helptext) in enumerate(items, start=1):
                row = 3 * r
                wide = kind in ("crop", "filters")
                ttk.Label(g, text=label).grid(row=row, column=0, sticky="nw", pady=(6, 0),
                                              padx=(0, 10))
                w = self._make_widget(g, name, kind, extra)
                if wide:  # own row, indented under the label
                    w.grid(row=row + 1, column=0, columnspan=2, sticky="w", padx=(16, 0))
                else:
                    w.grid(row=row, column=1, sticky="w", pady=(6, 0))
                if helptext:
                    hl = ttk.Label(g, text=helptext, foreground="gray40", justify="left",
                                   font=("TkDefaultFont", 9))
                    hl.grid(row=row + 2, column=0, columnspan=2, sticky="w")
                    self._help_labels.append(hl)
            ttk.Button(g, text="Reset this section to defaults",
                       command=lambda it=items: self.reset([i[0] for i in it])
                       ).grid(row=3 * len(items) + 3, column=0, sticky="w", pady=(12, 0))

        self.nav.bind("<<ListboxSelect>>", lambda e: self.show_page())
        self.nav.selection_set(0)
        self.show_page()

        bb = ttk.Frame(right)
        bb.pack(fill="x", pady=(8, 0))
        ttk.Button(bb, text="Load settings…", command=lambda: self.load_settings()).pack(side="left")
        ttk.Button(bb, text="Save settings…", command=self.save_settings).pack(side="left", padx=4)
        ttk.Button(bb, text="Reset all", command=lambda: self.reset(None)).pack(side="left")
        self.cancel_btn = ttk.Button(bb, text="Cancel", command=self.cancel, state="disabled")
        self.cancel_btn.pack(side="right")
        self.run_btn = ttk.Button(bb, text="Run", command=self.run)
        self.run_btn.pack(side="right", padx=4)
        ttk.Button(bb, text="Output…", command=self.open_out).pack(side="right")
        ttk.Button(bb, text="View results…", command=self.view_results).pack(side="right", padx=4)

        pf = ttk.Frame(right)
        pf.pack(fill="x", pady=(6, 0))
        self.status = ttk.Label(pf, text="Ready")
        self.status.pack(side="left")
        self.progress = ttk.Progressbar(pf, mode="determinate", length=260)
        self.progress.pack(side="right")

        # --- log
        self.log = tk.Text(logf, height=10, wrap="none", font=("TkFixedFont", 9))
        lsb = ttk.Scrollbar(logf, orient="vertical", command=self.log.yview)
        self.log.configure(yscrollcommand=lsb.set, state="disabled")
        self.log.pack(side="left", fill="both", expand=True)
        lsb.pack(side="right", fill="y")
        self.log.tag_configure("err", foreground="red3")
        self.log.tag_configure("warn", foreground="darkorange3")
        self.log.tag_configure("step", foreground="blue4")

    def show_page(self):
        sel = self.nav.curselection()
        if sel:
            self.pages[sel[0]].tkraise()

    def _make_widget(self, parent, name, kind, extra):
        default = getattr(self.defaults, name)
        if kind == "bool":
            v = tk.BooleanVar(value=default)
            w = ttk.Checkbutton(parent, variable=v)
        elif kind == "choice":
            v = tk.StringVar(value=default)
            w = ttk.Combobox(parent, textvariable=v, values=extra, width=20, state="readonly")
        elif kind == "crop":
            w = v = CropEditor(parent)
            v.set(default)
        elif kind == "filters":
            w = v = FilterEditor(parent)
            v.set(default)
        elif kind == "file":
            v = tk.StringVar(value=default)
            w = ttk.Frame(parent)
            ttk.Entry(w, textvariable=v, width=26).pack(side="left")
            ttk.Button(w, text="…", width=2, command=lambda: v.set(filedialog.askopenfilename(
                filetypes=[("Inscopix image", "*.isxd"), ("All files", "*")]) or v.get())
                       ).pack(side="left", padx=2)
            ttk.Button(w, text="×", width=2, command=lambda: v.set("")).pack(side="left")
        else:  # int, float, opt_int, opt_float
            v = tk.StringVar(value="" if default is None else str(default))
            w = ttk.Entry(parent, textvariable=v, width=12)
        self.widgets[name] = (kind, v)
        return w

    # ---------------------------------------------------------------- params
    def get_params(self) -> Params:
        P = Params()
        labels = {n: l for _, items in TABS for n, l, *_ in items}
        for name, (kind, v) in self.widgets.items():
            lab = labels.get(name, name)
            if kind in ("crop", "filters"):
                val = v.get()
            elif kind in ("bool", "choice", "file"):
                val = v.get()
            else:
                txt = v.get().strip()
                if txt == "" and kind.startswith("opt_"):
                    val = None
                else:
                    try:
                        val = int(txt) if kind.endswith("int") else float(txt)
                    except ValueError:
                        raise ValueError(f"{lab}: '{txt}' is not a valid "
                                         f"{'whole number' if kind.endswith('int') else 'number'}")
                    if kind.endswith("int") and val < 0:
                        raise ValueError(f"{lab}: must not be negative")
            setattr(P, name, val)
        if P.spatial_downsample_factor < 1 or P.temporal_downsample_factor < 1:
            raise ValueError("Downsample factors must be at least 1")
        if not 0 <= P.sf_low_cutoff < P.sf_high_cutoff <= 0.5:
            raise ValueError("Spatial filter: need 0 <= low cutoff < high cutoff <= 0.5")
        if P.num_ics < 1:
            raise ValueError("PCA-ICA: number of ICs must be at least 1")
        if P.num_pcs is not None and P.num_pcs < P.num_ics:
            raise ValueError("PCA-ICA: number of PCs must be >= number of ICs")
        return P

    def set_params(self, d: dict):
        for name, (kind, v) in self.widgets.items():
            if name not in d:
                continue
            val = d[name]
            if kind == "crop":
                v.set(tuple(val) if val else None)
            elif kind == "filters":
                v.set(val)
            elif kind == "bool":
                v.set(bool(val))
            else:
                v.set("" if val is None else str(val))

    def reset(self, names):
        d = asdict(Params())
        self.set_params({k: d[k] for k in (names or d.keys())})

    def settings_dict(self):
        return {
            "params": asdict(self.get_params()),
            "run": {
                "files": list(self.files.get(0, "end")), "out_dir": self.out_dir.get(),
                "mode": self.mode.get(), "series_name": self.series_name.get(),
                "backend": self.backend.get(), "python": self.python.get(),
                "overwrite": self.overwrite.get(), "continue_on_error": self.cont.get(),
                "steps": [s for s in ALL_STEPS if self.steps[s].get()],
                "recursive": self.recursive.get(),
            },
        }

    def save_settings(self):
        try:
            d = self.settings_dict()
        except ValueError as e:
            return messagebox.showerror("Invalid parameter", str(e))
        path = filedialog.asksaveasfilename(defaultextension=".json",
                                            filetypes=[("Settings", "*.json")])
        if path:
            with open(path, "w") as f:
                json.dump(d, f, indent=2)
            self.status.configure(text=f"Settings saved to {os.path.basename(path)}")

    def load_settings(self, path=None, quiet=False):
        path = path or filedialog.askopenfilename(filetypes=[("Settings", "*.json")])
        if not path or not os.path.exists(path):
            return
        try:
            with open(path) as f:
                d = json.load(f)
            # Accept both GUI settings files and bare parameter files (--config).
            params = d.get("params", d)
            self.set_params(params)
            run = d.get("run", {})
            if "files" in run:
                self.files.delete(0, "end")
                for fpath in run["files"]:
                    self.files.insert("end", fpath)
            for key, var in (("out_dir", self.out_dir), ("mode", self.mode),
                             ("series_name", self.series_name), ("backend", self.backend),
                             ("python", self.python), ("overwrite", self.overwrite),
                             ("continue_on_error", self.cont), ("recursive", self.recursive)):
                if key in run:
                    var.set(run[key])
            if "steps" in run:
                for s in ALL_STEPS:
                    self.steps[s].set(s in run["steps"])
            self._count()
        except Exception as e:  # noqa: BLE001
            if not quiet:
                messagebox.showerror("Could not load settings", str(e))

    # ---------------------------------------------------------------- files
    def _count(self):
        n = self.files.size()
        self.count_lbl.configure(text=f"{n} file{'s' if n != 1 else ''}")

    def _add(self, paths):
        have = set(self.files.get(0, "end"))
        for p in paths:
            if p not in have:
                self.files.insert("end", p)
                have.add(p)
        self._count()

    def add_files(self):
        self._add(filedialog.askopenfilenames(
            filetypes=[("Inscopix movies", "*.isxd"), ("All files", "*")]))

    def add_folder(self):
        d = filedialog.askdirectory()
        if d:
            found = collect_inputs([d], self.recursive.get())
            if not found:
                messagebox.showinfo("No files", f"No .isxd movies found in {d}")
            self._add(found)

    def remove_files(self):
        for i in reversed(self.files.curselection()):
            self.files.delete(i)
        self._count()

    def move(self, delta):
        sel = list(self.files.curselection())
        if not sel:
            return
        items = list(self.files.get(0, "end"))
        order = sel if delta < 0 else sel[::-1]
        for i in order:
            j = i + delta
            if 0 <= j < len(items) and j not in sel:
                items[i], items[j] = items[j], items[i]
                sel[sel.index(i)] = j
        self.files.delete(0, "end")
        for it in items:
            self.files.insert("end", it)
        for i in sel:
            self.files.selection_set(i)

    def pick_out(self):
        d = filedialog.askdirectory()
        if d:
            self.out_dir.set(d)

    def view_results(self):
        """Open the results viewer on the output folder (or ask for a folder)."""
        from isx_viewer import ResultsViewer
        d = self.out_dir.get().strip()
        have = bool(d) and os.path.isdir(d)
        v = ResultsViewer(self, d if have else None)
        if not have:
            v.open_folder()

    def pick_python(self):
        p = filedialog.askopenfilename(title="Python interpreter (e.g. of the IDPS isx environment)")
        if p:
            self.python.set(p)

    def open_out(self):
        d = self.out_dir.get()
        if not d or not os.path.isdir(d):
            return
        if sys.platform.startswith("win"):
            os.startfile(d)  # noqa
        else:
            subprocess.Popen(["open" if sys.platform == "darwin" else "xdg-open", d])

    # ---------------------------------------------------------------- running
    def _write_log(self, line):
        tag = None
        if " ERROR " in line or "Traceback" in line or "Error" in line:
            tag = "err"
        elif " WARNING " in line:
            tag = "warn"
        elif "====" in line or "---- [" in line:
            tag = "step"
        self.log.configure(state="normal")
        self.log.insert("end", line + ("\n" if not line.endswith("\n") else ""), tag)
        self.log.see("end")
        self.log.configure(state="disabled")

    def run(self):
        files = list(self.files.get(0, "end"))
        out = self.out_dir.get().strip()
        steps = [s for s in ALL_STEPS if self.steps[s].get()]
        if not files:
            return messagebox.showerror("No input", "Add at least one .isxd file.")
        if not out:
            return messagebox.showerror("No output folder", "Choose an output folder.")
        if not steps:
            return messagebox.showerror("No steps", "Select at least one processing step.")
        if self.mode.get() == "series" and not re.fullmatch(r"[\w.-]+", self.series_name.get()):
            return messagebox.showerror("Series name", "Use letters, digits, '-', '_' or '.'.")
        try:
            P = self.get_params()
            settings = self.settings_dict()
        except ValueError as e:
            return messagebox.showerror("Invalid parameter", str(e))
        if P.method == "cnmfe" and self.backend.get() == "numpy" and "cells" in steps:
            return messagebox.showerror("CNMFe", "CNMFe needs the isx backend (IDPS API).")

        os.makedirs(out, exist_ok=True)
        stamp = time.strftime("%Y%m%d-%H%M%S")
        cfg = os.path.join(out, f"params_{stamp}.json")
        with open(cfg, "w") as f:
            json.dump(asdict(P), f, indent=2)
        with open(LAST_SETTINGS, "w") as f:
            json.dump(settings, f, indent=2)

        cmd = [self.python.get() or sys.executable, "-u", SCRIPT, *files, "-o", out,
               "--config", cfg, "--mode", self.mode.get(), "--series-name",
               self.series_name.get(), "--backend", self.backend.get(), "--steps", *steps]
        if self.overwrite.get():
            cmd.append("--overwrite")
        if self.cont.get():
            cmd.append("--continue-on-error")

        self.log.configure(state="normal")
        self.log.delete("1.0", "end")
        self.log.configure(state="disabled")
        self._write_log(f"Parameters saved to {cfg}")
        self._total = len(files) if self.mode.get() == "independent" else 1
        self.progress.configure(maximum=self._total * len(steps), value=0)
        self._done_steps = 0
        self._t0 = time.time()
        kw = {}
        if sys.platform.startswith("win"):
            kw["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
        try:
            self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                         text=True, bufsize=1, cwd=HERE, **kw)
        except OSError as e:
            return messagebox.showerror("Could not start", f"{e}\n\n{cmd[0]}")
        threading.Thread(target=self._reader, args=(self.proc,), daemon=True).start()
        self.run_btn.configure(state="disabled")
        self.cancel_btn.configure(state="normal")
        self.status.configure(text="Running…")

    def _reader(self, proc):
        for line in proc.stdout:
            self.q.put(("line", line.rstrip("\n")))
        proc.wait()
        self.q.put(("exit", proc.returncode))

    def _poll(self):
        try:
            while True:
                kind, val = self.q.get_nowait()
                if kind == "line":
                    self._write_log(val)
                    self._track(val)
                else:
                    self._finished(val)
        except queue.Empty:
            pass
        self.after(100, self._poll)

    def _track(self, line):
        m = re.search(r"---- \[(\d+)/(\d+)\] (.*)", line)
        if m:
            self._current = f"File {m.group(1)}/{m.group(2)}: {m.group(3)}"
        step_names = {"[preprocess]": "pp", "[spatial filter]": "bp", "[motion correct]": "mc",
                      "[dF/F]": "dff", "[PCA-ICA]": "cells", "[CNMFe]": "cells",
                      "[event detection]": "events", "[auto accept/reject]": "accept",
                      "[export cells]": "export"}
        for key, s in step_names.items():
            if key in line:
                self._done_steps += 1
                self.progress.configure(value=min(self._done_steps, self.progress["maximum"]))
                prefix = getattr(self, "_current", "")
                self.status.configure(text=f"{prefix}  –  {STEP_LABELS[s]}".strip(" –"))

    def _finished(self, code):
        self.proc = None
        self.run_btn.configure(state="normal")
        self.cancel_btn.configure(state="disabled")
        mins = (time.time() - self._t0) / 60
        if code == 0:
            self.progress.configure(value=self.progress["maximum"])
            self.status.configure(text=f"Finished in {mins:.1f} min")
            messagebox.showinfo("Done", f"Processing finished in {mins:.1f} min.")
        elif code is not None and code < 0 or getattr(self, "_cancelled", False):
            self.status.configure(text="Cancelled")
        else:
            self.status.configure(text=f"Finished with errors (exit code {code})")
            messagebox.showwarning("Errors", "Some files failed; see the log.")
        self._cancelled = False

    def cancel(self):
        if self.proc and self.proc.poll() is None:
            if messagebox.askyesno("Cancel", "Stop processing? The current step will be "
                                   "redone next time."):
                self._cancelled = True
                self.proc.terminate()
                self._write_log("Cancelled by user.")

    def on_close(self):
        if self.proc and self.proc.poll() is None:
            if not messagebox.askyesno("Quit", "Processing is running. Stop it and quit?"):
                return
            self.proc.terminate()
        try:
            with open(LAST_SETTINGS, "w") as f:
                json.dump(self.settings_dict(), f, indent=2)
        except Exception:  # noqa: BLE001 - never block closing
            pass
        self.destroy()


if __name__ == "__main__":
    App().mainloop()
