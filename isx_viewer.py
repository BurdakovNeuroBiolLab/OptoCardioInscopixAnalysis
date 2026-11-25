#!/usr/bin/env python
"""
Viewer for processed Inscopix results: cell map with ROI outlines, cell traces and
cell images.

    python isx_viewer.py [folder]          (or "View results…" in isx_batch_gui.py)

Reads, from a folder (searched recursively), the files written by isx_batch_process.py
or exported from IDPS:
    <name>_celltrace.csv          traces (IDPS layout: names row, status row, time + cells)
    <name>_celltrace-props.csv    props (Name, Status, CentroidX, CentroidY, Size, ...)
    <name>_cellmap.tif            cell map with the ROI outlines drawn in white
    <name>_cellimage_C00.tif ...  optional per-cell images (IDPS: <name>_C00.tif)

Window:
    left            the cell map, all ROI outlines in yellow
    right, bottom   one trace per cell (10 visible, scroll for more), in yellow, with two
                    check boxes: the left one selects the cell (trace and ROI get their
                    own colour, several cells can be selected); the right one shows that
                    cell's image above (one at a time)
    right, top      the image of the cell ticked on the right, or "Cell Images Not Found"

Outline pixels are matched to cells the way plot_cellmap_and_traces.m does it: the map's
scale relative to the cell set is estimated from the outlines and the props centroids,
then every outline pixel goes to the nearest centroid within 0.6 x Size.

Needs only NumPy and Tkinter; uses tifffile or Pillow to read TIFFs when installed,
otherwise a built-in reader for uncompressed TIFFs.
"""

import csv
import os
import queue
import re
import struct
import subprocess
import sys
import tempfile
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np

YELLOW = "#ffd400"
SELECT_COLOURS = ["#00e5ff", "#ff4081", "#76ff03", "#ff9100", "#d500f9", "#2979ff",
                  "#ff1744", "#1de9b6", "#c6ff00", "#f50057", "#00b0ff", "#ffab40"]
BG = "#111111"
CONTOUR_MAX_FRAC = 0.6     # as plot_cellmap_and_traces.m
VISIBLE_TRACES = 10


# --------------------------------------------------------------------------------------
# File reading
# --------------------------------------------------------------------------------------
def read_tiff(path):
    """First page of a TIFF as a NumPy array (H x W or H x W x C)."""
    try:
        import tifffile
        return tifffile.imread(path, key=0)
    except ImportError:
        pass
    except Exception:  # noqa: BLE001 - fall through to the other readers
        pass
    try:
        from PIL import Image
        with Image.open(path) as im:
            return np.array(im)
    except Exception:  # noqa: BLE001
        pass
    return _read_tiff_basic(path)


def _read_tiff_basic(path):
    """Minimal reader: uncompressed, strip-based, chunky TIFF (first page)."""
    with open(path, "rb") as f:
        data = f.read()
    bo = {b"II": "<", b"MM": ">"}.get(data[:2])
    if bo is None or struct.unpack(bo + "H", data[2:4])[0] != 42:
        raise ValueError(f"{path}: not a TIFF file")
    off = struct.unpack(bo + "I", data[4:8])[0]
    n = struct.unpack(bo + "H", data[off:off + 2])[0]
    sizes = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8, 16: 8}
    codes = {1: "B", 3: "H", 4: "I", 6: "b", 8: "h", 9: "i", 11: "f", 12: "d", 16: "Q"}
    tags = {}
    for i in range(n):
        e = off + 2 + 12 * i
        tag, typ, cnt = struct.unpack(bo + "HHI", data[e:e + 8])
        nbytes = sizes.get(typ, 1) * cnt
        vo = e + 8 if nbytes <= 4 else struct.unpack(bo + "I", data[e + 8:e + 12])[0]
        if typ in codes:
            tags[tag] = list(struct.unpack(bo + codes[typ] * cnt, data[vo:vo + nbytes]))
    w, h = tags[256][0], tags[257][0]
    bits = tags.get(258, [1])[0]
    spp = tags.get(277, [1])[0]
    if tags.get(259, [1])[0] != 1:
        raise ValueError(f"{path}: compressed TIFF; install tifffile or Pillow to read it")
    fmt = tags.get(339, [1])[0]
    kind = {1: "u", 2: "i", 3: "f"}[fmt]
    dt = np.dtype(f"{bo}{kind}{bits // 8}")
    raw = b"".join(data[o:o + c] for o, c in zip(tags[273], tags[279]))
    arr = np.frombuffer(raw, dtype=dt, count=w * h * spp)
    return arr.reshape(h, w, spp) if spp > 1 else arr.reshape(h, w)


def read_props(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows or not {"Name", "CentroidX", "CentroidY"} <= set(k.strip() for k in rows[0]):
        raise ValueError(f"{os.path.basename(path)} does not look like a props file")
    rows = [{k.strip(): (v or "").strip() for k, v in r.items()} for r in rows]

    def num(r, k, d):
        try:
            return float(r[k])
        except (KeyError, ValueError):
            return d
    return {
        "names": [r["Name"] for r in rows],
        "status": [r.get("Status", "undecided").lower() for r in rows],
        "cx": np.array([num(r, "CentroidX", np.nan) for r in rows]),
        "cy": np.array([num(r, "CentroidY", np.nan) for r in rows]),
        "size": np.array([num(r, "Size", 10.0) for r in rows]),
    }


def read_traces(path):
    """Returns (time, names, statuses or None, data[T x N])."""
    with open(path, newline="") as f:
        l1 = f.readline()
        l2 = f.readline()
    names = [x.strip().strip('"') for x in l1.rstrip("\r\n").split(",")][1:]
    row2 = [x.strip().strip('"') for x in l2.rstrip("\r\n").split(",")][1:]

    def isnum(x):
        try:
            float(x)
            return True
        except ValueError:
            return x == ""
    has_status = not all(isnum(x) for x in row2)
    skip = 2 if has_status else 1
    try:
        import pandas as pd
        d = pd.read_csv(path, skiprows=skip, header=None, skipinitialspace=True).to_numpy(float)
    except ImportError:
        d = np.genfromtxt(path, delimiter=",", skip_header=skip)
    d = np.atleast_2d(d)[:, :len(names) + 1]
    return d[:, 0], names, ([s.lower() for s in row2] if has_status else None), d[:, 1:]


def find_result_sets(folder):
    """All <name>_celltrace.csv under folder, each with its props, map and cell images."""
    sets = []
    for root, dirs, files in os.walk(folder):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        lower = {f.lower(): f for f in files if not f.startswith(".")}
        for f in sorted(files):
            if f.startswith(".") or not re.search(r"celltrace\.csv$", f, re.I):
                continue
            base = f[:-len("celltrace.csv")]
            props = lower.get((base + "celltrace-props.csv").lower())
            cmap = lower.get((base + "cellmap.tif").lower()) or lower.get((base + "cellmap.tiff").lower())
            if cmap is None:  # a single cell map in the folder will do
                maps = [v for k, v in lower.items() if re.search(r"cellmap\.tiff?$", k)]
                cmap = maps[0] if len(maps) == 1 else None
            sets.append({
                "label": os.path.relpath(os.path.join(root, f), folder),
                "dir": root, "base": base, "traces": os.path.join(root, f),
                "props": os.path.join(root, props) if props else None,
                "map": os.path.join(root, cmap) if cmap else None,
                "files": lower,
            })
    return sets


def find_cell_image(rs, name):
    """Per-cell image of `name` (e.g. C015), or None."""
    files = rs["files"]
    for cand in (f"{rs['base']}cellimage_{name}", f"{rs['base']}{name}",
                 f"{rs['base'].rstrip('_-')}_{name}"):
        for ext in (".tif", ".tiff"):
            hit = files.get((cand + ext).lower())
            if hit:
                return os.path.join(rs["dir"], hit)
    pat = re.compile(rf"[_-]{re.escape(name)}\.tiff?$", re.I)
    for k, v in files.items():
        if pat.search(k) and "cellmap" not in k:
            return os.path.join(rs["dir"], v)
    return None


# --------------------------------------------------------------------------------------
# Cell map: outlines and their owners
# --------------------------------------------------------------------------------------
def analyse_map(img, props):
    """Background (H x W float 0-1), outline pixel coordinates (y, x) and the index of
    the cell each outline pixel belongs to (-1 = none)."""
    img = np.asarray(img)
    maxval = np.iinfo(img.dtype).max if np.issubdtype(img.dtype, np.integer) else 1.0
    f = img.astype(np.float64)
    if f.ndim == 3:
        f = f[..., :3] if f.shape[2] >= 3 else f[..., :1]
        white = np.all(f == maxval, axis=2)
        gray = f.mean(2)
    else:
        white = f == maxval
        gray = f
    bg = gray.copy()
    if white.any() and (~white).any():
        bg[white] = np.median(gray[~white])
    lo, hi = np.percentile(bg, [1, 99.8])
    bg = np.clip((bg - lo) / max(hi - lo, 1e-12), 0, 1)

    wy, wx = np.nonzero(white)
    owner = np.full(wy.size, -1)
    ok = np.isfinite(props["cx"]) & np.isfinite(props["cy"])
    if wy.size and ok.any():
        cx, cy, size = props["cx"][ok], props["cy"][ok], props["size"][ok]
        X, Y = wx + 1.0, wy + 1.0          # 1-based, as MATLAB
        step = max(1, wy.size // 5000)
        sx, sy = X[::step], Y[::step]

        def score(s):
            d2 = ((cx[:, None] + .5) * s + .5 - sx[None]) ** 2 + \
                 ((cy[:, None] + .5) * s + .5 - sy[None]) ** 2
            return np.median(np.sqrt(d2.min(0)))
        coarse = np.arange(0.5, 12.01, 0.1)
        b = coarse[np.argmin([score(s) for s in coarse])]
        fine = np.arange(max(0.5, b - 0.1), b + 0.1001, 0.01)
        scale = fine[np.argmin([score(s) for s in fine])]
        cxI, cyI = (cx + .5) * scale + .5, (cy + .5) * scale + .5
        idx_ok = np.nonzero(ok)[0]
        for s in range(0, wy.size, 20000):          # memory-friendly chunks
            d2 = (cxI[:, None] - X[None, s:s + 20000]) ** 2 + (cyI[:, None] - Y[None, s:s + 20000]) ** 2
            near = d2.argmin(0)
            dmin = np.sqrt(d2[near, np.arange(near.size)])
            good = dmin < size[near] * scale * CONTOUR_MAX_FRAC
            o = np.full(near.size, -1)
            o[good] = idx_ok[near[good]]
            owner[s:s + 20000] = o
    return bg, wy, wx, owner


def to_photo(rgb):
    """uint8 H x W x 3 array -> tk.PhotoImage (via PPM, no extra packages)."""
    h, w = rgb.shape[:2]
    data = b"P6\n%d %d\n255\n" % (w, h) + np.ascontiguousarray(rgb, dtype=np.uint8).tobytes()
    return tk.PhotoImage(data=data, format="PPM")


def hex_rgb(c):
    return tuple(int(c[i:i + 2], 16) for i in (1, 3, 5))


# --------------------------------------------------------------------------------------
class ScrollRows(ttk.Frame):
    """Vertically scrolling frame for the trace rows."""

    def __init__(self, master):
        super().__init__(master)
        self.canvas = tk.Canvas(self, highlightthickness=0, borderwidth=0, background=BG)
        self.sb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = tk.Frame(self.canvas, background=BG)
        self.inner.bind("<Configure>",
                        lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.canvas.bind("<Configure>",
                         lambda e: self.canvas.itemconfigure(self._win, width=e.width))
        self.canvas.configure(yscrollcommand=self.sb.set)
        self.sb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        for w in (self.canvas, self.inner):
            w.bind("<Enter>", lambda e: self._wheel(True))
            w.bind("<Leave>", lambda e: self._wheel(False))

    def _wheel(self, on):
        if on:
            self.canvas.bind_all("<MouseWheel>", lambda e: self.canvas.yview_scroll(
                int(-e.delta / 120) if abs(e.delta) >= 120 else -int(e.delta), "units"))
            self.canvas.bind_all("<Button-4>", lambda e: self.canvas.yview_scroll(-1, "units"))
            self.canvas.bind_all("<Button-5>", lambda e: self.canvas.yview_scroll(1, "units"))
        else:
            for ev in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                self.canvas.unbind_all(ev)


class ResultsViewer(tk.Toplevel):
    def __init__(self, master=None, folder=None):
        super().__init__(master)
        self.title("Inscopix results viewer")
        self.geometry("1400x860")
        self.minsize(900, 600)
        self.configure(background=BG)
        self.sets, self.rs = [], None
        self.selected = {}          # cell index -> colour
        self.image_cell = None      # cell index shown in the image panel
        self.rows = []
        self._photo_map = self._photo_img = None
        self._after = {}

        self._build()
        if folder:
            self.after(50, lambda: self.open_folder(folder))

    # ------------------------------------------------------------------ layout
    def _build(self):
        bar = ttk.Frame(self, padding=4)
        bar.pack(fill="x")
        ttk.Button(bar, text="Open folder…", command=self.open_folder).pack(side="left")
        ttk.Label(bar, text="Result:").pack(side="left", padx=(10, 2))
        self.set_var = tk.StringVar()
        self.set_box = ttk.Combobox(bar, textvariable=self.set_var, state="readonly", width=60)
        self.set_box.pack(side="left")
        self.set_box.bind("<<ComboboxSelected>>", lambda e: self.load_set(self.set_box.current()))
        self.hide_rejected = tk.BooleanVar(value=False)
        ttk.Checkbutton(bar, text="Hide rejected cells", variable=self.hide_rejected,
                        command=self._build_rows).pack(side="left", padx=10)
        ttk.Button(bar, text="Clear selection", command=self.clear_selection).pack(side="left")
        self.status = ttk.Label(bar, text="Open a folder with processed results.")
        self.status.pack(side="left", padx=10)

        pw = tk.PanedWindow(self, orient="horizontal", background="#333", sashwidth=5)
        pw.pack(fill="both", expand=True)
        self.map_canvas = tk.Canvas(pw, background=BG, highlightthickness=0)
        pw.add(self.map_canvas, stretch="always", width=760)
        right = tk.PanedWindow(pw, orient="vertical", background="#333", sashwidth=5)
        pw.add(right, stretch="always")

        top = tk.Frame(right, background=BG)
        self.img_title = tk.Label(top, text="", fg="white", bg=BG, font=("TkDefaultFont", 10, "bold"))
        self.img_title.pack(fill="x")
        self.img_canvas = tk.Canvas(top, background=BG, highlightthickness=0)
        self.img_canvas.pack(fill="both", expand=True)
        right.add(top, stretch="always", height=380)

        bottom = tk.Frame(right, background=BG)
        hdr = tk.Frame(bottom, background=BG)
        hdr.pack(fill="x")
        tk.Label(hdr, text="Select  Image   Cell", fg="#bbb", bg=BG,
                 font=("TkDefaultFont", 8)).pack(side="left", padx=4)
        self.time_lbl = tk.Label(hdr, text="", fg="#bbb", bg=BG, font=("TkDefaultFont", 8))
        self.time_lbl.pack(side="right", padx=20)
        self.trace_area = ScrollRows(bottom)
        self.trace_area.pack(fill="both", expand=True)
        right.add(bottom, stretch="always")

        self.map_canvas.bind("<Configure>", lambda e: self._schedule(self.draw_map))
        self.img_canvas.bind("<Configure>", lambda e: self._schedule(self.draw_cell_image))
        self.trace_area.canvas.bind("<Configure>", self._on_trace_resize, add="+")

    def _schedule(self, fn, delay=60):
        """Run fn once after `delay` ms, coalescing repeated requests (resizes)."""
        key = fn.__name__
        if key in self._after:
            self.after_cancel(self._after[key])
        self._after[key] = self.after(delay, lambda: (self._after.pop(key, None), fn()))

    # ------------------------------------------------------------------ loading
    def open_folder(self, folder=None):
        folder = folder or filedialog.askdirectory(parent=self, title="Folder with results")
        if not folder:
            return
        self.sets = find_result_sets(folder)
        if not self.sets:
            messagebox.showinfo("No results", f"No *_celltrace.csv files found in\n{folder}",
                                parent=self)
            return
        self.set_box.configure(values=[s["label"] for s in self.sets])
        self.set_box.current(0)
        self.load_set(0)

    def load_set(self, i):
        rs = self.sets[i]
        missing = [k for k in ("props", "map") if not rs[k]]
        if missing:
            messagebox.showerror("Missing files", f"{rs['label']}: no matching "
                                 + " / ".join({"props": "*celltrace-props.csv",
                                               "map": "*cellmap.tif"}[m] for m in missing),
                                 parent=self)
            return
        self.status.configure(text="Loading…")
        self.config(cursor="watch")

        q = queue.Queue()

        def work():  # no Tk calls in this thread
            try:
                props = read_props(rs["props"])
                t, names, statuses, data = read_traces(rs["traces"])
                img = read_tiff(rs["map"])
                bg, wy, wx, owner = analyse_map(img, props)
                q.put(("ok", dict(props=props, t=t, tnames=names, data=data,
                                  bg=bg, wy=wy, wx=wx, owner=owner)))
            except Exception as e:  # noqa: BLE001
                q.put(("err", f"{type(e).__name__}: {e}"))

        def poll():
            try:
                kind, val = q.get_nowait()
            except queue.Empty:
                self.after(100, poll)
                return
            if kind == "ok":
                self._loaded(rs, val)
            else:
                self._load_failed(val)
        threading.Thread(target=work, daemon=True).start()
        self.after(100, poll)

    def _load_failed(self, err):
        self.config(cursor="")
        self.status.configure(text="Could not load")
        messagebox.showerror("Could not load results", err, parent=self)

    def _loaded(self, rs, res):
        self.config(cursor="")
        self.rs, self.d = rs, res
        props = res["props"]
        col = {n: j for j, n in enumerate(res["tnames"])}
        # cells in props order that have a trace
        self.cells = [i for i, n in enumerate(props["names"]) if n in col]
        self.trace_col = {i: col[props["names"][i]] for i in self.cells}
        X = res["data"]
        self.norm = {}
        for i in self.cells:
            x = X[:, self.trace_col[i]]
            fin = x[np.isfinite(x)]
            lo, hi = (np.percentile(fin, [0.5, 99.5]) if fin.size else (0, 1))
            self.norm[i] = (lo, hi if hi > lo else lo + 1)
        self.selected, self.image_cell = {}, None
        t = res["t"]
        self.time_lbl.configure(text=f"time {t[0]:.0f} – {t[-1]:.0f} s")
        n_acc = sum(s == "accepted" for s in props["status"])
        n_rej = sum(s == "rejected" for s in props["status"])
        self.status.configure(text=f"{len(self.cells)} cells ({n_acc} accepted, {n_rej} rejected), "
                                   f"{len(t)} frames")
        self.title(f"Inscopix results viewer – {rs['label']}")
        self._build_rows()
        self.draw_map()
        self.draw_cell_image()

    # ------------------------------------------------------------------ traces
    def _visible_cells(self):
        st = self.d["props"]["status"]
        return [i for i in self.cells if not (self.hide_rejected.get() and st[i] == "rejected")]

    def _row_height(self):
        h = self.trace_area.canvas.winfo_height()
        return max(28, (h if h > 50 else 400) // VISIBLE_TRACES)

    def _build_rows(self):
        if not self.rs:
            return
        for r in self.rows:
            r["frame"].destroy()
        self.rows = []
        names, st = self.d["props"]["names"], self.d["props"]["status"]
        rh = self._row_height()
        for i in self._visible_cells():
            fr = tk.Frame(self.trace_area.inner, background=BG, height=rh)
            fr.pack(fill="x")
            fr.pack_propagate(False)
            sv = tk.BooleanVar(value=i in self.selected)
            iv = tk.BooleanVar(value=i == self.image_cell)
            tk.Checkbutton(fr, variable=sv, bg=BG, activebackground=BG, selectcolor="white",
                           command=lambda i=i, v=sv: self.toggle_select(i, v.get())
                           ).pack(side="left")
            tk.Checkbutton(fr, variable=iv, bg=BG, activebackground=BG, selectcolor="white",
                           command=lambda i=i, v=iv: self.toggle_image(i, v.get())
                           ).pack(side="left")
            lab = tk.Label(fr, text=names[i] + (" rej" if st[i] == "rejected" else ""),
                           width=9, anchor="w", bg=BG, fg=self._colour(i),
                           font=("TkFixedFont", 9))
            lab.pack(side="left")
            cv = tk.Canvas(fr, background=BG, highlightthickness=0, height=rh)
            cv.pack(side="left", fill="both", expand=True)
            cv.bind("<Configure>", lambda e, i=i: self._schedule(self.draw_traces, 80))
            self.rows.append({"cell": i, "frame": fr, "sel": sv, "img": iv, "label": lab,
                              "canvas": cv, "line": None})
        self.trace_area.canvas.yview_moveto(0)
        self.after(30, self.draw_traces)

    def _on_trace_resize(self, event):
        rh = max(28, event.height // VISIBLE_TRACES)
        for r in self.rows:
            if int(r["frame"].cget("height")) != rh:
                r["frame"].configure(height=rh)
                r["canvas"].configure(height=rh)
        self._schedule(self.draw_traces, 80)

    def draw_traces(self):
        if not self.rows:
            return
        X = self.d["data"]
        w = max(50, self.rows[0]["canvas"].winfo_width())
        h = max(20, self.rows[0]["canvas"].winfo_height())
        T = X.shape[0]
        nb = min(w, T)
        edges = np.linspace(0, T, nb + 1).astype(int)[:-1]
        xs = np.linspace(2, w - 2, nb)
        for r in self.rows:
            i = r["cell"]
            x = X[:, self.trace_col[i]].astype(float)
            if np.isnan(x).any():
                x = np.where(np.isnan(x), np.nanmean(x) if np.isfinite(x).any() else 0, x)
            mn, mx = np.minimum.reduceat(x, edges), np.maximum.reduceat(x, edges)
            lo, hi = self.norm[i]

            def ypix(v):
                return np.clip(h - 3 - (v - lo) / (hi - lo) * (h - 6), 1, h - 1)
            pts = np.empty((2 * nb, 2))
            pts[0::2, 0] = xs
            pts[1::2, 0] = xs
            pts[0::2, 1] = ypix(mn)
            pts[1::2, 1] = ypix(mx)
            cv = r["canvas"]
            cv.delete("all")
            cv.create_line(0, h - 0.5, w, h - 0.5, fill="#262626")
            r["line"] = cv.create_line(*pts.ravel().tolist(), fill=self._colour(i), width=1)

    # ------------------------------------------------------------------ selection
    def _colour(self, i):
        return self.selected.get(i, YELLOW)

    def toggle_select(self, i, on):
        if on:
            used = set(self.selected.values())
            free = [c for c in SELECT_COLOURS if c not in used]
            self.selected[i] = free[0] if free else SELECT_COLOURS[len(self.selected) % len(SELECT_COLOURS)]
        else:
            self.selected.pop(i, None)
        for r in self.rows:
            if r["cell"] == i:
                if r["line"]:
                    r["canvas"].itemconfigure(r["line"], fill=self._colour(i))
                r["label"].configure(fg=self._colour(i))
        self.draw_map()
        if self.image_cell == i:
            self.draw_cell_image()

    def clear_selection(self):
        for i in list(self.selected):
            self.selected.pop(i)
        for r in self.rows:
            r["sel"].set(False)
            r["label"].configure(fg=YELLOW)
            if r["line"]:
                r["canvas"].itemconfigure(r["line"], fill=YELLOW)
        self.draw_map()
        self.draw_cell_image()

    def toggle_image(self, i, on):
        self.image_cell = i if on else None
        for r in self.rows:  # only one image box ticked at a time
            if r["cell"] != i:
                r["img"].set(False)
        self.draw_cell_image()

    # ------------------------------------------------------------------ drawing
    def draw_map(self):
        if not self.rs:
            return
        cv = self.map_canvas
        cw, ch = cv.winfo_width(), cv.winfo_height()
        if cw < 20 or ch < 20:
            return
        bg = self.d["bg"]
        H, W = bg.shape
        s = min(cw / W, ch / H)
        dw, dh = max(1, int(W * s)), max(1, int(H * s))
        ys = np.minimum((np.arange(dh) / s).astype(int), H - 1)
        xs = np.minimum((np.arange(dw) / s).astype(int), W - 1)
        g = (bg[np.ix_(ys, xs)] * 255).astype(np.uint8)
        rgb = np.repeat(g[:, :, None], 3, axis=2)
        wy, wx, owner = self.d["wy"], self.d["wx"], self.d["owner"]
        if wy.size:
            k = max(1, int(np.ceil(s)))                     # outline thickness when enlarging
            cols = np.tile(np.array(hex_rgb(YELLOW), np.uint8), (wy.size, 1))
            for i, c in self.selected.items():
                cols[owner == i] = hex_rgb(c)
            # draw unselected first, selected on top
            order = np.argsort(np.isin(owner, list(self.selected)), kind="stable")
            for dy in range(k):
                for dx in range(k):
                    py = np.minimum((wy[order] * s).astype(int) + dy, dh - 1)
                    px = np.minimum((wx[order] * s).astype(int) + dx, dw - 1)
                    rgb[py, px] = cols[order]
        self._photo_map = to_photo(rgb)
        cv.delete("all")
        cv.create_image(cw // 2, ch // 2, image=self._photo_map, anchor="center")
        cv.create_text(8, 8, anchor="nw", fill="#ddd", font=("TkDefaultFont", 9),
                       text=os.path.basename(self.rs["map"]))

    def draw_cell_image(self):
        cv = self.img_canvas
        cv.delete("all")
        cw, ch = cv.winfo_width(), cv.winfo_height()
        if not self.rs:
            return
        if self.image_cell is None:
            self.img_title.configure(text="")
            cv.create_text(cw // 2, ch // 2, fill="#888", font=("TkDefaultFont", 11),
                           text="Tick the right-hand box of a trace to show its cell image")
            return
        name = self.d["props"]["names"][self.image_cell]
        self.img_title.configure(text=name, fg=self._colour(self.image_cell))
        path = find_cell_image(self.rs, name)
        img = None
        if path:
            try:
                img = np.asarray(read_tiff(path), dtype=np.float64)
                if img.ndim == 3:
                    img = img[..., :3].mean(2)
            except Exception:  # noqa: BLE001
                img = None
        if img is None:
            cv.create_text(cw // 2, ch // 2, fill="#ff6e6e", font=("TkDefaultFont", 14, "bold"),
                           text="Cell Images Not Found")
            return
        img = np.nan_to_num(img)
        lo, hi = np.percentile(img, [1, 99.9])
        g = np.clip((img - lo) / max(hi - lo, 1e-12), 0, 1)
        H, W = g.shape
        s = min(max(cw, 20) / W, max(ch, 20) / H)
        dw, dh = max(1, int(W * s)), max(1, int(H * s))
        ys = np.minimum((np.arange(dh) / s).astype(int), H - 1)
        xs = np.minimum((np.arange(dw) / s).astype(int), W - 1)
        g = g[np.ix_(ys, xs)]
        c = np.array(hex_rgb(self._colour(self.image_cell)), float) / 255
        rgb = (g[:, :, None] * c[None, None, :] * 255).astype(np.uint8)   # tinted
        self._photo_img = to_photo(rgb)
        cv.create_image(cw // 2, ch // 2, image=self._photo_img, anchor="center")
        cv.create_text(8, ch - 8, anchor="sw", fill="#999", font=("TkDefaultFont", 8),
                       text=os.path.basename(path))


# --------------------------------------------------------------------------------------
# Starting from Spyder / Jupyter
# --------------------------------------------------------------------------------------
def in_ipython_kernel():
    """True when running inside an IPython kernel (Spyder's console, Jupyter), where a
    Tkinter window clashes with the console's own event loop."""
    return any(m in sys.modules for m in ("ipykernel", "spyder_kernels")) \
        or "SPYDER_ARGS" in os.environ


def relaunch_outside_ipython(script, args=()):
    """When started from Spyder's console (Run / F5) or Jupyter, start `script` again as
    a separate Python process, so its window runs on its own and the console stays free.
    Returns True if it relaunched (the caller should then stop), False otherwise.
    Nothing to configure in Spyder: no 'Execute in an external system terminal'."""
    if os.environ.get("ISX_NO_RELAUNCH") or not in_ipython_kernel():
        return False
    script = os.path.abspath(script)
    env = dict(os.environ, ISX_NO_RELAUNCH="1")
    for k in ("MPLBACKEND", "SPYDER_ARGS"):
        env.pop(k, None)
    log_path = os.path.join(tempfile.gettempdir(),
                            os.path.splitext(os.path.basename(script))[0] + "_window.log")
    kw = {}
    if sys.platform.startswith("win"):
        kw["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000)
    else:
        kw["start_new_session"] = True
    with open(log_path, "w") as log:
        proc = subprocess.Popen([sys.executable, script, *args], stdin=subprocess.DEVNULL,
                                stdout=log, stderr=subprocess.STDOUT, env=env,
                                cwd=os.path.dirname(script), **kw)
    try:
        code = proc.wait(timeout=3)          # a quick exit means it failed to start
    except subprocess.TimeoutExpired:
        print(f"{os.path.basename(script)} opened in its own window "
              f"(messages go to {log_path}).")
        return True
    with open(log_path) as f:
        msg = f.read().strip()
    print(f"{os.path.basename(script)} could not open its window (exit code {code}):\n{msg}")
    return True


def main():
    root = tk.Tk()
    root.withdraw()
    v = ResultsViewer(root, sys.argv[1] if len(sys.argv) > 1 else None)
    v.protocol("WM_DELETE_WINDOW", root.destroy)
    if len(sys.argv) <= 1:
        root.after(100, v.open_folder)
    root.mainloop()


if __name__ == "__main__":
    if not relaunch_outside_ipython(__file__, sys.argv[1:]):
        main()
