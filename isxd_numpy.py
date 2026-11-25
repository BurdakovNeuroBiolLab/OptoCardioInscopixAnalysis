"""
Pure NumPy/SciPy reader, writer and processing backend for Inscopix .isxd files.

Use it when the IDPS `isx` Python API is not available. It exposes functions with the
same names and arguments as `isx` 1.9.x (preprocess, spatial_filter, project_movie,
motion_correct, dff, pca_ica, event_detection, auto_accept_reject,
export_cell_set_to_csv_tiff, export_event_set_to_csv), so isx_batch_process.py can use
either backend.

.isxd binary layout (as written by IDPS):
    [ binary data ][ JSON header ][ 0x00 ][ uint64 little-endian: JSON length ]
  Movie:    frames one after another, row-major, height x width pixels of dataType
            (0 = uint16, 1 = float32, 2 = uint8). nVista 3 raw files
            ("hasFrameHeaderFooter": true) carry a header before and a footer after each
            frame (2 sensor rows each, i.e. 4 image rows when 2x binned to 640 px); their
            size is measured from the file size.
  Image:    a movie with a single frame.
  Cell set: per cell, its image (height*width float32) then its trace (numTimes float32).

The algorithms follow what IDPS documents for each step. They are re-implementations, not
Inscopix's code, so results are close to but not identical with IDPS (cell sets differ most,
because ICA has a random start). Event sets written here use this module's own layout and
are not readable by IDPS; export them to CSV instead. CNMFe is not implemented.

Reading a file yourself:
    import isxd_numpy as ix
    mov = ix.IsxdMovie("rec.isxd")
    print(mov.shape, mov.period)          # (frames, height, width), seconds per frame
    frames = mov.get_frames(0, 100)       # float32 array (100, height, width)
    cs = ix.IsxdCellSet("rec-PCA-ICA.isxd")
    trace, image = cs.get_trace(0), cs.get_image(0)
"""

import copy
import csv
import json
import logging
import os
import struct

import numpy as np
from scipy import ndimage, signal

log = logging.getLogger("isx_batch")

__version__ = "numpy-backend-1.1"  # bump when results change: forces recomputation

DTYPES = {0: np.uint16, 1: np.float32, 2: np.uint8}
TYPE_MOVIE, TYPE_CELLSET, TYPE_IMAGE, TYPE_EVENTS = 0, 1, 4, 5
STATUS = {0: "accepted", 1: "undecided", 2: "rejected"}
CHUNK = 500  # frames processed at once


# ======================================================================================
# Reading / writing
# ======================================================================================
def read_header(path):
    """Return (header dict, byte offset where the binary data ends)."""
    with open(path, "rb") as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(fsize - 8)
        n = struct.unpack("<Q", f.read(8))[0]
        for extra in (1, 0, 2):  # JSON is normally followed by one null byte
            start = fsize - 8 - n - extra
            if start < 0:
                continue
            f.seek(start)
            raw = f.read(n + extra).strip(b"\x00")
            try:
                return json.loads(raw.decode("utf-8")), start
            except (UnicodeDecodeError, json.JSONDecodeError):
                continue
    raise ValueError(f"{path}: could not parse the .isxd JSON header")


def write_header(f, header):
    raw = json.dumps(header, indent=4).encode("utf-8")
    f.write(raw + b"\x00" + struct.pack("<Q", len(raw)))


class _IsxdBase:
    def __init__(self, path):
        self.path = path
        self.header, self.data_end = read_header(path)
        sp = self.header["spacingInfo"]["numPixels"]
        self.width, self.height = int(sp["x"]), int(sp["y"])
        ti = self.header["timingInfo"]
        self.num_times = int(ti["numTimes"])
        self.period = ti["period"]["num"] / ti["period"]["den"]
        self.dropped = [int(d) for d in ti.get("dropped", [])]

    @property
    def start_time(self):
        """Recording start, seconds since the Unix epoch."""
        s = self.header["timingInfo"]["start"]["secsSinceEpoch"]
        return s["num"] / s["den"]

    @property
    def timestamps(self):
        return np.arange(self.num_times) * self.period


class IsxdMovie(_IsxdBase):
    """Lazy, memory-mapped .isxd movie (or image). Dropped frames are filled with the
    nearest earlier recorded frame so downstream steps see a continuous movie."""

    def __init__(self, path):
        super().__init__(path)
        h = self.header
        self.dtype = DTYPES[int(h.get("dataType", 0))]
        bpp = np.dtype(self.dtype).itemsize
        img = self.height * self.width
        nt, nd = self.num_times, len(self.dropped)
        has_hf = bool(h.get("hasFrameHeaderFooter", False))

        # Frame layout from the file size, not from assumptions: every stored frame is
        # [header][height x width image][footer]. nVista 3 headers/footers are a fixed
        # number of pixels of the full sensor width (e.g. 2 x 1280 px), which is a
        # different number of image rows for binned / cropped movies, so the size must
        # be measured. Dropped frames may or may not be stored.
        layouts = []
        for n in sorted({nt, nt - nd}, reverse=True):
            if n > 0 and self.data_end % (n * bpp) == 0:
                fpx = self.data_end // (n * bpp)
                if fpx >= img:
                    layouts.append((n, fpx))
        if not layouts:
            raise ValueError(
                f"{os.path.basename(path)}: cannot work out the frame layout: "
                f"{self.data_end} data bytes for {nt} frames ({nd} dropped) of "
                f"{self.width} x {self.height} {np.dtype(self.dtype).name}. "
                "Please report this file's header (isxd_numpy.read_header).")
        # Several layouts can fit the file size (dropped frames stored or not). Keep the
        # one whose frames show the same field of view from start to end: a wrong frame
        # size makes every frame start a little further off, so they drift apart.
        if len(layouts) > 1:
            layouts.sort(key=lambda L: -self._layout_score(*L))
        stored, fpx = layouts[0]
        extra = fpx - img
        self._fpx, self._img = fpx, img
        self._mm = np.memmap(path, dtype=self.dtype, mode="r", shape=(stored, fpx))
        self._hdr = self._find_header_size(self._mm, extra) if extra else 0
        if extra:
            log.info("  %s: frame header %d px + footer %d px (%.1f + %.1f image rows)",
                     os.path.basename(path), self._hdr, extra - self._hdr,
                     self._hdr / self.width, (extra - self._hdr) / self.width)

        if stored == nt:  # dropped frames stored as blanks: replace them too
            idx = np.arange(nt)
            idx[[d for d in self.dropped if 0 <= d < nt]] = -1
        else:  # dropped frames not stored
            keep = np.setdiff1d(np.arange(nt), self.dropped)
            idx = np.full(nt, -1)
            idx[keep] = np.arange(len(keep))
        # Fill dropped frames with the nearest earlier stored frame.
        for t in range(nt):
            if idx[t] < 0:
                idx[t] = idx[t - 1] if t > 0 else -1
        if idx[0] < 0:
            idx[idx < 0] = idx[idx >= 0][0]
        self._idx = idx

    def _sample(self, mm, k=8):
        idx = np.linspace(0, mm.shape[0] - 1, min(k, mm.shape[0])).astype(int)
        return np.asarray(mm[idx], dtype=np.float64)

    def _find_header_size(self, mm, extra):
        """Split the extra pixels per frame into header and footer. Try every whole-row
        split; a wrong split puts a metadata/image boundary inside the image, which
        shows as a row-to-row jump far larger than the image's own row changes. Keep the
        split whose image has the smallest such jump. Without whole rows, equal halves."""
        w, h, img = self.width, self.height, self.width * self.height
        if extra % w:
            return extra // 2
        buf = self._sample(mm).mean(0)
        best, best_score = extra // 2, np.inf
        for hs in range(0, extra + 1, w):
            im = buf[hs:hs + img].reshape(h, w)
            d = np.abs(np.diff(im, axis=0)).mean(1)
            score = d.max() / (np.median(d) + 1e-9) + (1e-6 if hs != extra // 2 else 0)
            if score < best_score:
                best, best_score = hs, score
        return best

    def _layout_score(self, n, fpx):
        """Mean correlation of frames across the movie with its first frame."""
        mm = np.memmap(self.path, dtype=self.dtype, mode="r", shape=(n, fpx))
        extra = fpx - self.height * self.width
        hs = self._find_header_size(mm, extra) if extra else 0
        fr = self._sample(mm)[:, hs:hs + self.height * self.width]
        fr -= fr.mean(1, keepdims=True)
        ref = fr[0]
        c = [(f @ ref) / (np.linalg.norm(f) * np.linalg.norm(ref) + 1e-12) for f in fr[1:]]
        return float(np.mean(c)) if c else 0.0

    @property
    def shape(self):
        return (self.num_times, self.height, self.width)

    def get_frames(self, start=0, stop=None):
        stop = self.num_times if stop is None else min(stop, self.num_times)
        sel = self._mm[self._idx[start:stop], self._hdr:self._hdr + self._img]
        return np.asarray(sel, dtype=np.float32).reshape(-1, self.height, self.width)

    def get_frame(self, t):
        return self.get_frames(t, t + 1)[0]

    def chunks(self, size=CHUNK):
        for s in range(0, self.num_times, size):
            yield s, self.get_frames(s, s + size)


class IsxdWriter:
    """Writes float32 frames, then the JSON header on close()."""

    def __init__(self, path, header):
        self.path, self.header = path, header
        self.f = open(path, "wb")
        self.n = 0

    def write(self, frames):
        frames = np.ascontiguousarray(frames, dtype=np.float32)
        self.f.write(frames.tobytes())
        self.n += frames.shape[0] if frames.ndim == 3 else 1

    def close(self):
        write_header(self.f, self.header)
        self.f.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if exc[0] is None:
            self.close()
        else:
            self.f.close()
            os.remove(self.path)


class IsxdCellSet(_IsxdBase):
    def __init__(self, path):
        super().__init__(path)
        self.names = list(self.header.get("CellNames", []))
        self.num_cells = len(self.names)
        block = self.height * self.width + self.num_times
        self._mm = np.memmap(path, dtype=np.float32, mode="r",
                             shape=(self.num_cells, block))

    def get_image(self, i):
        return np.array(self._mm[i, :self.height * self.width]).reshape(self.height, self.width)

    def get_trace(self, i):
        return np.array(self._mm[i, self.height * self.width:])

    def get_status(self, i):
        return STATUS[int(self.header.get("CellStatuses", [1] * self.num_cells)[i])]

    def set_statuses(self, statuses):
        """Rewrite the header with new statuses ('accepted'/'undecided'/'rejected')."""
        inv = {v: k for k, v in STATUS.items()}
        self.header["CellStatuses"] = [inv[s] for s in statuses]
        del self._mm
        with open(self.path, "r+b") as f:
            f.seek(self.data_end)
            f.truncate()
            write_header(f, self.header)
        self.__init__(self.path)


class IsxdEventSet(_IsxdBase):
    """Event set in this module's own layout: float64 (time_s, amplitude) pairs per cell."""

    def __init__(self, path):
        super().__init__(path)
        self.names = self.header["CellNames"]
        counts = self.header["EventCounts"]
        data = np.fromfile(path, dtype=np.float64, count=2 * sum(counts)).reshape(-1, 2)
        offs = np.concatenate([[0], np.cumsum(counts)])
        self.events = [data[offs[i]:offs[i + 1]] for i in range(len(counts))]


# ======================================================================================
# Header helpers
# ======================================================================================
def _frac(v):
    return v["num"] / v["den"] if isinstance(v, dict) else float(v)


def _derive_header(h, dtype_code=1, width=None, height=None, spatial_factor=1,
                   crop_tl=(0, 0), num_times=None, temporal_factor=1, ds_type=None):
    h = copy.deepcopy(h)
    h["dataType"] = dtype_code
    h["hasFrameHeaderFooter"] = False
    if ds_type is not None:
        h["type"] = ds_type
    sp = h["spacingInfo"]
    if width is not None:
        sp["numPixels"] = {"x": int(width), "y": int(height)}
    psz = sp.get("pixelSize")
    tl = sp.get("topLeft")
    if psz and tl and crop_tl != (0, 0):
        for ax, off in (("y", crop_tl[0]), ("x", crop_tl[1])):
            tl[ax] = {"num": int(round((_frac(tl[ax]) + off * _frac(psz[ax])) * 1000)),
                      "den": 1000}
    if psz and spatial_factor != 1:
        for ax in ("x", "y"):
            v = psz[ax]
            psz[ax] = ({"num": v["num"] * spatial_factor, "den": v["den"]}
                       if isinstance(v, dict) else v * spatial_factor)
    ti = h["timingInfo"]
    if temporal_factor != 1:
        ti["period"] = {"num": ti["period"]["num"] * temporal_factor, "den": ti["period"]["den"]}
        ti["dropped"] = sorted({d // temporal_factor for d in ti.get("dropped", [])})
        ti["cropped"], ti["blank"] = [], []
    if num_times is not None:
        ti["numTimes"] = int(num_times)
        for k in ("dropped", "cropped", "blank"):
            if k in ti:
                ti[k] = [d for d in ti[k] if (d[0] if isinstance(d, list) else d) < num_times]
    h["producer"] = {"name": "isxd_numpy", "version": __version__}
    return h


def _as_list(x):
    return [x] if isinstance(x, str) else list(x)


def _check_out(paths):
    for p in paths:
        if os.path.exists(p):
            raise FileExistsError(f"Output file already exists: {p}")


# ======================================================================================
# Image filters
# ======================================================================================
def _bandpass_kernel(h, w, low, high):
    """Frequency response of the IDPS spatial band-pass: a Gaussian low-pass at `high`
    (removes pixel noise) minus a Gaussian low-pass at `low` (removes background).
    Cutoffs are in cycles/pixel."""
    fy = np.fft.fftfreq(h)[:, None]
    fx = np.fft.rfftfreq(w)[None, :]
    f2 = fx ** 2 + fy ** 2
    k = np.ones_like(f2) if not high else np.exp(-f2 / (2 * high ** 2))
    if low:
        k = k - np.exp(-f2 / (2 * low ** 2))
    return k


PAD = 16  # reflect-padding so the FFT filter does not wrap image edges around


def _bandpass(frames, low, high):
    h, w = frames.shape[-2:]
    p = min(PAD, h - 1, w - 1)
    x = np.pad(frames, ((0, 0), (p, p), (p, p)), mode="reflect")
    k = _bandpass_kernel(x.shape[1], x.shape[2], low, high)
    y = np.fft.irfft2(np.fft.rfft2(x) * k, s=x.shape[-2:])
    return y[:, p:p + h, p:p + w].astype(np.float32)


# ======================================================================================
# Processing steps (isx 1.9.x compatible signatures)
# ======================================================================================
def preprocess(input_movie_files, output_movie_files, temporal_downsample_factor=1,
               spatial_downsample_factor=1, crop_rect=None, fix_defective_pixels=True,
               trim_early_frames=True):
    ins, outs = _as_list(input_movie_files), _as_list(output_movie_files)
    _check_out(outs)
    tf, sf = int(temporal_downsample_factor), int(spatial_downsample_factor)
    if trim_early_frames:
        log.debug("trim_early_frames has no effect in the NumPy backend")
    for src, dst in zip(ins, outs):
        mov = IsxdMovie(src)
        top, left, bottom, right = crop_rect if crop_rect else (0, 0, mov.height - 1, mov.width - 1)
        ch, cw = bottom - top + 1, right - left + 1
        oh, ow = ch // sf, cw // sf

        bad = None
        if fix_defective_pixels:
            # Defective (hot/dead) pixels: deviate from their 3x3 median by > 10 noise SDs
            # in (the median over) the first frames, noise measured within single frames.
            fr = mov.get_frames(0, min(50, mov.num_times))
            devs = fr - ndimage.median_filter(fr, size=(1, 3, 3))
            noise = 1.4826 * np.median(np.abs(devs - np.median(devs))) + 1e-6
            bad = np.abs(np.median(devs, axis=0)) > 10 * noise
            log.info("  %s: %d defective pixel(s)", os.path.basename(src), int(bad.sum()))
            if not bad.any():
                bad = None

        n_out = mov.num_times // tf
        hdr = _derive_header(mov.header, width=ow, height=oh, spatial_factor=sf,
                             crop_tl=(top, left), num_times=n_out, temporal_factor=tf)
        with IsxdWriter(dst, hdr) as w:
            step = max(tf, (CHUNK // tf) * tf)
            for s in range(0, n_out * tf, step):
                fr = mov.get_frames(s, min(s + step, n_out * tf))
                if bad is not None:
                    med = ndimage.median_filter(fr, size=(1, 3, 3))
                    fr[:, bad] = med[:, bad]
                fr = fr[:, top:top + oh * sf, left:left + ow * sf]
                fr = fr.reshape(-1, oh, sf, ow, sf).mean(axis=(2, 4))
                fr = fr.reshape(-1, tf, oh, ow).mean(axis=1)
                w.write(fr)


def spatial_filter(input_movie_files, output_movie_files, low_cutoff=0.005, high_cutoff=0.500,
                   retain_mean=False, subtract_global_minimum=True):
    ins, outs = _as_list(input_movie_files), _as_list(output_movie_files)
    _check_out(outs)
    for src, dst in zip(ins, outs):
        mov = IsxdMovie(src)
        gmin = np.inf
        with IsxdWriter(dst, _derive_header(mov.header)) as w:
            for _, fr in mov.chunks():
                out = _bandpass(fr, low_cutoff, high_cutoff)
                if retain_mean:
                    out += fr.mean(axis=(1, 2), keepdims=True)
                gmin = min(gmin, float(out.min()))
                w.write(out)
        if subtract_global_minimum:
            out = IsxdMovie(dst)
            mm = np.memmap(dst, dtype=np.float32, mode="r+", shape=out.shape)
            for s in range(0, out.num_times, CHUNK):
                mm[s:s + CHUNK] -= gmin
            mm.flush()
            del mm


def project_movie(input_movie_files, output_image_file, stat_type="mean"):
    ins = _as_list(input_movie_files)
    _check_out([output_image_file])
    n, acc, acc2 = 0, None, None
    for src in ins:
        mov = IsxdMovie(src)
        for _, fr in mov.chunks():
            fr = fr.astype(np.float64)
            if acc is None:
                acc = np.zeros(fr.shape[1:]) if stat_type in ("mean", "standard_deviation") \
                    else fr[0].copy()
                acc2 = np.zeros(fr.shape[1:])
            if stat_type == "mean":
                acc += fr.sum(0)
            elif stat_type == "max":
                acc = np.maximum(acc, fr.max(0))
            elif stat_type == "min":
                acc = np.minimum(acc, fr.min(0))
            elif stat_type == "standard_deviation":
                acc += fr.sum(0)
                acc2 += (fr ** 2).sum(0)
            else:
                raise ValueError(f"stat_type {stat_type!r}")
            n += fr.shape[0]
    if stat_type == "mean":
        img = acc / n
    elif stat_type == "standard_deviation":
        img = np.sqrt(np.maximum(acc2 / n - (acc / n) ** 2, 0))
    else:
        img = acc
    hdr = _derive_header(IsxdMovie(ins[0]).header, num_times=1, ds_type=TYPE_IMAGE)
    with IsxdWriter(output_image_file, hdr) as w:
        w.write(img[None].astype(np.float32))


def _xcorr_shift(ref_f, frame, max_t):
    """(dy, dx) to apply to `frame` to align it to the reference: cross-correlation of the
    mean-subtracted images (like IDPS), parabolic sub-pixel refinement. `ref_f` is the FFT
    of the mean-subtracted reference."""
    F = np.fft.fft2(frame - frame.mean())
    r = np.fft.ifft2(ref_f * np.conj(F)).real
    h, w = r.shape
    # restrict to +/- max_t
    ys = np.r_[0:max_t + 1, h - max_t:h]
    xs = np.r_[0:max_t + 1, w - max_t:w]
    sub = r[np.ix_(ys, xs)]
    iy, ix = np.unravel_index(np.argmax(sub), sub.shape)
    py, px = ys[iy], xs[ix]

    def refine(c_m, c_0, c_p):
        d = c_m - 2 * c_0 + c_p
        return 0.0 if d == 0 else 0.5 * (c_m - c_p) / d

    dy = py + refine(r[(py - 1) % h, px], r[py, px], r[(py + 1) % h, px])
    dx = px + refine(r[py, (px - 1) % w], r[py, px], r[py, (px + 1) % w])
    dy = dy - h if dy > h / 2 else dy
    dx = dx - w if dx > w / 2 else dx
    return float(np.clip(dy, -max_t, max_t)), float(np.clip(dx, -max_t, max_t))


def motion_correct(input_movie_files, output_movie_files, max_translation=20,
                   low_bandpass_cutoff=0.004, high_bandpass_cutoff=0.016, roi=None,
                   reference_segment_index=0, reference_frame_index=0, reference_file_name="",
                   global_registration_weight=1.0, output_translation_files=None,
                   output_crop_rect_file=None):
    ins, outs = _as_list(input_movie_files), _as_list(output_movie_files)
    _check_out(outs)
    if roi is not None:
        log.warning("motion_correct: roi is ignored by the NumPy backend")
    first = IsxdMovie(ins[0])
    if reference_file_name:
        ref = IsxdMovie(reference_file_name).get_frame(0)
    else:
        ref = IsxdMovie(ins[reference_segment_index]).get_frame(reference_frame_index)
    h, w = ref.shape
    use_bp = bool(low_bandpass_cutoff) or bool(high_bandpass_cutoff)
    prep = ((lambda x: _bandpass(x, low_bandpass_cutoff, high_bandpass_cutoff))
            if use_bp else (lambda x: x))
    rp = prep(ref[None])[0]
    ref_f = np.fft.fft2(rp - rp.mean())
    max_t = int(max_translation)

    all_shifts = []
    for i, (src, dst) in enumerate(zip(ins, outs)):
        mov = IsxdMovie(src)
        shifts = []
        with IsxdWriter(dst, _derive_header(mov.header)) as wr:
            for _, fr in mov.chunks():
                pf = prep(fr)
                out = np.empty_like(fr)
                for j in range(fr.shape[0]):
                    dy, dx = _xcorr_shift(ref_f, pf[j], max_t)
                    shifts.append((dy, dx))
                    out[j] = ndimage.shift(fr[j], (dy, dx), order=1, mode="constant", cval=0)
                wr.write(out)
        shifts = np.array(shifts)
        all_shifts.append(shifts)
        at_limit = np.mean(np.any(np.abs(shifts) >= max_t - 0.5, axis=1))
        log.info("  %s: shift x %.2f +/- %.2f px, y %.2f +/- %.2f px",
                 os.path.basename(src), shifts[:, 1].mean(), shifts[:, 1].std(),
                 shifts[:, 0].mean(), shifts[:, 0].std())
        if at_limit > 0.01:
            log.warning("  %s: %.0f%% of frames hit the maximum translation (%d px). "
                        "Motion correction is probably failing; check the movie and "
                        "the translations file.", os.path.basename(src), 100 * at_limit, max_t)
        if output_translation_files:
            with open(output_translation_files[i], "w", newline="") as f:
                cw = csv.writer(f)
                cw.writerow(["translationX", "translationY"])
                cw.writerows(np.round(shifts[:, ::-1], 4).tolist())
    if output_crop_rect_file:
        s = np.concatenate(all_shifts)
        top = int(np.ceil(max(s[:, 0].max(), 0)))
        left = int(np.ceil(max(s[:, 1].max(), 0)))
        bottom = first.height - 1 - int(np.ceil(max(-s[:, 0].min(), 0)))
        right = first.width - 1 - int(np.ceil(max(-s[:, 1].min(), 0)))
        with open(output_crop_rect_file, "w") as f:
            f.write("top,left,bottom,right\n")
            f.write(f"{top},{left},{bottom},{right}\n")


def dff(input_movie_files, output_movie_files, f0_type="mean"):
    """(F - F0) / F0 with F0 the per-pixel mean (or min) over all input segments."""
    ins, outs = _as_list(input_movie_files), _as_list(output_movie_files)
    _check_out(outs)
    tmp = outs[0] + ".f0.tmp.isxd"
    if os.path.exists(tmp):
        os.remove(tmp)
    project_movie(ins, tmp, stat_type="mean" if f0_type == "mean" else "min")
    f0 = IsxdMovie(tmp).get_frame(0)
    os.remove(tmp)
    safe = np.where(np.abs(f0) > 1e-6, f0, np.nan)
    for src, dst in zip(ins, outs):
        mov = IsxdMovie(src)
        with IsxdWriter(dst, _derive_header(mov.header)) as w:
            for _, fr in mov.chunks():
                w.write(np.nan_to_num((fr - safe) / safe))


def _fastica(X, n_comp, max_iter, tol, seed=0):
    """Symmetric FastICA on rows of X (already decorrelated, unit variance), using the
    skewness contrast g(u)=u^2 as in Mukamel et al. 2009 (CellSort)."""
    rng = np.random.default_rng(seed)
    W = np.linalg.qr(rng.standard_normal((n_comp, n_comp)))[0]
    N = X.shape[1]
    for it in range(max_iter):
        U = W @ X
        W_new = (U ** 2) @ X.T / N - np.diag((2 * U).mean(1)) @ W
        s, V = np.linalg.eigh(W_new @ W_new.T)
        W_new = V @ np.diag(1 / np.sqrt(np.maximum(s, 1e-12))) @ V.T @ W_new
        conv = np.max(np.abs(np.abs(np.diag(W_new @ W.T)) - 1))
        W = W_new
        if conv < tol:
            break
    return W, it + 1


def pca_ica(input_movie_files, output_cell_set_files, num_pcs, num_ics=120,
            unmix_type="spatial", ica_temporal_weight=0, max_iterations=100,
            convergence_threshold=1e-5, block_size=1000, auto_estimate_num_ics=False,
            average_cell_diameter=10):
    ins, outs = _as_list(input_movie_files), _as_list(output_cell_set_files)
    _check_out(outs)
    movs = [IsxdMovie(p) for p in ins]
    h, w = movs[0].height, movs[0].width
    lens = [m.num_times for m in movs]
    M = np.concatenate([m.get_frames().reshape(m.num_times, -1) for m in movs])  # T x P
    T, P = M.shape
    M -= M.mean(0)

    num_pcs = int(min(num_pcs, T - 1, P - 1))
    # PCA through the smaller covariance matrix
    if T <= P:
        C = (M @ M.T).astype(np.float64)
        s, U = np.linalg.eigh(C)
        o = np.argsort(s)[::-1][:num_pcs]
        s, U = np.maximum(s[o], 1e-12), U[:, o]           # temporal PCs (T x k)
        V = (U.T @ M) / np.sqrt(s)[:, None]                # spatial PCs (k x P)
    else:
        C = (M.T @ M).astype(np.float64)
        s, Vt = np.linalg.eigh(C)
        o = np.argsort(s)[::-1][:num_pcs]
        s, V = np.maximum(s[o], 1e-12), Vt[:, o].T
        U = (M @ V.T) / np.sqrt(s)[None, :]

    if auto_estimate_num_ics:
        # Keep PCs above the Marchenko-Pastur noise edge.
        noise = np.median(s)
        num_ics = int(max(1, np.sum(s > noise * (1 + np.sqrt(P / T)) ** 2)))
        log.info("  auto-estimated number of ICs: %d", num_ics)
    num_ics = int(min(num_ics, num_pcs))

    sig = (U.T / U.T.std(1, keepdims=True))                      # k x T
    filt = (V / V.std(1, keepdims=True))                          # k x P
    mu = float(ica_temporal_weight) if unmix_type == "both" else (
        1.0 if unmix_type == "temporal" else 0.0)
    X = np.concatenate([mu * sig, (1 - mu) * filt], axis=1) if 0 < mu < 1 else (
        sig if mu == 1 else filt)
    X = X - X.mean(1, keepdims=True)
    # whiten (rows already ~orthogonal; this makes it exact) and keep num_ics dims
    cov = X @ X.T / X.shape[1]
    e, E = np.linalg.eigh(cov)
    o = np.argsort(e)[::-1][:num_ics]
    Wh = (E[:, o] / np.sqrt(np.maximum(e[o], 1e-12))).T          # num_ics x k
    Wica, n_it = _fastica(Wh @ X, num_ics, max_iterations, convergence_threshold)
    log.info("  ICA: %d components, %d iterations", num_ics, n_it)
    A = Wica @ Wh                                                 # num_ics x k

    images = A @ filt                                             # num_ics x P
    # make cells positive (positive skew), normalise images
    sk = ((images - images.mean(1, keepdims=True)) ** 3).mean(1)
    images *= np.sign(sk)[:, None]
    images /= np.abs(images).max(1, keepdims=True) + 1e-12
    # traces: image-weighted average of the movie, in the input's units (e.g. dF/F)
    wts = np.clip(images, 0, None)
    wts /= wts.sum(1, keepdims=True) + 1e-12

    traces = np.concatenate([m.get_frames().reshape(m.num_times, -1) for m in movs]) @ wts.T

    order = np.argsort(-traces.std(0))   # most active first
    images, traces = images[order], traces[:, order]
    n = images.shape[0]
    names = [f"C{i:0{max(2, len(str(n - 1)))}d}" for i in range(n)]
    off = 0
    for mov, dst, L in zip(movs, outs, lens):
        hdr = _derive_header(mov.header, ds_type=TYPE_CELLSET)
        hdr.update({"CellNames": names, "CellStatuses": [1] * n,
                    "CellColors": [[255, 255, 255]] * n, "CellActivity": [[True]] * n})
        with IsxdWriter(dst, hdr) as wr:
            for c in range(n):
                wr.write(images[c].reshape(1, h, w))
                wr.f.write(traces[off:off + L, c].astype(np.float32).tobytes())
        off += L


def cnmfe(*args, **kwargs):
    raise NotImplementedError("CNMFe is not available in the NumPy backend; use "
                              "--method pca-ica, or install the IDPS isx API.")


def _detect_events(trace, period, threshold, tau, time_ref, ignore_negative):
    med = np.median(trace)
    mad = np.median(np.abs(trace - med)) + 1e-12
    x = trace - med
    dist = max(1, int(round(tau / period)))
    peaks, _ = signal.find_peaks(x, height=threshold * mad, prominence=threshold * mad,
                                 distance=dist)
    if ignore_negative:
        # a transient must decay for at least tau after its peak
        peaks = [p for p in peaks if p + dist >= len(x) or x[p:p + dist + 1].min() < x[p]]
    ev = []
    for p in peaks:
        t = p
        if time_ref in ("beginning", "mid_rise"):
            while t > 0 and x[t - 1] < x[t] and x[t - 1] > mad:
                t -= 1
            if time_ref == "mid_rise":  # first frame above half the rise
                half = x[t] + 0.5 * (x[p] - x[t])
                t = t + int(np.argmax(x[t:p + 1] >= half))
        ev.append((t * period, float(x[p])))
    return np.array(ev, dtype=np.float64).reshape(-1, 2), mad


def event_detection(input_cell_set_files, output_event_set_files, threshold=5, tau=0.2,
                    event_time_ref="beginning", ignore_negative_transients=True,
                    accepted_cells_only=False):
    ins, outs = _as_list(input_cell_set_files), _as_list(output_event_set_files)
    _check_out(outs)
    for src, dst in zip(ins, outs):
        cs = IsxdCellSet(src)
        evs, snrs, names = [], [], []
        for i in range(cs.num_cells):
            if accepted_cells_only and cs.get_status(i) != "accepted":
                continue
            tr = cs.get_trace(i)
            ev, mad = _detect_events(tr, cs.period, threshold, tau, event_time_ref,
                                     ignore_negative_transients)
            sd = 1.4826 * mad
            snrs.append(float(ev[:, 1].mean() / sd) if len(ev) else 0.0)
            evs.append(ev)
            names.append(cs.names[i])
        hdr = copy.deepcopy(cs.header)
        for k in ("CellStatuses", "CellColors", "CellActivity"):
            hdr.pop(k, None)
        hdr.update({"type": TYPE_EVENTS, "CellNames": names,
                    "EventCounts": [len(e) for e in evs], "SNR": snrs,
                    "producer": {"name": "isxd_numpy", "version": __version__}})
        with open(dst, "wb") as f:
            for e in evs:
                f.write(e.astype(np.float64).tobytes())
            write_header(f, hdr)


def _cell_metrics(img):
    """# Comps, Cell Size (pixels of the largest component), centroid."""
    mask = img >= 0.5 * img.max()
    lab, n = ndimage.label(mask)
    if n == 0:
        return 0, 0, (np.nan, np.nan)
    sizes = ndimage.sum(mask, lab, range(1, n + 1))
    big = np.argmax(sizes) + 1
    cy, cx = ndimage.center_of_mass(lab == big)
    return n, int(sizes.max()), (cx, cy)


def auto_accept_reject(input_cell_set_files, input_event_set_files, filters=None):
    """Cells are judged once over all segments (a series shares one status per cell)."""
    import operator
    ops = {"<": operator.lt, ">": operator.gt, "=": operator.eq, "==": operator.eq,
           "<=": operator.le, ">=": operator.ge}
    css = [IsxdCellSet(p) for p in _as_list(input_cell_set_files)]
    ess = [IsxdEventSet(p) for p in _as_list(input_event_set_files)]
    dur = sum(cs.num_times * cs.period for cs in css)
    snr, count = {}, {}
    for es in ess:
        for n, e, v in zip(es.names, es.events, es.header["SNR"]):
            snr.setdefault(n, []).append((v, len(e)))
            count[n] = count.get(n, 0) + len(e)
    statuses = []
    cs0 = css[0]
    for i, name in enumerate(cs0.names):
        comps, size, _ = _cell_metrics(cs0.get_image(i))
        pairs = snr.get(name, [])
        n_ev = sum(c for _, c in pairs)
        m = {"# Comps": comps, "Cell Size": size,
             "SNR": sum(v * c for v, c in pairs) / n_ev if n_ev else 0.0,
             "Event Rate": count.get(name, 0) / dur}
        ok = all(ops[op](m[key], val) for key, op, val in (filters or []))
        statuses.append("accepted" if ok else "rejected")
    for cs in css:
        cs.set_statuses(statuses)
    log.info("  %d/%d cells accepted", statuses.count("accepted"), len(statuses))


def write_tiff(path, img):
    """Minimal uncompressed single-page grayscale TIFF (uint8, uint16 or float32);
    no extra packages needed, readable by MATLAB imread, ImageJ, Python."""
    img = np.ascontiguousarray(img)
    if img.dtype not in (np.uint8, np.uint16, np.float32):
        img = img.astype(np.float32)
    h, w = img.shape
    bits = img.dtype.itemsize * 8
    fmt = 3 if img.dtype == np.float32 else 1           # SampleFormat: 3 float, 1 uint
    data = img.astype(img.dtype.newbyteorder("<")).tobytes()
    tags = [(256, 4, w), (257, 4, h), (258, 3, bits), (259, 3, 1), (262, 3, 1),
            (273, 4, 0), (277, 3, 1), (278, 4, h), (279, 4, len(data)), (284, 3, 1),
            (339, 3, fmt)]
    ifd_off = 8
    data_off = ifd_off + 2 + 12 * len(tags) + 4
    out = bytearray(b"II*\x00" + struct.pack("<I", ifd_off))
    out += struct.pack("<H", len(tags))
    for tag, typ, val in tags:
        if tag == 273:
            val = data_off
        if typ == 3:
            out += struct.pack("<HHIHH", tag, typ, 1, val, 0)
        else:
            out += struct.pack("<HHII", tag, typ, 1, val)
    out += struct.pack("<I", 0)
    with open(path, "wb") as f:
        f.write(bytes(out) + data)


def cell_outline_masks(images):
    """Per cell: boolean mask of the pixels >= 50% of the cell image maximum."""
    return [(img >= 0.5 * img.max()) if img.max() > 0 else np.zeros(img.shape, bool)
            for img in images]


def write_cell_map(path, images, max_width=800):
    """IDPS-style cell map: maximum projection of the (non-negative) cell images,
    enlarged by an integer factor, with every cell's outline drawn in white
    (= 65535, the image maximum) on a uint16 background kept below that value."""
    images = [np.nan_to_num(np.asarray(im, dtype=np.float64)) for im in images]
    h, w = images[0].shape
    scale = int(max(1, min(8, max_width // max(w, 1))))
    bg = np.max(np.clip(np.stack(images), 0, None), axis=0) if images else np.zeros((h, w))
    bg = bg / (bg.max() or 1) * 60000
    bg = np.kron(bg, np.ones((scale, scale)))
    out = bg.astype(np.uint16)
    for m in cell_outline_masks(images):
        big = np.kron(m, np.ones((scale, scale), bool))
        edge = big & ~ndimage.binary_erosion(big, border_value=0)
        out[edge] = 65535
    write_tiff(path, out)
    return scale


def _times(objs, time_ref):
    if time_ref == "tsc":
        raise ValueError("time_ref='tsc' is not supported by the NumPy backend")
    t0 = objs[0].start_time
    out = []
    for o in objs:
        base = o.start_time if time_ref == "unix" else o.start_time - t0
        out.append(base + o.timestamps)
    return out


def export_cell_set_to_csv_tiff(input_cell_set_files, output_csv_file, output_tiff_file,
                                time_ref="start", output_props_file=""):
    css = [IsxdCellSet(p) for p in _as_list(input_cell_set_files)]
    cs0 = css[0]
    with open(output_csv_file, "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow([" "] + [f" {n}" for n in cs0.names])
        cw.writerow(["Time(s)/Cell Status"] + [f" {cs0.get_status(i)}"
                                               for i in range(cs0.num_cells)])
        for cs, times in zip(css, _times(css, time_ref)):
            traces = np.stack([cs.get_trace(i) for i in range(cs.num_cells)], 1) \
                if cs.num_cells else np.zeros((cs.num_times, 0))
            for t, row in zip(times, traces):
                cw.writerow([f"{t:.6f}"] + [f"{v:.6g}" for v in row])
    if output_tiff_file:
        base = os.path.splitext(output_tiff_file)[0]
        for i, n in enumerate(cs0.names):
            write_tiff(f"{base}_{n}.tif", cs0.get_image(i).astype(np.float32))
    if output_props_file:
        # Same columns as the IDPS props export. Size = equivalent diameter (pixels)
        # of the largest component of the cell image (>= 50% of its maximum).
        with open(output_props_file, "w", newline="") as f:
            cw = csv.writer(f)
            cw.writerow(["Name", "Status", "ColorR", "ColorG", "ColorB", "CentroidX",
                         "CentroidY", "NumComponents", "Size"])
            for i, n in enumerate(cs0.names):
                comps, area, (cx, cy) = _cell_metrics(cs0.get_image(i))
                diam = 2 * np.sqrt(area / np.pi)
                cw.writerow([n, cs0.get_status(i), 255, 255, 255, f"{cx:.4f}", f"{cy:.4f}",
                             comps, f"{diam:.4f}"])


def export_event_set_to_csv(input_event_set_files, output_csv_file, time_ref="start",
                            output_props_file="", sparse_output=True, write_amplitude=True):
    ess = [IsxdEventSet(p) for p in _as_list(input_event_set_files)]
    starts = _times(ess, time_ref)
    with open(output_csv_file, "w", newline="") as f:
        cw = csv.writer(f)
        if sparse_output:
            cw.writerow(["Time (s)", " Cell Name", " Value"])
            rows = []
            for es, ts in zip(ess, starts):
                for name, ev in zip(es.names, es.events):
                    for t, a in ev:
                        rows.append((ts[0] + t, name, a if write_amplitude else 1))
            for t, n, a in sorted(rows):
                cw.writerow([f"{t:.6f}", f" {n}", f" {a:.6g}"])
        else:  # dense: one row per frame, one column per cell
            names = ess[0].names
            cw.writerow(["Time (s)"] + [f" {n}" for n in names])
            for es, ts in zip(ess, starts):
                table = np.zeros((es.num_times, len(names)))
                for c, ev in enumerate(es.events):
                    for t, a in ev:
                        i = min(int(round(t / es.period)), es.num_times - 1)
                        table[i, c] = a if write_amplitude else 1
                for t, row in zip(ts, table):
                    cw.writerow([f"{t:.6f}"] + [f"{v:.6g}" for v in row])
