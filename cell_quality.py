"""
Cell likelihood score and duplicate-ROI rejection, shared by the viewer (isx_viewer.py)
and the batch pipeline (isx_batch_process.py). NumPy only.

Score (higher = more likely a real cell):
    log(1 + SNR)                     SNR = (99.9th percentile - median) of the trace,
                                     slow drift removed, divided by its frame-to-frame
                                     noise (robust SD)
    - 0.5 x (NumComponents - 1)      ROIs made of several pieces are less cell-like
    - 0.5 x |log(Size / median Size)| unusually small or large ROIs are less cell-like
    + 1 if accepted, - 3 if rejected  the software's earlier decision

Duplicates: two ROIs whose centroids are closer than DUP_DIST_FRAC x the average cell
diameter (mean props Size) are the same cell unless their traces differ substantially
(Pearson correlation below DUP_MIN_CORR). Of each such group the highest-scoring ROI is
kept and the others are rejected as duplicates of it. Cells already rejected never
suppress others.
"""

import csv
import os

import numpy as np

DUP_DIST_FRAC = 0.25     # fraction of the average cell diameter
DUP_MIN_CORR = 0.5       # traces less correlated than this are different cells


def _detrend(x):
    """Remove slow drift (e.g. bleaching) with a running mean over ~1% of the
    recording (50-1000 frames); calcium transients are much shorter and survive."""
    w = int(min(max(x.size // 100, 50), 1000))
    if x.size <= 2 * w:
        return x - np.median(x)
    c = np.cumsum(np.concatenate([[0.0], x]))
    base = (c[w:] - c[:-w]) / w
    pad = w // 2
    base = np.concatenate([np.full(pad, base[0]), base, np.full(x.size - base.size - pad, base[-1])])
    return x - base


def trace_snr(x):
    """Transient size over noise: (99.9th percentile - median) of the drift-free trace
    divided by its frame-to-frame noise (robust SD)."""
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size < 10:
        return 0.0
    x = _detrend(x)
    d = np.diff(x)
    noise = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)
    signal = np.percentile(x, 99.9) - np.median(x)
    return float(max(signal, 0) / max(noise, 1e-12))


def cell_scores(traces, status, size=None, ncomp=None):
    """traces: T x N array; status: N strings; size, ncomp: N numbers or None."""
    n = traces.shape[1]
    snr = np.array([trace_snr(traces[:, j]) for j in range(n)])
    score = np.log1p(snr)
    if ncomp is not None:
        nc = np.nan_to_num(np.asarray(ncomp, float), nan=1.0)
        score -= 0.5 * np.clip(nc - 1, 0, None)
    if size is not None:
        sz = np.asarray(size, float)
        ok = np.isfinite(sz) & (sz > 0)
        if ok.any():
            med = np.median(sz[ok])
            score[ok] -= 0.5 * np.abs(np.log(sz[ok] / med))
    st = [s.lower() for s in status]
    score += np.array([1.0 if s == "accepted" else -3.0 if s == "rejected" else 0.0 for s in st])
    return score, snr


def _z(x):
    x = np.nan_to_num(np.asarray(x, np.float64) - np.nanmean(x))
    sd = x.std()
    return x / sd if sd > 0 else x


def find_duplicates(cx, cy, size, traces, score, status,
                    dist_frac=DUP_DIST_FRAC, min_corr=DUP_MIN_CORR):
    """Returns {duplicate index: index of the kept cell}."""
    cx, cy, size = (np.asarray(a, float) for a in (cx, cy, size))
    ok_size = size[np.isfinite(size) & (size > 0)]
    if not ok_size.size:
        return {}
    thr = dist_frac * ok_size.mean()
    kept, dup, zc = [], {}, {}

    def z(i):
        if i not in zc:
            zc[i] = _z(traces[:, i])
        return zc[i]
    for i in np.argsort(-np.asarray(score), kind="stable"):
        i = int(i)
        if status[i].lower() == "rejected" or not (np.isfinite(cx[i]) and np.isfinite(cy[i])):
            continue
        for k in kept:
            if np.hypot(cx[i] - cx[k], cy[i] - cy[k]) <= thr:
                r = float(z(i) @ z(k)) / len(z(i))
                if r >= min_corr:
                    dup[i] = k
                    break
        else:
            kept.append(i)
    return dup


def reject_duplicates_in_exports(props_csv, traces_csv, dist_frac=DUP_DIST_FRAC,
                                 min_corr=DUP_MIN_CORR):
    """Apply duplicate rejection to an exported props/traces CSV pair (IDPS layout):
    sets Status to 'rejected' in both files and adds a RejectReason column to props.
    Returns the list of (duplicate, kept) cell names."""
    with open(props_csv, newline="") as f:
        rd = csv.DictReader(f)
        fields = [h.strip() for h in rd.fieldnames]
        rows = [{k.strip(): (v or "").strip() for k, v in r.items()} for r in rd]
    with open(traces_csv, newline="") as f:
        lines = f.read().splitlines()
    names = [x.strip().strip('"') for x in lines[0].split(",")][1:]
    col = {n: j for j, n in enumerate(names)}

    def isnum(x):
        try:
            float(x)
            return True
        except ValueError:
            return x.strip() == ""
    has_status = len(lines) > 1 and not all(isnum(x) for x in lines[1].split(",")[1:])
    skip = 2 if has_status else 1
    try:
        data = np.loadtxt(traces_csv, delimiter=",", skiprows=skip)
    except ValueError:
        data = np.genfromtxt(traces_csv, delimiter=",", skip_header=skip)
    data = np.atleast_2d(data)[:, 1:len(names) + 1]
    idx = [col.get(r["Name"]) for r in rows]
    use = [i for i, c in enumerate(idx) if c is not None]
    X = data[:, [idx[i] for i in use]]

    def num(r, k, d=np.nan):
        try:
            return float(r[k])
        except (KeyError, ValueError):
            return d
    status = [rows[i].get("Status", "undecided") for i in use]
    size = [num(rows[i], "Size") for i in use]
    score, _ = cell_scores(X, status, size, [num(rows[i], "NumComponents", 1) for i in use])
    dup = find_duplicates([num(rows[i], "CentroidX") for i in use],
                          [num(rows[i], "CentroidY") for i in use],
                          size, X, score, status, dist_frac, min_corr)
    pairs = []
    if "RejectReason" not in fields:
        fields.append("RejectReason")
    for r in rows:
        r.setdefault("RejectReason", "")
    for j, k in dup.items():
        r = rows[use[j]]
        r["Status"] = "rejected"
        r["RejectReason"] = f"duplicate of {rows[use[k]]['Name']}"
        pairs.append((r["Name"], rows[use[k]]["Name"]))
    tmp = props_csv + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, props_csv)

    # traces: status row (line 2) follows the props
    if has_status and pairs:
        st = {r["Name"]: r.get("Status", "undecided") for r in rows}
        head = lines[1].split(",")
        lines[1] = ",".join([head[0]] + [f" {st.get(n, h.strip())}"
                                         for n, h in zip(names, head[1:])])
        tmp = traces_csv + ".tmp"
        with open(tmp, "w", newline="") as f:
            f.write("\n".join(lines) + "\n")
        os.replace(tmp, traces_csv)
    return pairs
