#!/usr/bin/env python
"""07_register_he.py -- bring an H&E image into the CosMx coordinate frame
(microns of 01: x_um / y_um), so that pathology regions drawn on the H&E can
be assigned to cells in 09_he_regions.R.

The CosMx flat files contain no slide image to register against, so the
transform is fitted from LANDMARKS: points you can find in both the H&E and
the cell map. Two ways to give them (--landmarks, a CSV):

  core centres   columns core_id, moving_x, moving_y
                 moving_x / moving_y = centre of that core in the H&E, in
                 full-resolution pixels (QuPath shows them in the status bar).
                 The matching CosMx position is taken from the dearraying
                 (tables/02_cores_all.csv). Easiest on a TMA: 6 or more cores
                 spread over the slide, not all in one row.
                 NOTE: on CosMx a "core" is the part covered by FOVs, so its
                 centre in 02 can be off by a few hundred um for half-imaged
                 cores -- prefer cores that were imaged whole, or use points.
  free points    columns moving_x, moving_y, fixed_x, fixed_y
                 fixed_x / fixed_y in microns, read off the axes of
                 figures/01_counts_spatial_by_slide.png or 00_fov_map_<slide>.png
                 (tissue corners, gaps, single cores).

With >= 6 landmarks an affine transform is fitted, otherwise a similarity
transform (rotation, uniform scale, shift; a mirror image is allowed and
reported). Landmarks further than 3 MAD from the fit are dropped and listed.

H&E on the CosMx slide itself (stained after the run) is the same section: a
good landmark fit is then accurate to a few cells. H&E of a neighbouring
section aligns tissue structure only.

Output: data/registration/transforms_he_<slide>.json   (H&E px -> CosMx um; read by 09)
        tables/07_he_landmark_residuals_<slide>.csv, figures/07_he_landmarks_<slide>.png

Example
-------
  python 07_register_he.py --project-dir /work/Spatial_TMA_WTX --slide-id TMA_ORGAN \
      --landmarks /work/Spatial_TMA_WTX/data/he/TMA_ORGAN_landmarks.csv
"""

import argparse
import json
import os

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# transforms (3x3 homogeneous, applied to column vectors [x, y, 1])
# ---------------------------------------------------------------------------
def apply_tf(M, xy):
    xy = np.asarray(xy, dtype=float)
    h = np.c_[xy, np.ones(len(xy))] @ M.T
    return h[:, :2]


def fit_affine(src, dst):
    A = np.c_[src, np.ones(len(src))]
    P, *_ = np.linalg.lstsq(A, dst, rcond=None)       # 3x2
    M = np.eye(3)
    M[:2, :] = P.T
    return M


def fit_similarity(src, dst):
    """Umeyama: rotation + uniform scale + translation (a mirror image if it fits better)."""
    mu_s, mu_d = src.mean(0), dst.mean(0)
    s, d = src - mu_s, dst - mu_d
    U, S, Vt = np.linalg.svd(d.T @ s / len(src))
    R = U @ Vt
    scale = S.sum() / s.var(0).sum()
    M = np.eye(3)
    M[:2, :2] = scale * R
    M[:2, 2] = mu_d - scale * R @ mu_s
    return M


def robust_fit(src, dst, model="auto", n_iter=5):
    """Least squares with iterative rejection of pairs > 3 MAD (min 50 um)."""
    src, dst = np.asarray(src, float), np.asarray(dst, float)
    if len(src) < 3:
        raise ValueError(f"need >= 3 landmark pairs, got {len(src)}")
    if model == "auto":
        model = "affine" if len(src) >= 6 else "similarity"
    fit = fit_affine if model == "affine" else fit_similarity
    keep = np.ones(len(src), bool)
    for _ in range(n_iter):
        M = fit(src[keep], dst[keep])
        res = np.linalg.norm(apply_tf(M, src) - dst, axis=1)
        thr = max(50.0, np.median(res[keep]) + 3 * 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep]))))
        new_keep = res <= thr
        if new_keep.sum() < 3 or (new_keep == keep).all():
            break
        keep = new_keep
    M = fit(src[keep], dst[keep])
    res = np.linalg.norm(apply_tf(M, src) - dst, axis=1)
    return M, res, keep, model


def rot_deg(M):
    return float(np.degrees(np.arctan2(M[1, 0], M[0, 0])))


# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-dir", required=True)
    p.add_argument("--slide-id", required=True, help="slide_id (slides.csv)")
    p.add_argument("--landmarks", required=True,
                   help="CSV: core_id,moving_x,moving_y (core centres in H&E px) or "
                        "moving_x,moving_y,fixed_x,fixed_y (H&E px -> CosMx um)")
    p.add_argument("--global-model", choices=["auto", "similarity", "affine"], default="auto")
    return p.parse_args()


def main():
    a = parse_args()
    reg_dir = os.path.join(a.project_dir, "data", "registration")
    tab_dir = os.path.join(a.project_dir, "tables")
    fig_dir = os.path.join(a.project_dir, "figures")
    for d in (reg_dir, fig_dir):
        os.makedirs(d, exist_ok=True)

    # ---- 1. landmarks: H&E pixels -> CosMx microns
    lm = pd.read_csv(a.landmarks)
    if {"fixed_x", "fixed_y"} <= set(lm.columns):
        names = lm["core_id"].astype(str).tolist() if "core_id" in lm else [f"landmark{i + 1}" for i in range(len(lm))]
    elif "core_id" in lm.columns:
        cores = pd.read_csv(os.path.join(tab_dir, "02_cores_all.csv"))
        cores = cores[(cores.slide_id == a.slide_id) & cores.core_id.notna()].set_index("core_id")
        missing = sorted(set(lm.core_id) - set(cores.index))
        if missing:
            print(f"core_id not in 02_cores_all.csv for {a.slide_id} (ignored): {missing}")
        lm = lm[lm.core_id.isin(cores.index)].copy()
        lm["fixed_x"] = cores.loc[lm.core_id, "x_center"].to_numpy()
        lm["fixed_y"] = cores.loc[lm.core_id, "y_center"].to_numpy()
        names = lm["core_id"].astype(str).tolist()
    else:
        raise SystemExit("--landmarks needs core_id,moving_x,moving_y or moving_x,moving_y,fixed_x,fixed_y")
    src = lm[["moving_x", "moving_y"]].to_numpy(float)
    dst = lm[["fixed_x", "fixed_y"]].to_numpy(float)

    # ---- 2. fit
    G, res, keep, model = robust_fit(src, dst, a.global_model)
    det = np.linalg.det(G[:2, :2])
    print(f"H&E -> CosMx, {model} from {len(src)} landmarks: median residual {np.median(res[keep]):.1f} um, "
          f"rotation {rot_deg(G):.1f} deg, {np.sqrt(abs(det)):.4f} um per H&E pixel"
          + ("  ** MIRRORED image (det<0) **" if det < 0 else "")
          + (f"; dropped as outliers: {[n for n, k in zip(names, keep) if not k]}" if (~keep).any() else ""))
    if np.median(res[keep]) > 100:
        print("WARNING: residuals above 100 um -- check the landmark table (swapped x/y? pixels of a downsampled image? "
              "core labels?).")

    # ---- 3. outputs
    qc = pd.DataFrame({"landmark": names, "moving_x": src[:, 0], "moving_y": src[:, 1], "fixed_x": dst[:, 0],
                       "fixed_y": dst[:, 1], "residual_um": np.round(res, 1), "used": keep})
    qc.to_csv(os.path.join(tab_dir, f"07_he_landmark_residuals_{a.slide_id}.csv"), index=False)
    with open(os.path.join(reg_dir, f"transforms_he_{a.slide_id}.json"), "w") as fh:
        json.dump({"moving": "he_px_level0", "fixed": "cosmx_um", "global_model": model, "global": G.tolist(),
                   "cores": {}, "landmark_residuals_um": dict(zip(names, np.round(res, 1).tolist()))}, fh, indent=1)
    print(f"Wrote transforms_he_{a.slide_id}.json")

    # the check: landmarks after the transform (x) on their CosMx positions (o), over the cells if 04 has run
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 9))
    cells_csv = os.path.join(reg_dir, f"cosmx_cells_{a.slide_id}.csv.gz")
    if os.path.exists(cells_csv):
        cells = pd.read_csv(cells_csv, usecols=["x_um", "y_um"])
        cells = cells.sample(min(len(cells), 200_000), random_state=0)
        ax.scatter(cells.x_um, cells.y_um, s=0.2, c="#c8c8c8", linewidths=0)
    moved = apply_tf(G, src)
    ax.scatter(dst[:, 0], dst[:, 1], s=60, facecolors="none", edgecolors="#2a78d6", label="CosMx position")
    ax.scatter(moved[:, 0], moved[:, 1], s=40, marker="x", c=np.where(keep, "#e34948", "#888888"), label="H&E landmark, transformed")
    for n, (x, y) in zip(names, moved):
        ax.annotate(n, (x, y), fontsize=6, xytext=(4, 4), textcoords="offset points")
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.legend(loc="best", fontsize=8)
    ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
    ax.set_title(f"{a.slide_id}: H&E landmarks in the CosMx frame (median residual {np.median(res[keep]):.0f} um)")
    fig.savefig(os.path.join(fig_dir, f"07_he_landmarks_{a.slide_id}.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
