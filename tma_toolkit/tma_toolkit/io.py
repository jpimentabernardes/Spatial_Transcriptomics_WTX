"""Load CosMx flat files (or Xenium outs / any segmentation in the 10x layout) into AnnData.

One AnnData for all slides:
  * X / layers["counts"]  genes only (raw counts)
  * obs                   slide_id, cell_id, x_um, y_um, cell_area, nucleus_area, n_counts, n_genes,
                          n_<control class> (negative-control probe / codeword counts per cell);
                          CosMx: also fov, x_global_px, y_global_px
  * obsm["spatial"]       x_um, y_um (microns, y grows downwards as in an image)
  * uns["tma"]["n_control_features"]  number of features per control class and slide, needed to
                          turn control counts into a per-probe background rate
  * obs_names             "<slide_id>_<cell_id>" (cell ids repeat across slides; CosMx cell_id = <fov>_<cell_ID>)

CosMx: negative probes ("Negative*") are stored as n_ctrl_probe and system controls
("SystemControl*") as n_ctrl_codeword -- the same columns as Xenium's negative control probes /
codewords, so everything downstream treats both platforms alike.

Negative controls are kept (as per-cell counts) because they are the background model for
probe-specificity questions; they are never mixed into X.
"""

from __future__ import annotations

import os
import re

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

CONTROL_TYPES = {
    "Negative Control Probe": "ctrl_probe",
    "Negative Control Codeword": "ctrl_codeword",
    "Unassigned Codeword": "ctrl_unassigned",
    "Blank Codeword": "ctrl_unassigned",
    "Genomic Control": "ctrl_genomic",
    "Deprecated Codeword": "ctrl_deprecated",
}


def _read_10x_mtx(mtx_dir: str):
    m = scipy.io.mmread(os.path.join(mtx_dir, "matrix.mtx.gz")).T.tocsr()   # cells x features
    feats = pd.read_csv(os.path.join(mtx_dir, "features.tsv.gz"), sep="\t", header=None)
    feats = feats.iloc[:, :3]
    feats.columns = ["id", "name", "type"][: feats.shape[1]]
    if "type" not in feats:
        feats["type"] = "Gene Expression"
    bc = pd.read_csv(os.path.join(mtx_dir, "barcodes.tsv.gz"), header=None)[0].astype(str).to_numpy()
    return m, feats, bc


def _read_10x_h5(path: str):
    import h5py
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        shape = tuple(g["shape"][:])
        m = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).T.tocsr()
        dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])  # noqa: E731
        feats = pd.DataFrame({"id": dec(g["features/id"][:]), "name": dec(g["features/name"][:]),
                              "type": dec(g["features/feature_type"][:])})
        bc = dec(g["barcodes"][:])
    return m, feats, bc


def read_feature_matrix(seg_dir: str):
    """cells x features CSR, features DataFrame (id, name, type), barcodes."""
    mtx_dir = os.path.join(seg_dir, "cell_feature_matrix")
    h5 = os.path.join(seg_dir, "cell_feature_matrix.h5")
    if os.path.isdir(mtx_dir):
        return _read_10x_mtx(mtx_dir)
    if os.path.exists(h5):
        return _read_10x_h5(h5)
    raise FileNotFoundError(f"No cell_feature_matrix/ or cell_feature_matrix.h5 in {seg_dir}")


def read_cells_table(seg_dir: str) -> pd.DataFrame:
    pq = os.path.join(seg_dir, "cells.parquet")
    csv = os.path.join(seg_dir, "cells.csv.gz")
    if os.path.exists(pq):
        cells = pd.read_parquet(pq)
    elif os.path.exists(csv):
        cells = pd.read_csv(csv)
    else:
        raise FileNotFoundError(f"No cells.parquet or cells.csv.gz in {seg_dir}")
    cells["cell_id"] = cells["cell_id"].map(lambda v: v.decode() if isinstance(v, bytes) else str(v))
    for col in ("cell_area", "nucleus_area"):
        if col not in cells:
            cells[col] = np.nan
    return cells


def load_xenium(seg_dir: str, slide_id: str, segmentation: str = "xenium") -> ad.AnnData:
    """One slide -> AnnData (see module docstring for the layout)."""
    m, feats, bc = read_feature_matrix(seg_dir)
    cells = read_cells_table(seg_dir).set_index("cell_id")
    keep = np.isin(bc, cells.index)
    if (~keep).any():
        print(f"{slide_id}: {(~keep).sum()} matrix cells without centroid dropped")
    m, bc = m[keep], bc[keep]
    cells = cells.loc[bc]

    is_gene = (feats["type"] == "Gene Expression").to_numpy()
    obs = pd.DataFrame(index=[f"{slide_id}_{b}" for b in bc])
    obs["slide_id"] = slide_id
    obs["cell_id"] = bc
    obs["x_um"] = cells["x_centroid"].to_numpy(float)
    obs["y_um"] = cells["y_centroid"].to_numpy(float)
    obs["cell_area"] = cells["cell_area"].to_numpy(float)
    obs["nucleus_area"] = cells["nucleus_area"].to_numpy(float)
    obs["segmentation"] = segmentation
    n_ctrl = {}
    for ftype, short in CONTROL_TYPES.items():
        sel = (feats["type"] == ftype).to_numpy()
        if sel.any():
            prev = obs[f"n_{short}"].to_numpy() if f"n_{short}" in obs else 0
            obs[f"n_{short}"] = prev + np.asarray(m[:, sel].sum(1)).ravel()
            n_ctrl[short] = n_ctrl.get(short, 0) + int(sel.sum())

    X = m[:, is_gene].astype(np.float32)
    var = pd.DataFrame(index=feats.loc[is_gene, "name"].to_numpy())
    var["feature_id"] = feats.loc[is_gene, "id"].to_numpy()
    obs["n_counts"] = np.asarray(X.sum(1)).ravel()
    obs["n_genes"] = np.asarray((X > 0).sum(1)).ravel()

    a = ad.AnnData(X=X, obs=obs, var=var)
    a.var_names_make_unique()
    a.layers["counts"] = a.X.copy()
    a.obsm["spatial"] = obs[["x_um", "y_um"]].to_numpy()
    a.uns["tma"] = {"n_control_features": {slide_id: n_ctrl}, "n_genes_panel": {slide_id: int(is_gene.sum())}}
    return a


# ---------------------------------------------------------------------------
# CosMx flat files (AtoMx export): *_exprMat_file.csv(.gz) + *_metadata_file.csv(.gz)
# ---------------------------------------------------------------------------
COSMX_PIXEL_SIZE_UM = 0.12028
COSMX_NEG_REGEX = r"^(Negative|NegPrb|NegProbe)"
COSMX_SYS_REGEX = r"^(SystemControl|FalseCode)"


def find_cosmx_file(directory: str, pattern: str, prefix: str | None = None) -> str | None:
    """The one file below `directory` whose name matches `pattern` (AtoMx nests flatFiles/<slide>/)."""
    hits = []
    for root, _, files in os.walk(directory):
        for f in files:
            if re.search(pattern, f, flags=re.I) and (not prefix or f.startswith(prefix)):
                hits.append(os.path.join(root, f))
    hits = [h for h in hits if h + ".gz" not in hits]
    if len(hits) > 1:
        raise ValueError(f"Several files match {pattern} below {directory}: {[os.path.basename(h) for h in hits]}. "
                         "Point cosmx_dir at one slide's folder or set cosmx_prefix.")
    return hits[0] if hits else None


def load_cosmx(directory: str, slide_id: str, prefix: str | None = None, flip_y: bool = True,
               pixel_size_um: float = COSMX_PIXEL_SIZE_UM, chunk_rows: int = 2000) -> ad.AnnData:
    """One CosMx slide (flat files) -> AnnData.

    The expression matrix is a dense csv with one column per target (~19,000 for the whole
    transcriptome), so it is read in chunks of rows and made sparse chunk by chunk. Rows with
    cell_ID 0 (transcripts outside every cell) are dropped. Coordinates: global pixels x pixel
    size; y is mirrored (flip_y) so that it grows downwards like every other plot here.
    """
    expr = find_cosmx_file(directory, r"exprMat_file\.csv(\.gz)?$", prefix)
    meta = find_cosmx_file(directory, r"metadata_file\.csv(\.gz)?$", prefix)
    if expr is None or meta is None:
        raise FileNotFoundError(f"No *exprMat_file.csv(.gz) / *metadata_file.csv(.gz) below {directory}")
    id_cols = ["fov", "cell_ID", "cellID", "cell", "cell_id"]
    blocks, fov, cid, feats = [], [], [], None
    for chunk in pd.read_csv(expr, chunksize=chunk_rows):
        if feats is None:
            feats = [c for c in chunk.columns if c not in id_cols]
        fov.append(chunk["fov"].to_numpy()); cid.append(chunk["cell_ID"].to_numpy())
        blocks.append(sp.csr_matrix(chunk[feats].to_numpy(dtype=np.float32)))
    m = sp.vstack(blocks).tocsr()
    fov, cid = np.concatenate(fov), np.concatenate(cid)
    m, fov, cid = m[cid != 0], fov[cid != 0], cid[cid != 0]
    bc = np.array([f"{f}_{c}" for f, c in zip(fov, cid)])

    md = pd.read_csv(meta)
    md.index = [f"{f}_{c}" for f, c in zip(md["fov"], md["cell_ID"])]
    keep = np.isin(bc, md.index)
    if (~keep).any():
        print(f"{slide_id}: {(~keep).sum()} matrix cells without a metadata row dropped")
    m, bc, fov = m[keep], bc[keep], fov[keep]
    md = md.loc[bc]
    y_px = md["CenterY_global_px"].to_numpy(float)
    y_ref = np.nanmax(pd.read_csv(meta, usecols=["CenterY_global_px"])["CenterY_global_px"].to_numpy(float)) if flip_y else 0.0

    feats = np.array(feats)
    is_neg = np.array([bool(re.match(COSMX_NEG_REGEX, f)) for f in feats])
    is_sys = np.array([bool(re.match(COSMX_SYS_REGEX, f)) for f in feats])
    is_gene = ~(is_neg | is_sys)
    obs = pd.DataFrame(index=[f"{slide_id}_{b}" for b in bc])
    obs["slide_id"] = slide_id
    obs["cell_id"] = bc
    obs["fov"] = fov
    obs["x_um"] = md["CenterX_global_px"].to_numpy(float) * pixel_size_um
    obs["y_um"] = ((y_ref - y_px) if flip_y else y_px) * pixel_size_um
    obs["x_global_px"] = md["CenterX_global_px"].to_numpy(float)
    obs["y_global_px"] = y_px
    obs["cell_area"] = (md["Area.um2"] if "Area.um2" in md else md["Area"] * pixel_size_um ** 2).to_numpy(float)
    obs["nucleus_area"] = (md["NucArea"] * pixel_size_um ** 2).to_numpy(float) if "NucArea" in md else np.nan
    obs["segmentation"] = "cosmx"
    n_ctrl = {}
    for sel, short in ((is_neg, "ctrl_probe"), (is_sys, "ctrl_codeword")):
        if sel.any():
            obs[f"n_{short}"] = np.asarray(m[:, sel].sum(1)).ravel()
            n_ctrl[short] = int(sel.sum())
    X = m[:, is_gene].astype(np.float32)
    var = pd.DataFrame(index=feats[is_gene])
    var["feature_id"] = feats[is_gene]
    obs["n_counts"] = np.asarray(X.sum(1)).ravel()
    obs["n_genes"] = np.asarray((X > 0).sum(1)).ravel()
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.var_names_make_unique()
    a.layers["counts"] = a.X.copy()
    a.obsm["spatial"] = obs[["x_um", "y_um"]].to_numpy()
    a.uns["tma"] = {"n_control_features": {slide_id: n_ctrl}, "n_genes_panel": {slide_id: int(is_gene.sum())}}
    return a


def load_slides(slides: pd.DataFrame, base_dir: str = ".", segmentation: str | None = None,
                seg_root: str | None = None) -> ad.AnnData:
    """All slides of slides.csv -> one AnnData.

    slides.csv has slide_id and cosmx_dir (CosMx flat files) or xenium_dir (Xenium outs).
    segmentation: None = the platform's own ("cosmx" if there is a cosmx_dir column, else
    "xenium"); any other name reads <seg_root>/<segmentation>/<slide_id>/ in the 10x layout.
    """
    is_cosmx = "cosmx_dir" in slides.columns
    if segmentation is None:
        segmentation = "cosmx" if is_cosmx else "xenium"
    parts, n_ctrl, n_genes = [], {}, {}
    for row in slides.itertuples(index=False):
        if segmentation in ("cosmx", "xenium"):
            raw = row.cosmx_dir if segmentation == "cosmx" else row.xenium_dir
            d = raw if os.path.isabs(raw) else os.path.join(base_dir, raw)
        else:
            d = os.path.join(seg_root or os.path.join(base_dir, "data", "segmentation"), segmentation, row.slide_id)
        print(f"Loading {row.slide_id} [{segmentation}] from {d}")
        if segmentation == "cosmx":
            prefix = getattr(row, "cosmx_prefix", None)
            a = load_cosmx(d, row.slide_id, prefix if isinstance(prefix, str) and prefix else None)
        else:
            a = load_xenium(d, row.slide_id, segmentation)
        n_ctrl.update(a.uns["tma"]["n_control_features"])
        n_genes.update(a.uns["tma"]["n_genes_panel"])
        parts.append(a)
    adata = ad.concat(parts, join="outer", merge="first") if len(parts) > 1 else parts[0]
    for col in [c for c in adata.obs if c.startswith("n_ctrl_")]:
        adata.obs[col] = adata.obs[col].fillna(0)
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr()
    adata.uns["tma"] = {"n_control_features": n_ctrl, "n_genes_panel": n_genes}
    return adata
