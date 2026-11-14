# Spatial transcriptomics on TMAs: CosMx whole transcriptome (WTX), and Xenium vs CosMx

Analysis code for tissue microarrays (`TMA_ORGAN`, `TMA_TUMOR`) measured with the **CosMx SMI
whole-transcriptome panel**. It mirrors the Xenium pipeline in
[spatial_transcriptomics_proteomics](https://github.com/jpimentabernardes/spatial_transcriptomics_proteomics)
step by step, and ends with a comparison of the two platforms on the same TMAs.

| Folder | What it is |
|---|---|
| [`experiment1_tma/`](experiment1_tma/README.md) | The pipeline, R (+ one Python step), run step by step: `00` what the CosMx export contains, `01` object and QC (cells and fields of view), `02` TMA dearraying, `03` segmentation QC, `04` cell-type annotation, `05` first TMA insights, `07` / `09` H&E regions, **`10` Xenium (probe panel) vs CosMx (WTX): transcript quantification, accuracy, cell annotation, probe crowding**. Settings live in `experiment1_tma/config/`. |
| [`tma_toolkit/`](tma_toolkit/README.md) | Python package for large TMAs (e.g. 200 cores): core-metadata sheet, automatic dearraying, selecting and saving sets of cores, core/group comparisons on AnnData. Reads CosMx flat files and Xenium outs. |
| `tools/make_notebooks.py` | Builds the Jupyter notebooks in `*/notebooks/` from the Python scripts (same code, one cell per step). |

Start with [`experiment1_tma/README.md`](experiment1_tma/README.md): it explains what a CosMx
run delivers and which files are needed, lists every step with its inputs and outputs, and
describes how the platform comparison is done and how to read it.
