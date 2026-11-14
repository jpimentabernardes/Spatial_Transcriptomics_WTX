# Experiment 1 on CosMx: multi-organ and multi-tumour TMAs, whole transcriptome (WTX)

The same TMAs (`TMA_ORGAN`, `TMA_TUMOR`) and the same steps as the Xenium pipeline
([spatial_transcriptomics_proteomics](https://github.com/jpimentabernardes/spatial_transcriptomics_proteomics)),
with **CosMx SMI whole-transcriptome** data as input, and at the end a
**Xenium (probe panel) vs CosMx (WTX) comparison** (`10`).

Conventions are unchanged: editable parameters in `config/config.R`, `save_fig()` /
`save_table()`, explicit sanity checks, one saved `.rds` per stage, and
**the statistical unit is the core (or patient), never the cell.**

## What CosMx delivers, and what this pipeline reads

A CosMx run is exported from AtoMx. You will get one or more of these:

| Export | What it is | Used here |
|---|---|---|
| **Flat files** (csv, one set per slide) | `*_exprMat_file.csv.gz`: counts, one row per cell, one column per target (~19,000 genes + `Negative*` probes + `SystemControl*` codes); rows with `cell_ID` 0 hold the transcripts of a field of view that lie outside every cell | **yes, required** |
| | `*_metadata_file.csv.gz`: one row per cell (position in pixels, area, mean stain intensities, AtoMx QC flags) | **yes, required** |
| | `*_fov_positions_file.csv.gz`: where each field of view (FOV) sits on the slide | optional |
| | `*-polygons.csv.gz`: cell outlines | optional (zoom figures in `01`) |
| | `*_tx_file.csv.gz`: one row per transcript with its position and compartment (nucleus / cytoplasm / membrane / outside). Very large | optional (`03`) |
| Raw files | `CellStatsDir/` (morphology images and cell masks per FOV), `AnalysisResults/`, `RunSummary/` | no |
| Seurat `.RDS` | the same data as an R object | yes: set `cosmx_dir` to the `.RDS` (flat files preferred) |
| TileDB | array store used by AtoMx | no: export flat files |

**Ask for the flat-file export of each slide** (including the transcripts file if you want `03`
section 5). Then run `00_cosmx_inspect_input.R`: it tells you which files it found, how many
genes / negative probes / system controls the matrix has, whether the pixel size matches, and
draws a map of the FOV numbers. Nothing else needs to be known in advance.

How CosMx differs from Xenium, as far as this code is concerned:

| | Xenium | CosMx |
|---|---|---|
| unit of imaging | whole slide region | **fields of view** (0.51 x 0.51 mm) that you place; only they are imaged, and each is segmented on its own |
| coordinates | microns | **pixels** (0.12028 um), y upwards. `01` converts to microns, y downwards (`tables/01_cosmx_frame.csv`) |
| cell ids | unique per slide | numbered **within each FOV** -> cells are named `<slide>_<fov>_<cell_ID>` |
| counts | sparse matrix | dense csv, read in chunks and cached (first run: minutes per slide) |
| background controls | negative control probes / codewords | negative probes / system controls, stored under the same names (`nCount_ControlProbe`, `nCount_ControlCodeword`) |
| extra per cell | - | mean intensity of the morphology stains (`IF_PanCK`, `IF_CD45`, ...) |
| genes | a few hundred designed probes | ~19,000: clustering uses variable genes (`N_HVG`), QC thresholds are higher |

## Aims → where they are answered

| Question | Script / output |
|---|---|
| **Aim 1.** Which genes are measured above background in which organ; do immune targets generalise? | `05` (§1.1–1.4) → `tables/05_aim1_*` |
| – is a gene above background in organ X? | `05` §1.1: signal/background vs negative probes, **per core** (Xenium-panel genes marked) |
| – organ-specific, or just more immune cells in that organ? | `05` §1.2: tau on all cells **vs** tau within immune cells |
| – does the immune signal land in immune cells in every organ? | `05` §1.3: on-target enrichment + off-target/background ratio (read with MECR from `03`) |
| – do T cells/macrophages look the same across organs? | `05` §1.4 |
| **Aim 2.** How heterogeneous is the data? Can we compare organs, locations, patients? | `05` §2.1–2.5, `09` |
| **Platform comparison.** Transcript quantification, accuracy, cell annotation, probe crowding | `10` → `tables/10_*`, `figures/10_*` |

## Pipeline

```
CosMx flat files ─► 00 inspect ─► 01 object + QC ─► 02 TMA dearray + metadata ─► 04 annotation ─► 05 Aim 1 / Aim 2
                                        (cells, FOVs)      │  (cores, FOV-vs-core QC)      │               │
                                                           └─► 03 segmentation QC          │               │
Xenium project: xenium_<seg>_annotated.rds ────────────────────────────────────────────────┴─► 10 Xenium vs CosMx
H&E ─► 07 landmarks → transform ─► QuPath regions ─► 09 regions → cells
```

| # | Script | Lang | What it does |
|---|---|---|---|
| 00 | `00_cosmx_inspect_input.R` | R | Reads only headers and the cell metadata: which files exist, panel size, controls, pixel size, FOV map, template for `fov_to_core.csv`. **Run first.** |
| 01 | `01_cosmx_object_prep.R` | R | Loads all slides into one Seurat v5 object (assay `CosMx`, one centroid image per slide). Run metrics, per-cell QC flags (counts, area, FOV-border fragments, control fraction), **FOV QC**, negative-probe background, zoom figures with cell outlines. |
| 02a | `02a_make_pseudo_tma.R` | R | **Testing only.** Virtual TMA cores on whole sections. |
| 02 | `02_tma_dearray.R` | R | Cells → cores, three ways: polygons (`config/core_selections/`), **FOV → core table** (`config/fov_to_core.csv`, CosMx only) or automatic (DBSCAN + grid). Adds patient/organ/location, core QC, and compares every FOV with the other FOVs of its core. |
| 03 | `03_segmentation_qc.R` | R | Quality of the CosMx segmentation per core (transcripts inside cells, counts, area, **MECR**, control rate), where transcripts sit (nucleus / cytoplasm / membrane / outside; from the tx file), and a scorecard against a second segmentation if you provide one. |
| 04 | `04_annotation_wtx.R` | R | Clustering on variable + marker genes, FindAllMarkers, dot plot, clusters in the tissue, stains per cluster. You name the clusters in `CLUSTER_TO_CELLTYPE` / `CLUSTER_TO_LINEAGE` (config.R), starting from the template 04 writes to `annotation/`. |
| 05 | `05_tma_first_insights.R` | R | Aim 1 and Aim 2 analyses plus cellular niches. |
| 05c, 05d | `05c_checkpoint_genes_by_compartment.R`, `05d_plot_gene_core.R` | R | Checkpoint genes per compartment; one gene in one core. |
| 07 | `07_register_he.py` | Py | H&E pixels → CosMx microns from landmarks (core centres or free points). |
| 09 | `09_he_regions.R` | R | Pathology regions drawn on H&E (QuPath GeoJSON) → cells → composition and immune targets per region. |
| 10 | `10_xenium_vs_cosmx.R` | R | The platform comparison, see below. |

### Not carried over from the Xenium pipeline

The re-segmentation steps (`03a` cellpose, `03b` / `03d` / `03e` segger) are not in this
repository: they read the Xenium morphology image and Xenium-format transcript tables. On CosMx,
re-segmenting from images needs the **raw** export (`CellStatsDir/Morphology2D`, one image per
FOV), which is not implemented here. What works instead: a re-segmentation made in AtoMx (export
its flat files into `data/segmentation/<name>/<slide_id>/`) or any tool's output in the common
10x-style format; `03` then compares it with the original on the same cores and `SEGMENTATION`
in config.R selects it. The CellScape steps were already parked in the Xenium repository.

## 10: Xenium (probe panel) vs CosMx (WTX)

Both pipelines dearray with the same TMA map, so **a `core_id` present in both objects is the
same piece of tissue**. Consecutive sections do not contain the same cells: nothing is compared
cell to cell. Comparisons are per core, per core x lineage, or — after registering each CosMx
core onto its Xenium core from the cell centroids — per 50-um bin. CosMx images only the FOVs you
placed, Xenium the whole core; for registered cores the core-level numbers use the region
covered by both.

| Topic | Sections | What is measured |
|---|---|---|
| **Transcript quantification** | 1a–1e | counts per cell and per um2 of tissue on the shared genes; correlation of the gene profile per core and per core x lineage; per-gene sensitivity ratio; and per gene, whether it varies **across cores** the same way on both platforms (independent of each platform's efficiency for that gene) |
| **Accuracy** | 2a–2c | there is no ground truth, so: genes above the negative-probe background per core (same rule on both), negative-probe rate, MECR (lineage markers that must not share a cell), on-target enrichment of the genes in `probe_expectations.csv`; plus how many genes the whole transcriptome detects at all |
| **Cell annotation** | 3a–3d | lineage composition per core; label transfer on the shared genes in both directions (which Xenium cell type each CosMx cell resembles, and vice versa: shows what the whole transcriptome splits or merges); lineage profiles platform against platform; lineage agreement of neighbouring cells across the registered sections |
| **Probe / optical crowding** | 4a–4d | transcripts per um2 of cell; log-log slope of gene abundance (compression of abundant genes); the densest cells of a core, platform against platform (no registration needed); and the CosMx/Xenium ratio of the shared genes against the local transcript density in registered bins |

`tables/10_summary.csv` has one line per measure. Requirements: `04` has run in both projects
(lineage is needed; cell-type names may differ), `XENIUM_PROJECT_DIR` in config.R points at the
Xenium output folder, and the `core_id`s match (or give `config/core_id_map.csv`). Genes are
matched by name ignoring case; CosMx multi-gene probes (`"A/B"`) are listed but not used unless
you say so in `config/gene_name_map.csv` (`tables/10_shared_genes.csv` shows every decision).

How to read the crowding section. Optical crowding = more fluorescent spots per area than the
optics can resolve, so transcripts are lost where RNA is dense. With ~19,000 genes imaged at
once this is the expected weak point of a whole-transcriptome panel; a probe panel can saturate
too, for its most abundant genes. 4c and 4d ask the same question two ways and should agree in
sign: negative = CosMx loses relative to Xenium in dense cells / dense tissue, positive = the
reverse, about zero = neither. 4b (gene level) also contains probe design: 10x deliberately
reduces the efficiency of probes for very abundant genes.

## Key design decisions (and why)

- **Same object structure as the Xenium project.** One Seurat v5 object for all slides, cells
  named `<slide>_<fov>_<cell_ID>`, controls as per-cell metadata, centroids only. `05` and `10`
  therefore treat both platforms with the same code.
- **Microns, y downwards.** `01` converts CosMx global pixels once; the original pixels stay in
  the metadata (`x_global_px`, `y_global_px`) for going back to AtoMx / napari-cosmx. Check the
  orientation on `figures/00_fov_map_<slide>.png` (`COSMX_FLIP_Y`).
- **The TMA map is the source of truth** for sample identity (`config/tma_map.csv`), with the
  **same `core_id` as in the Xenium project**. Dearraying is checked visually
  (`figures/02_dearray_<slide>.png`).
- **FOVs are a QC unit.** A single FOV can fail on CosMx. `01` flags FOVs on absolute numbers,
  `02` against the other FOVs of the same core (same tissue, so a difference is technical).
- **Lenient cell filter.** `MIN_COUNTS = 50` for WTX is on the permissive side; `01` prints what
  a data-driven rule would give. Whole cores that are dim are flagged, not silently emptied.
- **You name the clusters.** For `10`, use the Xenium names where the cell type is the same.

## Getting started

1. Environments: `Rscript env/install_r.R`; for `07` / the toolkit `conda env create -f env/py_imaging.yml`.
2. Fill in `config/slides.csv` (`cosmx_dir` = folder with the flat files of that slide;
   `xenium_slide_id` = its name in the Xenium project) and copy `config/tma_map.csv` and
   `config/probe_expectations.csv` from the Xenium project.
3. Set `EXP1_PROJECT_DIR` (a folder of its own, not the Xenium one) and run from inside `experiment1_tma/`:
   ```
   export EXP1_PROJECT_DIR=/work_beegfs/sukmb430/Spatial_TMA_WTX
   export EXP1_XENIUM_PROJECT_DIR=/work_beegfs/sukmb430/Spatial_TMA      # for 10
   Rscript 00_cosmx_inspect_input.R                                      # read the report, look at figures/00_fov_map_*.png
   Rscript 01_cosmx_object_prep.R && Rscript 02_tma_dearray.R           # then CHECK figures/02_dearray_*.png
   Rscript 03_segmentation_qc.R
   Rscript 04_annotation_wtx.R                                           # fill in CLUSTER_TO_CELLTYPE / _LINEAGE, re-run
   Rscript 10_xenium_vs_cosmx.R                                          # also writes config/xenium_panel_genes.csv
   Rscript 05_tma_first_insights.R
   python 07_register_he.py --project-dir $EXP1_PROJECT_DIR --slide-id <slide> --landmarks <csv>
   (QuPath annotations on the H&E) ; Rscript 09_he_regions.R
   ```
   Memory: a WTX object of two TMAs with ~500,000 cells needs roughly 30–60 GB in `04`
   (see `N_HVG`); `10` holds both platforms' objects for a moment.

### Trying it without data

`tests/make_synthetic_data.R` writes a small synthetic data set in the real file formats of both
platforms (two TMAs, consecutive sections, a rotated slide, half-imaged cores, crowding on
CosMx only). `tests/run_synthetic.sh <out_dir> <xenium repo>/experiment1_tma` runs the Xenium
pipeline on the Xenium half, this pipeline on the CosMx half and then `10`. With
`SYN_CROWD_K=Inf` no crowding is simulated: section 4 of `10` should then report about zero.

### Python step as a Jupyter notebook

`notebooks/07_register_he.ipynb` has the same code as `07_register_he.py`, one cell per step
(run it in the `exp1-imaging` environment). The `.py` script is the source; after changing it,
regenerate with `python tools/make_notebooks.py` from the repo root (this overwrites the notebook).

## Things to verify with the real data (flagged in the scripts as ADAPT / CHECK)

- **The export itself**: run `00`. Column names changed between AtoMx versions; the reader copes
  with the documented variants and stops with a clear message otherwise.
- **Control names.** `COSMX_NEG_REGEX` / `COSMX_SYS_REGEX` must catch all control columns; `00`
  lists what was classified as what. Without negative probes `05` and `10` have no background model.
- **Pixel size** (`COSMX_PIXEL_SIZE_UM`): `00` compares it with `Area.um2 / Area`.
- **Orientation** (`COSMX_FLIP_Y`, `row_flip` / `col_flip` in slides.csv): a CosMx slide may be
  mounted differently from the Xenium one. Check `00_fov_map`, `02_dearray` and, in `10`,
  `10_registration_check.png` (a mirrored or rotated section is reported there).
- **QC thresholds** (`MIN_COUNTS`, areas, `MAX_CONTROL_FRAC`): look at the `01` figures first.
- **Species / gene symbols.** The marker files use human symbols; multi-gene probes and aliases
  show up as "missing" in `tables/04_panel_coverage.csv`.
- **Seurat `.RDS` input** was written from the documented structure and tested on a mock object only.
- **Confounding in the TMA design.** If all cores of one organ sit on one slide, organ and
  slide effects cannot be separated (`05` §2.2 says so).
- **`05c`** calls a helper (`add_site_and_tumour`) that is not in `R/helpers.R` — it was already
  missing in the Xenium repository this one was copied from, so `05c` does not run yet.
