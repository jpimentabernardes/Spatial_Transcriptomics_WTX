## =============================================================================
## config/config.R
##
## Single place for paths and parameters. Every R script in experiment1_tma/
## starts with source("config/config.R"), so run scripts with the working
## directory set to experiment1_tma/ (or set EXP1_CODE_DIR).
##
## Anything marked ADAPT is a judgement call that depends on your data --
## check it the first time real data lands, don't just accept the default.
## =============================================================================

## -----------------------------------------------------------------------
## Paths
## -----------------------------------------------------------------------
## ADAPT: where the raw data and outputs live on your cluster. Code stays in
## the git repo; data never goes in the repo. Use a DIFFERENT folder from the
## Xenium project: both pipelines write tables/02_cores_all.csv, figures/04_*, ...
PROJECT_DIR   <- Sys.getenv("EXP1_PROJECT_DIR", "/work_beegfs/sukmb430/Spatial_TMA_WTX")
CODE_DIR      <- Sys.getenv("EXP1_CODE_DIR", getwd())

CONFIG_DIR    <- file.path(CODE_DIR, "config")
PROCESSED_DIR <- file.path(PROJECT_DIR, "data", "processed")
FIGURES_DIR   <- file.path(PROJECT_DIR, "figures")
TABLES_DIR    <- file.path(PROJECT_DIR, "tables")
ANNOT_DIR     <- file.path(PROJECT_DIR, "annotation")   # editable cluster->label maps
SEG_DIR       <- file.path(PROJECT_DIR, "data", "segmentation")
REG_DIR       <- file.path(PROJECT_DIR, "data", "registration")

for (d in c(PROCESSED_DIR, FIGURES_DIR, TABLES_DIR, ANNOT_DIR, SEG_DIR, REG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

SLIDES_CSV        <- file.path(CONFIG_DIR, "slides.csv")
TMA_MAP_CSV       <- file.path(CONFIG_DIR, "tma_map.csv")
MARKERS_RNA_CSV   <- file.path(CONFIG_DIR, "marker_panel_rna.csv")                 # human TMAs
## MARKERS_RNA_CSV <- file.path(CONFIG_DIR, "marker_panel_rna_mouse_colon.csv")  # mouse
PROBE_EXPECT_CSV  <- file.path(CONFIG_DIR, "probe_expectations.csv")
FOV_TO_CORE_CSV   <- file.path(CONFIG_DIR, "fov_to_core.csv")        # optional, see 02 (mode D)
GENE_NAME_MAP_CSV <- file.path(CONFIG_DIR, "gene_name_map.csv")      # optional, see 10
XENIUM_PANEL_CSV  <- file.path(CONFIG_DIR, "xenium_panel_genes.csv") # optional; 10 writes it

## Paths inside slides.csv are relative to PROJECT_DIR unless absolute.
resolve_path <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p) || !nzchar(p)) return(NA_character_)
  if (grepl("^(/|[A-Za-z]:)", p)) p else file.path(PROJECT_DIR, p)
}

## -----------------------------------------------------------------------
## Platform. The object has ONE assay, named ASSAY, with the genes; Seurat then
## names the per-cell totals nCount_CosMx / nFeature_CosMx.
## -----------------------------------------------------------------------
ASSAY <- "CosMx"

## -----------------------------------------------------------------------
## CosMx input (00, 01). See README "What CosMx delivers".
## -----------------------------------------------------------------------
## Pixel edge of the CosMx camera. 0.12028 um on all commercial instruments;
## the 2022 NSCLC public data set (beta instrument) used 0.18. Only used where the
## export has no um column (Area.um2 is used when present).
COSMX_PIXEL_SIZE_UM <- 0.12028
## Columns of the expression matrix that are NOT genes, by name:
##   negative probes  = real probes against sequences that do not exist in the
##                      transcriptome -> off-target hybridisation + background
##   system controls  = barcodes that no probe carries ("false codes") -> decoding errors
## They play the role of Xenium's "negative control probe" / "negative control
## codeword" and are stored under the same names (nCount_ControlProbe /
## nCount_ControlCodeword), so 05 and 10 treat both platforms the same way.
COSMX_NEG_REGEX <- "^(Negative|NegPrb|NegProbe)"
COSMX_SYS_REGEX <- "^(SystemControl|FalseCode)"
## CosMx global coordinates have y growing UPWARDS (stage frame); every plot here
## uses the image convention (y downwards, scale_y_reverse). TRUE = y is mirrored
## on loading, so the tissue looks as in AtoMx. The original pixel coordinates
## stay in the metadata (x_global_px, y_global_px) to go back to AtoMx / napari.
COSMX_FLIP_Y <- TRUE
## Rows of the expression matrix read per chunk. A WTX matrix has ~19,000
## columns: one chunk of 2,000 cells is ~300 MB while it is being parsed.
COSMX_CHUNK_ROWS <- 2000

## -----------------------------------------------------------------------
## Segmentation used for everything downstream.
##   "cosmx" = the segmentation in the flat files (AtoMx; default, always exists)
##   other   = a folder data/segmentation/<name>/<slide_id>/ holding either another
##             flat-file export (e.g. a re-segmentation done in AtoMx) or the common
##             format (cell_feature_matrix/ + cells.csv.gz). Compare with 03, then
##             set the name here and re-run 01 -> 02.
## -----------------------------------------------------------------------
SEGMENTATION <- Sys.getenv("EXP1_SEGMENTATION", "cosmx")

## -----------------------------------------------------------------------
## Cell QC (01). ADAPT. These defaults are for the whole-transcriptome panel
## (~19,000 genes; typically several hundred to >1,000 transcripts per cell) and
## are deliberately lenient: 01 prints what a "2.5 SD below the mean of
## log10(counts)" rule would give on YOUR data, and 02 shows the result per core.
## -----------------------------------------------------------------------
MIN_COUNTS       <- 50
MIN_FEATURES     <- 30
MIN_CELL_AREA    <- 10     # um^2 -- below this is almost always debris
MAX_CELL_AREA    <- 1000   # um^2 -- hard ceiling; a per-slide MAD ceiling is also applied
AREA_MAD_UPPER   <- 5      # flag cells > median + 5 MAD of log(area)
MAX_CONTROL_FRAC <- 0.1    # flag cells where >10% of counts are negative probes / system controls (AtoMx default)
## Cells cut by a field-of-view (FOV) border: CosMx segments each FOV on its own,
## so a cell on the border is two half cells. SplitRatioToLocal (metadata) is the
## area of a border cell relative to the mean cell area of its FOV (0 = not on a
## border). Border cells below this ratio are flagged. NA = do not flag.
MIN_SPLIT_RATIO  <- 0.5
## FOV QC. An FOV whose cells have few transcripts is flagged (not removed) in 01
## (absolute) and in 02 (relative to the other FOVs of the SAME core, which is the
## fair comparison on a multi-organ TMA).
FOV_MIN_CELLS           <- 50
FOV_MIN_MEDIAN_COUNTS   <- 50
FOV_MIN_RATIO_TO_CORE   <- 0.5   # FOV median counts / median of its core's other FOVs

## QC maps in 01 (whole-slide plots)
QC_MAP_POINT_SIZE <- 0.3     # bigger = denser-looking tissue
QC_MAP_MAX_CELLS  <- 1e6     # cells drawn per figure (random subsample above this)
## Zoomed QC inspection (01, section 7): cell outlines from the polygons file
INSPECT_WINDOW_UM   <- 150     # side of each zoom window
INSPECT_N_PER_SLIDE <- 3       # windows chosen automatically per slide (mix of pass + fail)
## Your own windows instead (centre in um, as read off 01_qc_fail_spatial_by_slide.png):
## INSPECT_CENTERS <- data.frame(slide_id = "TMA_ORGAN", x = 5200, y = 3100)
INSPECT_CENTERS     <- NULL

## -----------------------------------------------------------------------
## TMA dearraying (02). Cores are usually 0.6-2 mm diameter with >= 200 um
## gaps; if neighbouring cores get merged, set a fixed DBSCAN_EPS_UM below the gap.
## On CosMx only the FOVs you placed are imaged (0.51 x 0.51 mm each), so a
## "core" is the part of the core covered by FOVs.
## -----------------------------------------------------------------------
DBSCAN_EPS_UM        <- NA    # NA = adaptive (3 x median kNN distance); or a fixed value in um
DBSCAN_MIN_PTS       <- 10
MIN_CELLS_PER_CORE   <- 200   # smaller DBSCAN clusters = debris / tissue fragments
CORE_MERGE_DIST_UM   <- 400   # fragments of one torn core closer than this get merged
CORE_MAX_DIST_FACTOR <- 1.25  # assign cells up to 1.25 x core radius from centre

## Virtual TMA for TESTING on whole sections (02a_make_pseudo_tma.R only).
PSEUDO_TMA_SLIDES       <- NULL       # NULL = all slides; or c("TMA_ORGAN")
PSEUDO_CORE_DIAMETER_UM <- 600        # typical TMA needle: 0.6-1.5 mm
PSEUDO_CORE_PITCH_UM    <- 1000       # centre-to-centre distance of the punch grid
PSEUDO_MIN_CELLS        <- 500        # a punch must contain this many QC-pass cells
PSEUDO_PATIENT_ID       <- NULL       # NULL = slide_id
PSEUDO_ORGAN            <- "colon"
PSEUDO_CENTERS          <- NULL       # roll centre per slide, e.g. data.frame(slide_id = "TMA_ORGAN", x = 5000, y = 4800)

## Core-level QC: cores failing these are flagged (not removed) in 02.
CORE_MIN_CELLS         <- 500
CORE_MIN_MEDIAN_COUNTS <- 100

## -----------------------------------------------------------------------
## Clustering / annotation (04)
## -----------------------------------------------------------------------
## Whole transcriptome: PCA runs on the N_HVG most variable genes PLUS every gene
## of marker_panel_rna.csv, CANONICAL_MARKERS and probe_expectations.csv (so the
## genes you annotate with are never dropped). ScaleData makes a dense
## (N_HVG x cells) matrix: 2,000 genes x 500,000 cells = 8 GB.
N_HVG            <- 2000
N_PCS            <- 30
USE_HARMONY      <- TRUE     # correct slide effects, NOT organ/patient effects
RES_CLUSTERS     <- 0.8      # 04: clustering of all cells (ST_clusters); higher = more, finer clusters
TOP_N_PER_CLUSTER <- 5       # 04: top FindAllMarkers genes per cluster shown in the dot plot
MARKERS_MAX_CELLS <- 2000    # 04: cells per cluster used by FindAllMarkers (speed; Inf = all)
## 04: canonical markers always shown in the dot plot, next to the top markers of
## each cluster. Matched to the panel ignoring case (Epcam = EPCAM); genes not on
## the panel are listed and dropped. ADAPT to tissue.
CANONICAL_MARKERS <- c(
  "EPCAM", "KRT8", "KRT19", "KRT20", "MKI67", "LGR5", "OLFM4", "MUC2", "CHGA",
  "PTPRC", "CD3E", "CD4", "CD8A", "FOXP3", "NKG7", "MS4A1", "CD79A", "JCHAIN", "MZB1",
  "CD68", "CD163", "C1QA", "CD14", "S100A8", "TPSAB1", "KIT", "CLEC9A", "LILRA4",
  "COL1A1", "DCN", "PDGFRA", "ACTA2", "MYH11", "DES", "RGS5", "PECAM1", "VWF", "LYVE1",
  "PLP1", "S100B", "ALB", "SFTPC", "PMEL"
)
## 04: cluster -> cell type, YOUR call, from figures/04_top_markers_dotplot_transcriptomics.png
## and 04_ST_clusters_spatial_<slide>.png. 04 writes a ready-to-fill template with
## each cluster's top markers: annotation/cluster_to_celltype_template_<seg>.R --
## copy it here, fill in the names, re-run 04. NULL = cells are labelled cluster_<n>.
## Cluster numbers change with the data, SEGMENTATION and RES_CLUSTERS: re-check
## the mapping whenever those change (04 stops if a cluster has no name).
## For the Xenium <-> CosMx comparison (10) use the SAME names as in the Xenium
## project wherever the cell type is the same; cell types that only the whole
## transcriptome resolves get new names.
CLUSTER_TO_CELLTYPE <- NULL
## Lineage per cluster, one of Epithelial / Immune / Fibroblast / Endothelial /
## Smooth_muscle_pericyte / Neural_glia (level-1 names of marker_panel_rna.csv):
## 05 uses lineage == "Immune" and matches probe_expectations.csv to these names,
## 10 compares the two platforms lineage by lineage.
## NULL = the marker-signature suggestion written in the template.
CLUSTER_TO_LINEAGE <- NULL
USE_UCELL        <- TRUE     # rank-based scores. Falls back to AddModuleScore

## -----------------------------------------------------------------------
## Aim 1 / Aim 2 (05)
## -----------------------------------------------------------------------
SBR_LOG2_MIN      <- 1       # gene counted "detected" in a core if >= 2x neg-probe background...
SBR_PADJ_MAX      <- 0.01    # ...and Poisson-test BH-adjusted p below this
MIN_CELLS_PSEUDOBULK <- 20   # min cells per core x cell type for pseudobulk
NICHE_K_NEIGHBOURS <- 20
N_NICHES           <- 8

## -----------------------------------------------------------------------
## Xenium (probe panel) <-> CosMx (whole transcriptome) comparison (10)
## -----------------------------------------------------------------------
## ADAPT: output folder of the Xenium pipeline (spatial_transcriptomics_proteomics,
## its EXP1_PROJECT_DIR) and the segmentation used there. 10 reads
## <XENIUM_PROJECT_DIR>/data/processed/xenium_<XENIUM_SEGMENTATION>_annotated.rds
## (or _insights.rds). Cores are matched by core_id: use the same core_id in both
## projects' tma_map.csv (or give config/core_id_map.csv: cosmx_core_id,xenium_core_id).
XENIUM_PROJECT_DIR  <- Sys.getenv("EXP1_XENIUM_PROJECT_DIR", "/work_beegfs/sukmb430/Spatial_TMA")
XENIUM_SEGMENTATION <- Sys.getenv("EXP1_XENIUM_SEGMENTATION", "xenium")
CORE_ID_MAP_CSV     <- file.path(CONFIG_DIR, "core_id_map.csv")       # optional
CMP_MIN_CELLS_CORE  <- 200    # a core is compared only with this many cells on BOTH platforms
CMP_MIN_CELLS_GROUP <- 50     # same, per core x lineage
CMP_REGISTER        <- TRUE   # align each CosMx core onto the Xenium core (needed for 10 section 4d)
CMP_BIN_UM          <- 50     # spatial bins for the registered comparison (serial sections: not the same cells)
CMP_MIN_NCC         <- 0.3    # a core counts as registered if the two cell-density maps correlate at least this much

set.seed(42)
