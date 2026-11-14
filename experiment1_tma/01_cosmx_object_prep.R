## =============================================================================
## 01_cosmx_object_prep.R
##
## QC and object preparation for all CosMx whole-transcriptome (WTX) TMA slides.
## Same structure as the Xenium project's 01, so 02 - 05 run unchanged and the
## two platforms can be compared object against object in 10.
##
## INPUT FORMAT -- the AtoMx "flat files" of each slide (run 00 first to see
## what your export contains): the expression matrix (cells x ~19,000 targets)
## and the cell metadata are read; transcripts and polygons are not needed here.
##
## OBJECT TYPE -- why this structure (decided once here, used everywhere):
##   * One Seurat v5 object for ALL slides, one centroid FOV per slide
##     (Seurat's "FOV" = one image = one slide, named slide_id).
##   * Cells renamed <slide_id>_<fov>_<cell_ID>: CosMx numbers cells within
##     each field of view. The field of view is kept in metadata (`fov`).
##   * Genes in assay "CosMx"; negative probes and system controls kept as
##     per-cell metadata (nCount_ControlProbe, nCount_ControlCodeword). They are
##     the background model for Aim 1 and for the comparison with Xenium, so
##     they are never thrown away.
##   * Coordinates in MICRONS (x_um, y_um), y growing downwards like an image
##     (COSMX_FLIP_Y); the original global pixels stay as x_global_px / y_global_px.
##   * Centroids only, no cell polygons (many GB on multi-slide TMAs; not needed
##     for any statistic). The polygons file is only used for the zoom figures.
##   * Stain intensities per cell (Mean.PanCK, Mean.CD45, ...) are kept as
##     IF_<stain>: a protein-level check of the lineage calls in 04.
##   * TMA metadata (core, patient, organ, tissue type, location) is added per
##     cell in 02 from config/tma_map.csv -- the TMA map, not the file name, is
##     the source of truth for sample identity.
##
## CosMx-specific QC (not in the Xenium 01):
##   * FOV level: CosMx images one field of view (0.51 mm) at a time, and a
##     single FOV can fail (focus, bubbles, registration). 01 flags FOVs by
##     absolute numbers; 02 compares each FOV with the other FOVs of its core.
##   * FOV-border cells: cells cut by an FOV edge are segmented as fragments.
##   * Transcripts outside cells: the cell_ID 0 rows of the expression matrix.
##
## INPUT : config/slides.csv, CosMx flat files (or data/segmentation/<name>/<slide>/)
## OUTPUT: data/processed/cosmx_<seg>_raw.rds  (all cells, QC flags, not filtered)
##         tables/01_run_metrics.csv, 01_cosmx_frame.csv, 01_fov_qc.csv,
##         01_cell_qc_by_slide.csv, 01_slide_background.csv
## =============================================================================

source("config/config.R")
source("R/helpers.R")

slides <- read_slides()
print(slides[, c("slide_id", "tma_type", "cosmx_dir")])

## -----------------------------------------------------------------------
## 1. Load every slide with an identical structure and merge
##    The expression matrix is parsed in chunks and cached
##    (data/processed/cosmx_cache/): the first run takes minutes per slide.
## -----------------------------------------------------------------------
loaded <- lapply(seq_len(nrow(slides)), function(i) load_slide(slides[i, ], SEGMENTATION))
names(loaded) <- slides$slide_id

frames <- dplyr::bind_rows(lapply(names(loaded), function(sid) {
  f <- loaded[[sid]]$frame
  if (is.null(f)) return(NULL)
  data.frame(slide_id = sid, pixel_size_um = f$pixel_size_um, flip_y = f$flip_y, y_ref_px = f$y_ref_px)
}))
if (nrow(frames) > 0) {
  ## x_um = x_global_px * pixel_size_um ; y_um = (y_ref_px - y_global_px) * pixel_size_um  (flip_y = TRUE)
  save_table(frames, "01_cosmx_frame")
}
outside <- dplyr::bind_rows(lapply(names(loaded), function(sid) {
  u <- loaded[[sid]]$unassigned
  if (is.null(u)) NULL else cbind(slide_id = sid, u)
}))

cmx <- merge_slides(lapply(loaded, `[[`, "obj"))
polygon_files <- vapply(loaded, function(l) l$files$polygons %||% NA_character_, "")
frame_list <- lapply(loaded, `[[`, "frame")
rm(loaded); invisible(gc())

cmx$tma_type <- slides$tma_type[match(cmx$slide_id, slides$slide_id)]
print(cmx)
print(Images(cmx))   # one per slide
print(table(cmx$slide_id))

## -----------------------------------------------------------------------
## 2. Per-cell QC metrics
## -----------------------------------------------------------------------
meta <- cmx[[]]
ctrl_cols <- grep("^nCount_(ControlProbe|ControlCodeword|BlankCodeword)$", colnames(meta), value = TRUE)
meta$nCount_controls <- if (length(ctrl_cols) > 0) rowSums(meta[, ctrl_cols, drop = FALSE]) else 0
meta$control_frac <- meta$nCount_controls / pmax(meta$nCount_CosMx + meta$nCount_controls, 1)
meta$counts_per_um2 <- meta$nCount_CosMx / meta$cell_area
if (!"fov" %in% colnames(meta)) meta$fov <- NA_integer_
if (!"split_ratio" %in% colnames(meta)) meta$split_ratio <- NA_real_

## What a data-driven rule would give (CosMx WTX practice: 2.5 SD below the mean
## of log10 counts / genes). Printed for orientation; the thresholds USED are
## MIN_COUNTS / MIN_FEATURES in config.R.
suggest <- meta |>
  group_by(slide_id) |>
  summarise(counts_minus_2.5sd = round(10^(mean(log10(nCount_CosMx + 1)) - 2.5 * stats::sd(log10(nCount_CosMx + 1)))),
            genes_minus_2.5sd = round(10^(mean(log10(nFeature_CosMx + 1)) - 2.5 * stats::sd(log10(nFeature_CosMx + 1)))),
            counts_p01 = stats::quantile(nCount_CosMx, 0.01), counts_median = stats::median(nCount_CosMx), .groups = "drop")
message("Data-driven low-count thresholds per slide (config.R uses MIN_COUNTS = ", MIN_COUNTS,
        ", MIN_FEATURES = ", MIN_FEATURES, "):")
print(as.data.frame(suggest))

## Per-slide upper area cut: merged cells/doublets from segmentation. MAD on
## log(area) because area is right-skewed.
meta <- meta |>
  group_by(slide_id) |>
  mutate(
    area_upper = exp(median(log(cell_area), na.rm = TRUE) +
                       AREA_MAD_UPPER * mad(log(cell_area), na.rm = TRUE)),
    area_upper = pmin(area_upper, MAX_CELL_AREA)
  ) |>
  ungroup() |>
  as.data.frame()
rownames(meta) <- colnames(cmx)

## Flags, not deletions: each reason kept separately so you can see WHY cells
## go, per core, in 02. Cells with missing area (some segmentations) are not
## penalised on area.
meta$qc_low_counts   <- meta$nCount_CosMx < MIN_COUNTS | meta$nFeature_CosMx < MIN_FEATURES
meta$qc_small        <- !is.na(meta$cell_area) & meta$cell_area < MIN_CELL_AREA
meta$qc_large        <- !is.na(meta$cell_area) & meta$cell_area > meta$area_upper
meta$qc_high_control <- meta$control_frac > MAX_CONTROL_FRAC
## fragment of a cell cut by the FOV edge (split_ratio = 0: not on an edge)
meta$qc_fov_border   <- !is.na(MIN_SPLIT_RATIO) & !is.na(meta$split_ratio) &
  meta$split_ratio > 0 & meta$split_ratio < MIN_SPLIT_RATIO
meta$qc_pass <- !(meta$qc_low_counts | meta$qc_small | meta$qc_large | meta$qc_high_control | meta$qc_fov_border)
meta$qc_reason <- qc_reason(meta)          # first failing criterion, or "pass"

for (col in c("nCount_controls", "control_frac", "counts_per_um2", "area_upper", "fov", "split_ratio",
              "qc_low_counts", "qc_small", "qc_large", "qc_high_control", "qc_fov_border", "qc_pass", "qc_reason")) {
  cmx[[col]] <- meta[[col]]
}

qc_summary <- meta |>
  group_by(slide_id) |>
  summarise(
    n_cells = n(),
    n_fovs = n_distinct(fov),
    median_counts = median(nCount_CosMx),
    median_genes = median(nFeature_CosMx),
    median_area = median(cell_area, na.rm = TRUE),
    median_counts_per_um2 = median(counts_per_um2, na.rm = TRUE),
    pct_low_counts = 100 * mean(qc_low_counts),
    pct_small = 100 * mean(qc_small),
    pct_large = 100 * mean(qc_large),
    pct_fov_border = 100 * mean(qc_fov_border),
    pct_high_control = 100 * mean(qc_high_control),
    pct_pass = 100 * mean(qc_pass),
    .groups = "drop"
  )
print(as.data.frame(qc_summary))
save_table(qc_summary, "01_cell_qc_by_slide")

## -----------------------------------------------------------------------
## 3. Run-level metrics (CosMx has no metrics_summary.csv; this is the
##    equivalent). Look at these BEFORE any cell-level QC: few transcripts per
##    cell or a high negative-probe rate over a whole slide is a run / tissue
##    problem, not something cell filtering can fix.
## -----------------------------------------------------------------------
run_metrics <- meta |>
  group_by(slide_id) |>
  summarise(
    num_cells = n(), num_fovs = n_distinct(fov),
    genes_on_panel = first(nFeatTotal_Gene),
    transcripts_in_cells = sum(nCount_CosMx),
    mean_transcripts_per_cell = mean(nCount_CosMx),
    median_transcripts_per_cell = median(nCount_CosMx),
    median_genes_per_cell = median(nFeature_CosMx),
    transcripts_per_um2_cell_area = sum(nCount_CosMx) / sum(cell_area, na.rm = TRUE),
    negative_probe_counts_per_probe_per_cell =
      if ("nCount_ControlProbe" %in% colnames(meta)) mean(nCount_ControlProbe / pmax(nFeatTotal_ControlProbe, 1)) else NA,
    system_control_counts_per_code_per_cell =
      if ("nCount_ControlCodeword" %in% colnames(meta)) mean(nCount_ControlCodeword / pmax(nFeatTotal_ControlCodeword, 1)) else NA,
    .groups = "drop"
  ) |>
  mutate(mean_counts_per_gene_per_cell = mean_transcripts_per_cell / genes_on_panel,
         ## how many times above a negative probe the AVERAGE gene is
         signal_to_background = mean_counts_per_gene_per_cell / negative_probe_counts_per_probe_per_cell)
if (nrow(outside) > 0) {
  out_s <- outside |> group_by(slide_id) |> summarise(transcripts_outside_cells = sum(n_outside_genes), .groups = "drop")
  run_metrics <- run_metrics |> left_join(out_s, by = "slide_id") |>
    mutate(fraction_transcripts_in_cells = transcripts_in_cells / (transcripts_in_cells + transcripts_outside_cells))
}
print(as.data.frame(run_metrics))
save_table(run_metrics, "01_run_metrics")

## -----------------------------------------------------------------------
## 4. FOV-level QC. One row per slide x field of view. Flagged here on absolute
##    numbers only: on a multi-organ TMA a "dim" FOV may simply be a dim organ,
##    so 02 repeats this relative to the other FOVs of the same core.
## -----------------------------------------------------------------------
fov_qc <- meta |>
  group_by(slide_id, fov) |>
  summarise(
    n_cells = n(), x_um = mean(x_um), y_um = mean(y_um),
    median_counts = median(nCount_CosMx), median_genes = median(nFeature_CosMx),
    counts_per_um2 = sum(nCount_CosMx) / sum(cell_area, na.rm = TRUE),
    control_frac = sum(nCount_controls) / max(sum(nCount_CosMx + nCount_controls), 1),
    pct_cells_pass = 100 * mean(qc_pass),
    atomx_fov_flag = if ("atomx_qcFlagsFOV" %in% colnames(meta)) paste(unique(atomx_qcFlagsFOV), collapse = ";") else NA_character_,
    .groups = "drop"
  )
if (nrow(outside) > 0) {
  fov_qc <- fov_qc |> left_join(outside[, c("slide_id", "fov", "n_outside_genes")], by = c("slide_id", "fov"))
  in_cells <- meta |> group_by(slide_id, fov) |> summarise(n_in = sum(nCount_CosMx), .groups = "drop")
  fov_qc <- fov_qc |> left_join(in_cells, by = c("slide_id", "fov")) |>
    mutate(pct_transcripts_in_cells = 100 * n_in / pmax(n_in + n_outside_genes, 1)) |>
    dplyr::select(-n_in, -n_outside_genes)
}
fov_qc$fov_flag <- with(fov_qc, ifelse(n_cells < FOV_MIN_CELLS, "few_cells",
                                       ifelse(median_counts < FOV_MIN_MEDIAN_COUNTS, "low_counts", "ok")))
save_table(fov_qc, "01_fov_qc")
message("FOVs flagged: ", paste(names(table(fov_qc$fov_flag)), table(fov_qc$fov_flag), sep = " = ", collapse = ", "))

## -----------------------------------------------------------------------
## 5. Slide-level background: counts per negative probe vs counts per gene.
##    On the whole transcriptome most genes are NOT expressed in a given tissue,
##    so counts_per_gene is an average over many silent genes: read the ratio as
##    "how much of an average gene's signal is background", and use 05 (per core,
##    per gene) for decisions.
## -----------------------------------------------------------------------
bg <- meta |>
  group_by(slide_id) |>
  summarise(
    gene_counts = sum(nCount_CosMx),
    n_genes = first(nFeatTotal_Gene),
    ctrl_probe_counts = if ("nCount_ControlProbe" %in% colnames(meta)) sum(nCount_ControlProbe) else NA,
    n_ctrl_probes = if ("nFeatTotal_ControlProbe" %in% colnames(meta)) first(nFeatTotal_ControlProbe) else NA,
    ctrl_codeword_counts = if ("nCount_ControlCodeword" %in% colnames(meta)) sum(nCount_ControlCodeword) else NA,
    n_ctrl_codewords = if ("nFeatTotal_ControlCodeword" %in% colnames(meta)) first(nFeatTotal_ControlCodeword) else NA,
    .groups = "drop"
  ) |>
  mutate(
    counts_per_gene = gene_counts / n_genes,
    counts_per_ctrl_probe = ctrl_probe_counts / n_ctrl_probes,
    counts_per_ctrl_codeword = ctrl_codeword_counts / n_ctrl_codewords,
    background_ratio_probe = counts_per_ctrl_probe / counts_per_gene,
    background_ratio_codeword = counts_per_ctrl_codeword / counts_per_gene
  )
print(as.data.frame(bg))
save_table(bg, "01_slide_background")

## -----------------------------------------------------------------------
## 6. Figures
## -----------------------------------------------------------------------
## Violins on a log scale (a handful of huge or very bright cells otherwise
## squash everything), with the QC thresholds as dashed lines and the % of
## cells each threshold removes.
slide_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")
qc_violin <- function(col, label, trans, breaks, lo = NULL, hi = NULL) {
  pct_lo <- if (!is.null(lo)) 100 * mean(meta[[col]] < lo, na.rm = TRUE) else NA
  pct_hi <- if (!is.null(hi)) 100 * mean(meta[[col]] > hi, na.rm = TRUE) else NA
  sub <- paste(c(if (!is.na(pct_lo)) sprintf("%.1f%% below %g", pct_lo, lo),
                 if (!is.na(pct_hi)) sprintf("%.1f%% above %g", pct_hi, hi)), collapse = ", ")
  ggplot(meta, aes(slide_id, .data[[col]], fill = slide_id)) +
    geom_violin(scale = "width", linewidth = 0.2, colour = "grey30") +
    geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    { if (!is.null(lo)) geom_hline(yintercept = lo, linetype = 2) } +
    { if (!is.null(hi)) geom_hline(yintercept = hi, linetype = 2) } +
    scale_y_continuous(trans = trans, breaks = breaks) +
    scale_fill_manual(values = rep_len(slide_cols, nrow(slides)), guide = "none") +
    theme_bw(base_size = 10) + labs(x = NULL, y = NULL, title = label, subtitle = sub)
}
p_vln <- qc_violin("nCount_CosMx", "transcripts per cell", "log1p", c(0, 10, 30, 100, 300, 1000, 3000, 10000), lo = MIN_COUNTS) |
  qc_violin("nFeature_CosMx", "genes per cell", "log1p", c(0, 10, 30, 100, 300, 1000, 3000), lo = MIN_FEATURES) |
  qc_violin("cell_area", "cell area (um2)", "log1p", c(0, 5, 10, 30, 100, 300, 1000, 3000),
            lo = MIN_CELL_AREA, hi = MAX_CELL_AREA) |
  qc_violin("counts_per_um2", "transcripts per um2 of cell", "log1p", c(0, 0.3, 1, 3, 10, 30, 100)) |
  qc_violin("control_frac", "negative probe + system control fraction", "sqrt", c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1),
            hi = MAX_CONTROL_FRAC)
save_fig(p_vln, "01_qc_violin_by_slide", width = 20, height = 5)

## Spatial maps with micron axes: read coordinates here for INSPECT_CENTERS
## (section 7) and for core selections. QC_MAP_POINT_SIZE / QC_MAP_MAX_CELLS are set in config.R.
p_spatial <- plot_cores_spatial(
  meta |> mutate(log10_counts = log10(nCount_CosMx + 1)),
  color_by = "log10_counts", facet = FALSE, axes = TRUE,
  size = QC_MAP_POINT_SIZE, max_cells = QC_MAP_MAX_CELLS
) + facet_wrap(~ slide_id) +
  ggtitle("log10 transcripts per cell (look for whole cores or single FOVs that are dim = tissue/run problem)")
save_fig(p_spatial, "01_counts_spatial_by_slide", width = 8 * nrow(slides), height = 8)

## Coloured by WHY a cell fails -- tells you which threshold to look at.
p_pass <- plot_cores_spatial(meta[order(meta$qc_reason != "pass"), ],      # failing cells drawn on top
                             color_by = "qc_reason", facet = FALSE, axes = TRUE,
                             size = QC_MAP_POINT_SIZE, max_cells = QC_MAP_MAX_CELLS) +
  facet_wrap(~ slide_id) +
  scale_color_manual(values = c(QC_REASON_COLORS[-1], pass = "grey80"), drop = FALSE, name = "cell QC") +
  ggtitle("Cell QC: pass (grey) and the first criterion each failing cell fails")
save_fig(p_pass, "01_qc_fail_spatial_by_slide", width = 8 * nrow(slides), height = 8)

## FOV map: every field of view as a tile at its position, coloured by the
## median transcripts per cell, labelled with its number; flagged FOVs outlined.
fov_ext <- meta |> group_by(slide_id, fov) |>
  summarise(x0 = min(x_um), x1 = max(x_um), y0 = min(y_um), y1 = max(y_um), .groups = "drop") |>
  left_join(fov_qc[, c("slide_id", "fov", "median_counts", "fov_flag")], by = c("slide_id", "fov"))
p_fov <- ggplot(fov_ext) +
  geom_rect(aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = median_counts), colour = NA) +
  geom_rect(data = fov_ext[fov_ext$fov_flag != "ok", ], aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
            fill = NA, colour = "#e34948", linewidth = 0.8) +
  geom_text(aes((x0 + x1) / 2, (y0 + y1) / 2, label = fov), size = 2.2, colour = "white") +
  scale_fill_viridis_c(trans = "log10", name = "median\ntranscripts\nper cell") +
  scale_y_reverse() + coord_fixed() + facet_wrap(~ slide_id) + theme_bw() +
  labs(x = "x (um)", y = "y (um)", title = "FOV QC: median transcripts per cell per field of view (red outline = flagged)")
save_fig(p_fov, "01_fov_qc_map", width = 8 * nrow(slides), height = 8)

## >>> CHECK <<<
## If failing cells cluster at core EDGES or along FOV borders that is expected
## (cut cells). A single FOV that is much dimmer than its neighbours on the SAME
## core is a run problem (02 flags it). If a WHOLE core fails, that core is a
## tissue-quality issue (necrosis, fixation, detachment) -- keep it flagged in 02
## rather than silently letting the cell filter delete it; for Aim 2 a missing
## core is a missing sample.

## -----------------------------------------------------------------------
## 7. Zoomed inspection: QC pass / fail cells with their outlines
##    For each window: where it is on the slide | the cells of the window,
##    outlines from the polygons file coloured by QC result. Look for: failing
##    cells in tissue gaps (debris = correctly removed), very large outlines
##    (merged cells), fragments along the dashed FOV borders, or normal-looking
##    cells that fail only on counts (a threshold too strict for this tissue).
##    Windows are picked automatically where pass and fail cells mix; set
##    INSPECT_CENTERS in config.R to look at a specific place instead.
##    The stains themselves are not shown: they are only in the raw export
##    (CellStatsDir/Morphology2D) -- look at those places in AtoMx / napari-cosmx.
## -----------------------------------------------------------------------
xm <- cmx[[]]
xm$cell_id <- as.character(xm$cell_id)
windows <- if (!is.null(INSPECT_CENTERS)) {
  h <- INSPECT_WINDOW_UM / 2
  with(INSPECT_CENTERS, data.frame(slide_id = slide_id, x0 = x - h, y0 = y - h, x1 = x + h, y1 = y + h))
} else pick_inspection_windows(xm, INSPECT_WINDOW_UM, INSPECT_N_PER_SLIDE)
print(windows)
save_table(windows, "01_inspection_windows")
for (sid in unique(windows$slide_id)) {
  polys <- NULL
  if (SEGMENTATION == "cosmx" && !is.na(polygon_files[[sid]]) && !is.null(frame_list[[sid]])) {
    polys <- tryCatch(read_cosmx_polygons(polygon_files[[sid]], frame_list[[sid]]),
                      error = function(e) { message("polygons of ", sid, ": ", conditionMessage(e)); NULL })
  }
  for (k in which(windows$slide_id == sid)) {
    w <- windows[k, ]
    tag <- sprintf("%s_x%d_y%d", sid, round(w$x0), round(w$y0))
    p_zoom <- plot_inspection(xm[xm$slide_id == sid, ], polys, w,
                              title = sprintf("%s  x %.0f-%.0f, y %.0f-%.0f um", sid, w$x0, w$x1, w$y0, w$y1))
    save_fig(p_zoom, paste0("01_qc_zoom_", tag), width = 13, height = 6.5, dpi = 200)
  }
  rm(polys); invisible(gc())
}

saveRDS(cmx, obj_path("raw"))
message("Saved -> ", obj_path("raw"))
