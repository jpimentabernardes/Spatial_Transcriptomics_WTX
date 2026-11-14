## =============================================================================
## 10_xenium_vs_cosmx.R
##
## Xenium (targeted probe panel) vs CosMx (whole transcriptome, WTX) on the
## same TMAs (consecutive sections of TMA_ORGAN and TMA_TUMOR).
##
## The two platforms are compared CORE BY CORE: both pipelines dearray with the
## same TMA map, so a core_id present in both objects is the same piece of
## tissue. Consecutive sections do not contain the same cells, so nothing here
## is cell-to-cell: comparisons are per core, per core x lineage, or (after
## registration) per 50-um bin. The statistical unit is the core.
##
##   0. Shared genes, shared cores, registration of every CosMx core onto its
##      Xenium core (cell centroids only), region covered by both.
##   1. TRANSCRIPT QUANTIFICATION (shared genes)
##      1a per core: cells, tissue area, counts per cell and per um2 of tissue
##      1b per core: correlation of the gene profile between platforms
##      1c per gene: sensitivity ratio CosMx / Xenium
##      1d per gene: does it go up and down ACROSS CORES the same way on both
##         platforms? (independent of each platform's efficiency for the gene)
##      1e per core x lineage: 1b within one lineage
##   2. ACCURACY -- there is no ground truth; "accuracy" is read from controls
##      and from things that must not happen:
##      2a detection above the negative-probe background, per gene x core
##      2b MECR: marker pairs of different lineages in one cell (spill-over)
##      2c on-target enrichment of lineage genes (config/probe_expectations.csv)
##   3. CELL ANNOTATION
##      3a lineage composition per core
##      3b label transfer on the shared genes, both directions: which Xenium
##         cell type does each CosMx cell look like, and vice versa (= what the
##         whole transcriptome splits or merges compared with the panel)
##      3c lineage expression profiles, platform against platform
##      3d registered cores: do neighbouring cells of the two sections agree?
##   4. PROBE / OPTICAL CROWDING -- an imaging platform can only resolve so many
##      spots per um2. With ~19,000 genes imaged at once (WTX) the spot density
##      is far higher than with a few hundred, so dense cells may lose
##      transcripts; a probe panel can saturate too, for its top genes.
##      4a transcripts per um2 of cell, per platform and lineage
##      4b gene level: are highly expressed genes compressed on one platform?
##      4c cell level (no registration needed): the densest cells of a core,
##         platform against platform
##      4d bin level (registered cores): CosMx / Xenium ratio of the shared
##         genes against the local transcript density
##   5. Summary table.
##
## INPUT : cosmx_<seg>_annotated.rds / _insights.rds (this pipeline, 04 / 05)
##         <XENIUM_PROJECT_DIR>/data/processed/xenium_<seg>_annotated.rds / _insights.rds
##         (Xenium pipeline, 04 / 05). Paths: config.R, section "Xenium <-> CosMx".
## OUTPUT: tables/10_*.csv, figures/10_*.png, config/xenium_panel_genes.csv (if absent)
##
## READ WITH CARE
##   * Use the same lineage names in both projects (level-1 names of
##     marker_panel_rna.csv). Cell TYPE names may differ; 3b maps them.
##   * Sections 1c / 4b compare absolute counts: they include probe design
##     (number of probes per gene, deliberately "detuned" Xenium probes for very
##     abundant genes), chemistry and segmentation. 1d and 3 do not.
##   * Segmentation differs (Xenium: nucleus expansion or multimodal; CosMx:
##     membrane-stain based). Per-um2-of-tissue numbers avoid that; per-cell
##     numbers include it.
## =============================================================================

source("config/config.R")
source("R/helpers.R")
source("R/compare_helpers.R")

PLATFORM_COLOURS <- c(Xenium = "#2a78d6", CosMx = "#eb6834")
first_existing <- function(paths, what) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) stop("No ", what, " object found. Looked for:\n  ", paste(paths, collapse = "\n  "))
  ok[1]
}

## =============================================================================
## 0. Load and harmonise
## =============================================================================
cos_path <- first_existing(c(obj_path("insights"), obj_path("annotated")), "CosMx")
xen_path <- first_existing(file.path(XENIUM_PROJECT_DIR, "data", "processed",
                                     sprintf("xenium_%s_%s.rds", XENIUM_SEGMENTATION, c("insights", "annotated"))), "Xenium")
message("CosMx : ", cos_path, "\nXenium: ", xen_path)
slides <- read_slides()

## one tidy table per platform + the count matrix; the Seurat objects are dropped
tidy <- function(obj, assay, platform) {
  m <- obj[[]]
  ctrl <- if ("nCount_ControlProbe" %in% colnames(m)) "ControlProbe" else "ControlCodeword"
  d <- data.frame(
    cell = colnames(obj), platform = platform, slide_id = m$slide_id, core_id = m$core_id,
    x = m$x_um, y = m$y_um, cell_area = m$cell_area,
    n_counts = m[[paste0("nCount_", assay)]], n_genes = m[[paste0("nFeature_", assay)]],
    n_ctrl = m[[paste0("nCount_", ctrl)]], n_ctrl_feat = m[[paste0("nFeatTotal_", ctrl)]],
    lineage = if ("lineage" %in% colnames(m)) as.character(m$lineage) else NA_character_,
    cell_type = if ("cell_type" %in% colnames(m)) as.character(m$cell_type) else NA_character_,
    organ = if ("organ" %in% colnames(m)) as.character(m$organ) else NA_character_,
    tissue_type = if ("tissue_type" %in% colnames(m)) as.character(m$tissue_type) else NA_character_,
    stringsAsFactors = FALSE)
  attr(d, "control") <- ctrl
  list(meta = d, counts = GetAssayData(obj, assay = assay, layer = "counts"))
}
tmp <- readRDS(cos_path); cos <- tidy(tmp, "CosMx", "CosMx")
tmp <- readRDS(xen_path); xen <- tidy(tmp, "Xenium", "Xenium")
rm(tmp); invisible(gc())
message("Background model: CosMx ", attr(cos$meta, "control"), ", Xenium ", attr(xen$meta, "control"),
        " (ControlProbe = negative probes on both platforms)")
if (all(is.na(cos$meta$lineage)) || all(is.na(xen$meta$lineage))) {
  stop("Both objects need `lineage` (run 04 in both projects; lineage is suggested even before you name the clusters).")
}

## ---- core ids: the same id = the same tissue. Optional translation table. -------
if (file.exists(CORE_ID_MAP_CSV)) {
  cm <- read_config_csv(CORE_ID_MAP_CSV, stringsAsFactors = FALSE)
  stopifnot("core_id_map.csv needs cosmx_core_id and xenium_core_id" = all(c("cosmx_core_id", "xenium_core_id") %in% colnames(cm)))
  k <- match(xen$meta$core_id, cm$xenium_core_id)
  xen$meta$core_id <- ifelse(is.na(k), xen$meta$core_id, cm$cosmx_core_id[k])
}
n_c <- table(cos$meta$core_id); n_x <- table(xen$meta$core_id)
shared_cores <- intersect(names(n_c)[n_c >= CMP_MIN_CELLS_CORE], names(n_x)[n_x >= CMP_MIN_CELLS_CORE])
core_tab <- data.frame(core_id = union(names(n_c), names(n_x)))
core_tab$n_cells_cosmx <- as.integer(n_c[core_tab$core_id]); core_tab$n_cells_xenium <- as.integer(n_x[core_tab$core_id])
core_tab[is.na(core_tab)] <- 0L
core_tab$compared <- core_tab$core_id %in% shared_cores
message(length(shared_cores), " cores with >= ", CMP_MIN_CELLS_CORE, " cells on both platforms; ",
        sum(!core_tab$compared), " cores on one platform only or too small (tables/10_cores.csv).")
if (length(shared_cores) == 0) {
  stop("No shared core. The core_id values of the two projects' tma_map.csv must match (or give config/core_id_map.csv).\n",
       "  CosMx : ", paste(utils::head(names(n_c), 5), collapse = ", "), " ...\n  Xenium: ", paste(utils::head(names(n_x), 5), collapse = ", "), " ...")
}
core_info <- cos$meta[!duplicated(cos$meta$core_id), c("core_id", "slide_id", "organ", "tissue_type")]
rownames(core_info) <- core_info$core_id
core_info <- core_info[shared_cores, ]
## Xenium slide of every CosMx slide: slides.csv xenium_slide_id, else read off the shared cores
xen_slide_of_core <- xen$meta$slide_id[match(shared_cores, xen$meta$core_id)]
core_info$xenium_slide_id <- xen_slide_of_core

## ---- shared genes --------------------------------------------------------------
genes <- match_genes(rownames(xen$counts), rownames(cos$counts))
save_table(genes, "10_shared_genes")
message(sum(!is.na(genes$cosmx_gene)), " of ", nrow(genes), " Xenium panel genes found on the CosMx panel (",
        paste(names(table(genes$match)), table(genes$match), sep = " = ", collapse = ", "), "); ",
        sum(genes$use), " used for 1:1 comparisons.")
not_found <- genes$xenium_gene[is.na(genes$cosmx_gene)]
if (length(not_found) > 0) message("  Xenium genes without a CosMx partner (custom probes? other symbol? -> config/gene_name_map.csv): ",
                                   paste(utils::head(not_found, 30), collapse = ", "), if (length(not_found) > 30) " ...")
if (!file.exists(XENIUM_PANEL_CSV)) {
  utils::write.csv(data.frame(gene = genes$xenium_gene, cosmx_gene = genes$cosmx_gene), XENIUM_PANEL_CSV, row.names = FALSE)
  message("Wrote ", XENIUM_PANEL_CSV, " (05 uses it to mark the Xenium-panel genes; re-run 05 to get those columns).")
}
g <- genes[genes$use, ]
stopifnot("Fewer than 20 shared genes" = nrow(g) >= 20)

## cells of the shared cores; shared-gene matrices with IDENTICAL rows (Xenium names)
xen$meta <- xen$meta[xen$meta$core_id %in% shared_cores, ]; cos$meta <- cos$meta[cos$meta$core_id %in% shared_cores, ]
X <- xen$counts[g$xenium_gene, xen$meta$cell, drop = FALSE]
C <- cos$counts[g$cosmx_gene, cos$meta$cell, drop = FALSE]; rownames(C) <- g$xenium_gene
xen$meta$n_shared <- Matrix::colSums(X); cos$meta$n_shared <- Matrix::colSums(C)
n_cos_genes <- nrow(cos$counts); n_xen_genes <- nrow(xen$counts)
cos_all <- cos$counts[, cos$meta$cell, drop = FALSE]        # whole transcriptome, for 2a
xen$counts <- NULL; cos$counts <- NULL; invisible(gc())
xm <- xen$meta; cm <- cos$meta
rownames(xm) <- xm$cell; rownames(cm) <- cm$cell

## ---- registration: CosMx cells into the Xenium frame, core by core ------------------
## Slide level: similarity transform (rotation, scale, shift, mirror image if the
## section was flipped) from the centres of the shared cores. Core level: see
## register_core() in R/compare_helpers.R. A core counts as registered when the
## two cell-density maps correlate (CMP_MIN_NCC); only those enter 3d and 4d,
## and for them sections 1-2 use the region covered by both platforms.
cm$x_reg <- NA_real_; cm$y_reg <- NA_real_
reg_qc <- NULL
if (isTRUE(CMP_REGISTER)) {
  centre <- function(d) d |> group_by(core_id) |> summarise(x = stats::median(x), y = stats::median(y), .groups = "drop") |> as.data.frame()
  cen_c <- centre(cm); cen_x <- centre(xm)
  reg_rows <- list()
  for (sid in unique(cm$slide_id)) {
    cores_s <- intersect(unique(cm$core_id[cm$slide_id == sid]), shared_cores)
    if (length(cores_s) < 3) {
      message(sid, ": fewer than 3 shared cores -> no slide-level transform, registration skipped for this slide.")
      next
    }
    src <- as.matrix(cen_c[match(cores_s, cen_c$core_id), c("x", "y")]); dst <- as.matrix(cen_x[match(cores_s, cen_x$core_id), c("x", "y")])
    gfit <- robust_similarity(src, dst)
    G <- gfit$M
    message(sprintf("%s: slide transform from %d cores -- rotation %.1f deg, scale %.3f, median residual %.0f um%s%s",
                    sid, length(cores_s), rot_deg(G), sqrt(abs(det(G[1:2, 1:2]))), stats::median(gfit$residual[gfit$keep]),
                    if (det(G[1:2, 1:2]) < 0) "  ** MIRRORED section **" else "",
                    if (any(!gfit$keep)) paste0("; outlier cores: ", paste(cores_s[!gfit$keep], collapse = ", ")) else ""))
    if (stats::median(gfit$residual[gfit$keep]) > 400) {
      message("  WARNING: large residuals -- the core labels of the two projects probably disagree (compare both 02_dearray plots).")
    }
    for (cid in cores_s) {
      ci <- which(cm$core_id == cid); xi <- which(xm$core_id == cid)
      mov <- apply_tf(G, cm[ci, c("x", "y")]); fix <- as.matrix(xm[xi, c("x", "y")])
      r <- register_core(mov, fix, cm$lineage[ci], xm$lineage[xi])
      ok <- !is.na(r$ncc) && r$ncc >= CMP_MIN_NCC
      reg <- apply_tf(r$M, mov)
      lift <- lineage_lift(reg, cm$lineage[ci], fix, xm$lineage[xi])
      if (ok) { cm$x_reg[ci] <- reg[, 1]; cm$y_reg[ci] <- reg[, 2] }
      reg_rows[[cid]] <- data.frame(core_id = cid, slide_id = sid, n_cosmx = length(ci), n_xenium = length(xi),
                                    ncc_slide_transform = round(r$ncc_start, 3), ncc_registered = round(r$ncc, 3),
                                    rotation_deg = round(r$rotation_deg, 1), shift_um = round(r$shift_um),
                                    lineage_agreement = round(lift[["agreement"]], 3), lineage_lift = round(lift[["lift"]], 2),
                                    registered = ok)
    }
  }
  reg_qc <- dplyr::bind_rows(reg_rows)
  if (nrow(reg_qc) > 0) {
    save_table(reg_qc, "10_registration_qc")
    message(sum(reg_qc$registered), " of ", nrow(reg_qc), " cores registered (density correlation >= ", CMP_MIN_NCC, ").")
    ## the plot to look at: Xenium cells grey, registered CosMx cells coloured by lineage
    show <- utils::head(reg_qc$core_id[order(reg_qc$ncc_registered)], 30)       # worst first
    dd <- rbind(data.frame(core_id = xm$core_id, x = xm$x, y = xm$y, lineage = "Xenium cells")[xm$core_id %in% show, ],
                data.frame(core_id = cm$core_id, x = cm$x_reg, y = cm$y_reg, lineage = cm$lineage)[cm$core_id %in% show & !is.na(cm$x_reg), ])
    dd$core_id <- factor(dd$core_id, levels = show)
    lab <- stats::setNames(sprintf("%s  r=%.2f", reg_qc$core_id, reg_qc$ncc_registered), reg_qc$core_id)
    p_reg <- ggplot(dd, aes(x, y, colour = lineage)) +
      geom_point(data = dd[dd$lineage == "Xenium cells", ], size = 0.5, stroke = 0, colour = "grey75") +
      geom_point(data = dd[dd$lineage != "Xenium cells", ], size = 0.35, stroke = 0) +
      facet_wrap(~ core_id, scales = "free", labeller = labeller(core_id = lab)) + scale_y_reverse() +
      guides(colour = guide_legend(override.aes = list(size = 3))) + theme_void() +
      theme(aspect.ratio = 1, strip.text = element_text(size = 7), legend.position = "bottom") +
      ggtitle("Registration: Xenium cells (grey) and registered CosMx cells (by lineage); least well aligned cores first")
    save_fig(p_reg, "10_registration_check", width = 16, height = 3 * ceiling(length(show) / 6) + 1.5)
  }
}
registered_cores <- if (!is.null(reg_qc) && nrow(reg_qc) > 0) reg_qc$core_id[reg_qc$registered] else character(0)

## Region covered by BOTH platforms (registered cores): bins of CMP_BIN_UM in the
## Xenium frame holding cells of both. On CosMx only the FOVs are imaged, on
## Xenium the whole core; core-level numbers are only fair on the common region.
xm$bin <- paste(xm$core_id, floor(xm$x / CMP_BIN_UM), floor(xm$y / CMP_BIN_UM))
cm$bin <- ifelse(is.na(cm$x_reg), NA_character_, paste(cm$core_id, floor(cm$x_reg / CMP_BIN_UM), floor(cm$y_reg / CMP_BIN_UM)))
both_bins <- intersect(xm$bin, cm$bin[!is.na(cm$bin)])
xm$use <- !(xm$core_id %in% registered_cores) | xm$bin %in% both_bins
cm$use <- !(cm$core_id %in% registered_cores) | cm$bin %in% both_bins
core_info$coverage_matched <- core_info$core_id %in% registered_cores
message("Cells used for core-level comparisons: Xenium ", sum(xm$use), " of ", nrow(xm), ", CosMx ", sum(cm$use), " of ", nrow(cm),
        " (registered cores: common region only; other cores: all cells, coverage may differ).")
core_tab$registered <- core_tab$core_id %in% registered_cores
save_table(core_tab, "10_cores")

## =============================================================================
## 1. TRANSCRIPT QUANTIFICATION
## =============================================================================
xu <- xm[xm$use, ]; cu <- cm[cm$use, ]
pbX <- as.matrix(pseudobulk(X[, xu$cell, drop = FALSE], factor(xu$core_id, levels = shared_cores)))   # genes x cores
pbC <- as.matrix(pseudobulk(C[, cu$cell, drop = FALSE], factor(cu$core_id, levels = shared_cores)))

## ---- 1a. per core -----------------------------------------------------------------
per_core <- function(d, platform) {
  d |> group_by(core_id) |>
    summarise(n_cells = n(), tissue_mm2 = occupied_area(x, y) / 1e6,
              cells_per_mm2 = n_cells / tissue_mm2,
              median_cell_area = stats::median(cell_area, na.rm = TRUE),
              shared_counts_per_cell = mean(n_shared), shared_genes_counts_median = stats::median(n_shared),
              all_counts_per_cell = mean(n_counts), genes_per_cell = mean(n_genes),
              shared_counts_per_100um2_tissue = sum(n_shared) / (tissue_mm2 * 1e4),
              all_counts_per_100um2_tissue = sum(n_counts) / (tissue_mm2 * 1e4),
              all_counts_per_um2_cell = sum(n_counts) / sum(cell_area, na.rm = TRUE),
              .groups = "drop") |>
    mutate(platform = platform)
}
core_long <- rbind(per_core(xu, "Xenium"), per_core(cu, "CosMx"))
core_wide <- tidyr::pivot_wider(core_long, id_cols = core_id, names_from = platform,
                                values_from = setdiff(colnames(core_long), c("core_id", "platform"))) |>
  left_join(core_info, by = "core_id")
save_table(as.data.frame(core_wide), "10_quant_per_core")
p_1a <- (paired_core_plot(core_wide, "cells_per_mm2_Xenium", "cells_per_mm2_CosMx", "cells per mm2 of tissue") |
           paired_core_plot(core_wide, "median_cell_area_Xenium", "median_cell_area_CosMx", "median cell area (um2)") |
           paired_core_plot(core_wide, "shared_counts_per_cell_Xenium", "shared_counts_per_cell_CosMx",
                            "shared-gene transcripts per cell", log = TRUE) |
           paired_core_plot(core_wide, "shared_counts_per_100um2_tissue_Xenium", "shared_counts_per_100um2_tissue_CosMx",
                            "shared-gene transcripts per 100 um2 of tissue", log = TRUE)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = sprintf("1a. Per core, %d shared genes (dashed = equal)", nrow(g)),
                  subtitle = "CosMx also measures the rest of the transcriptome: see all_counts_* in tables/10_quant_per_core.csv")
save_fig(p_1a, "10_1a_quant_per_core", width = 20, height = 5.5)

## ---- 1b. per core: gene profile, platform against platform -------------------------
## mean counts per cell of every shared gene; correlation across genes
ncell_x <- as.numeric(table(factor(xu$core_id, levels = shared_cores))); ncell_c <- as.numeric(table(factor(cu$core_id, levels = shared_cores)))
mX <- sweep(pbX, 2, ncell_x, "/"); mC <- sweep(pbC, 2, ncell_c, "/")
pc <- 0.01
prof_cor <- data.frame(core_id = shared_cores,
                       spearman = vapply(shared_cores, function(k) stats::cor(mX[, k], mC[, k], method = "spearman"), numeric(1)),
                       pearson_log = vapply(shared_cores, function(k) stats::cor(log10(mX[, k] + pc), log10(mC[, k] + pc)), numeric(1)),
                       median_log2_ratio = vapply(shared_cores, function(k) stats::median(log2((mC[, k] + pc) / (mX[, k] + pc))), numeric(1))) |>
  left_join(core_info, by = "core_id")
save_table(prof_cor, "10_quant_profile_correlation_per_core")
pool <- data.frame(gene = rownames(mX), xenium = rowMeans(mX), cosmx = rowMeans(mC))     # every core weighs the same
p_1b <- (ggplot(pool, aes(xenium + pc, cosmx + pc)) +
           geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") + geom_point(alpha = 0.5, size = 1.2) +
           scale_x_log10() + scale_y_log10() + coord_fixed() + theme_bw() +
           labs(x = "Xenium: mean transcripts per cell", y = "CosMx: mean transcripts per cell",
                title = sprintf("all cores pooled: Spearman %.2f", stats::cor(pool$xenium, pool$cosmx, method = "spearman")))) |
  (ggplot(prof_cor, aes(organ, spearman, colour = slide_id)) + geom_boxplot(outlier.shape = NA, colour = "grey60") +
     geom_jitter(width = 0.15, size = 2) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     labs(x = NULL, y = "Spearman r across shared genes", title = "one dot per core"))
save_fig(p_1b + plot_annotation(title = "1b. Gene expression profile, CosMx against Xenium"), "10_1b_profile_correlation", width = 14, height = 6)

## ---- 1c. per gene: sensitivity ratio ---------------------------------------------------
## log2(CosMx / Xenium) of mean counts per cell, per core; cores where the gene is
## essentially absent on both platforms say nothing and are left out.
lr <- log2((mC + pc) / (mX + pc)); lr[(pbX + pbC) < 20] <- NA
gene_tab <- data.frame(
  gene = rownames(mX), cosmx_gene = g$cosmx_gene, match = g$match,
  xenium_counts_per_cell = rowMeans(mX), cosmx_counts_per_cell = rowMeans(mC),
  log2_ratio_cosmx_vs_xenium = apply(lr, 1, stats::median, na.rm = TRUE),
  log2_ratio_iqr = apply(lr, 1, stats::IQR, na.rm = TRUE),
  n_cores_informative = rowSums(!is.na(lr)))
p_1c <- (ggplot(gene_tab, aes(log2_ratio_cosmx_vs_xenium)) + geom_histogram(bins = 50, fill = "grey40") +
           geom_vline(xintercept = 0, linetype = 2) +
           geom_vline(xintercept = stats::median(gene_tab$log2_ratio_cosmx_vs_xenium, na.rm = TRUE), colour = "#e34948") +
           theme_bw() + labs(x = "log2 (CosMx / Xenium) transcripts per cell, median over cores", y = "genes",
                             title = sprintf("median gene: CosMx detects %.2f x the Xenium count",
                                             2^stats::median(gene_tab$log2_ratio_cosmx_vs_xenium, na.rm = TRUE)))) |
  (ggplot(gene_tab, aes(xenium_counts_per_cell + pc, log2_ratio_cosmx_vs_xenium)) +
     geom_hline(yintercept = 0, linetype = 2) + geom_point(alpha = 0.5, size = 1.2) +
     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "#e34948", linewidth = 0.7) +
     ggrepel::geom_text_repel(data = gene_tab[rank(-abs(gene_tab$log2_ratio_cosmx_vs_xenium -
                                                           stats::median(gene_tab$log2_ratio_cosmx_vs_xenium, na.rm = TRUE))) <= 15, ],
                              aes(label = gene), size = 2.6, max.overlaps = Inf) +
     scale_x_log10() + theme_bw() +
     labs(x = "Xenium: mean transcripts per cell", y = "log2 (CosMx / Xenium)", title = "ratio against abundance (labels: 15 most extreme genes)"))
save_fig(p_1c + plot_annotation(title = "1c. Sensitivity per gene"), "10_1c_gene_sensitivity", width = 15, height = 6)

## ---- 1d. per gene: agreement ACROSS CORES ----------------------------------------------
## Each platform's own scale is removed (CPM within the shared genes, per core),
## then: does gene g vary over the cores the same way on both? This is the
## quantification question that matters for biology, and it is independent of
## how many probes or how much efficiency each platform has for the gene.
cpmX <- sweep(pbX, 2, pmax(colSums(pbX), 1), "/") * 1e6; cpmC <- sweep(pbC, 2, pmax(colSums(pbC), 1), "/") * 1e6
if (length(shared_cores) >= 6) {
  gene_tab$across_core_spearman <- vapply(seq_len(nrow(cpmX)), function(i) {
    if (stats::sd(cpmX[i, ]) == 0 || stats::sd(cpmC[i, ]) == 0) NA_real_ else stats::cor(cpmX[i, ], cpmC[i, ], method = "spearman")
  }, numeric(1))
  ## what to expect from noise alone is not 1: a gene with few counts cannot agree
  gene_tab$abundance_bin <- cut(rank(gene_tab$xenium_counts_per_cell) / nrow(gene_tab), c(0, 0.25, 0.5, 0.75, 1),
                                labels = c("lowest 25%", "25-50%", "50-75%", "top 25%"))
  p_1d <- ggplot(gene_tab, aes(abundance_bin, across_core_spearman)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.size = 0.5, fill = "grey90") + theme_bw() +
    labs(x = "abundance of the gene on Xenium", y = "Spearman r across cores (CosMx vs Xenium)",
         title = sprintf("1d. Does a gene vary across the %d cores the same way on both platforms? median r = %.2f",
                         length(shared_cores), stats::median(gene_tab$across_core_spearman, na.rm = TRUE)),
         subtitle = "low r for rare genes is noise; low r for abundant genes is a real disagreement -- look at those genes")
  save_fig(p_1d, "10_1d_across_core_concordance", width = 10, height = 6)
} else message("1d skipped: fewer than 6 shared cores.")
save_table(gene_tab[order(gene_tab$log2_ratio_cosmx_vs_xenium), ], "10_quant_per_gene")

## ---- 1e. per core x lineage ---------------------------------------------------------------
grp_x <- paste(xu$core_id, xu$lineage, sep = "||"); grp_c <- paste(cu$core_id, cu$lineage, sep = "||")
ok_grp <- intersect(names(which(table(grp_x) >= CMP_MIN_CELLS_GROUP)), names(which(table(grp_c) >= CMP_MIN_CELLS_GROUP)))
ok_grp <- ok_grp[!grepl("\\|\\|NA$", ok_grp)]
lin_cor <- NULL
if (length(ok_grp) > 0) {
  gX <- as.matrix(pseudobulk(X[, xu$cell[grp_x %in% ok_grp], drop = FALSE], factor(grp_x[grp_x %in% ok_grp], levels = ok_grp)))
  gC <- as.matrix(pseudobulk(C[, cu$cell[grp_c %in% ok_grp], drop = FALSE], factor(grp_c[grp_c %in% ok_grp], levels = ok_grp)))
  gX <- sweep(gX, 2, as.numeric(table(factor(grp_x, levels = ok_grp))), "/"); gC <- sweep(gC, 2, as.numeric(table(factor(grp_c, levels = ok_grp))), "/")
  lin_cor <- data.frame(group = ok_grp, core_id = sub("\\|\\|.*$", "", ok_grp), lineage = sub("^.*\\|\\|", "", ok_grp),
                        spearman = vapply(ok_grp, function(k) stats::cor(gX[, k], gC[, k], method = "spearman"), numeric(1)),
                        median_log2_ratio = vapply(ok_grp, function(k) stats::median(log2((gC[, k] + pc) / (gX[, k] + pc))), numeric(1)))
  save_table(lin_cor, "10_quant_profile_correlation_per_core_lineage")
  p_1e <- ggplot(lin_cor, aes(lineage, spearman)) + geom_boxplot(outlier.shape = NA, fill = "grey90") +
    geom_jitter(width = 0.15, size = 1.3, alpha = 0.7) + theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = NULL, y = "Spearman r across shared genes", title = "1e. Profile agreement within a lineage (one dot = one core x lineage)",
         subtitle = "uses each platform's own lineage calls: a lineage that agrees poorly is either measured or annotated differently")
  save_fig(p_1e, "10_1e_profile_correlation_by_lineage", width = 10, height = 6)
}

## =============================================================================
## 2. ACCURACY (specificity / background)
## =============================================================================
## ---- 2a. detection above the negative-probe background ------------------------------------
## Same rule as 05: lambda = counts one negative probe collects in the core;
## a gene is detected if it is >= 2^SBR_LOG2_MIN x lambda and significantly above it.
detect <- function(pb, d) {
  bg <- d |> group_by(core_id) |> summarise(ctrl = sum(n_ctrl), n_feat = first(n_ctrl_feat), .groups = "drop")
  lambda <- stats::setNames(bg$ctrl / pmax(bg$n_feat, 1), bg$core_id)[colnames(pb)]
  L <- matrix(lambda, nrow(pb), ncol(pb), byrow = TRUE)
  p <- stats::ppois(pb - 1, L, lower.tail = FALSE)
  padj <- matrix(stats::p.adjust(p, "BH"), nrow(pb))
  list(detected = log2((pb + 0.5) / (L + 0.5)) >= SBR_LOG2_MIN & padj < SBR_PADJ_MAX, lambda = lambda)
}
dX <- detect(pbX, xu); dC <- detect(pbC, cu)
dimnames(dX$detected) <- dimnames(pbX); dimnames(dC$detected) <- dimnames(pbC)
## whole transcriptome on CosMx: how many genes are above background at all?
pbC_all <- as.matrix(pseudobulk(cos_all[, cu$cell, drop = FALSE], factor(cu$core_id, levels = shared_cores)))
dC_all <- detect(pbC_all, cu)
rm(cos_all); invisible(gc())
det_core <- data.frame(
  core_id = shared_cores,
  pct_shared_detected_xenium = 100 * colMeans(dX$detected), pct_shared_detected_cosmx = 100 * colMeans(dC$detected),
  neg_per_probe_per_cell_xenium = dX$lambda / ncell_x, neg_per_probe_per_cell_cosmx = dC$lambda / ncell_c,
  ## mean gene / negative probe: how far the average shared gene is above background
  signal_to_background_xenium = colMeans(pbX) / pmax(dX$lambda, 0.5), signal_to_background_cosmx = colMeans(pbC) / pmax(dC$lambda, 0.5),
  n_genes_detected_cosmx_whole_transcriptome = colSums(dC_all$detected), n_genes_cosmx_panel = n_cos_genes,
  n_genes_detected_xenium_panel = colSums(dX$detected), n_genes_xenium_shared = nrow(g)) |>
  left_join(core_info, by = "core_id")
save_table(det_core, "10_accuracy_detection_per_core")
gene_det <- data.frame(gene = rownames(pbX), frac_cores_detected_xenium = rowMeans(dX$detected), frac_cores_detected_cosmx = rowMeans(dC$detected))
gene_det$class <- with(gene_det, ifelse(frac_cores_detected_xenium >= 0.5 & frac_cores_detected_cosmx >= 0.5, "both",
                                 ifelse(frac_cores_detected_xenium >= 0.5, "Xenium only",
                                 ifelse(frac_cores_detected_cosmx >= 0.5, "CosMx only", "neither"))))
save_table(gene_det, "10_accuracy_detection_per_gene")
message("2a. Shared genes detected in >= half the cores: ", paste(names(table(gene_det$class)), table(gene_det$class), sep = " = ", collapse = ", "))
p_2a <- (paired_core_plot(det_core, "pct_shared_detected_xenium", "pct_shared_detected_cosmx", "% of shared genes above background") |
           paired_core_plot(det_core, "signal_to_background_xenium", "signal_to_background_cosmx",
                            "mean shared gene / one negative probe", log = TRUE) |
           paired_core_plot(det_core, "neg_per_probe_per_cell_xenium", "neg_per_probe_per_cell_cosmx",
                            "negative-probe counts per probe per cell", log = TRUE)) +
  plot_layout(guides = "collect") + plot_annotation(title = "2a. Signal over background per core (dashed = equal)")
save_fig(p_2a, "10_2a_detection_background", width = 16, height = 5.5)

## ---- 2b. MECR on the shared lineage markers --------------------------------------------------
markers <- read_marker_panel(MARKERS_RNA_CSV)
mecr_sets <- lapply(markers$mecr, function(v) intersect(match_case(v, rownames(X)), rownames(X)))
mecr_sets <- mecr_sets[lengths(mecr_sets) > 0]
mecr_tab <- NULL
if (length(mecr_sets) >= 2) {
  mecr_tab <- data.frame(core_id = shared_cores,
                         mecr_xenium = vapply(shared_cores, function(k) compute_mecr(X[, xu$cell[xu$core_id == k], drop = FALSE], mecr_sets), numeric(1)),
                         mecr_cosmx = vapply(shared_cores, function(k) compute_mecr(C[, cu$cell[cu$core_id == k], drop = FALSE], mecr_sets), numeric(1))) |>
    left_join(core_info, by = "core_id")
  save_table(mecr_tab, "10_accuracy_mecr_per_core")
  p_2b <- paired_core_plot(mecr_tab, "mecr_xenium", "mecr_cosmx", "2b. MECR per core: marker pairs of different lineages in one cell") +
    labs(subtitle = paste0("lower = less spill-over between cells. Same ", length(unlist(mecr_sets)), " genes on both platforms.\n",
                           "More counts per cell raise MECR by themselves: read together with 1a."))
  save_fig(p_2b, "10_2b_mecr", width = 7.5, height = 6.5)
} else message("2b skipped: fewer than 2 lineages with MECR markers among the shared genes.")

## ---- 2c. on-target enrichment of lineage genes -------------------------------------------------
probes <- read_config_csv(PROBE_EXPECT_CSV, stringsAsFactors = FALSE)
probes$gene <- match_case(as_seurat_features(probes$gene), rownames(X))
probes <- probes[probes$gene %in% rownames(X) & !is.na(probes$expected_lineage) & nzchar(probes$expected_lineage), ]
spec_tab <- NULL
if (nrow(probes) > 0) {
  on_target <- function(M, d) {
    vapply(seq_len(nrow(probes)), function(j) {
      on <- d$lineage == probes$expected_lineage[j] & !is.na(d$lineage)
      v <- M[probes$gene[j], d$cell]
      (sum(v[on]) / max(sum(v), 1)) / (sum(d$n_shared[on]) / max(sum(d$n_shared), 1))
    }, numeric(1))
  }
  spec_tab <- data.frame(gene = probes$gene, expected_lineage = probes$expected_lineage,
                         enrichment_xenium = on_target(X, xu), enrichment_cosmx = on_target(C, cu))
  save_table(spec_tab, "10_accuracy_lineage_specificity")
  p_2c <- ggplot(spec_tab, aes(enrichment_xenium, enrichment_cosmx, colour = expected_lineage)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") + geom_point(size = 2.5) +
    ggrepel::geom_text_repel(aes(label = gene), size = 3, show.legend = FALSE) +
    scale_x_log10() + scale_y_log10() + coord_fixed() + theme_bw() +
    labs(x = "Xenium", y = "CosMx", title = "2c. On-target enrichment of lineage genes",
         subtitle = "share of the gene's counts in its expected lineage / that lineage's share of all shared-gene counts (1 = no specificity)")
  save_fig(p_2c, "10_2c_lineage_specificity", width = 8.5, height = 7)
}

## =============================================================================
## 3. CELL ANNOTATION
## =============================================================================
## ---- 3a. lineage composition per core -----------------------------------------------------------
comp <- rbind(data.frame(platform = "Xenium", core_id = xu$core_id, lineage = xu$lineage),
              data.frame(platform = "CosMx", core_id = cu$core_id, lineage = cu$lineage)) |>
  filter(!is.na(lineage)) |> count(platform, core_id, lineage) |>
  tidyr::complete(platform, core_id, lineage, fill = list(n = 0)) |>
  group_by(platform, core_id) |> mutate(frac = n / sum(n)) |> ungroup()
comp_w <- tidyr::pivot_wider(comp[, c("platform", "core_id", "lineage", "frac")], names_from = platform, values_from = frac) |>
  left_join(core_info, by = "core_id")
save_table(as.data.frame(comp_w), "10_annotation_lineage_composition")
js <- comp_w |> group_by(core_id) |> summarise(jensen_shannon = jensen_shannon(Xenium + 1e-9, CosMx + 1e-9), .groups = "drop")
lin_r <- comp_w |> group_by(lineage) |> summarise(r = suppressWarnings(stats::cor(Xenium, CosMx)), .groups = "drop")
p_3a <- ggplot(comp_w, aes(Xenium, CosMx, colour = organ)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") + geom_point(size = 2) +
  geom_text(data = lin_r, aes(x = -Inf, y = Inf, label = sprintf(" r = %.2f", r)), inherit.aes = FALSE, hjust = 0, vjust = 1.5, size = 3.3) +
  facet_wrap(~ lineage, scales = "free") + theme_bw() + theme(aspect.ratio = 1) +
  labs(x = "fraction of the core's cells, Xenium", y = "fraction of the core's cells, CosMx",
       title = "3a. Lineage composition per core (each platform's own annotation)",
       subtitle = sprintf("median Jensen-Shannon divergence between the two compositions of a core: %.3f (0 = identical, 1 = disjoint)",
                          stats::median(js$jensen_shannon)))
save_fig(p_3a, "10_3a_lineage_composition", width = 13, height = 9)

## ---- 3b. label transfer on the shared genes, both directions --------------------------------------
## All cells of the shared cores (not only the common region): more cells per label.
lin_of <- function(d) { t <- table(d$cell_type, d$lineage); stats::setNames(colnames(t)[max.col(t, ties.method = "first")], rownames(t)) }
tr_c <- transfer_labels(X[, xm$cell], xm$cell_type, C[, cm$cell])        # CosMx cells <- Xenium cell types
tr_x <- transfer_labels(C[, cm$cell], cm$cell_type, X[, xm$cell])        # Xenium cells <- CosMx cell types
cm$xenium_cell_type <- tr_c$label; cm$xenium_lineage <- unname(lin_of(xm)[tr_c$label]); cm$transfer_score <- tr_c$score
xm$cosmx_cell_type <- tr_x$label; xm$cosmx_lineage <- unname(lin_of(cm)[tr_x$label]); xm$transfer_score <- tr_x$score
data.table::fwrite(cm[, c("cell", "slide_id", "core_id", "cell_type", "lineage", "xenium_cell_type", "xenium_lineage", "transfer_score")],
                   file.path(TABLES_DIR, "10_annotation_label_transfer_cosmx_cells.csv.gz"))
data.table::fwrite(xm[, c("cell", "slide_id", "core_id", "cell_type", "lineage", "cosmx_cell_type", "cosmx_lineage", "transfer_score")],
                   file.path(TABLES_DIR, "10_annotation_label_transfer_xenium_cells.csv.gz"))
agree_core <- rbind(
  cm |> group_by(core_id) |> summarise(direction = "CosMx cells, Xenium labels", lineage_agreement = mean(lineage == xenium_lineage, na.rm = TRUE), .groups = "drop"),
  xm |> group_by(core_id) |> summarise(direction = "Xenium cells, CosMx labels", lineage_agreement = mean(lineage == cosmx_lineage, na.rm = TRUE), .groups = "drop")) |>
  left_join(core_info, by = "core_id")
save_table(as.data.frame(agree_core), "10_annotation_label_transfer_agreement_per_core")
message(sprintf("3b. Lineage of the transferred label = the cell's own lineage: %.1f%% of CosMx cells, %.1f%% of Xenium cells.",
                100 * mean(cm$lineage == cm$xenium_lineage, na.rm = TRUE), 100 * mean(xm$lineage == xm$cosmx_lineage, na.rm = TRUE)))
confusion <- function(own, other, own_name, other_name) {
  t <- as.data.frame(table(own = own, other = other), stringsAsFactors = FALSE)
  t <- t |> group_by(own) |> mutate(frac = Freq / sum(Freq)) |> ungroup()
  ## order rows / columns so that the best matches sit on a diagonal
  mat <- tidyr::pivot_wider(t[, c("own", "other", "frac")], names_from = other, values_from = frac)
  m <- as.matrix(mat[, -1]); rownames(m) <- mat$own
  ord_c <- colnames(m)[order(apply(m, 2, which.max))]
  p <- ggplot(t, aes(factor(other, levels = ord_c), factor(own, levels = rev(rownames(m)[order(max.col(m[, ord_c, drop = FALSE]))])), fill = frac)) +
    geom_tile() + geom_text(aes(label = ifelse(frac >= 0.1, sprintf("%.0f", 100 * frac), "")), size = 2.6) +
    scale_fill_gradient(low = "white", high = "#2a78d6", limits = c(0, 1), name = "fraction\nof row") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = other_name, y = own_name)
  list(table = t, plot = p)
}
cf_c <- confusion(cm$cell_type, cm$xenium_cell_type, "CosMx cell type (own annotation)", "most similar Xenium cell type")
cf_x <- confusion(xm$cell_type, xm$cosmx_cell_type, "Xenium cell type (own annotation)", "most similar CosMx cell type")
save_table(rbind(cbind(direction = "CosMx cells, Xenium labels", cf_c$table), cbind(direction = "Xenium cells, CosMx labels", cf_x$table)),
           "10_annotation_label_transfer_confusion")
n_types <- max(length(unique(cm$cell_type)), length(unique(xm$cell_type)))
save_fig((cf_c$plot | cf_x$plot) +
           plot_annotation(title = "3b. Label transfer on the shared genes (numbers = % of the row)",
                           subtitle = paste("One dominant cell per row = the cell type exists on both platforms. A row spread over several columns =",
                                            "the other platform splits it (or cannot tell it apart).")),
         "10_3b_label_transfer_confusion", width = max(14, 0.9 * n_types + 8), height = max(6, 0.35 * n_types + 3))

## ---- 3c. lineage profiles, platform against platform ------------------------------------------------
lins <- intersect(names(which(table(xm$lineage) >= CMP_MIN_CELLS_GROUP)), names(which(table(cm$lineage) >= CMP_MIN_CELLS_GROUP)))
if (length(lins) >= 2) {
  prof <- function(M, d) {
    pb <- as.matrix(pseudobulk(M[, d$cell[d$lineage %in% lins], drop = FALSE], factor(d$lineage[d$lineage %in% lins], levels = lins)))
    z <- log1p(sweep(pb, 2, colSums(pb), "/") * 1e4)
    t(scale(t(z)))                                    # per gene, across lineages, within the platform
  }
  cc <- stats::cor(prof(X, xm), prof(C, cm), use = "pairwise.complete.obs")
  cc_df <- as.data.frame(as.table(cc)); colnames(cc_df) <- c("xenium_lineage", "cosmx_lineage", "r")
  save_table(cc_df, "10_annotation_lineage_profile_correlation")
  p_3c <- ggplot(cc_df, aes(cosmx_lineage, xenium_lineage, fill = r)) + geom_tile() +
    geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) + theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = "CosMx lineage", y = "Xenium lineage", title = "3c. Lineage expression profiles (shared genes), platform against platform",
         subtitle = "a clear diagonal = the same lineage name means the same expression programme on both platforms")
  save_fig(p_3c, "10_3c_lineage_profile_correlation", width = 8.5, height = 7)
}

## ---- 3d. registered cores: neighbours across sections ------------------------------------------------
if (!is.null(reg_qc) && any(reg_qc$registered)) {
  p_3d <- ggplot(reg_qc[reg_qc$registered, ], aes(lineage_lift)) + geom_histogram(bins = 30, fill = "grey40") +
    geom_vline(xintercept = 1, linetype = 2) + theme_bw() +
    labs(x = "lineage agreement of neighbouring cells across sections / chance", y = "cores",
         title = "3d. Do the two annotations put the same lineage in the same place?",
         subtitle = "nearest Xenium cell within 20 um of every registered CosMx cell; 1 = no better than chance")
  save_fig(p_3d, "10_3d_spatial_lineage_agreement", width = 8, height = 5)
}

## =============================================================================
## 4. PROBE / OPTICAL CROWDING
## =============================================================================
## ---- 4a. transcript density per cell -----------------------------------------------------------------
dens <- rbind(data.frame(platform = "Xenium", lineage = xm$lineage, core_id = xm$core_id, density = xm$n_counts / xm$cell_area),
              data.frame(platform = "CosMx", lineage = cm$lineage, core_id = cm$core_id, density = cm$n_counts / cm$cell_area))
dens <- dens[is.finite(dens$density) & !is.na(dens$lineage), ]
dens_tab <- dens |> group_by(platform, lineage) |>
  summarise(n_cells = n(), median = stats::median(density), p90 = stats::quantile(density, 0.9), p99 = stats::quantile(density, 0.99),
            p99.9 = stats::quantile(density, 0.999), max = max(density), .groups = "drop")
save_table(as.data.frame(dens_tab), "10_crowding_density_per_cell")
p_4a <- ggplot(dens, aes(density, colour = platform)) + geom_density(linewidth = 0.7) +
  scale_x_log10() + scale_colour_manual(values = PLATFORM_COLOURS) + facet_wrap(~ lineage, scales = "free_y") + theme_bw() +
  labs(x = "detected transcripts per um2 of cell (all genes of the platform's panel)", y = "density of cells",
       title = "4a. How crowded are the cells?",
       subtitle = "a distribution that piles up against an upper limit suggests saturation; compare p99 / median in tables/10_crowding_density_per_cell.csv")
save_fig(p_4a, "10_4a_density_per_cell", width = 13, height = 8)

## ---- 4b. gene level: compression of highly expressed genes ------------------------------------------------
## log-log slope of CosMx against Xenium abundance. 1 = proportional. < 1: the more
## abundant a gene, the less of it CosMx sees relative to Xenium (CosMx compresses
## the top); > 1: Xenium compresses. Total least squares: both axes are noisy.
ok_g <- pool$xenium > 0 & pool$cosmx > 0
slope <- tls_slope(log10(pool$xenium[ok_g]), log10(pool$cosmx[ok_g]))
slope_ci <- boot_ci(log10(pool$xenium[ok_g]), log10(pool$cosmx[ok_g]))
slope_core <- vapply(shared_cores, function(k) { o <- mX[, k] > 0 & mC[, k] > 0; tls_slope(log10(mX[o, k]), log10(mC[o, k])) }, numeric(1))
save_table(data.frame(core_id = shared_cores, loglog_slope_cosmx_vs_xenium = slope_core) |> left_join(core_info, by = "core_id"),
           "10_crowding_gene_slope_per_core")
p_4b <- (ggplot(pool[ok_g, ], aes(xenium, cosmx)) + geom_point(alpha = 0.5, size = 1.2) +
           geom_abline(slope = 1, intercept = stats::median(log10(pool$cosmx[ok_g]) - log10(pool$xenium[ok_g])), linetype = 2, colour = "grey40") +
           geom_abline(slope = slope, intercept = mean(log10(pool$cosmx[ok_g])) - slope * mean(log10(pool$xenium[ok_g])), colour = "#e34948") +
           scale_x_log10() + scale_y_log10() + coord_fixed() + theme_bw() +
           labs(x = "Xenium: mean transcripts per cell", y = "CosMx: mean transcripts per cell",
                title = sprintf("slope %.2f (95%% CI %.2f - %.2f)", slope, slope_ci[1], slope_ci[2]),
                subtitle = "red = fitted, dashed = proportional (slope 1)")) |
  (ggplot(data.frame(slope = slope_core), aes(slope)) + geom_histogram(bins = 25, fill = "grey40") +
     geom_vline(xintercept = 1, linetype = 2) + theme_bw() + labs(x = "log-log slope per core", y = "cores", title = "the same, core by core"))
save_fig(p_4b + plot_annotation(title = "4b. Are highly expressed genes compressed on one platform?",
                                subtitle = "slope < 1: CosMx under-counts abundant genes relative to Xenium; > 1: the reverse. Probe design contributes (see header)."),
         "10_4b_gene_level_compression", width = 14, height = 6.5)

## ---- 4c. cell level, no registration: the densest cells of a core -------------------------------------------
## Within one core x lineage both platforms sample the same cell population.
## If nothing saturates, the distributions of shared-gene counts per cell differ
## by a constant factor: quantile for quantile, CosMx / Xenium is flat. If dense
## cells lose transcripts on one platform, the ratio bends at the high quantiles.
## Shown relative to the ratio at the median, on a log2 scale.
probs <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
qq <- dplyr::bind_rows(lapply(ok_grp, function(k) {
  a <- xm$n_shared[paste(xm$core_id, xm$lineage, sep = "||") == k]; b <- cm$n_shared[paste(cm$core_id, cm$lineage, sep = "||") == k]
  if (length(a) < CMP_MIN_CELLS_GROUP || length(b) < CMP_MIN_CELLS_GROUP) return(NULL)
  qa <- stats::quantile(a, probs); qb <- stats::quantile(b, probs)
  if (qa[1] < 5 || qb[1] < 5) return(NULL)                    # too few counts: quantiles are shot noise
  data.frame(group = k, core_id = sub("\\|\\|.*$", "", k), lineage = sub("^.*\\|\\|", "", k), quantile = probs,
             log2_rel_ratio = log2(qb / qa) - log2(qb[1] / qa[1]))
}))
crowd_cell <- NULL
if (!is.null(qq) && nrow(qq) > 0) {
  save_table(qq, "10_crowding_cell_quantiles")
  crowd_cell <- qq |> filter(quantile == 0.99) |> group_by(core_id) |> summarise(v = mean(log2_rel_ratio), .groups = "drop")
  qsum <- qq |> group_by(lineage, quantile) |> summarise(med = stats::median(log2_rel_ratio), lo = stats::quantile(log2_rel_ratio, 0.25),
                                                         hi = stats::quantile(log2_rel_ratio, 0.75), .groups = "drop")
  p_4c <- ggplot(qsum, aes(factor(quantile), med, group = lineage, colour = lineage)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_line() + geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3) +
    theme_bw() +
    labs(x = "quantile of shared-gene transcripts per cell (within core x lineage)",
         y = "log2 (CosMx / Xenium), relative to the median cell",
         title = "4c. Do the cells with most RNA lose transcripts on one platform?",
         subtitle = "flat = no; falling = CosMx compresses its densest cells (crowding on CosMx); rising = Xenium does. Median and IQR over cores.")
  save_fig(p_4c, "10_4c_cell_level_quantiles", width = 11, height = 6.5)
}

## ---- 4d. bin level, registered cores: ratio against local transcript density -------------------------------------
## CMP_BIN_UM bins of the common region. Per bin:
##   response = log2( observed CosMx shared-gene counts per cell
##                    / expected from the Xenium counts of the same bin ),
##              where expected = sum over genes of (Xenium counts x that gene's
##              overall CosMx/Xenium ratio) -- so a bin full of genes that CosMx
##              detects poorly everywhere does not look like crowding.
##              Centred per core.
##   load     = CosMx transcripts (ALL genes) per um2 of cell area in the
##              NEIGHBOURING bins: the spot density the optics have to resolve
##              around there. Taken from the neighbours, not the bin itself, on
##              purpose: the two sections hold different cells, and a bin that
##              happens to catch RNA-rich cells on the CosMx section would have a
##              high load AND a high ratio by chance. Neighbouring bins share the
##              tissue (density is similar over ~100 um) but not that chance.
## slope per core = change of the response per doubling of the load.
## Negative = CosMx loses transcripts relative to Xenium where the tissue is dense.
## The same with the Xenium density of the neighbours as load is the cross-check
## (both loads measure the same tissue density; the two slopes should agree in sign).
crowd_bin <- NULL
if (length(registered_cores) >= 3) {
  xb <- xm[xm$bin %in% both_bins & xm$core_id %in% registered_cores, ]; cb <- cm[!is.na(cm$bin) & cm$bin %in% both_bins, ]
  r_g <- (rowSums(pbC[, registered_cores, drop = FALSE]) + 1) / (rowSums(pbX[, registered_cores, drop = FALSE]) + 1)
  xb$expected <- as.numeric(Matrix::crossprod(X[, xb$cell, drop = FALSE], r_g))          # per Xenium cell: sum_g r_g x_g
  agg <- function(d, cols) { a <- rowsum(as.matrix(d[, cols]), d$bin); data.frame(bin = rownames(a), a, n = as.numeric(table(d$bin)[rownames(a)])) }
  bx <- agg(xb, c("expected", "n_counts", "cell_area")); bc <- agg(cb, c("n_shared", "n_counts", "cell_area"))
  bins <- merge(bx, bc, by = "bin", suffixes = c("_x", "_c"))
  key <- do.call(rbind, strsplit(bins$bin, " ", fixed = TRUE))
  bins$core_id <- key[, 1]; bins$ix <- as.integer(key[, 2]); bins$iy <- as.integer(key[, 3])
  ## the 8 neighbouring bins of the same core
  nb <- matrix(0, nrow(bins), 5, dimnames = list(NULL, c("cnt_c", "area_c", "cnt_x", "area_x", "k")))
  for (dx in -1:1) for (dy in -1:1) {
    if (dx == 0 && dy == 0) next
    j <- match(paste(bins$core_id, bins$ix + dx, bins$iy + dy), bins$bin)
    ok <- !is.na(j)
    nb[ok, ] <- nb[ok, ] + cbind(bins$n_counts_c[j[ok]], bins$cell_area_c[j[ok]], bins$n_counts_x[j[ok]], bins$cell_area_x[j[ok]], 1)
  }
  bins$load_cosmx <- nb[, "cnt_c"] / nb[, "area_c"]; bins$load_xenium <- nb[, "cnt_x"] / nb[, "area_x"]
  bins$response <- log2((bins$n_shared / bins$n_c + 1) / (bins$expected / bins$n_x + 1))
  bin_df <- bins[nb[, "k"] >= 3 & bins$n_x >= 3 & bins$n_c >= 3 & is.finite(bins$load_cosmx) & bins$load_cosmx > 0 &
                   is.finite(bins$load_xenium) & bins$load_xenium > 0, ]
  bin_df$response <- bin_df$response - stats::ave(bin_df$response, bin_df$core_id, FUN = stats::median)
  crowd_bin <- bin_df |> group_by(core_id) |> filter(n() >= 15) |>
    summarise(slope_vs_cosmx_load = stats::coef(stats::lm(response ~ log2(load_cosmx)))[2],
              slope_vs_xenium_load = stats::coef(stats::lm(response ~ log2(load_xenium)))[2],
              n_bins = n(), load_p10 = stats::quantile(load_cosmx, 0.1), load_p90 = stats::quantile(load_cosmx, 0.9), .groups = "drop") |>
    left_join(core_info, by = "core_id")
  if (nrow(crowd_bin) > 0) {
    save_table(as.data.frame(crowd_bin), "10_crowding_bins_per_core")
    wt <- if (nrow(crowd_bin) >= 6) stats::wilcox.test(crowd_bin$slope_vs_cosmx_load)$p.value else NA_real_
    bin_df$load_bin <- cut(rank(bin_df$load_cosmx) / nrow(bin_df), seq(0, 1, 0.1), include.lowest = TRUE)
    trend <- bin_df |> group_by(load_bin) |> summarise(load = stats::median(load_cosmx), med = stats::median(response),
                                                       lo = stats::quantile(response, 0.25), hi = stats::quantile(response, 0.75), .groups = "drop")
    p_4d <- (ggplot(bin_df, aes(load_cosmx, response)) + geom_point(alpha = 0.15, size = 0.6, stroke = 0) +
               geom_hline(yintercept = 0, linetype = 2) +
               geom_linerange(data = trend, aes(load, ymin = lo, ymax = hi), inherit.aes = FALSE, colour = "#e34948") +
               geom_line(data = trend, aes(load, med), inherit.aes = FALSE, colour = "#e34948") +
               geom_point(data = trend, aes(load, med), inherit.aes = FALSE, colour = "#e34948", size = 1.8) +
               scale_x_log10() + theme_bw() +
               labs(x = "local load: CosMx transcripts per um2 of cell, neighbouring bins",
                    y = "log2 (CosMx observed / expected from Xenium), centred per core",
                    title = sprintf("%d-um bins of %d registered cores (red: median and IQR per decile of load)", CMP_BIN_UM, nrow(crowd_bin)))) |
      (ggplot(crowd_bin, aes(slope_vs_cosmx_load, slope_vs_xenium_load, colour = organ)) +
         geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) + geom_point(size = 2.2) + theme_bw() +
         labs(x = "slope per core, load from CosMx", y = "slope per core, load from Xenium (cross-check)",
              title = sprintf("median slope %.3f per doubling of load%s", stats::median(crowd_bin$slope_vs_cosmx_load),
                              if (is.na(wt)) "" else sprintf(" (Wilcoxon across cores, p = %.2g)", wt))))
    save_fig(p_4d + plot_annotation(title = "4d. Does CosMx lose transcripts, relative to Xenium, where the tissue is dense?",
                                    subtitle = "negative slope = yes (consistent with optical crowding on the whole-transcriptome panel); positive = Xenium loses; ~0 = neither"),
             "10_4d_bin_level_crowding", width = 15, height = 6.5)
  }
} else message("4d skipped: fewer than 3 registered cores (CMP_REGISTER = FALSE, or see tables/10_registration_qc.csv).")

## =============================================================================
## 5. Summary
## =============================================================================
med <- function(v) if (is.null(v) || all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE)
summary_tab <- data.frame(
  section = c("0", "0", "0", "1a", "1a", "1a", "1b", "1c", "1d", "2a", "2a", "2a", "2a", "2b", "2b", "3a", "3b", "3b", "4b", "4c", "4d"),
  measure = c("shared cores compared", "cores registered", "shared genes used (1:1)",
              "median cells per core: CosMx / Xenium", "median shared-gene transcripts per cell: CosMx / Xenium",
              "median all transcripts per cell: CosMx (whole transcriptome) / Xenium (panel)",
              "median per-core Spearman of gene profiles", "median gene: CosMx / Xenium transcripts per cell",
              "median per-gene Spearman across cores",
              "median % shared genes above background: Xenium", "median % shared genes above background: CosMx",
              "median genes above background per core: CosMx whole transcriptome", "median negative-probe counts per probe per cell: CosMx / Xenium",
              "median MECR: Xenium", "median MECR: CosMx",
              "median Jensen-Shannon divergence of lineage composition", "lineage agreement, CosMx cells with Xenium labels (%)",
              "lineage agreement, Xenium cells with CosMx labels (%)",
              "log-log slope CosMx vs Xenium gene abundance (1 = proportional)",
              "log2 relative CosMx/Xenium ratio at the 99th percentile cell (0 = no crowding; < 0 = CosMx)",
              "bin-level slope of log2 CosMx/Xenium per doubling of transcript load (< 0 = CosMx)"),
  value = c(length(shared_cores), length(registered_cores), nrow(g),
            med(core_wide$n_cells_CosMx / core_wide$n_cells_Xenium), med(core_wide$shared_counts_per_cell_CosMx / core_wide$shared_counts_per_cell_Xenium),
            med(core_wide$all_counts_per_cell_CosMx / core_wide$all_counts_per_cell_Xenium),
            med(prof_cor$spearman), 2^med(gene_tab$log2_ratio_cosmx_vs_xenium),
            if ("across_core_spearman" %in% colnames(gene_tab)) med(gene_tab$across_core_spearman) else NA,
            med(det_core$pct_shared_detected_xenium), med(det_core$pct_shared_detected_cosmx),
            med(det_core$n_genes_detected_cosmx_whole_transcriptome), med(det_core$neg_per_probe_per_cell_cosmx / det_core$neg_per_probe_per_cell_xenium),
            med(mecr_tab$mecr_xenium), med(mecr_tab$mecr_cosmx),
            med(js$jensen_shannon), 100 * mean(cm$lineage == cm$xenium_lineage, na.rm = TRUE), 100 * mean(xm$lineage == xm$cosmx_lineage, na.rm = TRUE),
            slope, med(crowd_cell$v), med(crowd_bin$slope_vs_cosmx_load)))
summary_tab$value <- signif(summary_tab$value, 4)
cat(sprintf("  %-3s %-92s %s\n", summary_tab$section, summary_tab$measure, format(summary_tab$value, scientific = FALSE, drop0trailing = TRUE)), sep = "")
save_table(summary_tab, "10_summary")
message("Unit of every 'median' above: the core (or the gene, where it says so). Figures: figures/10_*.png")
