## =============================================================================
## 04_annotation_wtx.R
##
## Cell type calling from the whole transcriptome: cluster all cells, find the
## marker genes of each cluster, look at them next to canonical markers and at
## where each cluster sits in the tissue, then YOU name the clusters.
##
##   1. Coverage of the marker signatures (tables/04_panel_coverage.csv)
##   2. Normalise, pick variable genes, cluster -> ST_clusters (Harmony on slide only)
##   3. FindAllMarkers per cluster + dot plot of canonical markers (config.R) and
##      the top markers of each cluster   -> figures/04_top_markers_dotplot_transcriptomics.png
##   4. Where is each cluster?             -> figures/04_ST_clusters_spatial_<slide>.png
##   5. Template for your cluster names     -> annotation/cluster_to_celltype_template_<seg>.R
##      (top markers + a lineage suggestion from the marker signatures, per cluster)
##      + the CosMx stains per cluster (PanCK, CD45, ...) as a protein-level check
##   6. cell_type / lineage from CLUSTER_TO_CELLTYPE / CLUSTER_TO_LINEAGE (config.R);
##      before you have written them, cells are labelled cluster_<n>.
##   7. Cell types on the UMAP, in the tissue and per core.
##
## Workflow: run 04 -> look at the dot plot + spatial plots -> copy the template
## into config.R and fill in the names -> run 04 again.
##
## Differences from the Xenium (probe panel) version of this step:
##   * ~19,000 genes instead of a few hundred: clustering uses the N_HVG most
##     variable genes plus all marker genes, not every gene.
##   * nearly every marker of config/marker_panel_rna.csv is on the panel, so the
##     lineage suggestion is usually well supported, and clusters can be split
##     further than the probe panel allowed. For the comparison in 10 keep the
##     Xenium names where the cell type is the same.
##
## Gene names in config/ are matched to the panel ignoring case, so human-style
## (EPCAM) and mouse-style (Epcam) names both work.
##
## INPUT : cosmx_<seg>_cores.rds (02)
## OUTPUT: cosmx_<seg>_annotated.rds, annotation/cluster_to_celltype_template_<seg>.R,
##         tables/04_*.csv, data/registration/cosmx_cells_<slide>.csv.gz (for 07)
## =============================================================================

source("config/config.R")
source("R/helpers.R")

cmx <- readRDS(obj_path("cores"))
panel_genes <- rownames(cmx[["CosMx"]])
markers <- read_marker_panel(MARKERS_RNA_CSV)

## -----------------------------------------------------------------------
## 1. Coverage of the marker signatures (used for the lineage suggestion in
##    the template). On the whole transcriptome a missing marker is usually a
##    naming issue (alias, or a multi-gene probe such as "TPSAB1/B2"): the
##    table lists them.
## -----------------------------------------------------------------------
l1 <- restrict_to_panel(markers$level1, panel_genes, "Level 1")
save_table(l1$coverage, "04_panel_coverage")
message("Level-1 signature markers on the panel: ", sum(l1$coverage$n_present), " of ", sum(l1$coverage$n_markers))

## -----------------------------------------------------------------------
## 2. Normalise + cluster all cells
##    LogNormalize with scale.factor = median counts per cell (not 1e4: keeps
##    the numbers on the scale of the data). Variable genes: the N_HVG most
##    variable ones (vst), plus every gene you annotate with (marker panel,
##    CANONICAL_MARKERS, probe_expectations.csv), so no marker is left out of
##    the PCA because it is not "variable enough" across a multi-organ TMA.
##    ScaleData / PCA run on these genes only (a scaled matrix of all 19,000
##    genes x all cells would not fit in memory).
## -----------------------------------------------------------------------
DefaultAssay(cmx) <- "CosMx"
cmx <- NormalizeData(cmx, scale.factor = median(cmx$nCount_CosMx), verbose = FALSE)
cmx <- FindVariableFeatures(cmx, selection.method = "vst", nfeatures = min(N_HVG, length(panel_genes)), verbose = FALSE)
probe_genes <- if (file.exists(PROBE_EXPECT_CSV)) read_config_csv(PROBE_EXPECT_CSV, stringsAsFactors = FALSE)$gene else character(0)
annot_genes <- unique(c(unlist(markers$level1), unlist(markers$level2), CANONICAL_MARKERS, probe_genes))
annot_genes <- intersect(match_case(as_seurat_features(annot_genes), panel_genes), panel_genes)
VariableFeatures(cmx) <- union(VariableFeatures(cmx), annot_genes)
message(length(VariableFeatures(cmx)), " genes used for PCA (", N_HVG, " most variable + ",
        length(VariableFeatures(cmx)) - min(N_HVG, length(panel_genes)), " annotation genes that were not among them)")
cmx <- ScaleData(cmx, features = VariableFeatures(cmx), verbose = FALSE)
cmx <- RunPCA(cmx, features = VariableFeatures(cmx), npcs = N_PCS, verbose = FALSE)
red <- "pca"
## Harmony corrects SLIDE (technical batch). Never put organ/patient here:
## those are the biology Aim 1/2 asks about.
if (USE_HARMONY && length(unique(cmx$slide_id)) > 1 && requireNamespace("harmony", quietly = TRUE)) {
  cmx <- harmony::RunHarmony(cmx, group.by.vars = "slide_id", reduction.use = "pca",
                             reduction.save = "harmony", verbose = FALSE)
  red <- "harmony"
}
cmx <- FindNeighbors(cmx, reduction = red, dims = 1:N_PCS, verbose = FALSE)
cmx <- FindClusters(cmx, resolution = RES_CLUSTERS, cluster.name = "ST_clusters", verbose = FALSE)
cmx <- RunUMAP(cmx, reduction = red, dims = 1:N_PCS, reduction.name = "umap", verbose = FALSE)
clusters <- levels(cmx$ST_clusters)
message(length(clusters), " clusters at RES_CLUSTERS = ", RES_CLUSTERS)

p_umap_cl <- plot_umap_labelled(cmx, "ST_clusters", label_size = 6, legend = FALSE, title = "ST_clusters") |
  plot_umap_labelled(cmx, "slide_id", label = FALSE, title = "Slide",
                     pal = stats::setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")[seq_along(unique(cmx$slide_id))],
                                           sort(unique(cmx$slide_id))))
save_fig(p_umap_cl, "04_ST_clusters_umap", width = 22, height = 10, dpi = 200)

## -----------------------------------------------------------------------
## 3. Cell type calling, based on transcriptomics: marker genes per cluster,
##    shown next to canonical cell-type markers.
##    max.cells.per.ident subsamples big clusters (speed); install presto for
##    a much faster Wilcoxon test.
## -----------------------------------------------------------------------
Idents(cmx) <- "ST_clusters"
cl_markers <- FindAllMarkers(cmx, assay = "CosMx", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
                             max.cells.per.ident = MARKERS_MAX_CELLS, verbose = FALSE)
save_table(cl_markers, "04_cluster_markers_all")
top_markers <- cl_markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = TOP_N_PER_CLUSTER) |>
  ungroup()
save_table(top_markers, "04_cluster_markers_top")
print(top_markers |> dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2), n = 400)

canonical <- match_case(as_seurat_features(CANONICAL_MARKERS), panel_genes)
all_features <- unique(c(canonical, top_markers$gene))
present_features <- intersect(all_features, panel_genes)
missing_features <- setdiff(all_features, present_features)
if (length(missing_features) > 0) {
  message("Not on the CosMx panel under this name, dropped from plot: ", paste(missing_features, collapse = ", "))
}
if (length(present_features) == 0) stop("No canonical or marker genes on the panel -- check CANONICAL_MARKERS in config.R")

p_dotplot <- suppressMessages(marker_dotplot(cmx, features = present_features, group.by = "ST_clusters") +
  scale_colour_gradient2(low = "cornflowerblue", mid = "white", high = "red", midpoint = 0,
                         name = "Average\nexpression\n(scaled)")) +
  theme(axis.text.x = element_text(size = 7)) +
  ggtitle("Canonical markers + top markers per cluster (ST_clusters)")
save_fig(p_dotplot, "04_top_markers_dotplot_transcriptomics",
         width = max(16, 0.16 * length(present_features) + 4), height = max(6, 0.3 * length(clusters) + 3))

## -----------------------------------------------------------------------
## 4. Spatial location per cluster, one figure per slide
## -----------------------------------------------------------------------
for (sid in unique(cmx$slide_id)) {
  save_fig(plot_clusters_spatial(cmx[[]], "ST_clusters", slide = sid),
           paste0("04_ST_clusters_spatial_", sid), width = 20, height = 18)
}

## -----------------------------------------------------------------------
## 5. Template for the cluster names: top markers and a lineage suggestion
##    (level-1 signature of config/marker_panel_rna.csv with the highest
##    mean score; margin = best - second best, small = ambiguous).
## -----------------------------------------------------------------------
top_by_cluster <- top_markers |> group_by(cluster) |> summarise(genes = paste(gene, collapse = ", "), .groups = "drop")
top_by_cluster <- stats::setNames(top_by_cluster$genes, as.character(top_by_cluster$cluster))[clusters]
top_by_cluster[is.na(top_by_cluster)] <- "(no markers)"
suggested <- stats::setNames(rep("Unknown", length(clusters)), clusters)
margin <- stats::setNames(rep(NA_real_, length(clusters)), clusters)
usable <- l1$signatures[lengths(l1$signatures) >= 2]
if (length(usable) >= 2) {
  sc <- score_signatures(cmx, usable)[[paste0("score_", names(usable))]]
  agg <- stats::aggregate(sc, by = list(cluster = as.character(cmx$ST_clusters)), FUN = mean)
  m <- as.matrix(agg[, -1, drop = FALSE]); dimnames(m) <- list(agg$cluster, names(usable))
  m <- m[clusters, , drop = FALSE]
  suggested[] <- colnames(m)[apply(m, 1, which.max)]
  srt <- t(apply(m, 1, sort, decreasing = TRUE))
  margin[] <- round(srt[, 1] - srt[, 2], 3)
} else {
  message("Fewer than 2 level-1 signatures with >= 2 markers on the panel -- no lineage suggestion.")
}
summary_tab <- data.frame(cluster = clusters, n_cells = as.integer(table(cmx$ST_clusters)[clusters]),
                          top_markers = unname(top_by_cluster), suggested_lineage = unname(suggested),
                          lineage_margin = unname(margin))
save_table(summary_tab, "04_cluster_summary")

tmpl_path <- file.path(ANNOT_DIR, sprintf("cluster_to_celltype_template_%s.R", SEGMENTATION))
n <- length(clusters); comma <- ifelse(seq_len(n) < n, ",", " ")
writeLines(c(
  sprintf("## %d clusters (SEGMENTATION = %s, RES_CLUSTERS = %s). Copy into config/config.R,",
          n, SEGMENTATION, RES_CLUSTERS),
  "## fill in the cell type names (see figures/04_top_markers_dotplot_transcriptomics.png),",
  "## check the lineages, and re-run 04. Comments: top markers of the cluster | lineage suggestion margin.",
  "CLUSTER_TO_CELLTYPE <- c(",
  sprintf('  "%s" = ""%s   # %s', clusters, comma, top_by_cluster),
  ")",
  "CLUSTER_TO_LINEAGE <- c(",
  sprintf('  "%s" = "%s"%s   # margin %s', clusters, suggested, comma, margin),
  ")"), tmpl_path)
message("Template for your cluster names -> ", tmpl_path)

## CosMx stains per cluster: the morphology markers (PanCK, CD45, CD68, membrane,
## DAPI; names depend on the segmentation kit) are measured per cell and are
## independent of the RNA. An "Epithelial" cluster should be PanCK-high, an
## "Immune" one CD45-high -- if not, look at that cluster again.
if_cols <- grep("^IF_", colnames(cmx[[]]), value = TRUE)
if_cols <- if_cols[vapply(if_cols, function(k) any(!is.na(cmx[[k, drop = TRUE]])), logical(1))]
if (length(if_cols) > 0) {
  if_tab <- cmx[[]][, c("ST_clusters", if_cols)] |>
    group_by(ST_clusters) |>
    summarise(across(all_of(if_cols), \(v) mean(log1p(v), na.rm = TRUE)), .groups = "drop")
  save_table(as.data.frame(if_tab), "04_IF_stains_by_cluster")
  if_long <- tidyr::pivot_longer(if_tab, all_of(if_cols), names_to = "stain", values_to = "mean_log") |>
    group_by(stain) |> mutate(z = as.numeric(scale(mean_log))) |> ungroup() |>
    mutate(stain = sub("^IF_", "", stain), z = ifelse(is.na(z), 0, z),
           lineage_suggestion = suggested[as.character(ST_clusters)])
  p_if <- ggplot(if_long, aes(stain, ST_clusters, fill = z)) +
    geom_tile() + scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", name = "z-score\nacross\nclusters") +
    facet_grid(lineage_suggestion ~ ., scales = "free_y", space = "free_y") + theme_bw() +
    theme(strip.text.y = element_text(angle = 0)) +
    labs(x = "CosMx stain (mean intensity per cell)", y = "ST_clusters",
         title = "Protein-level check: stains per cluster, grouped by the RNA lineage suggestion")
  save_fig(p_if, "04_IF_stains_by_cluster", width = 8, height = max(5, 0.28 * length(clusters) + 2))
}

## -----------------------------------------------------------------------
## 6. cell_type and lineage from your mapping (config.R)
## -----------------------------------------------------------------------
apply_map <- function(map, what) {
  unlabelled <- setdiff(clusters, names(map)[nzchar(map)])
  if (length(unlabelled) > 0) {
    stop(what, " has no (or an empty) name for cluster(s): ", paste(unlabelled, collapse = ", "),
         "\n  The clusters changed (data / SEGMENTATION / RES_CLUSTERS?) or the template is not filled in yet.",
         "\n  See ", tmpl_path)
  }
  unname(map[as.character(cmx$ST_clusters)])
}
if (is.null(CLUSTER_TO_CELLTYPE)) {
  cmx$cell_type <- paste0("cluster_", cmx$ST_clusters)
  message("CLUSTER_TO_CELLTYPE not set yet -- cells labelled cluster_<n>. Fill in ", tmpl_path)
} else {
  cmx$cell_type <- apply_map(CLUSTER_TO_CELLTYPE, "CLUSTER_TO_CELLTYPE")
}
cmx$lineage <- if (is.null(CLUSTER_TO_LINEAGE)) unname(suggested[as.character(cmx$ST_clusters)]) else
  apply_map(CLUSTER_TO_LINEAGE, "CLUSTER_TO_LINEAGE")
bad_lin <- setdiff(unique(cmx$lineage), c(names(markers$level1), "Unknown", "Low_quality"))
if (length(bad_lin) > 0) {
  message("NOTE: lineage values not among the level-1 names of ", basename(MARKERS_RNA_CSV), ": ",
          paste(bad_lin, collapse = ", "), " -- 05 selects immune cells with lineage == \"Immune\".")
}
cmx$celltype <- factor(cmx$cell_type)     # same labels; `cell_type` is what 05-09 use
print(table(cmx$celltype))
print(table(cell_type = cmx$cell_type, lineage = cmx$lineage))

## -----------------------------------------------------------------------
## 7. Cell types on the UMAP, in the tissue and per core
## -----------------------------------------------------------------------
## colours: one hue family per lineage (Epithelial blue-purple, Immune orange,
## ...), the same in every figure below (celltype_palette() in R/helpers.R)
pal_ct <- celltype_palette(cmx$cell_type, cmx$lineage)
save_fig(plot_celltype_umap(cmx), "04_celltype_UMAP", width = 16, height = 12, dpi = 200)
n_slides <- length(unique(cmx$slide_id))
save_fig(plot_celltype_sections(cmx[[]], panel_in = 11), "04_celltype_sections",
         width = 11 * n_slides + 4, height = 12, dpi = 250)
for (sid in unique(cmx$slide_id)) {
  save_fig(plot_clusters_spatial(cmx[[]], "cell_type", slide = sid),
           paste0("04_celltype_spatial_facet_", sid), width = 20, height = 18)
}

## composition: one bar per core, cell types stacked in legend (lineage) order
meta_ct <- cmx[[]]
meta_ct$cell_type <- factor(meta_ct$cell_type, levels = names(pal_ct))
facet_col <- if ("organ" %in% colnames(meta_ct) && length(unique(meta_ct$organ)) > 1) "organ" else "slide_id"
comp <- meta_ct[!is.na(meta_ct$core_id), ] |>
  count(.data[[facet_col]], core_id, cell_type) |>
  group_by(core_id) |> mutate(frac = n / sum(n)) |> ungroup()
p_comp <- ggplot(comp, aes(core_id, frac, fill = cell_type)) +
  geom_col(width = 0.9) + facet_grid(stats::as.formula(paste("~", facet_col)), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_ct, name = "Cell type") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(), strip.text = element_text(face = "bold")) +
  labs(x = NULL, y = "fraction of cells", title = "Cell-type composition per core")
save_fig(p_comp, "04_composition_per_core", width = 20, height = 8, dpi = 200)

## Spatial gallery: the 3 cores with most cells of each slide, dots ~ one cell wide
gallery_cores <- meta_ct[!is.na(meta_ct$core_id), ] |> count(slide_id, core_id) |>
  group_by(slide_id) |> slice_max(n, n = 3, with_ties = FALSE) |> pull(core_id)
gal <- meta_ct[meta_ct$core_id %in% gallery_cores, ]
core_ext <- gal |> group_by(core_id) |>
  summarise(e = max(diff(range(x_um)), diff(range(y_um))), .groups = "drop") |> pull(e) |> max()
gal_cols <- 3
gal_rows <- ceiling(length(gallery_cores) / gal_cols)
p_gal <- plot_cores_spatial(gal, "cell_type", size = cell_point_size(core_ext, panel_in = 4.6)) +
  scale_colour_manual(values = pal_ct, name = "Cell type", drop = TRUE) +
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  theme(strip.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12),
        plot.title = element_text(size = 15)) +
  ggtitle("Cell types in the largest cores")
save_fig(p_gal, "04_spatial_gallery_celltype", width = 5 * gal_cols + 4, height = 5 * gal_rows + 0.8, dpi = 200)

## -----------------------------------------------------------------------
## 8. Export of the cell table (centroids + labels), used by
##    07_register_he.py (check figure) and handy for looking at cells outside R.
## -----------------------------------------------------------------------
for (sid in unique(cmx$slide_id)) {
  cells <- colnames(cmx)[cmx$slide_id == sid]
  out <- data.frame(
    cell = cells,
    slide_id = sid,
    core_id = cmx$core_id[cells],
    fov = cmx$fov[cells],
    x_um = cmx$x_um[cells],
    y_um = cmx$y_um[cells],
    lineage = cmx$lineage[cells],
    cell_type = cmx$cell_type[cells]
  )
  f <- file.path(REG_DIR, paste0("cosmx_cells_", sid, ".csv.gz"))
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(out, f) else utils::write.csv(out, gzfile(f), row.names = FALSE)
  message("Exported ", nrow(out), " cells -> ", f)
}

saveRDS(cmx, obj_path("annotated"))
message("Saved -> ", obj_path("annotated"))
