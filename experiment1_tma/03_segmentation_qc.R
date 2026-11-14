## =============================================================================
## 03_segmentation_qc.R
##
## How good is the cell segmentation the CosMx data came with -- and, if you
## have a second segmentation, which one is better?
##
## CosMx segments on the instrument / in AtoMx (Cellpose-based, on the membrane
## + nuclear stains), one field of view at a time. There is no ground truth, so
## no single number decides. What matters for THIS project is that immune genes
## end up in immune cells and not smeared into neighbouring epithelium/stroma
## -- that is exactly what would fake "immune target expressed in organ X" in
## Aim 1, and it differs between platforms (10 compares it with Xenium). Per core:
##   * MECR (mutually exclusive co-expression rate: e.g. PTPRC+EPCAM+ cells) --
##     marker pairs of DIFFERENT lineages that should never share a cell. Lower
##     = less spill-over between neighbouring cells.
##   * transcripts inside cells (%), counts and genes per cell, cell area,
##     negative-probe rate.
##   * from the transcripts file (optional, large): where the transcripts are --
##     nucleus / cytoplasm / membrane / outside every cell -- per FOV and per gene.
##
## A SECOND SEGMENTATION is compared on the same cores when it exists in
##   data/segmentation/<name>/<slide_id>/
## either as another flat-file export (a re-segmentation done in AtoMx with other
## parameters: export it and put the exprMat + metadata files there) or in the
## common format (cell_feature_matrix/ + cells.csv.gz, e.g. converted from
## proseg / Baysor; coordinates in the microns of tables/01_cosmx_frame.csv).
## List the names in METHODS. With one method the script still writes the per-core
## segmentation QC; the method-vs-method scorecard needs two.
##
## Not in this repository (Xenium project: 03a cellpose, 03b / 03d segger):
## re-segmenting from the images needs the RAW export (CellStatsDir/Morphology2D)
## and is not implemented here.
##
## INPUT : flat files (+ data/segmentation/<name>/<slide>/), tables/02_cores_all.csv (02)
## OUTPUT: tables/03_segmentation_metrics_by_core.csv, 03_segmentation_summary.csv,
##         03_scorecard.csv (>= 2 methods), 03_tx_compartment_by_fov.csv, 03_tx_compartment_by_gene.csv
##         figures/03_segmentation_qc_by_core.png, 03_cell_numbers.png, 03_scorecard.png, 03_tx_compartments.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")

## methods without output are skipped. "cosmx" = the flat files of slides.csv.
METHODS <- c("cosmx", "atomx_reseg", "proseg")
METHOD_COLOURS <- c(cosmx = "#00a6c7", atomx_reseg = "#e0a800", proseg = "#8e6bbf")
## ADAPT: benchmark cores; NULL = all cores present in every method's output.
BENCH_CORES <- NULL
## Transcripts file: rows read at most (Inf = the whole file; a WTX slide has
## hundreds of millions of rows -- the file is streamed, memory stays small,
## but it takes a while). FALSE skips section 5.
USE_TX_FILE <- TRUE
TX_MAX_ROWS <- 5e7
TX_CHUNK_ROWS <- 2e6

slides <- read_slides()
tma_map <- read_tma_map()
cores_all <- utils::read.csv(file.path(TABLES_DIR, "02_cores_all.csv"), stringsAsFactors = FALSE)
cores_all <- cores_all[!is.na(cores_all$core_id), ]
markers <- read_marker_panel(MARKERS_RNA_CSV)

## -----------------------------------------------------------------------
## 1. Load each method, keep benchmark cores only
## -----------------------------------------------------------------------
outside_fov <- list()                          # transcripts outside cells per FOV (cell_ID 0 rows)
load_method <- function(method) {
  objs <- list()
  for (i in seq_len(nrow(slides))) {
    sid <- slides$slide_id[i]
    d <- seg_dir_for(slides[i, ], method)
    if (is.na(d) || !file.exists(d)) { message("  no ", method, " output for ", sid); next }
    l <- load_slide(slides[i, ], method)
    o <- l$obj
    if (!is.null(l$unassigned)) outside_fov[[paste(method, sid)]] <<- cbind(method = method, slide_id = sid, l$unassigned)
    cores_s <- cores_all[cores_all$slide_id == sid, ]
    ## same core definition for every method: nearest core centre within its radius
    o$core_id <- assign_cells_to_cores(o$x_um, o$y_um, cores_s, factor = 1.05)
    keep <- !is.na(o$core_id)
    if (!is.null(BENCH_CORES)) keep <- keep & o$core_id %in% BENCH_CORES
    if (!any(keep)) {
      message("  ", method, " / ", sid, ": no cells inside the cores of 02_cores_all.csv -- skipped")
      next
    }
    objs[[sid]] <- subset(o, cells = colnames(o)[keep])
  }
  if (length(objs) == 0) return(NULL)
  merge_slides(objs)
}

seg <- lapply(stats::setNames(METHODS, METHODS), load_method)
seg <- seg[!vapply(seg, is.null, logical(1))]
message("Methods with output: ", paste(names(seg), collapse = ", "))
stopifnot("No segmentation could be loaded" = length(seg) >= 1)

## MECR uses the use_for_mecr genes of marker_panel_rna.csv: pairs of genes from
## DIFFERENT lineages should never be in the same cell
mecr_genes <- lapply(markers$mecr, function(g) intersect(match_case(g, rownames(seg[[1]])), rownames(seg[[1]])))
mecr_genes <- mecr_genes[lengths(mecr_genes) > 0]
message("MECR marker genes on the panel: ",
        paste(sprintf("%s = %s", names(mecr_genes), vapply(mecr_genes, paste, "", collapse = "/")), collapse = "; "))
if (length(mecr_genes) < 2) {
  message("  fewer than 2 lineages with MECR markers on the panel -> MECR cannot be computed; ",
          "check use_for_mecr in ", basename(MARKERS_RNA_CSV))
}

common_cores <- Reduce(intersect, lapply(seg, function(o) unique(o$core_id)))
message(length(common_cores), " cores present in all methods.")
if (length(common_cores) == 0) stop("No core was segmented by all methods. Run them on the same cores, or set BENCH_CORES.")

## -----------------------------------------------------------------------
## 2. Metrics per core x method
## -----------------------------------------------------------------------
metrics <- dplyr::bind_rows(lapply(names(seg), function(method) {
  o <- seg[[method]]
  counts <- GetAssayData(o, assay = "CosMx", layer = "counts")
  m <- o[[]]
  ctrl_cols <- grep("^nCount_(ControlProbe|ControlCodeword|BlankCodeword)$", colnames(m), value = TRUE)
  m$ctrl <- if (length(ctrl_cols)) rowSums(m[, ctrl_cols, drop = FALSE]) else 0
  dplyr::bind_rows(lapply(common_cores, function(cid) {
    cells <- rownames(m)[m$core_id == cid]
    mm <- m[cells, ]
    data.frame(
      method = method,
      core_id = cid,
      n_cells = length(cells),
      median_counts = median(mm$nCount_CosMx),
      median_genes = median(mm$nFeature_CosMx),
      median_area = median(mm$cell_area, na.rm = TRUE),
      counts_per_um2 = sum(mm$nCount_CosMx) / sum(mm$cell_area, na.rm = TRUE),
      pct_cells_min_counts = 100 * mean(mm$nCount_CosMx >= MIN_COUNTS),
      pct_fov_border_fragments = if ("split_ratio" %in% colnames(mm) && !is.na(MIN_SPLIT_RATIO))
        100 * mean(!is.na(mm$split_ratio) & mm$split_ratio > 0 & mm$split_ratio < MIN_SPLIT_RATIO) else NA_real_,
      ctrl_per_1000 = 1000 * sum(mm$ctrl) / max(sum(mm$nCount_CosMx), 1),
      mecr = compute_mecr(counts[, cells, drop = FALSE], markers$mecr)
    )
  }))
}))

## Transcripts inside cells: (counts in cells of the core) / (those + the
## cell_ID-0 counts of the core's FOVs). Needs the flat-file exprMat of the method.
if (length(outside_fov) > 0) {
  outside_fov <- dplyr::bind_rows(outside_fov)
  inside <- dplyr::bind_rows(lapply(names(seg), function(method) {
    m <- seg[[method]][[]]
    if (!"fov" %in% colnames(m)) return(NULL)
    m |> group_by(slide_id, core_id, fov) |> summarise(n_in = sum(nCount_CosMx), .groups = "drop") |> mutate(method = method)
  }))
  tx_in <- inside |> left_join(outside_fov, by = c("method", "slide_id", "fov")) |>
    filter(!is.na(n_outside_genes)) |>
    group_by(method, core_id) |>
    summarise(pct_tx_assigned = 100 * sum(n_in) / max(sum(n_in) + sum(n_outside_genes), 1), .groups = "drop")
  metrics <- metrics |> left_join(tx_in, by = c("method", "core_id"))
}
## ... or a summary written by the tool that made the segmentation (region_id, n_tx_genes, n_assigned_<method>)
for (mth in setdiff(names(seg), "cosmx")) {
  fs <- file.path(SEG_DIR, mth, slides$slide_id, "assignment_summary.csv")
  fs <- fs[file.exists(fs)]
  if (length(fs) == 0) next
  s <- dplyr::bind_rows(lapply(fs, utils::read.csv))
  n_col <- paste0("n_assigned_", mth)
  if (!all(c("region_id", "n_tx_genes", n_col) %in% colnames(s))) next
  k <- match(paste(mth, s$region_id), paste(metrics$method, metrics$core_id))
  if (!"pct_tx_assigned" %in% colnames(metrics)) metrics$pct_tx_assigned <- NA_real_
  metrics$pct_tx_assigned[k[!is.na(k)]] <- (100 * s[[n_col]] / pmax(s$n_tx_genes, 1))[!is.na(k)]
}

metrics <- metrics |> left_join(tma_map[, c("core_id", "organ", "tissue_type")], by = "core_id")
metrics$slide_id <- cores_all$slide_id[match(metrics$core_id, cores_all$core_id)]
metrics$method <- factor(metrics$method, levels = intersect(METHODS, unique(metrics$method)))
save_table(metrics, "03_segmentation_metrics_by_core")

## -----------------------------------------------------------------------
## 3. Per-core segmentation QC of every method, by organ. Segmentation
##    quality depends on the tissue: dense lymphoid tissue, large hepatocytes
##    and thin-walled lung behave very differently. A core far from the other
##    cores of its organ deserves a look in AtoMx.
## -----------------------------------------------------------------------
show <- intersect(c("pct_tx_assigned", "median_counts", "median_area", "counts_per_um2", "mecr", "ctrl_per_1000",
                    "pct_fov_border_fragments"), colnames(metrics))
show <- show[vapply(show, function(k) any(is.finite(metrics[[k]])), logical(1))]
long <- tidyr::pivot_longer(metrics, all_of(show), names_to = "measure", values_to = "value")
long$measure <- factor(long$measure, levels = show)
p_qc <- ggplot(long, aes(organ, value, colour = method)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), size = 1.8) +
  scale_colour_manual(values = METHOD_COLOURS, na.value = "grey40") +
  facet_wrap(~ measure, scales = "free_y", ncol = 4) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "Segmentation QC per core (one dot = one core)",
       subtitle = "mecr = marker pairs of different lineages in the same cell (lower = cleaner); ctrl_per_1000 = control counts per 1000 gene counts")
save_fig(p_qc, "03_segmentation_qc_by_core", width = 16, height = 4 * ceiling(length(show) / 4) + 1)

metric_cols <- intersect(c("n_cells", "median_counts", "median_genes", "pct_tx_assigned", "pct_cells_min_counts", "mecr",
                           "median_area", "counts_per_um2", "ctrl_per_1000", "pct_fov_border_fragments"), colnames(metrics))
decision <- metrics |>
  group_by(method, organ) |>
  summarise(n_cores = n(), across(all_of(metric_cols), \(v) round(median(v, na.rm = TRUE), 3)), .groups = "drop") |>
  arrange(organ, method)
print(as.data.frame(decision))
save_table(decision, "03_segmentation_summary")

## -----------------------------------------------------------------------
## 4. Two or more methods: cell numbers and a scorecard, paired per core.
##    Higher is better for transcripts assigned to cells, counts per cell and
##    % cells passing the count threshold; lower is better for MECR and the
##    control rate. More cells is NOT better by itself: split cells and debris
##    also add cells.
## -----------------------------------------------------------------------
others <- setdiff(levels(metrics$method), "cosmx")
if (length(seg) >= 2 && "cosmx" %in% names(seg)) {
  tot <- metrics |> group_by(slide_id, method) |> summarise(n_cells = sum(n_cells), .groups = "drop")
  p_tot <- ggplot(tot, aes(slide_id, n_cells, fill = method)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(aes(label = format(n_cells, big.mark = ",")), position = position_dodge(0.8), vjust = -0.3, size = 3) +
    scale_fill_manual(values = METHOD_COLOURS) +
    theme_bw() + labs(x = NULL, y = "cells in the compared cores", title = "Cells called per method")
  wide <- tidyr::pivot_wider(metrics[, c("core_id", "slide_id", "method", "n_cells")],
                             names_from = method, values_from = n_cells)
  p_cores <- lapply(others, function(mth) {
    lims <- range(c(0, wide$cosmx, wide[[mth]]), na.rm = TRUE)      # same scale on both axes: diagonal = equal
    ggplot(wide, aes(cosmx, .data[[mth]], colour = slide_id)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(size = 2) + coord_equal(xlim = lims, ylim = lims) + theme_bw() +
      labs(x = "cells per core: CosMx (AtoMx)", y = paste("cells per core:", mth),
           title = sprintf("%s / CosMx cells: median ratio %.2f", mth, stats::median(wide[[mth]] / wide$cosmx, na.rm = TRUE)))
  })
  save_fig(patchwork::wrap_plots(c(list(p_tot), p_cores), nrow = 1), "03_cell_numbers",
           width = 6 * (1 + length(others)), height = 5.5)

  better <- c(pct_tx_assigned = "higher", median_counts = "higher", pct_cells_min_counts = "higher",
              mecr = "lower", ctrl_per_1000 = "lower")
  better <- better[names(better) %in% colnames(metrics)]
  score <- dplyr::bind_rows(lapply(names(better), function(k) {
    w <- tidyr::pivot_wider(metrics[, c("core_id", "method", k)], names_from = method, values_from = all_of(k))
    dplyr::bind_rows(lapply(others, function(mth) {
      d <- w[[mth]] - w$cosmx
      ok <- !is.na(d)
      wins <- if (better[[k]] == "higher") sum(d[ok] > 0) else sum(d[ok] < 0)
      losses <- sum(ok) - wins - sum(d[ok] == 0)
      data.frame(measure = k, better_if = better[[k]], method = mth, n_cores = sum(ok),
                 cosmx_median = round(stats::median(w$cosmx, na.rm = TRUE), 3),
                 method_median = round(stats::median(w[[mth]], na.rm = TRUE), 3),
                 method_wins = wins, cosmx_wins = losses,
                 winner = ifelse(wins > sum(ok) / 2, mth, ifelse(losses > sum(ok) / 2, "cosmx", "tie")))
    }))
  }))
  print(score)
  save_table(score, "03_scorecard")

  long <- tidyr::pivot_longer(metrics, all_of(c("n_cells", names(better))), names_to = "measure", values_to = "value")
  lab <- score |> group_by(measure) |>
    summarise(lab = paste0(measure[1], " (", better_if[1], " = better): ",
                           paste(sprintf("cosmx better in %d, %s in %d of %d cores", cosmx_wins, method,
                                         method_wins, n_cores), collapse = "; ")), .groups = "drop")
  lab <- rbind(lab, data.frame(measure = "n_cells", lab = "n_cells (descriptive: more is not better by itself)"))
  long$panel <- lab$lab[match(long$measure, lab$measure)]
  p_score <- ggplot(long, aes(method, value, group = core_id)) +
    geom_line(colour = "grey70") + geom_point(aes(colour = method), size = 2) +
    scale_colour_manual(values = METHOD_COLOURS, guide = "none") +
    facet_wrap(~ panel, scales = "free_y", ncol = 2, labeller = label_wrap_gen(55)) +
    theme_bw() + labs(x = NULL, y = NULL, title = "Segmentation scorecard -- each line is one core")
  save_fig(p_score, "03_scorecard", width = 12, height = 4 * ceiling((length(better) + 1) / 2))
  ## >>> DECISION <<<
  ## Pick ONE segmentation for all slides (mixing methods between organs makes
  ## Aim 1/2 comparisons uninterpretable). Set SEGMENTATION in config/config.R,
  ## re-run 01 and 02, and state the choice + these numbers in your methods.
} else {
  message("One segmentation only: the method-vs-method scorecard is skipped (see the header for how to add one).")
}

## -----------------------------------------------------------------------
## 5. Transcripts file (optional): where are the transcripts?
##    CellComp of every transcript: Nuclear / Cytoplasm / Membrane, or None
##    (legacy: 0) = outside every cell. Per FOV: a low share inside cells means
##    under-segmentation or much extracellular / damaged tissue. Per gene: a
##    gene found mostly OUTSIDE cells behaves like background; negative probes
##    are the reference for that.
##    The file is streamed in chunks of TX_CHUNK_ROWS; only counts are kept.
## -----------------------------------------------------------------------
if (isTRUE(USE_TX_FILE)) {
  by_fov <- list(); by_gene <- list()
  for (i in seq_len(nrow(slides))) {
    sid <- slides$slide_id[i]
    files <- tryCatch(find_cosmx_files(resolve_path(slides$cosmx_dir[i]), slides$cosmx_prefix[i]), error = function(e) NULL)
    if (is.null(files) || is.na(files$tx)) { message(sid, ": no transcripts file -- section 5 skipped for this slide"); next }
    message(sid, ": streaming ", basename(files$tx), " (up to ", format(TX_MAX_ROWS, big.mark = ","), " rows)")
    con <- open_text(files$tx)
    header <- gsub('^"|"$', "", strsplit(readLines(con, n = 1), ",", fixed = TRUE)[[1]])
    need <- c("fov", "target", "CellComp")
    if (!all(need %in% header)) { close(con); message("  columns ", paste(setdiff(need, header), collapse = ", "), " missing -- skipped"); next }
    n_read <- 0; acc_f <- NULL; acc_g <- NULL
    repeat {
      lines <- readLines(con, n = min(TX_CHUNK_ROWS, TX_MAX_ROWS - n_read))
      if (length(lines) == 0) break
      idx <- sort(match(need, header))
      d <- data.table::fread(text = lines, header = FALSE, sep = ",", select = idx, showProgress = FALSE)
      data.table::setnames(d, header[idx])
      n_read <- n_read + nrow(d)
      d[, comp := ifelse(CellComp %in% c("None", "0", ""), "outside", as.character(CellComp))]
      acc_f <- data.table::rbindlist(list(acc_f, d[, .N, by = .(fov, comp)]))[, .(N = sum(N)), by = .(fov, comp)]
      acc_g <- data.table::rbindlist(list(acc_g, d[, .N, by = .(target, comp)]))[, .(N = sum(N)), by = .(target, comp)]
      if (n_read >= TX_MAX_ROWS) { message("  stopped at TX_MAX_ROWS (the file has more rows: per-FOV numbers cover the first FOVs only)"); break }
    }
    close(con)
    by_fov[[sid]] <- cbind(slide_id = sid, as.data.frame(acc_f))
    by_gene[[sid]] <- cbind(slide_id = sid, as.data.frame(acc_g))
  }
  if (length(by_fov) > 0) {
    tx_fov <- dplyr::bind_rows(by_fov) |> group_by(slide_id, fov) |> mutate(pct = 100 * N / sum(N)) |> ungroup()
    save_table(tx_fov, "03_tx_compartment_by_fov")
    tx_gene <- dplyr::bind_rows(by_gene) |> group_by(target, comp) |> summarise(N = sum(N), .groups = "drop") |>
      group_by(target) |> mutate(n_total = sum(N), pct = 100 * N / n_total) |> ungroup()
    tx_gene$feature_type <- cosmx_feature_type(tx_gene$target)
    save_table(tx_gene, "03_tx_compartment_by_gene")
    p_tf <- ggplot(tx_fov, aes(factor(fov), pct, fill = comp)) + geom_col(width = 0.9) +
      facet_grid(~ slide_id, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = c(Nuclear = "#2a78d6", Cytoplasm = "#1baf7a", Membrane = "#eda100", outside = "grey70"),
                        na.value = "grey40", name = NULL) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
      labs(x = "FOV", y = "% of transcripts", title = "Where the transcripts are, per field of view")
    out_gene <- tx_gene |> filter(comp == "outside", n_total >= 50)
    p_tg <- ggplot(out_gene, aes(pct, fill = feature_type)) +
      geom_histogram(bins = 50, position = "identity", alpha = 0.6) +
      scale_fill_manual(values = c(`Gene Expression` = "#2a78d6", `Negative Control Probe` = "#e34948",
                                   `Negative Control Codeword` = "#eda100"), name = NULL) +
      scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100, 1000, 10000)) + theme_bw() +
      labs(x = "% of the target's transcripts outside every cell", y = "targets",
           title = "Genes vs negative probes: share outside cells",
           subtitle = "a gene sitting with the negative probes is mostly background")
    save_fig(p_tf / p_tg, "03_tx_compartments", width = 16, height = 10)
  }
}
