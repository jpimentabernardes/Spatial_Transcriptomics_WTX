## =============================================================================
## R/helpers.R
##
## Shared functions for all experiment1_tma R scripts. Sourced right after
## config/config.R. Nothing in here reads data on its own -- only functions.
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## -----------------------------------------------------------------------
## Settings added to config/config.R after its first version. A config.R
## edited before they existed (your own paths) still works: any missing
## setting gets its default here, with a message -- copy it into config.R
## to change it.
## -----------------------------------------------------------------------
local({
  defaults <- list(
    ASSAY = "CosMx", COSMX_PIXEL_SIZE_UM = 0.12028, COSMX_FLIP_Y = TRUE, COSMX_CHUNK_ROWS = 2000,
    COSMX_NEG_REGEX = "^(Negative|NegPrb|NegProbe)", COSMX_SYS_REGEX = "^(SystemControl|FalseCode)",
    MIN_SPLIT_RATIO = 0.5, FOV_MIN_CELLS = 50, FOV_MIN_MEDIAN_COUNTS = 50, FOV_MIN_RATIO_TO_CORE = 0.5,
    QC_MAP_POINT_SIZE = 0.3, QC_MAP_MAX_CELLS = 1e6,
    INSPECT_WINDOW_UM = 150, INSPECT_N_PER_SLIDE = 3, INSPECT_CENTERS = NULL,
    N_HVG = 2000, RES_CLUSTERS = 0.8, TOP_N_PER_CLUSTER = 5, MARKERS_MAX_CELLS = 2000,
    CLUSTER_TO_CELLTYPE = NULL, CLUSTER_TO_LINEAGE = NULL,
    CMP_MIN_CELLS_CORE = 200, CMP_MIN_CELLS_GROUP = 50, CMP_REGISTER = TRUE, CMP_BIN_UM = 50,
    CMP_MIN_NCC = 0.3
  )
  added <- character(0)
  for (k in names(defaults)) {
    if (!exists(k, envir = globalenv(), inherits = FALSE)) {
      assign(k, defaults[[k]], envir = globalenv())
      added <- c(added, k)
    }
  }
  if (length(added) > 0) {
    message("Not in your config/config.R (older version) -> defaults used: ", paste(added, collapse = ", "),
            "\n  Copy them from the repository's config/config.R to change them.")
  }
})

## Per-cell totals as Seurat names them for the assay: nCount_CosMx / nFeature_CosMx
COUNT_COL   <- paste0("nCount_", ASSAY)
FEATURE_COL <- paste0("nFeature_", ASSAY)

## -----------------------------------------------------------------------
## Output helpers (same pattern as the class scripts)
## -----------------------------------------------------------------------
save_fig <- function(plot_obj, name, width = 9, height = 7, dpi = 150) {
  path <- file.path(FIGURES_DIR, paste0(name, ".png"))
  grDevices::png(path, width = width, height = height, units = "in", res = dpi)
  print(plot_obj)
  grDevices::dev.off()
  message("Saved figure -> ", path)
  invisible(path)
}

save_table <- function(df, name) {
  path <- file.path(TABLES_DIR, paste0(name, ".csv"))
  utils::write.csv(df, path, row.names = FALSE)
  message("Saved table  -> ", path)
  invisible(path)
}

obj_path <- function(stage, segmentation = SEGMENTATION) {
  file.path(PROCESSED_DIR, sprintf("cosmx_%s_%s.rds", segmentation, stage))
}

## -----------------------------------------------------------------------
## Config readers
## -----------------------------------------------------------------------
## read.csv() for sheets you edit by hand. Excel and most editors save the last
## row without a trailing newline; R then warns "incomplete final line", which is
## harmless (the row is read) -- so that one warning is muffled here.
read_config_csv <- function(path, ...) {
  withCallingHandlers(
    utils::read.csv(path, ...),
    warning = function(w) {
      if (grepl("incomplete final line", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

read_slides <- function(path = SLIDES_CSV) {
  s <- read_config_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot("slides.csv needs slide_id and cosmx_dir" =
              all(c("slide_id", "cosmx_dir") %in% colnames(s)))
  bad <- s$slide_id[s$slide_id != make.names(s$slide_id)]
  if (length(bad) > 0) {
    stop("slide_id must be a syntactically valid R name (used as FOV name): ",
         paste(bad, collapse = ", "), ". Use letters, digits, _ and . only.")
  }
  for (col in c("row_flip", "col_flip")) {
    if (col %in% colnames(s)) s[[col]] <- as.logical(s[[col]]) %in% TRUE
  }
  for (col in c("cosmx_prefix", "xenium_slide_id")) {
    if (!col %in% colnames(s)) s[[col]] <- NA_character_
    s[[col]] <- as.character(s[[col]])
    s[[col]][!is.na(s[[col]]) & !nzchar(s[[col]])] <- NA_character_
  }
  s
}

read_tma_map <- function(path = TMA_MAP_CSV) {
  m <- read_config_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("slide_id", "core_row", "core_col", "core_id", "patient_id", "organ")
  stopifnot("tma_map.csv is missing required columns" = all(need %in% colnames(m)))
  for (opt in c("tissue_type", "diagnosis", "location")) if (!opt %in% colnames(m)) m[[opt]] <- NA_character_
  if (!"include" %in% colnames(m)) m$include <- TRUE
  m$include <- as.logical(m$include) %in% TRUE
  m$core_row <- as.character(m$core_row)
  m$core_col <- as.character(m$core_col)
  dup <- m$core_id[duplicated(m$core_id)]
  if (length(dup) > 0) stop("core_id must be unique across slides: ", paste(unique(dup), collapse = ", "))
  m
}

## Seurat turns "_" into "-" in feature names. Apply the same to every gene /
## protein name read from config, or custom probes like "IL23_A" silently
## never match.
as_seurat_features <- function(x) gsub("_", "-", x, fixed = TRUE)

## Config gene names -> the spelling used in the object, ignoring upper / lower case,
## so the human-style names in config/ (EPCAM, PTPRC, CD3E) also find mouse genes
## (Epcam, Ptprc, Cd3e). Exact matches win; names not in `features` come back
## unchanged (and are then reported as missing by the caller).
match_case <- function(x, features) {
  hit <- match(toupper(x), toupper(features))
  ifelse(x %in% features | is.na(hit), x, features[hit])
}

split_markers <- function(x) {
  x <- trimws(unlist(strsplit(ifelse(is.na(x), "", x), ";")))
  as_seurat_features(x[nzchar(x)])
}

fov_key <- function(id) paste0("fov", tolower(gsub("[^[:alnum:]]", "", id)), "_")

## Returns list(level1 = named list, level2 = list(parent -> named list), mecr = named list)
read_marker_panel <- function(path) {
  p <- read_config_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  to_list <- function(df) stats::setNames(lapply(df$markers, split_markers), df$cell_type)
  out <- list(level1 = to_list(p[p$level == 1, ]), level2 = list(), mecr = list())
  for (par in unique(p$parent[p$level == 2])) {
    out$level2[[par]] <- to_list(p[p$level == 2 & p$parent == par, ])
  }
  if ("use_for_mecr" %in% colnames(p)) {
    l1 <- p[p$level == 1, ]
    out$mecr <- stats::setNames(lapply(l1$use_for_mecr, split_markers), l1$cell_type)
    out$mecr <- out$mecr[lengths(out$mecr) > 0]
  }
  out
}

## Keep only markers present in the panel; report what is missing. On the whole
## transcriptome nearly everything is present; a missing gene is then usually a
## different symbol (alias) or a multi-gene probe such as "TPSAB1/B2".
restrict_to_panel <- function(sig_list, panel_genes, label = "", min_reliable = 2) {
  sig_list <- lapply(sig_list, match_case, features = panel_genes)
  cov <- data.frame(
    cell_type = names(sig_list),
    n_markers = lengths(sig_list),
    n_present = vapply(sig_list, function(g) sum(g %in% panel_genes), integer(1)),
    missing = vapply(sig_list, function(g) paste(setdiff(g, panel_genes), collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
  weak <- cov$cell_type[cov$n_present < min_reliable]
  if (length(weak) > 0) {
    message(label, " signatures with <", min_reliable, " markers in panel (unreliable): ", paste(weak, collapse = ", "))
  }
  sig <- lapply(sig_list, intersect, panel_genes)
  list(signatures = sig[lengths(sig) > 0], coverage = cov)
}

## -----------------------------------------------------------------------
## Feature types. Next to the genes, both platforms ship control features;
## they are the basis of every background / specificity estimate in this
## project, so keep them (as per-cell metadata), don't discard them.
##   CosMx negative probes  ("Negative1", ...)      -> ControlProbe
##   CosMx system controls  ("SystemControl1", ...) -> ControlCodeword
## (the 10x names are kept so that a segmentation written in the common
## 10x-style format, and the Xenium object in 10, use the same columns)
## -----------------------------------------------------------------------
CONTROL_SHORT <- c(
  "Negative Control Probe"    = "ControlProbe",
  "Negative Control Codeword" = "ControlCodeword",
  "Unassigned Codeword"       = "BlankCodeword",
  "Blank Codeword"            = "BlankCodeword",
  "Genomic Control"           = "GenomicControl",
  "Deprecated Codeword"       = "DeprecatedCodeword"
)

control_short_name <- function(type) {
  ifelse(type %in% names(CONTROL_SHORT), CONTROL_SHORT[type], make.names(type))
}

cosmx_feature_type <- function(names) {
  ifelse(grepl(COSMX_NEG_REGEX, names), "Negative Control Probe",
         ifelse(grepl(COSMX_SYS_REGEX, names), "Negative Control Codeword", "Gene Expression"))
}

## -----------------------------------------------------------------------
## CosMx flat files (the "Flat Files" export of AtoMx). One set per slide:
##   <prefix>_exprMat_file.csv.gz        cells x targets counts (fov, cell_ID, genes, Negative*, SystemControl*)
##   <prefix>_metadata_file.csv.gz       one row per cell: position, area, stain intensities, AtoMx QC
##   <prefix>_fov_positions_file.csv.gz  position of every field of view on the slide
##   <prefix>_tx_file.csv.gz             one row per transcript (large; not needed for 01-05, used by 03)
##   <prefix>-polygons.csv.gz            cell outlines (only used for the zoom figures of 01)
## -----------------------------------------------------------------------
COSMX_FILE_PATTERNS <- c(
  expr     = "exprMat_file\\.csv(\\.gz)?$",
  metadata = "metadata_file\\.csv(\\.gz)?$",
  fov      = "fov_positions_file\\.csv(\\.gz)?$",
  tx       = "tx_file\\.csv(\\.gz)?$",
  polygons = "polygons\\.csv(\\.gz)?$"
)

## Finds the flat files of one slide anywhere below `dir` (AtoMx nests them in
## flatFiles/<slide>/). Returns a named list of paths, NA where a file is absent.
find_cosmx_files <- function(dir, prefix = NA) {
  if (is.na(dir) || !dir.exists(dir)) stop("CosMx folder not found: ", dir, " (cosmx_dir in config/slides.csv)")
  all <- list.files(dir, recursive = TRUE, full.names = TRUE)
  out <- lapply(COSMX_FILE_PATTERNS, function(pat) {
    f <- all[grepl(pat, basename(all), ignore.case = TRUE)]
    if (!is.na(prefix) && nzchar(prefix)) f <- f[startsWith(basename(f), prefix)]
    ## the same file zipped and unzipped: keep the zipped one
    f[!(paste0(f, ".gz") %in% f)]
  })
  multi <- names(out)[lengths(out) > 1]
  if (length(multi) > 0) {
    stop("Several ", paste(multi, collapse = " / "), " files below ", dir, ":\n  ",
         paste(basename(unlist(out[multi])), collapse = "\n  "),
         "\n  Point cosmx_dir at the folder of ONE slide, or set cosmx_prefix in slides.csv ",
         "(the start of this slide's file names).")
  }
  lapply(out, function(f) if (length(f) == 0) NA_character_ else f)
}

open_text <- function(path) if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")

csv_header <- function(path) {
  con <- open_text(path)
  on.exit(close(con))
  gsub('^"|"$', "", strsplit(readLines(con, n = 1), ",", fixed = TRUE)[[1]])
}

## fov + per-FOV cell number from whatever id columns a table has. `cell_ID`
## (capital ID) is the number of the cell WITHIN its FOV; `cell` / `cell_id`
## is the study-wide name c_<slide>_<fov>_<cell_ID>.
cosmx_cell_keys <- function(d) {
  cn <- colnames(d)
  num <- intersect(c("cell_ID", "cellID", "CellId"), cn)[1]
  if ("fov" %in% cn && !is.na(num)) {
    return(data.frame(fov = as.integer(d$fov), cell_ID = as.integer(d[[num]])))
  }
  name <- intersect(c("cell", "cell_id"), cn)[1]
  if (is.na(name)) stop("No cell identifier found (need fov + cell_ID, or cell = c_<slide>_<fov>_<cell_ID>).")
  parts <- regmatches(as.character(d[[name]]), regexec("_([0-9]+)_([0-9]+)$", as.character(d[[name]])))
  ok <- lengths(parts) == 3
  if (!all(ok)) stop("Cell names in column `", name, "` are not of the form c_<slide>_<fov>_<cell_ID>.")
  data.frame(fov = as.integer(vapply(parts, `[`, "", 2)), cell_ID = as.integer(vapply(parts, `[`, "", 3)))
}

## The expression matrix is a DENSE csv: one row per cell, one column per target.
## For the whole transcriptome that is ~19,000 columns, so reading it in one go
## needs (cells x 19,000 x 4 bytes) -- 38 GB for 500,000 cells. It is read in
## chunks of rows instead; each chunk becomes a sparse matrix straight away.
## Rows with cell_ID 0 are not cells: they hold the transcripts of one FOV that
## lie outside every cell. They are returned separately (`unassigned`).
read_cosmx_expr <- function(path, chunk_rows = COSMX_CHUNK_ROWS) {
  con <- open_text(path)
  on.exit(close(con))
  header <- gsub('^"|"$', "", strsplit(readLines(con, n = 1), ",", fixed = TRUE)[[1]])
  is_id <- header %in% c("fov", "cell_ID", "cellID", "cell", "cell_id", "V1", "")
  feats <- header[!is_id]
  if (length(feats) < 10) stop("Fewer than 10 target columns in ", path, " -- is this the exprMat file?")
  pieces <- list(); keys <- list(); n_done <- 0L
  repeat {
    lines <- readLines(con, n = chunk_rows)
    if (length(lines) == 0) break
    d <- data.table::fread(text = lines, header = FALSE, sep = ",", showProgress = FALSE, data.table = FALSE)
    if (ncol(d) != length(header)) stop("Row with ", ncol(d), " fields, header has ", length(header), ": ", path)
    colnames(d) <- header
    keys[[length(keys) + 1L]] <- cosmx_cell_keys(d[, is_id, drop = FALSE])
    m <- as.matrix(d[, !is_id, drop = FALSE])
    storage.mode(m) <- "double"
    pieces[[length(pieces) + 1L]] <- methods::as(methods::as(t(m), "CsparseMatrix"), "generalMatrix")
    n_done <- n_done + nrow(d)
    if (length(pieces) %% 25 == 0) message("  ... ", format(n_done, big.mark = ","), " rows")
  }
  if (length(pieces) == 0) stop("No rows in ", path)
  counts <- do.call(cbind, pieces)
  rm(pieces)
  keys <- do.call(rbind, keys)
  dimnames(counts) <- list(feats, paste(keys$fov, keys$cell_ID, sep = "_"))
  type <- cosmx_feature_type(feats)

  outside <- keys$cell_ID == 0
  unassigned <- NULL
  if (any(outside)) {
    u <- counts[, outside, drop = FALSE]
    unassigned <- data.frame(fov = keys$fov[outside],
                             n_outside_genes = Matrix::colSums(u[type == "Gene Expression", , drop = FALSE]),
                             n_outside_controls = Matrix::colSums(u[type != "Gene Expression", , drop = FALSE]))
    counts <- counts[, !outside, drop = FALSE]
  }
  dup <- duplicated(colnames(counts))
  if (any(dup)) stop(sum(dup), " duplicated fov + cell_ID in ", path, " -- several slides in one file? ",
                     "Export one flat-file set per slide.")
  message("  ", format(ncol(counts), big.mark = ","), " cells x ", sum(type == "Gene Expression"), " genes, ",
          sum(type == "Negative Control Probe"), " negative probes, ",
          sum(type == "Negative Control Codeword"), " system controls")
  list(counts = lapply(split(seq_along(feats), factor(type, levels = unique(type))),
                       function(i) counts[i, , drop = FALSE]),
       unassigned = unassigned)
}

## Slide frame: how CosMx global pixels become the microns used everywhere here.
cosmx_frame <- function(y_px, flip = COSMX_FLIP_Y, pixel_size = COSMX_PIXEL_SIZE_UM) {
  list(pixel_size_um = pixel_size, flip_y = isTRUE(flip), y_ref_px = if (isTRUE(flip)) max(y_px, na.rm = TRUE) else 0)
}
cosmx_px_to_um <- function(x_px, y_px, frame) {
  data.frame(x = x_px * frame$pixel_size_um,
             y = (if (frame$flip_y) frame$y_ref_px - y_px else y_px) * frame$pixel_size_um)
}

## Cell metadata (flat file or the meta.data of an AtoMx Seurat export) -> the
## cells table of the common format: cell_id, x_centroid, y_centroid (um),
## cell_area, nucleus_area (um2) + everything worth keeping from CosMx:
##   fov, cosmx_cell (AtoMx name), x_global_px / y_global_px (original frame),
##   split_ratio (FOV-border cells), IF_<stain> (mean stain intensity per cell,
##   e.g. IF_PanCK, IF_CD45, IF_DAPI) and AtoMx's own QC flags (atomx_*).
standardise_cosmx_metadata <- function(md, fov_positions = NULL) {
  cn <- colnames(md)
  keys <- cosmx_cell_keys(md)
  px <- COSMX_PIXEL_SIZE_UM
  if (all(c("CenterX_global_px", "CenterY_global_px") %in% cn)) {
    x_px <- md$CenterX_global_px; y_px <- md$CenterY_global_px
  } else if (all(c("x_slide_mm", "y_slide_mm") %in% cn)) {
    x_px <- md$x_slide_mm * 1000 / px; y_px <- md$y_slide_mm * 1000 / px
  } else if (all(c("CenterX_local_px", "CenterY_local_px") %in% cn) && !is.null(fov_positions)) {
    ## legacy exports without global centres: FOV origin + position in the FOV
    k <- match(keys$fov, fov_positions$fov)
    if (anyNA(k)) stop("FOVs in the metadata that are not in the fov_positions file.")
    x_px <- md$CenterX_local_px + fov_positions$x_global_px[k]
    y_px <- md$CenterY_local_px + fov_positions$y_global_px[k]
    message("  no CenterX_global_px in the metadata: global = FOV position + local position (legacy rule).")
  } else {
    stop("The metadata has no cell positions (CenterX_global_px / CenterY_global_px).")
  }
  frame <- cosmx_frame(y_px)
  xy <- cosmx_px_to_um(x_px, y_px, frame)
  area <- if ("Area.um2" %in% cn) md$Area.um2 else if ("Area" %in% cn) md$Area * px^2 else NA_real_
  nuc <- if ("NucArea" %in% cn) ifelse(md$NucArea > 0, md$NucArea * px^2, NA_real_) else NA_real_
  cells <- data.frame(
    cell_id = paste(keys$fov, keys$cell_ID, sep = "_"),
    x_centroid = xy$x, y_centroid = xy$y, cell_area = area, nucleus_area = nuc,
    segmentation_method = "cosmx",
    fov = keys$fov,
    cosmx_cell = if ("cell" %in% cn) as.character(md$cell) else if ("cell_id" %in% cn) as.character(md$cell_id) else NA_character_,
    x_global_px = x_px, y_global_px = y_px,
    split_ratio = if ("SplitRatioToLocal" %in% cn) md$SplitRatioToLocal else NA_real_,
    stringsAsFactors = FALSE
  )
  stains <- grep("^Mean\\.", cn, value = TRUE)
  for (s in stains) cells[[paste0("IF_", make.names(sub("^Mean\\.", "", s)))]] <- md[[s]]
  keep <- intersect(c("qcFlagsFOV", "qcCellsFlagged", "qcFlagsRNACounts", "qcFlagsCellCounts", "qcFlagsCellPropNeg",
                      "qcFlagsCellComplex", "qcFlagsCellArea", "Panel", "version", "Run_Tissue_name",
                      "cellSegmentationSetName", "unassignedTranscripts"), cn)
  for (k in keep) cells[[paste0("atomx_", k)]] <- md[[k]]
  cells <- cells[keys$cell_ID != 0, , drop = FALSE]
  attr(cells, "frame") <- frame
  cells
}

read_cosmx_fov_positions <- function(path) {
  if (is.na(path)) return(NULL)
  f <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  colnames(f)[tolower(colnames(f)) == "fov"] <- "fov"
  ## pixel columns where present; some exports only have mm
  for (ax in c("x", "y")) {
    px <- paste0(ax, "_global_px"); mm <- intersect(c(paste0(ax, "_global_mm"), paste0(toupper(ax), "_mm"), paste0(ax, "_mm")), colnames(f))[1]
    if (!px %in% colnames(f) && !is.na(mm)) f[[px]] <- f[[mm]] * 1000 / COSMX_PIXEL_SIZE_UM
  }
  f
}

read_cosmx_metadata <- function(path, fov_path = NA) {
  md <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  standardise_cosmx_metadata(md, read_cosmx_fov_positions(fov_path))
}

## Cell outlines of the polygons file, in the same microns as the centroids.
read_cosmx_polygons <- function(path, frame) {
  p <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  keys <- cosmx_cell_keys(p)
  xy <- cosmx_px_to_um(p$x_global_px, p$y_global_px, frame)
  data.frame(cell_id = paste(keys$fov, keys$cell_ID, sep = "_"), vertex_x = xy$x, vertex_y = xy$y)
}

## AtoMx can also export a Seurat object (.RDS) instead of flat files: counts in
## assay "RNA", negative probes in "negprobes", system controls in "falsecode",
## the flat-file metadata columns in meta.data. Written from the documented
## structure and tested on a mock object only -- check 00's report on a real one.
read_cosmx_seurat <- function(path) {
  so <- readRDS(path)
  get_counts <- function(assay) {
    if (!assay %in% Assays(so)) return(NULL)
    tryCatch(GetAssayData(so, assay = assay, layer = "counts"),
             error = function(e) GetAssayData(so, assay = assay, slot = "counts"))
  }
  gex <- get_counts("RNA") %||% get_counts(DefaultAssay(so))
  cells <- standardise_cosmx_metadata(so[[]])
  frame <- attr(cells, "frame")
  keys <- cosmx_cell_keys(so[[]])
  nm <- paste(keys$fov, keys$cell_ID, sep = "_")
  counts <- list("Gene Expression" = gex)
  type <- cosmx_feature_type(rownames(gex))
  if (any(type != "Gene Expression")) {                 # controls left inside the RNA assay
    for (t in setdiff(unique(type), "Gene Expression")) counts[[t]] <- gex[type == t, , drop = FALSE]
    counts[["Gene Expression"]] <- gex[type == "Gene Expression", , drop = FALSE]
  }
  if (!is.null(neg <- get_counts("negprobes"))) counts[["Negative Control Probe"]] <- neg
  if (!is.null(sys <- get_counts("falsecode"))) counts[["Negative Control Codeword"]] <- sys
  counts <- lapply(counts, function(m) { colnames(m) <- nm; m[, keys$cell_ID != 0, drop = FALSE] })
  attr(cells, "frame") <- frame
  list(counts = counts, cells = cells, unassigned = NULL, frame = frame, files = list(seurat = path))
}

## One slide from flat files. The parsed expression matrix is cached next to
## the processed objects (parsing 19,000 columns takes minutes); the cache is
## rebuilt when the exprMat file changes.
read_cosmx_flat <- function(dir, prefix = NA, slide_id = basename(dir), cache = TRUE) {
  files <- find_cosmx_files(dir, prefix)
  if (is.na(files$expr) || is.na(files$metadata)) {
    stop("No *exprMat_file.csv(.gz) and/or *metadata_file.csv(.gz) below ", dir,
         "\n  Run 00_cosmx_inspect_input.R to see what the folder contains.")
  }
  cache_file <- file.path(PROCESSED_DIR, "cosmx_cache", paste0(make.names(slide_id), "_exprMat.rds"))
  stamp <- paste(normalizePath(files$expr), file.size(files$expr), as.numeric(file.mtime(files$expr)), COSMX_NEG_REGEX,
                 COSMX_SYS_REGEX)
  expr <- NULL
  if (cache && file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (identical(cached$stamp, stamp)) { expr <- cached$expr; message("  expression matrix from cache: ", cache_file) }
  }
  if (is.null(expr)) {
    message("  reading ", basename(files$expr), " in chunks of ", COSMX_CHUNK_ROWS, " rows ...")
    expr <- read_cosmx_expr(files$expr)
    if (cache) {
      dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
      saveRDS(list(stamp = stamp, expr = expr), cache_file)
    }
  }
  cells <- read_cosmx_metadata(files$metadata, files$fov)
  list(counts = expr$counts, cells = cells, unassigned = expr$unassigned, frame = attr(cells, "frame"), files = files)
}

## -----------------------------------------------------------------------
## Common segmentation format (10x-style), for a segmentation made outside
## AtoMx: <dir>/cell_feature_matrix/{matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz}
## and <dir>/cells.csv.gz (cell_id, x_centroid, y_centroid, cell_area, nucleus_area;
## microns in the frame of 01, see tables/01_cosmx_frame.csv).
## -----------------------------------------------------------------------
read_feature_matrix <- function(seg_dir) {
  mtx_dir <- file.path(seg_dir, "cell_feature_matrix")
  h5 <- file.path(seg_dir, "cell_feature_matrix.h5")
  counts <- if (dir.exists(mtx_dir)) {
    Seurat::Read10X(mtx_dir)
  } else if (file.exists(h5)) {
    Seurat::Read10X_h5(h5)
  } else {
    stop("No cell_feature_matrix/ or cell_feature_matrix.h5 in ", seg_dir)
  }
  if (!is.list(counts)) counts <- list("Gene Expression" = counts)
  ## a writer that did not label the controls: split them off by name
  gex <- counts[["Gene Expression"]]
  type <- cosmx_feature_type(rownames(gex))
  for (t in setdiff(unique(type), "Gene Expression")) {
    counts[[t]] <- if (is.null(counts[[t]])) gex[type == t, , drop = FALSE] else rbind(counts[[t]], gex[type == t, , drop = FALSE])
  }
  counts[["Gene Expression"]] <- gex[type == "Gene Expression", , drop = FALSE]
  counts
}

read_cells_table <- function(seg_dir) {
  pq <- file.path(seg_dir, "cells.parquet")
  cells <- if (file.exists(pq)) {
    as.data.frame(arrow::read_parquet(pq))
  } else if (file.exists(file.path(seg_dir, "cells.csv.gz"))) {
    utils::read.csv(file.path(seg_dir, "cells.csv.gz"), stringsAsFactors = FALSE)
  } else {
    stop("No cells.parquet or cells.csv.gz in ", seg_dir)
  }
  cells$cell_id <- as.character(cells$cell_id)
  for (col in c("cell_area", "nucleus_area")) if (!col %in% colnames(cells)) cells[[col]] <- NA_real_
  if (!"segmentation_method" %in% colnames(cells)) cells$segmentation_method <- NA_character_
  cells
}

## Builds one Seurat object per slide with an identical structure for every
## segmentation, so downstream code never branches on it:
##   - assay ASSAY ("CosMx") = genes only
##   - controls summed into metadata: nCount_ControlProbe (negative probes),
##     nCount_ControlCodeword (system controls) plus nFeatTotal_<type> (number of
##     control features on the panel), needed to turn control counts into a
##     per-probe background rate
##   - x_um / y_um centroids in metadata AND as a centroid FOV named slide_id
##     ("FOV" here is Seurat's word for one image = one slide; the CosMx field of
##     view of a cell is the metadata column `fov`)
##   - cells named <slide_id>_<fov>_<cell_ID> (CosMx cell ids repeat across FOVs and slides)
## Only centroids are loaded (no polygons): polygons for multi-slide TMAs cost
## many GB and are not needed for any statistic here.
build_slide_object <- function(slide_id, counts_list, cells, segmentation) {
  gex <- counts_list[["Gene Expression"]]
  stopifnot("No 'Gene Expression' features found" = !is.null(gex))

  common <- intersect(colnames(gex), cells$cell_id)
  if (length(common) < ncol(gex)) {
    message(slide_id, ": ", ncol(gex) - length(common),
            " cells in the matrix have no metadata row (no centroid) and are dropped.")
  }
  if (length(common) < nrow(cells)) {
    message(slide_id, ": ", nrow(cells) - length(common), " cells in the metadata are not in the matrix and are dropped.")
  }
  if (length(common) == 0) stop(slide_id, ": no cell is in both the expression matrix and the metadata.")
  gex <- gex[, common, drop = FALSE]
  cells <- cells[match(common, cells$cell_id), , drop = FALSE]
  new_names <- paste0(slide_id, "_", common)
  colnames(gex) <- new_names

  meta <- data.frame(
    slide_id = slide_id,
    cell_id = common,
    x_um = cells$x_centroid,
    y_um = cells$y_centroid,
    cell_area = cells$cell_area,
    nucleus_area = cells$nucleus_area,
    seg_boundary_source = cells$segmentation_method,
    segmentation = segmentation,
    row.names = new_names,
    stringsAsFactors = FALSE
  )
  extra <- setdiff(colnames(cells), c("cell_id", "x_centroid", "y_centroid", "cell_area", "nucleus_area",
                                      "segmentation_method", "core_id"))
  for (col in extra) meta[[col]] <- cells[[col]]
  for (type in setdiff(names(counts_list), "Gene Expression")) {
    short <- control_short_name(type)
    m <- counts_list[[type]][, common, drop = FALSE]
    prev <- if (paste0("nCount_", short) %in% colnames(meta)) meta[[paste0("nCount_", short)]] else 0
    meta[[paste0("nCount_", short)]] <- prev + Matrix::colSums(m)
    prev_n <- if (paste0("nFeatTotal_", short) %in% colnames(meta)) meta[[paste0("nFeatTotal_", short)]] else 0
    meta[[paste0("nFeatTotal_", short)]] <- prev_n + nrow(m)
  }
  meta$nFeatTotal_Gene <- nrow(gex)

  obj <- CreateSeuratObject(counts = gex, assay = ASSAY, meta.data = meta)
  fov <- CreateFOV(
    coords = data.frame(x = meta$x_um, y = meta$y_um, cell = rownames(meta)),
    type = "centroids",
    assay = ASSAY,
    key = fov_key(slide_id)
  )
  obj[[slide_id]] <- fov
  obj
}

seg_dir_for <- function(slide_row, segmentation) {
  if (segmentation == "cosmx") return(resolve_path(slide_row$cosmx_dir))
  file.path(SEG_DIR, segmentation, slide_row$slide_id)
}

## Returns list(obj, frame, unassigned, files). The input kind is recognised
## from what is there: an .rds file (AtoMx Seurat export), a folder with
## cell_feature_matrix/ (common format) or a folder with flat files.
load_slide <- function(slide_row, segmentation = SEGMENTATION) {
  d <- seg_dir_for(slide_row, segmentation)
  sid <- slide_row$slide_id
  message("Loading ", sid, " [", segmentation, "] from ", d)
  if (!is.na(d) && grepl("\\.rds$", d, ignore.case = TRUE) && file.exists(d)) {
    x <- read_cosmx_seurat(d)
  } else if (!is.na(d) && dir.exists(d) &&
             (dir.exists(file.path(d, "cell_feature_matrix")) || file.exists(file.path(d, "cell_feature_matrix.h5")))) {
    x <- list(counts = read_feature_matrix(d), cells = read_cells_table(d), unassigned = NULL, frame = NULL,
              files = list(common_format = d))
  } else {
    x <- read_cosmx_flat(d, slide_row$cosmx_prefix %||% NA,
                         slide_id = if (segmentation == "cosmx") sid else paste0(sid, "_", segmentation))
  }
  list(obj = build_slide_object(sid, x$counts, x$cells, segmentation),
       frame = x$frame, unassigned = x$unassigned, files = x$files)
}

load_segmentation <- function(slide_row, segmentation = SEGMENTATION) load_slide(slide_row, segmentation)$obj

## merge() keeps one FOV per slide (named by slide_id). JoinLayers() because
## v5 merge splits counts into counts.1, counts.2, ... per input object.
merge_slides <- function(objs) {
  if (length(objs) == 1) return(objs[[1]])
  obj <- merge(objs[[1]], y = objs[-1])
  obj <- JoinLayers(obj)
  ## control feature totals are per-slide constants; merge may leave NA for
  ## slides whose panel lacked a control class -- set those to 0
  for (col in grep("^(nCount|nFeatTotal)_", colnames(obj[[]]), value = TRUE)) {
    v <- obj[[col, drop = TRUE]]
    v[is.na(v)] <- 0
    obj[[col]] <- v
  }
  obj
}

## -----------------------------------------------------------------------
## TMA dearraying
## -----------------------------------------------------------------------

## Grid-based clustering: bins of `cell` um, bins with >= 2 cells are tissue,
## 8-connected tissue bins form one cluster. Same result as DBSCAN for
## well-separated TMA cores, with no compiled dependency (fallback when the
## dbscan package cannot be loaded). Returns a cluster id per point, 0 = noise.
grid_components <- function(xy, cell) {
  gx <- floor((xy[, 1] - min(xy[, 1])) / cell) + 2      # +2: one empty border bin
  gy <- floor((xy[, 2] - min(xy[, 2])) / cell) + 2
  nx <- max(gx) + 1; ny <- max(gy) + 1
  lin <- gx + (gy - 1) * nx
  occ <- matrix(tabulate(lin, nbins = nx * ny) >= 2, nx, ny)
  lab <- matrix(NA_real_, nx, ny)
  lab[occ] <- seq_len(sum(occ))
  shift <- function(m, dx, dy) {                        # m[i + dx, j + dy], NA outside
    out <- matrix(NA_real_, nx, ny)
    xs <- max(1, 1 - dx):min(nx, nx - dx); ys <- max(1, 1 - dy):min(ny, ny - dy)
    out[xs, ys] <- m[xs + dx, ys + dy]
    out
  }
  repeat {                                              # propagate the smallest label
    old <- lab
    for (d in list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))) {
      nb <- shift(lab, d[1], d[2])
      lab[occ] <- pmin(lab[occ], nb[occ], na.rm = TRUE)
    }
    if (identical(old, lab)) break
  }
  cl <- lab[lin]
  out <- integer(length(cl))                           # 0 = noise
  ok <- !is.na(cl)
  out[ok] <- as.integer(factor(cl[ok]))
  out
}

## DBSCAN on cell centroids -> core seeds -> merge fragments -> core table.
## Returns data.frame(core_idx, x_center, y_center, radius, n_cells_dbscan, detected_core)
detect_cores <- function(x, y,
                         eps = DBSCAN_EPS_UM, min_pts = DBSCAN_MIN_PTS,
                         min_cells = MIN_CELLS_PER_CORE, merge_dist = CORE_MERGE_DIST_UM) {
  xy <- cbind(x, y)
  if (is.na(eps)) {
    ## Adaptive: 3 x the median distance to the min_pts-th neighbour. Works for
    ## dense lymphoid cores and sparse lung/adipose alike, while staying far
    ## below the >= 200 um gaps between cores.
    eps <- 3 * stats::median(RANN::nn2(xy, k = min_pts + 1)$nn.dists[, min_pts + 1])
    message("Clustering radius eps (adaptive): ", round(eps, 1), " um")
  }
  cl <- if (requireNamespace("dbscan", quietly = TRUE)) {
    tryCatch(dbscan::dbscan(xy, eps = eps, minPts = min_pts)$cluster, error = function(e) NULL)
  }
  if (is.null(cl)) {
    message("dbscan package not usable -- using grid-based clustering instead")
    cl <- grid_components(xy, eps)
  }
  cl[cl == 0] <- NA
  if (all(is.na(cl))) stop("DBSCAN found no clusters -- increase DBSCAN_EPS_UM.")

  cent <- data.frame(cl = cl, x = x, y = y) |>
    filter(!is.na(cl)) |>
    group_by(cl) |>
    summarise(x = median(x), y = median(y), n = n(), .groups = "drop")

  ## Seeds = big clusters. Seeds closer than merge_dist are one torn core.
  seeds <- cent[cent$n >= min_cells, ]
  if (nrow(seeds) == 0) stop("No cluster reaches MIN_CELLS_PER_CORE cells.")
  seed_group <- if (nrow(seeds) > 1) {
    stats::cutree(stats::hclust(stats::dist(seeds[, c("x", "y")]), method = "single"), h = merge_dist)
  } else 1L
  lut <- stats::setNames(seed_group, seeds$cl)

  ## Small fragments attach to the nearest seed if within merge_dist.
  small <- cent[cent$n < min_cells, ]
  if (nrow(small) > 0) {
    nn <- RANN::nn2(seeds[, c("x", "y")], small[, c("x", "y")], k = 1)
    ok <- nn$nn.dists[, 1] <= merge_dist
    lut <- c(lut, stats::setNames(seed_group[nn$nn.idx[ok, 1]], small$cl[ok]))
  }
  core_idx <- unname(lut[as.character(cl)])

  cores <- data.frame(core_idx = core_idx, x = x, y = y) |>
    filter(!is.na(core_idx)) |>
    group_by(core_idx) |>
    summarise(
      x_center = median(x), y_center = median(y),
      radius = stats::quantile(sqrt((x - median(x))^2 + (y - median(y))^2), 0.95),
      n_cells_dbscan = n(),
      .groups = "drop"
    ) |>
    arrange(y_center, x_center) |>
    mutate(detected_core = seq_len(n())) |>
    as.data.frame()
  cores
}

## Nearest core centre, but only within CORE_MAX_DIST_FACTOR x radius. Works
## for any segmentation's centroids once cores are defined.
assign_cells_to_cores <- function(x, y, cores, factor = CORE_MAX_DIST_FACTOR) {
  if (nrow(cores) == 0) return(rep(NA_character_, length(x)))     # slide without cores
  nn <- RANN::nn2(cores[, c("x_center", "y_center")], cbind(x, y), k = 1)
  idx <- nn$nn.idx[, 1]
  ok <- nn$nn.dists[, 1] <= cores$radius[idx] * factor
  out <- cores$core_id[idx]
  out[!ok] <- NA
  out
}

## Rows = clusters of y centres, columns = clusters of x centres. Row A is at
## the TOP of the image (smallest y) unless row_flip -- set *_flip in
## slides.csv if the TMA was mounted rotated. ALWAYS check the labelled plot.
assign_grid <- function(cores, row_levels, col_levels, row_flip = FALSE, col_flip = FALSE) {
  cluster_axis <- function(v, k) {
    k <- min(k, length(unique(round(v))))
    if (k <= 1) return(rep(1L, length(v)))
    g <- stats::cutree(stats::hclust(stats::dist(v), method = "average"), k = k)
    rank_of_group <- rank(tapply(v, g, mean))
    as.integer(rank_of_group[as.character(g)])
  }
  r <- cluster_axis(cores$y_center, length(row_levels))
  c <- cluster_axis(cores$x_center, length(col_levels))
  if (max(r) < length(row_levels) || max(c) < length(col_levels)) {
    warning("Fewer detected rows/columns than in tma_map.csv (missing cores?). ",
            "Row/column labels may be shifted -- verify the dearray plot and use ",
            "config/core_overrides.csv if needed.")
  }
  if (row_flip) r <- max(r) + 1L - r
  if (col_flip) c <- max(c) + 1L - c
  cores$core_row <- row_levels[r]
  cores$core_col <- col_levels[c]
  cores
}

## TMA row names: 1 -> A, 26 -> Z, 27 -> AA, ...
LETTERS_EXT <- function(i) {
  vapply(i, function(k) {
    s <- ""
    while (k > 0) { s <- paste0(LETTERS[(k - 1) %% 26 + 1], s); k <- (k - 1) %/% 26 }
    s
  }, character(1))
}

sort_levels <- function(v) {
  v <- unique(as.character(v))
  num <- suppressWarnings(as.numeric(v))
  if (!anyNA(num)) v[order(num)] else sort(v)
}

## Optional polygon mode: one CSV per core in config/core_selections/<slide_id>/,
## file name = core_id, with X/Y columns in microns of the 01 frame (x_um / y_um,
## as on the axes of the 01 maps; lines starting with # are skipped).
read_core_selections <- function(sel_dir) {
  files <- list.files(sel_dir, pattern = "\\.csv$", full.names = TRUE)
  lapply(stats::setNames(files, tools::file_path_sans_ext(basename(files))), function(f) {
    d <- read_config_csv(f, comment.char = "#", check.names = FALSE)
    xcol <- grep("^x$", colnames(d), ignore.case = TRUE, value = TRUE)[1]
    ycol <- grep("^y$", colnames(d), ignore.case = TRUE, value = TRUE)[1]
    stopifnot("Selection CSV needs X and Y columns" = !is.na(xcol) && !is.na(ycol))
    data.frame(x = d[[xcol]], y = d[[ycol]])
  })
}

assign_cells_to_polygons <- function(x, y, polygons) {
  out <- rep(NA_character_, length(x))
  for (nm in names(polygons)) {
    p <- polygons[[nm]]
    inside <- sp::point.in.polygon(x, y, p$x, p$y) > 0
    out[inside & is.na(out)] <- nm
  }
  out
}

## -----------------------------------------------------------------------
## Scoring and annotation
## -----------------------------------------------------------------------

## Adds one score column per signature named score_<signature>.
score_signatures <- function(obj, sigs, use_ucell = USE_UCELL, assay = ASSAY) {
  sigs <- sigs[lengths(sigs) > 0]
  if (length(sigs) == 0) stop("score_signatures: no signature has any gene on the panel -- see the panel coverage table")
  DefaultAssay(obj) <- assay
  if (use_ucell && requireNamespace("UCell", quietly = TRUE)) {
    ## serial (ncores = 1): no BiocParallel workers; if UCell still fails (version /
    ## object-specific problems), say why and fall back to AddModuleScore below
    res <- tryCatch(
      UCell::AddModuleScore_UCell(obj, features = sigs, name = "_UCell", assay = assay, ncores = 1),
      error = function(e) {
        message("UCell failed (", conditionMessage(e), ") -- falling back to Seurat::AddModuleScore.")
        NULL
      })
    if (!is.null(res)) {
      for (s in names(sigs)) {
        res[[paste0("score_", s)]] <- res[[paste0(s, "_UCell"), drop = TRUE]]
        res[[paste0(s, "_UCell")]] <- NULL
      }
      return(res)
    }
  } else if (use_ucell) {
    message("UCell not installed -- falling back to AddModuleScore.")
  }
  ## control genes are sampled per expression bin, so ctrl must stay well
  ## below the bin size (a few hundred genes / nbin) on a targeted panel
  n_genes <- nrow(obj[[assay]])
  nbin <- 12
  obj <- AddModuleScore(obj, features = sigs, name = "tmpscore", assay = assay,
                        ctrl = max(2, min(50, floor(n_genes / (2 * nbin)))), nbin = nbin)
  for (i in seq_along(sigs)) {
    obj[[paste0("score_", names(sigs)[i])]] <- obj[[paste0("tmpscore", i), drop = TRUE]]
    obj[[paste0("tmpscore", i)]] <- NULL
  }
  obj
}

## Dot plot of marker genes per group (colour = scaled mean expression, size = % of
## cells with counts), computed straight from the assay's data layer. Replaces
## Seurat::DotPlot, which goes through FetchData() and fails ("non-numeric argument
## to mathematical function") when a metadata column is returned instead of a gene.
marker_dotplot <- function(obj, features, group.by, assay = ASSAY, layer = "data", log_data = TRUE) {
  features <- intersect(unique(features), rownames(obj[[assay]]))
  if (length(features) == 0) stop("marker_dotplot: none of the features are in assay ", assay)
  m <- GetAssayData(obj, assay = assay, layer = layer)[features, , drop = FALSE]
  if (!methods::is(m, "CsparseMatrix")) m <- methods::as(as.matrix(m), "CsparseMatrix")
  grp <- droplevels(as.factor(obj[[group.by, drop = TRUE]]))
  keep <- !is.na(grp)
  m <- m[, keep, drop = FALSE]; grp <- grp[keep]
  G <- Matrix::sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1,
                            dims = c(length(grp), nlevels(grp)), dimnames = list(NULL, levels(grp)))
  n <- Matrix::colSums(G)
  lin <- m
  if (log_data) lin@x <- expm1(lin@x)                       # log-normalised -> linear, like DotPlot
  avg <- as.matrix(lin %*% G) / rep(n, each = nrow(m))
  pct <- as.matrix((m > 0) %*% G) / rep(n, each = nrow(m)) * 100
  sc <- t(scale(t(avg)))
  sc[is.na(sc)] <- 0
  sc <- pmin(pmax(sc, -2.5), 2.5)
  d <- data.frame(feature = factor(rep(features, times = ncol(avg)), levels = features),
                  group = factor(rep(colnames(avg), each = nrow(avg)), levels = levels(grp)),
                  avg_scaled = as.vector(sc), pct = as.vector(pct))
  ggplot(d, aes(feature, group)) +
    geom_point(aes(size = pct, colour = avg_scaled)) +
    scale_colour_gradient(low = "lightgrey", high = "blue", name = "Average\nexpression\n(scaled)") +
    scale_size(range = c(0, 6), limits = c(0, 100), name = "% expressed") +
    labs(x = "Features", y = group.by) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

## Semi-manual annotation. Writes an editable CSV (cluster, suggested label,
## top markers, label, reviewed). You edit `label`, set reviewed = TRUE, and
## re-run: the file is then used as-is. This keeps the manual decision in a
## versionable file instead of a hard-coded vector in the script.
resolve_cluster_map <- function(clusters, score_df, map_path, top_markers = NULL) {
  clusters <- as.character(clusters)
  mean_scores <- stats::aggregate(score_df, by = list(cluster = clusters), FUN = mean)
  sig_names <- sub("^score_", "", colnames(score_df))
  m <- as.matrix(mean_scores[, -1, drop = FALSE])
  best <- apply(m, 1, which.max)
  srt <- t(apply(m, 1, sort, decreasing = TRUE))
  margin <- if (ncol(m) > 1) srt[, 1] - srt[, 2] else srt[, 1]
  template <- data.frame(
    cluster = mean_scores$cluster,
    n_cells = as.integer(table(clusters)[mean_scores$cluster]),
    suggested = sig_names[best],
    margin = round(margin, 3),
    top_markers = if (is.null(top_markers)) "" else unname(top_markers[mean_scores$cluster]),
    label = sig_names[best],
    reviewed = FALSE,
    stringsAsFactors = FALSE
  )
  template$top_markers[is.na(template$top_markers)] <- ""

  if (file.exists(map_path)) {
    existing <- read_config_csv(map_path, stringsAsFactors = FALSE, colClasses = c(cluster = "character"))
    if (setequal(existing$cluster, template$cluster)) {
      if (!any(as.logical(existing$reviewed) %in% TRUE)) {
        message("Using UNREVIEWED labels from ", map_path, " -- edit `label`, set reviewed=TRUE.")
      }
      return(stats::setNames(existing$label, existing$cluster))
    }
    bak <- paste0(map_path, ".", format(Sys.time(), "%Y%m%d%H%M%S"), ".bak")
    file.copy(map_path, bak)
    warning("Clusters changed since ", map_path, " was written (backup: ", bak, "). Rewriting template.")
  }
  utils::write.csv(template, map_path, row.names = FALSE)
  message("Wrote annotation template -> ", map_path,
          "\n  Low `margin` = ambiguous cluster; check it on the dotplot before accepting.")
  stats::setNames(template$label, template$cluster)
}

## Top-N positive markers per cluster as "G1;G2;G3" strings, on downsampled
## cells (FindAllMarkers on 1M cells is not a good use of a day).
top_markers_per_cluster <- function(obj, ident, n = 5, max_cells = 500) {
  Idents(obj) <- ident
  mk <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
                       max.cells.per.ident = max_cells, verbose = FALSE)
  if (nrow(mk) == 0) return(NULL)
  mk |>
    group_by(cluster) |>
    slice_max(avg_log2FC, n = n) |>
    summarise(genes = paste(gene, collapse = ";"), .groups = "drop") |>
    (\(d) stats::setNames(d$genes, as.character(d$cluster)))()
}

## -----------------------------------------------------------------------
## Quantitative helpers
## -----------------------------------------------------------------------

## Sum counts by group: returns features x groups (sparse).
pseudobulk <- function(counts, group) {
  g <- factor(group)
  counts %*% Matrix::t(Matrix::fac2sparse(g))
}

## Mutually exclusive co-expression rate (Hartman & Satija 2024 idea): for
## marker pairs from DIFFERENT lineages, fraction of cells expressing both
## among cells expressing either. Lower = cleaner segmentation (less
## transcript spill-over between neighbouring cells).
compute_mecr <- function(counts, lineage_markers) {
  lineage_markers <- lapply(lineage_markers, function(g) intersect(match_case(g, rownames(counts)), rownames(counts)))
  lineage_markers <- lineage_markers[lengths(lineage_markers) > 0]
  genes <- unique(unlist(lineage_markers))
  if (length(genes) < 2) return(NA_real_)
  B <- counts[genes, , drop = FALSE] > 0
  lin_of <- stats::setNames(rep(names(lineage_markers), lengths(lineage_markers)), unlist(lineage_markers))
  pairs <- utils::combn(genes, 2)
  pairs <- pairs[, lin_of[pairs[1, ]] != lin_of[pairs[2, ]], drop = FALSE]
  if (ncol(pairs) == 0) return(NA_real_)
  rates <- apply(pairs, 2, function(p) {
    a <- B[p[1], ]; b <- B[p[2], ]
    either <- sum(a | b)
    if (either == 0) NA_real_ else sum(a & b) / either
  })
  mean(rates, na.rm = TRUE)
}

## Tissue-specificity index tau (Yanai 2005): 0 = uniform, 1 = one organ only.
tau_index <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2 || max(x) <= 0) return(NA_real_)
  sum(1 - x / max(x)) / (length(x) - 1)
}

## -----------------------------------------------------------------------
## Plotting
## -----------------------------------------------------------------------

## ggplot on metadata centroids instead of ImageDimPlot(): fast on millions of
## cells, never depends on @images surviving a subset, and can facet by core.
## axes = TRUE draws micron axes + grid, so you can read off coordinates (e.g.
## for INSPECT_CENTERS in config.R); theme_void otherwise.
## Cell-type colours: one hue family per LINEAGE (Epithelial blue-purple,
## Immune orange-red, Fibroblast green, ...). Within a lineage, types differ in
## hue AND lightness (dark / mid / light per hue), so neighbouring cells of two
## epithelial types are still told apart in the tissue. Same colours in every
## 04 figure (UMAP, sections, composition, gallery).
LINEAGE_HUES <- c(Epithelial = 245, Immune = 35, Fibroblast = 150,      # HCL hue centres (degrees)
                  Smooth_muscle_pericyte = 80, Endothelial = 345, Neural_glia = 115)
celltype_palette <- function(cell_type, lineage) {
  d <- unique(data.frame(cell_type = as.character(cell_type), lineage = as.character(lineage)))
  d <- d[!is.na(d$cell_type) & !duplicated(d$cell_type), ]
  d$lineage[is.na(d$lineage) | !d$lineage %in% names(LINEAGE_HUES)] <- "other"
  d <- d[order(match(d$lineage, c(names(LINEAGE_HUES), "other")), d$cell_type), ]
  cols <- unlist(lapply(split(d$cell_type, factor(d$lineage, levels = unique(d$lineage))), function(ct) {
    lin <- d$lineage[match(ct[1], d$cell_type)]
    n <- length(ct)
    if (lin == "other") return(stats::setNames(rep("#8a8984", n), ct))
    k <- ceiling(n / 3)                                    # hues, 3 lightness steps each
    span <- min(110, 37 * (k - 1))
    hues <- LINEAGE_HUES[[lin]] + seq(-span / 2, span / 2, length.out = k)
    i <- seq_len(n) - 1
    j <- if (n == 1) 2 else (i %% 3) + 1                   # 1 dark, 2 mid, 3 light
    stats::setNames(grDevices::hcl(h = hues[i %/% 3 + 1], c = c(50, 70, 45)[j], l = c(33, 57, 80)[j],
                                   fixup = TRUE), ct)
  }), use.names = FALSE)
  stats::setNames(cols, d$cell_type)                      # names in lineage order -> legend order
}

## ggplot point size (mm, stroke 0) for a dot about one cell (cell_um) wide when
## extent_um of tissue fills panel_in inches: the tissue looks filled instead of
## dusted with specks (100k cells) or covered by blobs (a few hundred cells).
## (A filled ggplot dot is ~0.75 x its size in mm across, hence / 0.75.)
cell_point_size <- function(extent_um, panel_in, cell_um = 12, min_size = 0.1) {
  max(min_size, cell_um * panel_in * 25.4 / extent_um / 0.75)
}

## UMAP coloured by a grouping column, with large readable names at each group's
## centre. pal: named colours (default: ggplot hues); label = FALSE for a legend only.
plot_umap_labelled <- function(obj, group, pal = NULL, reduction = "umap", label = TRUE, label_size = 5.5,
                               pt_size = 0.6, max_cells = 3e5, title = group, legend = TRUE) {
  emb <- Embeddings(obj, reduction)[, 1:2]
  d <- data.frame(u1 = emb[, 1], u2 = emb[, 2], grp = as.character(obj[[group, drop = TRUE]]))
  lev <- if (is.factor(obj[[group, drop = TRUE]])) levels(obj[[group, drop = TRUE]]) else sort(unique(d$grp))
  if (is.null(pal)) pal <- stats::setNames(scales::hue_pal()(length(lev)), lev)
  d$grp <- factor(d$grp, levels = names(pal))
  if (nrow(d) > max_cells) d <- d[sample(nrow(d), max_cells), ]
  p <- ggplot(d[sample(nrow(d)), ], aes(u1, u2, colour = grp)) +           # shuffled: no group paints over all
    geom_point(size = pt_size, stroke = 0, alpha = 0.8) +
    scale_colour_manual(values = pal, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), ncol = 1)) +
    coord_equal() + theme_classic(base_size = 14) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size = 12)) +
    labs(x = "UMAP 1", y = "UMAP 2", title = title)
  if (label) {
    lab <- d |> group_by(grp) |> summarise(u1 = stats::median(u1), u2 = stats::median(u2), .groups = "drop")
    p <- p + ggrepel::geom_label_repel(data = lab, aes(label = grp), colour = "#0b0b0b", fill = "#ffffffcc",
                                       label.size = 0, size = label_size, fontface = "bold", seed = 1,
                                       box.padding = 0.4, min.segment.length = 0, max.overlaps = Inf)
  }
  if (!legend) p <- p + theme(legend.position = "none")
  p
}

plot_celltype_umap <- function(obj, reduction = "umap", label_size = 5.5, pt_size = 0.6, max_cells = 3e5) {
  plot_umap_labelled(obj, "cell_type", pal = celltype_palette(obj$cell_type, obj$lineage), reduction = reduction,
                     label_size = label_size, pt_size = pt_size, max_cells = max_cells,
                     title = "Cell types (transcriptomics)")
}

## Both sections side by side, cells coloured by cell type (same colours as the
## UMAP). panel_in: size of one section panel in inches (figure width =
## panel_in * number of slides + ~4 for the legend); dots are sized to ~one cell.
plot_celltype_sections <- function(meta, panel_in = 11, pt_size = NULL, cell_um = 12) {
  meta <- meta[!is.na(meta$cell_type), ]
  pal <- celltype_palette(meta$cell_type, meta$lineage)
  meta$cell_type <- factor(meta$cell_type, levels = names(pal))
  ps <- lapply(sort(unique(meta$slide_id)), function(sid) {
    d <- meta[meta$slide_id == sid, ]
    sz <- pt_size %||% cell_point_size(max(diff(range(d$x_um)), diff(range(d$y_um))), panel_in, cell_um)
    ggplot(d[sample(nrow(d)), ], aes(x_um, y_um, colour = cell_type)) +
      geom_point(size = sz, stroke = 0) +
      scale_colour_manual(values = pal, name = "Cell type", drop = FALSE) +
      scale_y_reverse() + coord_fixed() + theme_void(base_size = 16) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text = element_text(size = 13)) +
      ggtitle(sprintf("%s (%s cells)", sid, format(nrow(d), big.mark = ",")))
  })
  patchwork::wrap_plots(ps, nrow = 1) + patchwork::plot_layout(guides = "collect") &
    guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1))
}

## Where is each cluster? One panel per cluster: that cluster's cells in colour on
## top of all cells of the slide in grey (a background layer without the facet
## column is drawn in every panel).
plot_clusters_spatial <- function(meta, cluster_col, slide = NULL, max_bg = 1e5, size = NULL, ncol = NULL,
                                  colour = "firebrick") {
  d <- meta
  if (!is.null(slide)) d <- d[d$slide_id == slide, , drop = FALSE]
  d$cluster <- factor(d[[cluster_col]])
  ## point size from cell number: visible on a few hundred cells, not a blob on 100k+
  if (is.null(size)) size <- max(0.35, min(1.5, 25 / sqrt(nrow(d))))
  bg <- d[sample(nrow(d), min(nrow(d), max_bg)), c("x_um", "y_um")]
  ggplot(d, aes(x_um, y_um)) +
    geom_point(data = bg, colour = "grey80", size = size, stroke = 0) +
    geom_point(colour = colour, size = size, stroke = 0) +
    facet_wrap(~ cluster, ncol = ncol) +
    scale_y_reverse() + coord_fixed() + theme_void() +
    theme(strip.text = element_text(size = 10, face = "bold"),
          panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3)) +
    ggtitle(paste0(cluster_col, if (!is.null(slide)) paste0(" -- ", slide) else ""))
}

plot_cores_spatial <- function(meta, color_by, cores = NULL, facet = TRUE,
                               max_cells = 3e5, size = 0.15, title = NULL, axes = FALSE) {
  d <- meta
  if (!is.null(cores)) d <- d[d$core_id %in% cores, , drop = FALSE]
  ## sort() keeps the caller's row order, so rows placed last are drawn on top
  if (nrow(d) > max_cells) d <- d[sort(sample(nrow(d), max_cells)), , drop = FALSE]
  p <- ggplot(d, aes(x_um, y_um, color = .data[[color_by]])) +
    geom_point(size = size, stroke = 0) +
    ggtitle(title %||% color_by)
  p <- if (axes) {
    p + scale_x_continuous(breaks = scales::breaks_pretty(10)) +
      scale_y_reverse(breaks = scales::breaks_pretty(10)) +   # image convention: y grows downwards
      theme_minimal(base_size = 9) +
      theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
            panel.grid.minor = element_blank()) +
      labs(x = "x (um)", y = "y (um)")
  } else {
    p + scale_y_reverse() + theme_void()
  }
  ## coord_fixed() cannot be combined with free facet scales; cores are ~round,
  ## so a square panel keeps them undistorted enough.
  do_facet <- facet && "core_id" %in% colnames(d)
  p <- p + if (do_facet) theme(aspect.ratio = 1) else coord_fixed()
  if (is.numeric(d[[color_by]])) {
    p <- p + scale_color_viridis_c()
  } else {
    p <- p + guides(color = guide_legend(override.aes = list(size = 3)))
  }
  if (do_facet) {
    p <- p + facet_wrap(~ core_id, scales = "free") + theme(strip.text = element_text(size = 7))
  }
  p
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## -----------------------------------------------------------------------
## Zoomed QC inspection: cells (pass / fail) with their outlines
## -----------------------------------------------------------------------

## First failing QC reason per cell ("pass" if none) -- for colouring.
qc_reason <- function(meta) {
  r <- rep("pass", nrow(meta))
  r[meta$qc_high_control] <- "high control"
  if ("qc_fov_border" %in% colnames(meta)) r[meta$qc_fov_border] <- "FOV border"
  r[meta$qc_large] <- "too large"
  r[meta$qc_small] <- "too small"
  r[meta$qc_low_counts] <- "low counts"
  factor(r, levels = c("pass", "low counts", "too small", "too large", "FOV border", "high control"))
}
QC_REASON_COLORS <- c(`pass` = "#9ec5f4", `low counts` = "#eda100", `too small` = "#e34948",
                      `too large` = "#eb6834", `FOV border` = "#8e6bbf", `high control` = "#e87ba4")

## Windows (square, window_um wide) that contain BOTH passing and failing cells:
## the slide is cut into a grid and the windows with the most balanced mix
## (largest min(n_pass, n_fail)) are returned, n_per_slide per slide.
pick_inspection_windows <- function(meta, window_um = 150, n_per_slide = 3, min_cells = 30) {
  d <- meta[!is.na(meta$x_um), c("slide_id", "x_um", "y_um", "qc_pass")]
  d$gx <- floor(d$x_um / window_um)
  d$gy <- floor(d$y_um / window_um)
  w <- d |>
    group_by(slide_id, gx, gy) |>
    summarise(n = n(), n_pass = sum(qc_pass), n_fail = sum(!qc_pass), .groups = "drop") |>
    filter(n >= min_cells) |>
    mutate(score = pmin(n_pass, n_fail)) |>
    group_by(slide_id) |>
    slice_max(score, n = n_per_slide, with_ties = FALSE) |>
    ungroup()
  data.frame(slide_id = w$slide_id,
             x0 = w$gx * window_um, y0 = w$gy * window_um,
             x1 = (w$gx + 1) * window_um, y1 = (w$gy + 1) * window_um,
             n_cells = w$n, n_pass = w$n_pass, n_fail = w$n_fail)
}

## Two panels: where the window is on the slide | the cells of the window,
## outlines (polygons file) and centroids coloured by QC result, FOV borders as
## dashed lines. y is plotted as -y (image convention: y grows downwards).
## polys = read_cosmx_polygons() of the slide, or NULL (centroids only).
plot_inspection <- function(meta_slide, polys, win, title = "", point_size = 1.2) {
  cells <- meta_slide[meta_slide$x_um >= win$x0 - 20 & meta_slide$x_um <= win$x1 + 20 &
                        meta_slide$y_um >= win$y0 - 20 & meta_slide$y_um <= win$y1 + 20, ]
  cells$qc <- qc_reason(cells)
  frame <- list(coord_fixed(xlim = c(win$x0, win$x1), ylim = c(-win$y1, -win$y0), expand = FALSE),
                theme_void(), theme(plot.title = element_text(size = 9),
                                    panel.background = element_rect(fill = "#101418", colour = NA)))
  p_cells <- ggplot()
  if (!is.null(polys)) {
    b <- polys[polys$cell_id %in% cells$cell_id, ]
    b$qc <- cells$qc[match(b$cell_id, cells$cell_id)]
    if (nrow(b) > 0) {
      p_cells <- p_cells + geom_polygon(data = b, aes(vertex_x, -vertex_y, group = cell_id, colour = qc),
                                        fill = NA, linewidth = 0.35)
    }
  }
  ## FOV borders inside the window: where one FOV's cells end and the next begin
  if ("fov" %in% colnames(meta_slide)) {
    ext <- meta_slide |> group_by(fov) |>
      summarise(x0 = min(x_um), x1 = max(x_um), y0 = min(y_um), y1 = max(y_um), .groups = "drop")
    p_cells <- p_cells + geom_rect(data = ext, aes(xmin = x0, xmax = x1, ymin = -y1, ymax = -y0),
                                   fill = NA, colour = "grey60", linetype = 2, linewidth = 0.2)
  }
  inw <- cells$x_um >= win$x0 & cells$x_um <= win$x1 & cells$y_um >= win$y0 & cells$y_um <= win$y1
  p_cells <- p_cells +
    geom_point(data = cells, aes(x_um, -y_um, colour = qc), size = point_size) +
    scale_colour_manual(values = QC_REASON_COLORS, drop = FALSE, name = "cell QC") +
    frame + ggtitle(sprintf("cells in window: %d pass, %d fail (%s; dashed = FOV extent)",
                            sum(cells$qc[inw] == "pass"), sum(cells$qc[inw] != "pass"),
                            if (is.null(polys)) "centroids only, no polygons file" else "outlines = segmentation"))

  ov <- meta_slide[sample(nrow(meta_slide), min(nrow(meta_slide), 1e5)), ]
  p_where <- ggplot(ov, aes(x_um, -y_um)) +
    geom_point(size = 0.05, colour = "grey60", stroke = 0) +
    annotate("rect", xmin = win$x0, xmax = win$x1, ymin = -win$y1, ymax = -win$y0,
             fill = NA, colour = "#e34948", linewidth = 0.8) +
    coord_fixed() + theme_void() + ggtitle("window position")
  (p_where | p_cells) + plot_layout(widths = c(1, 1.6)) + plot_annotation(title = title)
}
