## =============================================================================
## 00_cosmx_inspect_input.R
##
## "What did the CosMx run deliver, and can 01 read it?" -- run this FIRST, on
## every new export. It only reads file headers and the cell metadata (seconds
## to a minute per slide), never the big expression / transcript files.
##
## For every slide in config/slides.csv it reports
##   * which flat files exist below cosmx_dir (and their size)
##   * the expression matrix: number of genes, negative probes, system controls
##     (as classified by COSMX_NEG_REGEX / COSMX_SYS_REGEX in config.R)
##   * the metadata: cells, FOVs, which columns 01 will use, panel / assay
##   * the pixel size implied by the file (Area.um2 / Area) vs COSMX_PIXEL_SIZE_UM
##   * the size of the imaged area and a map of the FOVs with their numbers
##   * other export types lying around (raw images, TileDB, Seurat .RDS)
## and writes a template for config/fov_to_core.csv (02, mode D: FOV -> core).
##
## What CosMx delivers (AtoMx export types):
##   Flat files   csv tables: THIS is what the pipeline reads.
##                  *_exprMat_file.csv.gz        counts, cells x targets (dense)
##                  *_metadata_file.csv.gz       one row per cell
##                  *_fov_positions_file.csv.gz  where each field of view is
##                  *_tx_file.csv.gz             one row per transcript (very large)
##                  *-polygons.csv.gz            cell outlines
##   Raw files    CellStatsDir/ (morphology images, cell label masks per FOV),
##                AnalysisResults/, RunSummary/. Only needed to re-segment or to
##                look at the stains; not read here.
##   Seurat .RDS  01 can read it (cosmx_dir = path of the .RDS), flat files are preferred.
##   TileDB       not read; export flat files instead.
##
## INPUT : config/slides.csv
## OUTPUT: tables/00_cosmx_input_report.csv, tables/00_fov_to_core_template.csv,
##         figures/00_fov_map_<slide>.png
## =============================================================================

source("config/config.R")
source("R/helpers.R")

slides <- read_slides()
size_gb <- function(f) if (is.na(f)) NA_real_ else round(file.size(f) / 1e9, 3)
USED_BY_01 <- c("fov", "cell_ID", "CenterX_global_px", "CenterY_global_px", "Area", "Area.um2", "NucArea",
                "SplitRatioToLocal", "cell", "cell_id")

report <- list(); fov_template <- list()
add <- function(sid, item, value, note = "") {
  report[[length(report) + 1]] <<- data.frame(slide_id = sid, item = item, value = as.character(value), note = note)
}

for (i in seq_len(nrow(slides))) {
  sid <- slides$slide_id[i]
  d <- resolve_path(slides$cosmx_dir[i])
  message("\n================ ", sid, " ================\n", d)
  if (is.na(d) || !file.exists(d)) { add(sid, "cosmx_dir", d, "NOT FOUND -- fix cosmx_dir in config/slides.csv"); next }
  if (grepl("\\.rds$", d, ignore.case = TRUE)) {
    add(sid, "input type", "Seurat .RDS (AtoMx export)", "01 reads assays RNA / negprobes / falsecode and meta.data")
    next
  }

  ## ---- other export types in the folder -----------------------------------
  sub <- list.dirs(d, recursive = TRUE, full.names = FALSE)
  for (raw in c("CellStatsDir", "AnalysisResults", "RunSummary")) {
    hit <- sub[basename(sub) == raw]
    if (length(hit) > 0) add(sid, paste0("raw folder ", raw), hit[1], "raw export: images / masks per FOV, not read by this pipeline")
  }
  if (any(grepl("tiledb", sub, ignore.case = TRUE))) add(sid, "TileDB export", "present", "not read: export flat files from AtoMx")
  rds <- list.files(d, pattern = "\\.rds$", ignore.case = TRUE, recursive = TRUE)
  if (length(rds) > 0) add(sid, "Seurat .RDS", rds[1], "can be used instead of flat files: set cosmx_dir to this file")

  files <- tryCatch(find_cosmx_files(d, slides$cosmx_prefix[i]), error = function(e) { message(conditionMessage(e)); NULL })
  if (is.null(files)) { add(sid, "flat files", "ambiguous", "several slides in one folder: set cosmx_prefix in slides.csv"); next }
  for (k in names(files)) {
    add(sid, paste0("file: ", k), if (is.na(files[[k]])) "MISSING" else basename(files[[k]]),
        if (is.na(files[[k]])) c(expr = "REQUIRED by 01", metadata = "REQUIRED by 01", fov = "optional",
                                 tx = "optional (03: transcripts inside / outside cells, by compartment)",
                                 polygons = "optional (01: cell outlines in the zoom figures)")[[k]]
        else paste(size_gb(files[[k]]), "GB"))
  }

  ## ---- expression matrix: header only ----------------------------------------
  if (!is.na(files$expr)) {
    h <- csv_header(files$expr)
    feats <- h[!h %in% c("fov", "cell_ID", "cellID", "cell", "cell_id", "V1", "")]
    type <- cosmx_feature_type(feats)
    add(sid, "exprMat: id columns", paste(setdiff(h, feats), collapse = ", "))
    add(sid, "exprMat: genes", sum(type == "Gene Expression"),
        paste("e.g.", paste(utils::head(feats[type == "Gene Expression"], 6), collapse = ", ")))
    add(sid, "exprMat: negative probes", sum(type == "Negative Control Probe"),
        paste(utils::head(feats[type == "Negative Control Probe"], 3), collapse = ", "))
    add(sid, "exprMat: system controls", sum(type == "Negative Control Codeword"),
        paste(utils::head(feats[type == "Negative Control Codeword"], 3), collapse = ", "))
    multi <- grep("/", feats[type == "Gene Expression"], value = TRUE)
    if (length(multi) > 0) add(sid, "exprMat: multi-gene probes (name contains /)", length(multi),
                               paste(utils::head(multi, 6), collapse = ", "))
    odd <- feats[type == "Gene Expression" & grepl("^(Neg|Sys|False|Blank|Control|Custom)", feats, ignore.case = TRUE)]
    if (length(odd) > 0) add(sid, "exprMat: counted as GENES but look like controls", paste(odd, collapse = ", "),
                             "CHECK COSMX_NEG_REGEX / COSMX_SYS_REGEX in config.R")
    if (sum(type == "Negative Control Probe") == 0) {
      add(sid, "exprMat: WARNING", "no negative probes found",
          "05 and 10 need them as background model: check COSMX_NEG_REGEX against the column names")
    }
  }

  ## ---- metadata -----------------------------------------------------------------
  if (!is.na(files$metadata)) {
    md <- data.table::fread(files$metadata, data.table = FALSE, showProgress = FALSE)
    cn <- colnames(md)
    add(sid, "metadata: cells", nrow(md))
    add(sid, "metadata: FOVs", length(unique(md$fov)))
    add(sid, "metadata: columns used by 01", paste(intersect(USED_BY_01, cn), collapse = ", "),
        paste("absent:", paste(setdiff(USED_BY_01, cn), collapse = ", ")))
    add(sid, "metadata: stain intensity columns", paste(grep("^Mean\\.", cn, value = TRUE), collapse = ", "),
        "kept per cell as IF_<stain>")
    for (k in intersect(c("Panel", "assay_type", "version", "Run_Tissue_name", "cellSegmentationSetName"), cn)) {
      add(sid, paste0("metadata: ", k), paste(utils::head(unique(md[[k]]), 3), collapse = " | "),
          if (k == "assay_type" && any(md[[k]] != "RNA")) "NOT an RNA assay?" else "")
    }
    atomx <- grep("^qc|cellType|nbclust|RNA_", cn, value = TRUE)
    if (length(atomx) > 0) add(sid, "metadata: AtoMx analysis columns", paste(utils::head(atomx, 12), collapse = ", "),
                               "QC flags are kept as atomx_*; AtoMx cell types are not used")
    if (all(c("Area", "Area.um2") %in% cn)) {
      implied <- stats::median(sqrt(md$Area.um2 / md$Area), na.rm = TRUE)
      add(sid, "pixel size implied by Area.um2 / Area", round(implied, 5),
          if (abs(implied - COSMX_PIXEL_SIZE_UM) > 0.002) paste0("DIFFERS from COSMX_PIXEL_SIZE_UM = ", COSMX_PIXEL_SIZE_UM,
                                                                " -> set it in config.R") else "matches COSMX_PIXEL_SIZE_UM")
    } else {
      add(sid, "pixel size", COSMX_PIXEL_SIZE_UM, "no Area.um2 column: COSMX_PIXEL_SIZE_UM from config.R is trusted")
    }
    cells <- tryCatch(standardise_cosmx_metadata(md, read_cosmx_fov_positions(files$fov)),
                      error = function(e) { add(sid, "metadata: ERROR", conditionMessage(e)); NULL })
    if (!is.null(cells)) {
      add(sid, "imaged area (bounding box, mm)", sprintf("%.1f x %.1f", diff(range(cells$x_centroid)) / 1000,
                                                         diff(range(cells$y_centroid)) / 1000))
      add(sid, "median cell area (um2)", round(stats::median(cells$cell_area, na.rm = TRUE), 1),
          "whole-cell segmentation: typically 40-150")
      fv <- cells |> group_by(fov) |>
        summarise(n_cells = n(), x_um = round(mean(x_centroid)), y_um = round(mean(y_centroid)),
                  w_um = diff(range(x_centroid)), h_um = diff(range(y_centroid)), .groups = "drop")
      add(sid, "FOV size from the cells (um)", sprintf("%.0f x %.0f", stats::median(fv$w_um), stats::median(fv$h_um)),
          "CosMx FOV = 4256 px = 512 um")
      add(sid, "cells per FOV (median)", stats::median(fv$n_cells))
      fov_template[[sid]] <- data.frame(slide_id = sid, fov = fv$fov, n_cells = fv$n_cells, x_um = fv$x_um, y_um = fv$y_um,
                                        core_id = "")
      bg <- cells[sample(nrow(cells), min(nrow(cells), 2e5)), ]
      p <- ggplot(bg, aes(x_centroid, y_centroid)) +
        geom_point(size = 0.05, colour = "grey70", stroke = 0) +
        geom_text(data = fv, aes(x_um, y_um, label = fov), size = 2.6, fontface = "bold", colour = "#c0392b") +
        scale_y_reverse() + coord_fixed() + theme_bw() +
        labs(x = "x (um)", y = "y (um)", title = paste0(sid, ": cells and FOV numbers (", nrow(fv), " FOVs)"),
             subtitle = "Orientation as the pipeline uses it. Compare with the AtoMx slide overview: if it is mirrored, change COSMX_FLIP_Y")
      save_fig(p, paste0("00_fov_map_", sid), width = 12, height = 12)
    }
  }

  ## ---- transcripts / polygons: header + a few rows --------------------------------
  for (k in c("tx", "polygons")) {
    if (is.na(files[[k]])) next
    head_k <- data.table::fread(files[[k]], nrows = 1000, data.table = FALSE, showProgress = FALSE)
    add(sid, paste0(k, ": columns"), paste(colnames(head_k), collapse = ", "))
    if (k == "tx" && "CellComp" %in% colnames(head_k)) {
      add(sid, "tx: CellComp values (first rows)", paste(unique(head_k$CellComp), collapse = ", "),
          "None / 0 = outside every cell")
    }
  }
}

report <- dplyr::bind_rows(report)
show_report <- function(r) {
  for (sid in unique(r$slide_id)) {
    cat("\n", sid, "\n", sep = "")
    rs <- r[r$slide_id == sid, ]
    cat(sprintf("  %-46s %s%s\n", rs$item, rs$value, ifelse(nzchar(rs$note), paste0("   [", rs$note, "]"), "")), sep = "")
  }
}
show_report(report)
save_table(report, "00_cosmx_input_report")
if (length(fov_template) > 0) {
  save_table(dplyr::bind_rows(fov_template), "00_fov_to_core_template")
  message("FOV -> core by hand (optional, 02 mode D): fill in core_id in tables/00_fov_to_core_template.csv ",
          "using figures/00_fov_map_<slide>.png and save it as config/fov_to_core.csv")
}
problems <- report[grepl("MISSING|NOT FOUND|WARNING|ERROR|DIFFERS|CHECK", paste(report$value, report$note)), ]
problems <- problems[!grepl("optional", problems$note), ]
if (nrow(problems) > 0) {
  message("\n>>> Look at these before running 01:")
  show_report(problems)
} else message("\nNothing blocking: run 01_cosmx_object_prep.R")
