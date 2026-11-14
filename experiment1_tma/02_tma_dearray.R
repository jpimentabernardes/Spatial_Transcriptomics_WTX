## =============================================================================
## 02_tma_dearray.R
##
## Assigns every cell to a TMA core and attaches the sample metadata
## (patient, organ, tissue type, diagnosis, location) from config/tma_map.csv.
## Then does core-level QC and applies the cell filter from 01.
##
## No TMA (whole tissue section, e.g. a Swiss roll)? Give the slide exactly ONE
## row in tma_map.csv (core_row A, core_col 1, your sample id as core_id): the
## whole section then becomes one sample ("whole-section mode", C below).
##
## Three ways to define cores on a TMA (per slide, picked automatically):
##   A. Polygon mode -- if config/core_selections/<slide_id>/ exists: one CSV
##      per core (file name = core_id) with X,Y vertices in microns of the 01
##      frame (x_um / y_um, as on the axes of figures/01_*_spatial_by_slide.png).
##      Most accurate for torn/irregular tumour cores.
##   D. FOV mode (CosMx only) -- if config/fov_to_core.csv has rows for the
##      slide (slide_id, fov, core_id): every cell of a field of view belongs
##      to that core. You placed the FOVs on the cores yourself, so this is
##      often the most direct route; 00 writes a template and a map of the FOV
##      numbers. Only valid if no FOV covers two cores.
##   B. Automatic mode -- DBSCAN on cell centroids finds cores, rows/columns
##      are inferred from the core positions and matched to core_row/core_col
##      in tma_map.csv. Fast, but the grid inference MUST be checked on the
##      labelled plot (missing cores or a rotated TMA shift labels; set
##      row_flip / col_flip in slides.csv for a slide mounted upside down).
##      Fix any remaining mistakes with config/core_overrides.csv:
##         slide_id,detected_core,core_row,core_col
##
## Use the SAME core_id as in the Xenium project's tma_map.csv: 10 matches the
## two platforms core by core.
##
## Core definitions are saved once (tables/02_cores_<slide>.csv) and reused
## when you re-run with another segmentation, so core identity never depends
## on the segmentation method. DELETE that file to re-dearray (e.g. after
## editing core_overrides.csv or the DBSCAN parameters).
##
## CosMx-specific: after dearraying, every field of view is compared with the
## other FOVs of the SAME core (tables/02_fov_qc_by_core.csv). A FOV with far
## fewer transcripts per cell than its neighbours on the same tissue is a run
## problem (focus, bubble), flagged as `fov_flag`; cells are not removed for it.
##
## INPUT : cosmx_<seg>_raw.rds (01), config/tma_map.csv
## OUTPUT: cosmx_<seg>_cores.rds (cells in included cores, QC-passing)
##         tables/02_cores_<slide>.csv, 02_cores_all.csv, 02_core_qc.csv,
##         02_fov_qc_by_core.csv
## =============================================================================

source("config/config.R")
source("R/helpers.R")

cmx <- readRDS(obj_path("raw"))
slides <- read_slides()
tma_map <- read_tma_map()
meta <- cmx[[]]

overrides_path <- file.path(CONFIG_DIR, "core_overrides.csv")
overrides <- if (file.exists(overrides_path)) {
  read_config_csv(overrides_path, stringsAsFactors = FALSE, colClasses = "character")
} else NULL
fov_map <- if (file.exists(FOV_TO_CORE_CSV)) {
  f <- read_config_csv(FOV_TO_CORE_CSV, stringsAsFactors = FALSE)
  stopifnot("fov_to_core.csv needs slide_id, fov, core_id" = all(c("slide_id", "fov", "core_id") %in% colnames(f)))
  f[!is.na(f$core_id) & nzchar(f$core_id), ]
} else NULL

## centre / radius / bounding box of a core from the cells assigned to it
cores_from_cells <- function(x, y, core_id) {
  ok <- !is.na(core_id)
  data.frame(x = x[ok], y = y[ok], core_id = core_id[ok]) |>
    group_by(core_id) |>
    summarise(x_center = (min(x) + max(x)) / 2, y_center = (min(y) + max(y)) / 2,
              radius = max(sqrt((x - x_center)^2 + (y - y_center)^2)),
              xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y), .groups = "drop") |>
    as.data.frame()
}

meta$core_id <- NA_character_
all_cores <- list()

for (i in seq_len(nrow(slides))) {
  sid <- slides$slide_id[i]
  in_slide <- meta$slide_id == sid
  cores_file <- file.path(TABLES_DIR, paste0("02_cores_", sid, ".csv"))
  sel_dir <- file.path(CONFIG_DIR, "core_selections", sid)
  map_s <- tma_map[tma_map$slide_id == sid, ]
  if (nrow(map_s) == 0) stop("No rows for slide ", sid, " in tma_map.csv")

  if (nrow(map_s) == 1 && !dir.exists(sel_dir)) {
    ## ---- C. whole section: one sample per slide (no TMA) ---------------------
    ## A single tma_map row for this slide means "the whole section is one
    ## sample" (e.g. a Swiss roll or a whole-tissue section): every cell gets
    ## that core_id and no cores are searched for.
    cid <- map_s$core_id
    xs <- meta$x_um[in_slide]; ys <- meta$y_um[in_slide]
    cx <- stats::median(xs); cy <- stats::median(ys)
    meta$core_id[in_slide] <- cid
    cores <- data.frame(core_id = cid, x_center = cx, y_center = cy,
                        radius = max(sqrt((xs - cx)^2 + (ys - cy)^2)),
                        xmin = min(xs), xmax = max(xs), ymin = min(ys), ymax = max(ys))
    message(sid, ": one row in tma_map.csv -> whole-section mode, all cells = ", cid)

  } else if (dir.exists(sel_dir)) {
    ## ---- A. polygon mode --------------------------------------------------
    polys <- read_core_selections(sel_dir)
    unknown <- setdiff(names(polys), map_s$core_id)
    if (length(unknown) > 0) warning(sid, ": selection files not in tma_map: ", paste(unknown, collapse = ", "))
    meta$core_id[in_slide] <- assign_cells_to_polygons(meta$x_um[in_slide], meta$y_um[in_slide], polys)
    cores <- do.call(rbind, lapply(names(polys), function(nm) {
      p <- polys[[nm]]
      cx <- mean(p$x); cy <- mean(p$y)
      data.frame(core_id = nm, x_center = cx, y_center = cy,
                 radius = max(sqrt((p$x - cx)^2 + (p$y - cy)^2)),
                 xmin = min(p$x), xmax = max(p$x), ymin = min(p$y), ymax = max(p$y))
    }))
    message(sid, ": polygon mode, ", nrow(cores), " cores.")

  } else if (!is.null(fov_map) && any(fov_map$slide_id == sid)) {
    ## ---- D. FOV mode: fov -> core from config/fov_to_core.csv ---------------
    fm <- fov_map[fov_map$slide_id == sid, ]
    unknown <- setdiff(fm$core_id, map_s$core_id)
    if (length(unknown) > 0) warning(sid, ": core_id in fov_to_core.csv not in tma_map: ", paste(unknown, collapse = ", "))
    meta$core_id[in_slide] <- fm$core_id[match(meta$fov[in_slide], fm$fov)]
    cores <- cores_from_cells(meta$x_um[in_slide], meta$y_um[in_slide], meta$core_id[in_slide])
    not_mapped <- setdiff(unique(meta$fov[in_slide]), fm$fov)
    message(sid, ": FOV mode, ", nrow(cores), " cores from ", nrow(fm), " FOVs",
            if (length(not_mapped) > 0) paste0("; FOVs without a core (cells dropped): ", paste(sort(not_mapped), collapse = ", ")) else "")

  } else if (file.exists(cores_file)) {
    ## ---- previously defined cores (e.g. re-run with another segmentation) --
    cores <- utils::read.csv(cores_file, stringsAsFactors = FALSE)
    cores <- cores[!is.na(cores$core_id), ]
    meta$core_id[in_slide] <- assign_cells_to_cores(meta$x_um[in_slide], meta$y_um[in_slide], cores)
    message(sid, ": reusing core definitions from ", cores_file)

  } else {
    ## ---- B. automatic dearraying ------------------------------------------
    ## Use cells with some signal so empty debris segments don't bridge cores.
    use <- in_slide & meta$nCount_CosMx >= MIN_COUNTS
    cores <- detect_cores(meta$x_um[use], meta$y_um[use])
    cores <- assign_grid(
      cores,
      row_levels = sort_levels(map_s$core_row),
      col_levels = sort_levels(map_s$core_col),
      row_flip = isTRUE(slides$row_flip[i]),
      col_flip = isTRUE(slides$col_flip[i])
    )
    if (!is.null(overrides)) {
      ov <- overrides[overrides$slide_id == sid, ]
      for (j in seq_len(nrow(ov))) {
        k <- cores$detected_core == as.integer(ov$detected_core[j])
        cores$core_row[k] <- ov$core_row[j]
        cores$core_col[k] <- ov$core_col[j]
      }
    }
    cores <- cores |>
      left_join(map_s[, c("core_row", "core_col", "core_id")], by = c("core_row", "core_col"))
    dups <- cores$core_id[!is.na(cores$core_id) & duplicated(cores$core_id)]
    if (length(dups) > 0) {
      warning(sid, ": several detected cores map to the same grid position: ",
              paste(unique(dups), collapse = ", "), " -- check the plot and add overrides.")
    }
    ## Bounding boxes cover the same area as the cell assignment: radius x CORE_MAX_DIST_FACTOR.
    half <- cores$radius * CORE_MAX_DIST_FACTOR
    cores$xmin <- cores$x_center - half
    cores$xmax <- cores$x_center + half
    cores$ymin <- cores$y_center - half
    cores$ymax <- cores$y_center + half
    unmatched <- cores[is.na(cores$core_id), ]
    if (nrow(unmatched) > 0) {
      message(sid, ": ", nrow(unmatched), " detected cores have no tma_map entry (detected_core ",
              paste(unmatched$detected_core, collapse = ", "), ") -- cells there are dropped.")
    }
    utils::write.csv(cores, cores_file, row.names = FALSE)
    cores <- cores[!is.na(cores$core_id), ]
    meta$core_id[in_slide] <- assign_cells_to_cores(meta$x_um[in_slide], meta$y_um[in_slide], cores)
  }

  cores$slide_id <- sid
  all_cores[[sid]] <- cores

  ## ---- the plot you MUST look at: detected cores + labels vs the TMA map --
  d <- meta[in_slide, ]
  d$core_label <- ifelse(is.na(d$core_id), "unassigned", "core")
  if (nrow(d) > 3e5) d <- d[sample(nrow(d), 3e5), ]
  lab <- cores
  ## drop the "<slide_id>_" prefix from labels -- the slide is in the title
  short_id <- sub(paste0("^", sid, "_"), "", lab$core_id)
  lab$text <- if ("detected_core" %in% colnames(lab)) {
    paste0(short_id, "\n(", lab$core_row, lab$core_col, ", #", lab$detected_core, ")")
  } else short_id
  p_dearray <- ggplot(d, aes(x_um, y_um)) +
    geom_point(aes(color = core_label), size = 0.05, stroke = 0) +
    geom_text(data = lab, aes(x_center, y_center, label = text), size = 2.5, fontface = "bold") +
    scale_color_manual(values = c(core = "grey60", unassigned = "red")) +
    scale_y_reverse() + coord_fixed() + theme_bw() +
    ggtitle(paste0(sid, ": core assignment -- compare every label with your TMA map"))
  save_fig(p_dearray, paste0("02_dearray_", sid), width = 12, height = 12)
}

cores_all <- dplyr::bind_rows(all_cores)
save_table(cores_all, "02_cores_all")

## -----------------------------------------------------------------------
## Attach TMA metadata to cells
## -----------------------------------------------------------------------
map_cols <- setdiff(colnames(tma_map), c("slide_id", "core_row", "core_col", "notes"))
meta_join <- meta[, c("slide_id", "core_id")] |>
  left_join(tma_map[, map_cols], by = "core_id")
for (col in setdiff(map_cols, "core_id")) cmx[[col]] <- meta_join[[col]]
cmx$core_id <- meta$core_id

message("Cells outside any core: ", sum(is.na(cmx$core_id)),
        " (", round(100 * mean(is.na(cmx$core_id)), 1), "%)")

## -----------------------------------------------------------------------
## Core-level QC. Flag, don't delete: whether a poor core is excluded is an
## explicit decision (set include = FALSE in tma_map.csv), not a side effect.
## -----------------------------------------------------------------------
m <- cmx[[]]
core_qc <- m |>
  filter(!is.na(core_id)) |>
  group_by(slide_id, core_id, patient_id, organ, tissue_type) |>
  summarise(
    n_cells = n(),
    n_fovs = n_distinct(fov),
    n_pass = sum(qc_pass),
    pct_pass = 100 * mean(qc_pass),
    median_counts = median(nCount_CosMx[qc_pass]),
    median_genes = median(nFeature_CosMx[qc_pass]),
    median_area = median(cell_area[qc_pass], na.rm = TRUE),
    control_frac = sum(nCount_controls) / sum(nCount_CosMx + nCount_controls),
    .groups = "drop"
  ) |>
  mutate(core_flag = case_when(
    n_pass < CORE_MIN_CELLS ~ "few_cells",
    median_counts < CORE_MIN_MEDIAN_COUNTS ~ "low_counts",
    TRUE ~ "ok"
  ))
print(as.data.frame(core_qc))
save_table(core_qc, "02_core_qc")

## -----------------------------------------------------------------------
## FOV QC within cores (CosMx): each FOV against the other FOVs of its core.
## Same tissue, same organ -> a large difference is technical.
## -----------------------------------------------------------------------
fov_core <- m |>
  filter(!is.na(core_id), !is.na(fov)) |>
  group_by(slide_id, core_id, fov) |>
  summarise(n_cells = n(), median_counts = median(nCount_CosMx), pct_pass = 100 * mean(qc_pass), .groups = "drop") |>
  group_by(core_id) |>
  mutate(n_fovs_in_core = n(),
         ## median of the OTHER FOVs of the core (NA if the core has a single FOV)
         core_other_median = vapply(seq_along(median_counts), function(k)
           if (length(median_counts) > 1) stats::median(median_counts[-k]) else NA_real_, numeric(1)),
         ratio_to_core = median_counts / core_other_median) |>
  ungroup() |>
  mutate(fov_flag = case_when(
    n_cells < FOV_MIN_CELLS ~ "few_cells",
    !is.na(ratio_to_core) & ratio_to_core < FOV_MIN_RATIO_TO_CORE ~ "low_vs_core",
    TRUE ~ "ok"))
save_table(as.data.frame(fov_core), "02_fov_qc_by_core")
bad_fov <- fov_core[fov_core$fov_flag == "low_vs_core", ]
if (nrow(bad_fov) > 0) {
  message("FOVs with < ", FOV_MIN_RATIO_TO_CORE, " x the transcripts per cell of the other FOVs of their core ",
          "(look at them in 01_fov_qc_map.png; to drop one, remove its row from config/fov_to_core.csv or ",
          "exclude the core):")
  print(as.data.frame(bad_fov[, c("slide_id", "core_id", "fov", "median_counts", "core_other_median", "ratio_to_core")]))
}
cmx$fov_flag <- fov_core$fov_flag[match(paste(cmx$slide_id, cmx$fov), paste(fov_core$slide_id, fov_core$fov))]

missing_cores <- setdiff(tma_map$core_id[tma_map$slide_id %in% slides$slide_id], core_qc$core_id)
if (length(missing_cores) > 0) {
  message("Cores in tma_map.csv with NO cells (lost during sectioning?): ",
          paste(missing_cores, collapse = ", "))
}

p_core_qc <- ggplot(core_qc, aes(reorder(core_id, median_counts), median_counts, fill = organ)) +
  geom_col() +
  geom_hline(yintercept = CORE_MIN_MEDIAN_COUNTS, linetype = 2) +
  facet_grid(~ slide_id, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
  labs(x = "core", y = "median transcripts / cell (QC-pass cells)",
       title = "Core quality -- differences here confound every organ comparison in Aim 1/2")
save_fig(p_core_qc, "02_core_qc_median_counts", width = 16, height = 6)

## -----------------------------------------------------------------------
## Filter: QC-pass cells in cores that are in the map and included
## -----------------------------------------------------------------------
cmx$core_flag <- core_qc$core_flag[match(cmx$core_id, core_qc$core_id)]
keep <- colnames(cmx)[!is.na(cmx$core_id) & cmx$qc_pass & cmx$include %in% TRUE]
message("Keeping ", length(keep), " of ", ncol(cmx), " cells.")
cmx <- subset(cmx, cells = keep)

saveRDS(cmx, obj_path("cores"))
message("Saved -> ", obj_path("cores"))
