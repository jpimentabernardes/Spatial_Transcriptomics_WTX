## =============================================================================
## tests/make_synthetic_data.R
##
## Writes a small synthetic data set in the real file formats of both platforms,
## so the whole pipeline (00 - 05 and the Xenium <-> CosMx comparison, 10) can be
## run and checked before real data exists:
##
##   <out>/cosmx_project/data/cosmx/<slide>/   CosMx flat files (exprMat, metadata,
##                                             fov_positions, polygons, tx for a few FOVs)
##   <out>/xenium_project/data/xenium/<slide>/ Xenium-style cell_feature_matrix/ + cells.csv.gz
##   <out>/config_cosmx/, <out>/config_xenium/  slides.csv, tma_map.csv for each pipeline
##   <out>/truth_*.csv                          the simulated truth (cell type per cell, gene efficiencies)
##
## Two TMAs (TMA_ORGAN: 9 cores, TMA_TUMOR: 6 cores), consecutive sections:
## the two platforms see the same tissue structure but not the same cells.
## What is simulated ON PURPOSE, so that 10 has something to find:
##   * gene-specific detection efficiency, different on each platform
##   * CosMx sees ~1,500 genes ("whole transcriptome"), Xenium a 250-gene panel
##   * optical crowding on CosMx: efficiency falls with the RNA density of the
##     cell (CROWD_K), nothing of the kind on Xenium
##   * more transcript spill-over between neighbouring cells on CosMx
##   * the CosMx slide of TMA_TUMOR is mounted rotated by 180 degrees, every core
##     is shifted / rotated a little, and TMA_ORGAN has cores only half covered by FOVs
##
##   Rscript tests/make_synthetic_data.R /path/to/out
## =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1) args[1] else file.path(tempdir(), "synthetic_tma")
set.seed(1)

N_GENES   <- 1500      # "whole transcriptome"
N_PANEL   <- 250       # Xenium panel (subset of the above)
PX        <- 0.12028   # CosMx pixel (um)
FOV_PX    <- 4256      # CosMx FOV edge (px) = 512 um
CORE_D    <- 900       # core diameter (um)
PITCH     <- 1600      # core pitch (um)
CELLS_PER_CORE <- 1400
CROWD_K   <- as.numeric(Sys.getenv("SYN_CROWD_K", "60"))   # CosMx efficiency = 1 / (1 + density / CROWD_K), density in true
                       # transcripts / um2 of cell; SYN_CROWD_K=Inf = no crowding (negative control for 10, section 4)

## ---- genes -------------------------------------------------------------------
lineage_markers <- list(
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CDH1", "CLDN4", "ELF3"),
  Immune = c("PTPRC", "CD52", "CORO1A", "LCP1"),
  Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "PLVAP", "CLDN5"),
  Smooth_muscle_pericyte = c("ACTA2", "DES", "MYH11", "CNN1", "RGS5", "PDGFRB"))
type_markers <- list(
  Colonocyte = c("CA1", "CA2", "SLC26A3", "AQP8", "FABP1"), Goblet = c("MUC2", "TFF3", "SPINK4", "CLCA1"),
  T_cell = c("CD3E", "CD3D", "CD2", "TRAC", "CD8A", "CD4", "FOXP3"), B_cell = c("MS4A1", "CD19", "CD79A", "BANK1"),
  Macrophage = c("CD68", "CD163", "C1QA", "C1QB", "CSF1R"))
types <- data.frame(
  cell_type = c("Colonocyte", "Goblet", "T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Smooth_muscle"),
  lineage = c("Epithelial", "Epithelial", "Immune", "Immune", "Immune", "Fibroblast", "Endothelial", "Smooth_muscle_pericyte"),
  area = c(110, 100, 40, 40, 80, 90, 70, 120),          # um2
  content = c(5200, 4200, 1500, 1700, 2600, 1800, 1500, 2000))   # true transcripts per cell (all genes)
named <- unique(c(unlist(lineage_markers), unlist(type_markers)))
genes <- c(named, sprintf("GENE%04d", seq_len(N_GENES - length(named))))
## a multi-gene CosMx probe and a case difference, to exercise the name matching of 10
cosmx_names <- genes; cosmx_names[genes == "GENE0001"] <- "GENE0001/2"
panel <- c(named, sample(setdiff(genes, named), N_PANEL - length(named)))

## expression profile per cell type (fractions summing to 1)
base <- stats::rlnorm(N_GENES, 0, 1.6)
prof <- sapply(seq_len(nrow(types)), function(k) {
  p <- base * stats::rlnorm(N_GENES, 0, 0.5)
  spec <- sample(setdiff(seq_len(N_GENES), seq_along(named)), 40)        # type-specific filler genes
  p[spec] <- p[spec] * 25
  lin <- lineage_markers[[types$lineage[k]]]; typ <- type_markers[[types$cell_type[k]]]
  p[genes %in% named] <- p[genes %in% named] * 0.01                   # markers are off ...
  p[genes %in% lin] <- stats::quantile(base, 0.97)                    # ... except in their lineage
  if (!is.null(typ)) p[genes %in% typ] <- stats::quantile(base, 0.98)
  p / sum(p)
})
dimnames(prof) <- list(genes, types$cell_type)
eff_x <- stats::setNames(stats::rlnorm(N_GENES, log(0.30), 0.5), genes)   # Xenium detection efficiency per gene
eff_c <- stats::setNames(stats::rlnorm(N_GENES, log(0.25), 0.5), genes)   # CosMx
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fwrite(data.frame(gene = genes, cosmx_name = cosmx_names, on_xenium_panel = genes %in% panel, eff_xenium = eff_x,
                  eff_cosmx = eff_c), file.path(OUT, "truth_genes.csv"))

## ---- TMA layout --------------------------------------------------------------
tma <- rbind(
  data.frame(slide_id = "TMA_ORGAN", expand.grid(core_col = 1:3, core_row = c("A", "B", "C")),
             organ = rep(c("colon", "lung", "tonsil"), each = 3), patient_id = sprintf("P%02d", 1:9),
             tissue_type = "normal", diagnosis = "none", location = rep(c("proximal", "distal", "proximal"), 3)),
  data.frame(slide_id = "TMA_TUMOR", expand.grid(core_col = 1:3, core_row = c("A", "B")),
             organ = rep(c("colon", "lung"), each = 3), patient_id = sprintf("P%02d", c(1:3, 4:6)),
             tissue_type = "tumor", diagnosis = "adenocarcinoma",
             location = rep(c("tumor_center", "invasive_margin", "tumor_center"), 2)))
tma$core_id <- paste0(ifelse(tma$slide_id == "TMA_ORGAN", "ORG_", "TUM_"), tma$core_row, tma$core_col)
tma$include <- TRUE; tma$notes <- "synthetic"
tma <- tma[, c("slide_id", "core_row", "core_col", "core_id", "patient_id", "organ", "tissue_type", "diagnosis",
               "location", "include", "notes")]
## type mix per organ: stroma / epithelium / immune
organ_mix <- list(colon = c(0.30, 0.20, 0.08, 0.06, 0.08, 0.14, 0.06, 0.08),
                  lung = c(0.22, 0.02, 0.10, 0.04, 0.16, 0.22, 0.16, 0.08),
                  tonsil = c(0.08, 0.01, 0.34, 0.36, 0.08, 0.07, 0.04, 0.02))

## ---- truth: cells of one core, with tissue structure ---------------------------
make_core <- function(organ) {
  n <- round(CELLS_PER_CORE * stats::runif(1, 0.8, 1.2))
  r <- CORE_D / 2 * sqrt(stats::runif(n)); a <- stats::runif(n, 0, 2 * pi)
  x <- r * cos(a); y <- r * sin(a)
  mix <- organ_mix[[organ]]
  ## glands (epithelium) and one lymphoid aggregate (T / B) give the core a density pattern
  glands <- cbind(stats::runif(6, -300, 300), stats::runif(6, -300, 300))
  agg <- c(stats::runif(1, -250, 250), stats::runif(1, -250, 250))
  in_gland <- apply(sapply(seq_len(nrow(glands)), function(i) (x - glands[i, 1])^2 + (y - glands[i, 2])^2 < 90^2), 1, any)
  in_agg <- (x - agg[1])^2 + (y - agg[2])^2 < 110^2
  w <- matrix(mix, n, 8, byrow = TRUE)
  w[in_gland, 1:2] <- w[in_gland, 1:2] * 12
  w[!in_gland, 1:2] <- w[!in_gland, 1:2] * 0.15
  w[in_agg, 3:4] <- w[in_agg, 3:4] * 15
  ct <- apply(w, 1, function(p) sample(8, 1, prob = p))
  data.frame(x = x, y = y, type = ct, area = types$area[ct] * stats::rlnorm(n, 0, 0.25),
             content = types$content[ct] * stats::rlnorm(n, 0, 0.35))
}

## counts for one section. spill = share of a cell's signal that comes from its
## neighbourhood (segmentation spill-over); crowd = TRUE applies the density effect.
simulate_counts <- function(cells, gene_set, eff, spill, crowd) {
  g <- match(gene_set, genes)
  dens <- cells$content / cells$area                               # true transcripts / um2
  cell_eff <- if (crowd) 1 / (1 + dens / CROWD_K) else rep(1, nrow(cells))
  core_mean <- rowMeans(prof[g, cells$type, drop = FALSE])
  lam <- (prof[g, cells$type, drop = FALSE] * (1 - spill) + core_mean * spill) *
    rep(cells$content * cell_eff, each = length(g)) * eff[g]
  m <- matrix(stats::rpois(length(lam), lam), nrow = length(g), dimnames = list(gene_set, NULL))
  methods::as(methods::as(m, "CsparseMatrix"), "generalMatrix")
}

rot <- function(xy, deg) { a <- deg * pi / 180; cbind(xy[, 1] * cos(a) - xy[, 2] * sin(a), xy[, 1] * sin(a) + xy[, 2] * cos(a)) }

truth_x <- list(); truth_c <- list()
for (sid in unique(tma$slide_id)) {
  map <- tma[tma$slide_id == sid, ]
  xen_cells <- list(); cos_cells <- list()
  for (k in seq_len(nrow(map))) {
    core <- make_core(map$organ[k])
    cx <- 1000 + (map$core_col[k] - 1) * PITCH; cy <- 1000 + (match(map$core_row[k], LETTERS) - 1) * PITCH
    ## Xenium section
    sx <- core[sample(nrow(core), round(0.9 * nrow(core))), ]
    sx$x <- sx$x + stats::rnorm(nrow(sx), 0, 1) + cx; sx$y <- sx$y + stats::rnorm(nrow(sx), 0, 1) + cy
    sx$core_id <- map$core_id[k]
    xen_cells[[k]] <- sx
    ## CosMx section, 5 um further: other cells of the same structures; the core sits a little differently
    sc <- core[sample(nrow(core), round(0.9 * nrow(core))), ]
    xy <- rot(cbind(sc$x + stats::rnorm(nrow(sc), 0, 4), sc$y + stats::rnorm(nrow(sc), 0, 4)), stats::runif(1, -6, 6))
    sc$x <- xy[, 1] + cx + stats::runif(1, -40, 40); sc$y <- xy[, 2] + cy + stats::runif(1, -40, 40)
    sc$core_id <- map$core_id[k]
    sc$core_cx <- cx; sc$core_cy <- cy
    cos_cells[[k]] <- sc
  }
  xen_cells <- do.call(rbind, xen_cells); cos_cells <- do.call(rbind, cos_cells)

  ## ---------------- Xenium outs (10x-style matrix + cells table) ----------------
  xd <- file.path(OUT, "xenium_project", "data", "xenium", sid)
  dir.create(file.path(xd, "cell_feature_matrix"), recursive = TRUE, showWarnings = FALSE)
  cnt <- simulate_counts(xen_cells, panel, eff_x, spill = 0.03, crowd = FALSE)
  n_neg <- 20; n_cw <- 30
  neg <- Matrix(stats::rpois(n_neg * ncol(cnt), 0.004), nrow = n_neg, sparse = TRUE)
  cw <- Matrix(stats::rpois(n_cw * ncol(cnt), 0.001), nrow = n_cw, sparse = TRUE)
  full <- rbind(cnt, neg, cw)
  feat <- data.frame(id = c(panel, sprintf("NegControlProbe_%05d", 1:n_neg), sprintf("NegControlCodeword_%04d", 1:n_cw)))
  feat$name <- feat$id
  feat$type <- c(rep("Gene Expression", length(panel)), rep("Negative Control Probe", n_neg),
                 rep("Negative Control Codeword", n_cw))
  ids <- sprintf("cell%06d", seq_len(ncol(full)))
  writeMM(methods::as(full, "TsparseMatrix"), file.path(xd, "cell_feature_matrix", "matrix.mtx"))
  R.utils::gzip(file.path(xd, "cell_feature_matrix", "matrix.mtx"), overwrite = TRUE)
  fwrite(feat, file.path(xd, "cell_feature_matrix", "features.tsv.gz"), sep = "\t", col.names = FALSE)
  fwrite(data.frame(ids), file.path(xd, "cell_feature_matrix", "barcodes.tsv.gz"), col.names = FALSE)
  fwrite(data.frame(cell_id = ids, x_centroid = xen_cells$x, y_centroid = xen_cells$y,
                    cell_area = xen_cells$area * 0.9, nucleus_area = xen_cells$area * 0.35),
         file.path(xd, "cells.csv.gz"))
  truth_x[[sid]] <- data.frame(cell = paste0(sid, "_", ids), slide_id = sid, core_id = xen_cells$core_id,
                               cell_type = types$cell_type[xen_cells$type], lineage = types$lineage[xen_cells$type])

  ## ---------------- CosMx flat files ---------------------------------------------
  ## slide frame: TMA_TUMOR is mounted rotated by 180 degrees; both slides sit at
  ## another place on the stage; global pixels, y UPWARDS
  W <- max(cos_cells$core_cx) + 1000; H <- max(cos_cells$core_cy) + 1000
  flip180 <- sid == "TMA_TUMOR"
  to_slide <- function(x, y) {
    if (flip180) { x <- W - x; y <- H - y }
    list(x_px = (x + 3000) / PX, y_px = (H - y + 2000) / PX)
  }
  ## FOVs: 2 x 2 per core; on TMA_ORGAN column 3 only the upper 1 x 2 (half-covered cores)
  fov_tab <- list(); cos_cells$fov <- NA_integer_; f <- 0L
  half <- FOV_PX * PX
  for (k in seq_len(nrow(map))) {
    cx <- cos_cells$core_cx[cos_cells$core_id == map$core_id[k]][1]; cy <- cos_cells$core_cy[cos_cells$core_id == map$core_id[k]][1]
    rows <- if (sid == "TMA_ORGAN" && map$core_col[k] == 3) 1 else 1:2
    for (iy in rows) for (ix in 1:2) {
      f <- f + 1L
      x0 <- cx + (ix - 2) * half; y0 <- cy + (iy - 2) * half
      inside <- cos_cells$core_id == map$core_id[k] & cos_cells$x >= x0 & cos_cells$x < x0 + half &
        cos_cells$y >= y0 & cos_cells$y < y0 + half
      cos_cells$fov[inside] <- f
      c0 <- to_slide(if (flip180) x0 + half else x0, if (flip180) y0 else y0 + half)   # FOV corner with smallest px
      fov_tab[[f]] <- data.frame(FOV = f, x_global_px = c0$x_px, y_global_px = c0$y_px,
                                 x_global_mm = c0$x_px * PX / 1000, y_global_mm = c0$y_px * PX / 1000)
    }
  }
  cos_cells <- cos_cells[!is.na(cos_cells$fov), ]
  cos_cells <- cos_cells[order(cos_cells$fov), ]
  cos_cells$cell_ID <- stats::ave(cos_cells$fov, cos_cells$fov, FUN = seq_along)
  fov_tab <- do.call(rbind, fov_tab)
  gp <- to_slide(cos_cells$x, cos_cells$y)
  cnt <- simulate_counts(cos_cells, genes, eff_c, spill = 0.07, crowd = TRUE)
  rownames(cnt) <- cosmx_names
  n_neg <- 20; n_sys <- 50
  neg <- matrix(stats::rpois(n_neg * ncol(cnt), 0.012), nrow = n_neg, dimnames = list(paste0("Negative", 1:n_neg), NULL))
  sys <- matrix(stats::rpois(n_sys * ncol(cnt), 0.002), nrow = n_sys, dimnames = list(paste0("SystemControl", 1:n_sys), NULL))
  ex <- data.frame(fov = cos_cells$fov, cell_ID = cos_cells$cell_ID, as.matrix(t(cnt)), t(neg), t(sys), check.names = FALSE)
  ## cell_ID 0 = transcripts of the FOV outside every cell
  zero <- do.call(rbind, lapply(unique(cos_cells$fov), function(ff) {
    tot <- colSums(ex[ex$fov == ff, -(1:2), drop = FALSE])
    data.frame(fov = ff, cell_ID = 0L, t(stats::rpois(length(tot), tot * 0.08)), check.names = FALSE)
  }))
  colnames(zero) <- colnames(ex)
  ex <- rbind(zero, ex); ex <- ex[order(ex$fov, ex$cell_ID), ]
  cd <- file.path(OUT, "cosmx_project", "data", "cosmx", sid, "flatFiles", sid)   # nested as AtoMx exports it
  dir.create(cd, recursive = TRUE, showWarnings = FALSE)
  fwrite(ex, file.path(cd, paste0(sid, "_exprMat_file.csv.gz")))
  fov_px <- fov_tab[match(cos_cells$fov, fov_tab$FOV), ]
  area_px <- round(cos_cells$area * 1.1 / PX^2)
  border <- pmin(gp$x_px - fov_px$x_global_px, gp$y_px - fov_px$y_global_px,
                 fov_px$x_global_px + FOV_PX - gp$x_px, fov_px$y_global_px + FOV_PX - gp$y_px) < 40
  md <- data.frame(
    fov = cos_cells$fov, Area = area_px, AspectRatio = round(stats::runif(nrow(cos_cells), 0.6, 1.5), 2),
    CenterX_local_px = round(gp$x_px - fov_px$x_global_px), CenterY_local_px = round(gp$y_px - fov_px$y_global_px),
    Width = round(sqrt(area_px)), Height = round(sqrt(area_px)),
    Mean.PanCK = round(ifelse(types$lineage[cos_cells$type] == "Epithelial", 8000, 300) * stats::rlnorm(nrow(cos_cells), 0, 0.4)),
    Max.PanCK = 20000,
    Mean.CD45 = round(ifelse(types$lineage[cos_cells$type] == "Immune", 5000, 250) * stats::rlnorm(nrow(cos_cells), 0, 0.4)),
    Max.CD45 = 15000, Mean.DAPI = round(4000 * stats::rlnorm(nrow(cos_cells), 0, 0.3)), Max.DAPI = 30000,
    SplitRatioToLocal = ifelse(border, round(stats::runif(nrow(cos_cells), 0.1, 1.4), 2), 0),
    NucArea = round(area_px * 0.35), cell_id = sprintf("c_1_%d_%d", cos_cells$fov, cos_cells$cell_ID),
    Panel = "Synthetic WTX", slide_ID = 1L,
    CenterX_global_px = round(gp$x_px, 1), CenterY_global_px = round(gp$y_px, 1), cell_ID = cos_cells$cell_ID,
    cell = sprintf("c_1_%d_%d", cos_cells$fov, cos_cells$cell_ID))
  if (sid == "TMA_ORGAN") md$Area.um2 <- round(area_px * PX^2, 2)     # one slide with, one without Area.um2
  fwrite(md, file.path(cd, paste0(sid, "_metadata_file.csv.gz")))
  fwrite(fov_tab, file.path(cd, paste0(sid, "_fov_positions_file.csv.gz")))
  ## polygons: an octagon per cell
  ang <- seq(0, 2 * pi, length.out = 9)[-9]
  rr <- sqrt(area_px / pi)
  poly <- data.frame(fov = rep(cos_cells$fov, each = 8), cellID = rep(cos_cells$cell_ID, each = 8),
                     cell = rep(md$cell, each = 8),
                     x_global_px = round(rep(gp$x_px, each = 8) + rep(rr, each = 8) * cos(ang)),
                     y_global_px = round(rep(gp$y_px, each = 8) + rep(rr, each = 8) * sin(ang)))
  poly$x_local_px <- poly$x_global_px - rep(fov_px$x_global_px, each = 8)
  poly$y_local_px <- poly$y_global_px - rep(fov_px$y_global_px, each = 8)
  fwrite(poly, file.path(cd, paste0(sid, "-polygons.csv.gz")))
  ## transcripts of the first 3 FOVs (the real file has every transcript of the slide)
  sel <- which(cos_cells$fov <= 3)
  tr <- summary(cnt[, sel, drop = FALSE])
  tr <- tr[rep(seq_len(nrow(tr)), tr$x), ]
  cell_i <- sel[tr$j]
  tx <- data.frame(fov = cos_cells$fov[cell_i], cell_ID = cos_cells$cell_ID[cell_i], cell = md$cell[cell_i],
                   x_global_px = gp$x_px[cell_i] + stats::rnorm(nrow(tr), 0, 20),
                   y_global_px = gp$y_px[cell_i] + stats::rnorm(nrow(tr), 0, 20), z = sample(0:8, nrow(tr), TRUE),
                   target = cosmx_names[tr$i],
                   CellComp = sample(c("Nuclear", "Cytoplasm", "Membrane"), nrow(tr), TRUE, prob = c(0.35, 0.55, 0.10)))
  out <- tx[sample(nrow(tx), round(0.08 * nrow(tx))), ]                 # transcripts outside cells
  out$cell_ID <- 0L; out$cell <- sprintf("c_1_%d_0", out$fov); out$CellComp <- "None"
  tx <- rbind(tx, out)
  tx$x_local_px <- tx$x_global_px - fov_tab$x_global_px[tx$fov]; tx$y_local_px <- tx$y_global_px - fov_tab$y_global_px[tx$fov]
  fwrite(tx, file.path(cd, paste0(sid, "_tx_file.csv.gz")))
  truth_c[[sid]] <- data.frame(cell = sprintf("%s_%d_%d", sid, cos_cells$fov, cos_cells$cell_ID), slide_id = sid,
                               core_id = cos_cells$core_id, cell_type = types$cell_type[cos_cells$type],
                               lineage = types$lineage[cos_cells$type], true_density = cos_cells$content / cos_cells$area)
  message(sid, ": ", nrow(xen_cells), " Xenium cells, ", nrow(cos_cells), " CosMx cells in ", f, " FOVs")
}
fwrite(do.call(rbind, truth_x), file.path(OUT, "truth_cells_xenium.csv"))
fwrite(do.call(rbind, truth_c), file.path(OUT, "truth_cells_cosmx.csv"))

## ---- config for both pipelines -------------------------------------------------
for (k in c("config_cosmx", "config_xenium")) dir.create(file.path(OUT, k), showWarnings = FALSE)
fwrite(tma, file.path(OUT, "config_cosmx", "tma_map.csv")); fwrite(tma, file.path(OUT, "config_xenium", "tma_map.csv"))
fwrite(data.frame(slide_id = c("TMA_ORGAN", "TMA_TUMOR"), tma_type = c("multi_organ", "multi_tumor"),
                  cosmx_dir = file.path("data", "cosmx", c("TMA_ORGAN", "TMA_TUMOR")), cosmx_prefix = "", panel = "WTX",
                  row_flip = c(FALSE, TRUE), col_flip = c(FALSE, TRUE),          # TMA_TUMOR is rotated by 180 degrees
                  xenium_slide_id = c("TMA_ORGAN", "TMA_TUMOR"), he_image = "", he_pixel_size_um = "", notes = "synthetic"),
       file.path(OUT, "config_cosmx", "slides.csv"))
fwrite(data.frame(slide_id = c("TMA_ORGAN", "TMA_TUMOR"), tma_type = c("multi_organ", "multi_tumor"),
                  xenium_dir = file.path("data", "xenium", c("TMA_ORGAN", "TMA_TUMOR")), panel = "synthetic",
                  row_flip = FALSE, col_flip = FALSE, notes = "synthetic"),
       file.path(OUT, "config_xenium", "slides.csv"))
message("Synthetic data written to ", OUT)
