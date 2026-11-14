## =============================================================================
## R/compare_helpers.R
##
## Functions for 10_xenium_vs_cosmx.R: gene-name matching between the two
## panels, registration of consecutive sections core by core (cell centroids
## only, no images), label transfer on the shared genes, and small statistics.
## Sourced after config/config.R and R/helpers.R. Only functions.
## =============================================================================

## -----------------------------------------------------------------------
## Gene names
## -----------------------------------------------------------------------
## Names a multi-gene CosMx probe may stand for: "HLA-DQB1/2" -> HLA-DQB1,
## HLA-DQB2; "CCL3/L1" -> CCL3, CCL3L1; "XCL1/2" -> XCL1, XCL2. The short forms
## are ambiguous, so several spellings are tried; only names that exist on the
## Xenium panel are ever used, and such pairs are listed, not silently trusted.
expand_multi_gene <- function(name) {
  parts <- strsplit(name, "/", fixed = TRUE)[[1]]
  first <- parts[1]
  cand <- first
  for (p in parts[-1]) {
    cand <- c(cand, p, paste0(first, p),
              if (nchar(p) < nchar(first)) paste0(substr(first, 1, nchar(first) - nchar(p)), p))
  }
  unique(cand)
}

## One row per Xenium gene that has a CosMx partner.
##   match   exact | case (same name, other case) | multi_gene_probe | manual
##   use     FALSE for multi-gene probes: the CosMx probe counts several genes,
##           the Xenium probe one -- not comparable 1:1 unless you decide so in
##           config/gene_name_map.csv (xenium_gene, cosmx_gene), which always wins.
match_genes <- function(xen_genes, cos_genes, map_csv = GENE_NAME_MAP_CSV) {
  out <- data.frame(xenium_gene = xen_genes, cosmx_gene = NA_character_, match = NA_character_, stringsAsFactors = FALSE)
  ex <- xen_genes %in% cos_genes
  out$cosmx_gene[ex] <- xen_genes[ex]; out$match[ex] <- "exact"
  k <- match(toupper(xen_genes), toupper(cos_genes))
  ci <- !ex & !is.na(k)
  out$cosmx_gene[ci] <- cos_genes[k[ci]]; out$match[ci] <- "case"
  multi <- grep("/", cos_genes, fixed = TRUE, value = TRUE)
  if (length(multi) > 0) {
    alias <- do.call(rbind, lapply(multi, function(m) data.frame(alias = toupper(expand_multi_gene(m)), cosmx_gene = m)))
    k <- match(toupper(xen_genes), alias$alias)
    mg <- is.na(out$cosmx_gene) & !is.na(k)
    out$cosmx_gene[mg] <- alias$cosmx_gene[k[mg]]; out$match[mg] <- "multi_gene_probe"
  }
  if (file.exists(map_csv)) {
    man <- read_config_csv(map_csv, stringsAsFactors = FALSE)
    stopifnot("gene_name_map.csv needs xenium_gene and cosmx_gene" = all(c("xenium_gene", "cosmx_gene") %in% colnames(man)))
    man$xenium_gene <- match_case(as_seurat_features(man$xenium_gene), xen_genes)
    man$cosmx_gene <- match_case(as_seurat_features(man$cosmx_gene), cos_genes)
    bad <- man[!man$xenium_gene %in% xen_genes | !man$cosmx_gene %in% cos_genes, ]
    if (nrow(bad) > 0) message("gene_name_map.csv rows with a name on neither panel (ignored): ",
                               paste(bad$xenium_gene, bad$cosmx_gene, sep = "=", collapse = ", "))
    man <- man[man$xenium_gene %in% xen_genes & man$cosmx_gene %in% cos_genes, ]
    k <- match(man$xenium_gene, out$xenium_gene)
    out$cosmx_gene[k] <- man$cosmx_gene; out$match[k] <- "manual"
  }
  out$use <- !is.na(out$cosmx_gene) & out$match != "multi_gene_probe"
  ## a CosMx target claimed by two Xenium genes cannot be compared 1:1 with both
  dup <- out$cosmx_gene[out$use][duplicated(out$cosmx_gene[out$use])]
  out$use[out$cosmx_gene %in% dup & out$match != "exact"] <- FALSE
  out
}

## -----------------------------------------------------------------------
## 2-D transforms (3x3 homogeneous matrices acting on [x, y, 1])
## -----------------------------------------------------------------------
apply_tf <- function(M, xy) {
  xy <- as.matrix(xy)
  cbind(xy[, 1] * M[1, 1] + xy[, 2] * M[1, 2] + M[1, 3], xy[, 1] * M[2, 1] + xy[, 2] * M[2, 2] + M[2, 3])
}

## rotation + uniform scale + translation, mirror image allowed (Umeyama)
fit_similarity <- function(src, dst) {
  mu_s <- colMeans(src); mu_d <- colMeans(dst)
  s <- sweep(src, 2, mu_s); d <- sweep(dst, 2, mu_d)
  sv <- svd(crossprod(d, s) / nrow(src))
  R <- sv$u %*% t(sv$v)
  scale <- sum(sv$d) / sum(colMeans(s^2))
  M <- diag(3); M[1:2, 1:2] <- scale * R; M[1:2, 3] <- mu_d - scale * R %*% mu_s
  M
}

## rotation + translation only (no mirror, no scale)
fit_rigid <- function(src, dst) {
  mu_s <- colMeans(src); mu_d <- colMeans(dst)
  sv <- svd(crossprod(sweep(dst, 2, mu_d), sweep(src, 2, mu_s)))
  D <- diag(c(1, sign(det(sv$u %*% t(sv$v)))))
  R <- sv$u %*% D %*% t(sv$v)
  M <- diag(3); M[1:2, 1:2] <- R; M[1:2, 3] <- mu_d - R %*% mu_s
  M
}

## similarity fit with pairs further than 3 MAD (at least 150 um) rejected
robust_similarity <- function(src, dst, n_iter = 5) {
  keep <- rep(TRUE, nrow(src))
  for (i in seq_len(n_iter)) {
    M <- fit_similarity(src[keep, , drop = FALSE], dst[keep, , drop = FALSE])
    res <- sqrt(rowSums((apply_tf(M, src) - dst)^2))
    thr <- max(150, stats::median(res[keep]) + 3 * stats::mad(res[keep]))
    new_keep <- res <= thr
    if (sum(new_keep) < 3 || identical(new_keep, keep)) break
    keep <- new_keep
  }
  list(M = M, residual = res, keep = keep)
}

rot_deg <- function(M) atan2(M[2, 1], M[1, 1]) * 180 / pi
rot_about <- function(deg, cx, cy, dx = 0, dy = 0) {
  a <- deg * pi / 180
  R <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
  M <- diag(3); M[1:2, 1:2] <- R; M[1:2, 3] <- c(cx, cy) - R %*% c(cx, cy) + c(dx, dy)
  M
}

## -----------------------------------------------------------------------
## Registration of one core: CosMx cells (moving) onto Xenium cells (fixed)
## -----------------------------------------------------------------------
## Correlation of the two cell-density maps (bins of `bin` um), over the part of
## the grid the MOVING cells cover -- on CosMx the FOVs may cover only part of a
## core, and the uncovered part must not count against the alignment.
density_ncc <- function(mov, fix, bin = 40) {
  x0 <- min(fix[, 1]) - 4 * bin; y0 <- min(fix[, 2]) - 4 * bin
  nx <- ceiling((max(fix[, 1]) - x0) / bin) + 4; ny <- ceiling((max(fix[, 2]) - y0) / bin) + 4
  to_img <- function(p) {
    ix <- floor((p[, 1] - x0) / bin) + 1; iy <- floor((p[, 2] - y0) / bin) + 1
    ok <- ix >= 1 & ix <= nx & iy >= 1 & iy <= ny
    list(img = matrix(tabulate(ix[ok] + (iy[ok] - 1) * nx, nbins = nx * ny), nx, ny), ix = ix[ok], iy = iy[ok])
  }
  f <- to_img(fix); m <- to_img(mov)
  if (length(m$ix) < 20) return(NA_real_)
  rx <- range(m$ix); ry <- range(m$iy)
  a <- as.vector(m$img[rx[1]:rx[2], ry[1]:ry[2]]); b <- as.vector(f$img[rx[1]:rx[2], ry[1]:ry[2]])
  if (stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
  stats::cor(a, b)
}

## Coarse search (rotation x shift around the core centre, scored by density_ncc),
## then trimmed ICP on the centroids: every moving cell is pulled towards the
## nearest fixed cell -- of the same lineage where both have lineages. Consecutive
## sections do not contain the same cells, so the ICP aligns tissue structure
## (glands, aggregates, edges), with the 30% worst pairs ignored in every step.
## Returns the 3x3 correction (applied AFTER the slide-level transform).
register_core <- function(mov, fix, mov_lab = NULL, fix_lab = NULL, bin = 40,
                          rotations = seq(-20, 20, by = 2), shifts = seq(-200, 200, by = 25),
                          icp_dist = 40, icp_iter = 30) {
  cx <- stats::median(mov[, 1]); cy <- stats::median(mov[, 2])
  best <- list(score = density_ncc(mov, fix, bin), M = diag(3))
  ncc_start <- best$score
  for (r in rotations) {
    rotated <- apply_tf(rot_about(r, cx, cy), mov)
    for (dx in shifts) for (dy in shifts) {
      sc <- density_ncc(cbind(rotated[, 1] + dx, rotated[, 2] + dy), fix, bin)
      if (!is.na(sc) && (is.na(best$score) || sc > best$score)) best <- list(score = sc, M = rot_about(r, cx, cy, dx, dy))
    }
  }
  Tm <- best$M
  cur <- apply_tf(Tm, mov)
  use_lab <- !is.null(mov_lab) && !is.null(fix_lab)
  groups <- if (use_lab) intersect(unique(mov_lab), unique(fix_lab)) else "all"
  for (it in seq_len(icp_iter)) {
    src <- list(); dst <- list(); dd <- list()
    for (g in groups) {
      mi <- if (use_lab) which(mov_lab == g) else seq_len(nrow(cur))
      fi <- if (use_lab) which(fix_lab == g) else seq_len(nrow(fix))
      if (length(mi) < 5 || length(fi) < 5) next
      nn <- RANN::nn2(fix[fi, , drop = FALSE], cur[mi, , drop = FALSE], k = 1)
      ok <- nn$nn.dists[, 1] <= icp_dist
      src[[g]] <- cur[mi[ok], , drop = FALSE]; dst[[g]] <- fix[fi[nn$nn.idx[ok, 1]], , drop = FALSE]; dd[[g]] <- nn$nn.dists[ok, 1]
    }
    if (length(src) == 0) break
    src <- do.call(rbind, src); dst <- do.call(rbind, dst); dd <- unlist(dd)
    if (nrow(src) < 30) break
    keep <- dd <= stats::quantile(dd, 0.7)
    step <- fit_rigid(src[keep, , drop = FALSE], dst[keep, , drop = FALSE])
    cur <- apply_tf(step, cur)
    Tm <- step %*% Tm
    if (sqrt(sum(step[1:2, 3]^2)) < 0.05 && abs(rot_deg(step)) < 0.01) break
  }
  ncc_icp <- density_ncc(cur, fix, bin)
  ## the ICP is a refinement: keep it only if it does not make the density match worse
  if (is.na(ncc_icp) || (!is.na(best$score) && ncc_icp < best$score - 0.02)) { Tm <- best$M; ncc_icp <- best$score }
  centre <- matrix(c(cx, cy), 1)
  list(M = Tm, ncc_start = ncc_start, ncc = ncc_icp, rotation_deg = rot_deg(Tm),
       shift_um = sqrt(sum((apply_tf(Tm, centre) - centre)^2)))
}

## Do neighbouring cells of the two sections have the same lineage more often
## than chance? lift = observed agreement / agreement expected from the lineage
## frequencies. ~1 = not aligned (or the annotations disagree); >> 1 = aligned.
lineage_lift <- function(mov, mov_lab, fix, fix_lab, radius = 20) {
  nn <- RANN::nn2(fix, mov, k = 1)
  ok <- nn$nn.dists[, 1] <= radius & !is.na(mov_lab)
  if (sum(ok) < 20) return(c(agreement = NA_real_, lift = NA_real_))
  f <- fix_lab[nn$nn.idx[ok, 1]]
  agree <- mean(mov_lab[ok] == f, na.rm = TRUE)
  labs <- union(mov_lab, fix_lab)
  expected <- sum(vapply(labs, function(l) mean(mov_lab[ok] == l, na.rm = TRUE) * mean(fix_lab == l, na.rm = TRUE), numeric(1)))
  c(agreement = agree, lift = if (expected > 0) agree / expected else NA_real_)
}

## -----------------------------------------------------------------------
## Label transfer on shared genes
## -----------------------------------------------------------------------
## Reference = one platform with labels, query = the other. Each query cell
## gets the label whose reference profile it resembles most (cosine similarity),
## using the shared genes only. Gene-specific detection efficiency differs
## between the platforms, so every gene is standardised WITHIN its platform
## first (z-score across that platform's cells): what is compared is "high or
## low for this gene on this platform", not the raw counts.
## Assumes the two platforms saw a similar mix of cells (same cores) -- true here.
## ref / qry: genes x cells sparse counts with IDENTICAL row order.
log_norm <- function(counts) {
  tot <- Matrix::colSums(counts)
  sf <- stats::median(tot[tot > 0])
  m <- counts
  m@x <- log1p(m@x / rep(pmax(tot, 1), diff(m@p)) * sf)
  m
}
transfer_labels <- function(ref, ref_labels, qry, min_cells = 20, chunk = 50000) {
  stopifnot(identical(rownames(ref), rownames(qry)))
  ok <- !is.na(ref_labels)
  ref <- ref[, ok, drop = FALSE]; ref_labels <- as.character(ref_labels[ok])
  labs <- names(which(table(ref_labels) >= min_cells))
  R <- log_norm(methods::as(ref, "CsparseMatrix"))
  mu_r <- Matrix::rowMeans(R); sd_r <- sqrt(pmax(Matrix::rowMeans(R^2) - mu_r^2, 0))
  G <- Matrix::fac2sparse(factor(ref_labels[ref_labels %in% labs], levels = labs))            # labels x cells
  cent <- as.matrix(R[, ref_labels %in% labs, drop = FALSE] %*% Matrix::t(G)) / rep(Matrix::rowSums(G), each = nrow(R))
  good <- sd_r > 0
  Q <- log_norm(methods::as(qry, "CsparseMatrix"))
  mu_q <- Matrix::rowMeans(Q); sd_q <- sqrt(pmax(Matrix::rowMeans(Q^2) - mu_q^2, 0))
  good <- good & sd_q > 0
  Zc <- (cent[good, , drop = FALSE] - mu_r[good]) / sd_r[good]                   # genes x labels
  Zc <- sweep(Zc, 2, sqrt(colSums(Zc^2)), "/")
  Q <- Q[good, , drop = FALSE]; mu <- mu_q[good]; sdv <- sd_q[good]
  const <- as.vector(crossprod(mu / sdv, Zc))                                    # per label
  mu2 <- sum((mu / sdv)^2)
  pred <- character(ncol(Q)); score <- numeric(ncol(Q)); margin <- numeric(ncol(Q))
  for (s in seq(1, ncol(Q), by = chunk)) {
    idx <- s:min(ncol(Q), s + chunk - 1)
    q <- Q[, idx, drop = FALSE]
    dot <- as.matrix(Matrix::crossprod(q / sdv, Zc))                             # cells x labels: (x/sd) . z
    dot <- sweep(dot, 2, const)                                                  # ((x - mu)/sd) . z
    nrm <- sqrt(pmax(Matrix::colSums((q / sdv)^2) - 2 * as.vector(Matrix::crossprod(q, mu / sdv^2)) + mu2, 1e-12))
    cs <- dot / nrm
    top <- max.col(cs, ties.method = "first")
    pred[idx] <- labs[top]
    score[idx] <- cs[cbind(seq_along(idx), top)]
    cs[cbind(seq_along(idx), top)] <- -Inf
    margin[idx] <- score[idx] - apply(cs, 1, max)
  }
  data.frame(label = pred, score = score, margin = margin, row.names = colnames(qry))
}

## -----------------------------------------------------------------------
## Small statistics
## -----------------------------------------------------------------------
## Slope of y on x when BOTH are measured with noise (total least squares =
## direction of the first principal component). Ordinary regression would
## flatten the slope and fake "compression".
tls_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(NA_real_)
  v <- stats::prcomp(cbind(x[ok], y[ok]))$rotation[, 1]
  unname(v[2] / v[1])
}
boot_ci <- function(x, y, fun = tls_slope, n = 500) {
  ok <- which(is.finite(x) & is.finite(y))
  if (length(ok) < 10) return(c(NA_real_, NA_real_))
  b <- vapply(seq_len(n), function(i) { k <- sample(ok, replace = TRUE); fun(x[k], y[k]) }, numeric(1))
  unname(stats::quantile(b, c(0.025, 0.975), na.rm = TRUE))
}

## Tissue area actually covered by cells: occupied bins x bin area (um2).
## Independent of how each platform draws cell boundaries.
occupied_area <- function(x, y, bin = 25) length(unique(paste(floor(x / bin), floor(y / bin)))) * bin^2

jensen_shannon <- function(p, q) {
  p <- p / sum(p); q <- q / sum(q); m <- (p + q) / 2
  kl <- function(a, b) sum(ifelse(a > 0, a * log2(a / b), 0))
  (kl(p, m) + kl(q, m)) / 2
}

## paired scatter with identity line, one dot per core
paired_core_plot <- function(d, xcol, ycol, title, colour = "organ", log = FALSE) {
  lim <- range(c(d[[xcol]], d[[ycol]]), na.rm = TRUE)
  p <- ggplot(d, aes(.data[[xcol]], .data[[ycol]], colour = .data[[colour]])) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 2.2) + theme_bw(base_size = 10) +
    labs(x = "Xenium", y = "CosMx", title = title, colour = colour)
  if (log) p + scale_x_log10(limits = lim) + scale_y_log10(limits = lim) + coord_fixed() else
    p + coord_fixed(xlim = lim, ylim = lim)
}
