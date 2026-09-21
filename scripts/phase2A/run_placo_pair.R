#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

write_tsv_gz <- function(x, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", quote = FALSE)
  gz <- Sys.which("gzip")
  if (!nzchar(gz)) stop("gzip executable is required for compressed PLACO outputs")
  unlink(path)
  status <- system2(gz, c("-c", shQuote(tmp)), stdout = path)
  if (!identical(status, 0L) || !file.exists(path) || file.info(path)$size == 0) {
    stop("Failed to create compressed output: ", path)
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: run_placo_pair.R PAIR INPUT_TSV_GZ")
pair <- args[[1]]
input_path <- args[[2]]
root <- normalizePath(file.path(dirname(input_path), "../../.."), mustWork = TRUE)
results_dir <- file.path(root, "results", "phase2A")
dir.create(file.path(results_dir, "pair_runs"), recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "data", "phase2A", "PLACO_source", "PLACO_v0.2.0.R"))

threads <- as.integer(Sys.getenv("PHASE2A_THREADS", "4"))
if (is.na(threads) || threads < 1) threads <- 1L

dat <- fread(input_path, sep = "\t", showProgress = FALSE)
required <- c("SNP", "CHR", "BP", "A1", "A2", "BETA_1", "SE_1", "Z_1", "P_1", "BETA_2", "SE_2", "Z_2", "P_2")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ","))

dat[, Z_1_RECALC := BETA_1 / SE_1]
dat[, Z_2_RECALC := BETA_2 / SE_2]
z1_diff <- max(abs(dat$Z_1 - dat$Z_1_RECALC), na.rm = TRUE)
z2_diff <- max(abs(dat$Z_2 - dat$Z_2_RECALC), na.rm = TRUE)
if (!is.finite(z1_diff) || !is.finite(z2_diff) || z1_diff > 1e-8 || z2_diff > 1e-8) {
  stop("Signed Z consistency check failed")
}

zmat <- as.matrix(dat[, .(Z_1, Z_2)])
pmat <- as.matrix(dat[, .(P_1, P_2)])
VarZ <- var.placo(zmat, pmat, p.threshold = 1e-4)
CorZ <- cor.pearson(zmat, pmat, p.threshold = 1e-4, returnMatrix = FALSE)

prune_path <- file.path(results_dir, "reference_freq", "eur_ld_pruned_r2_0.1_1000kb.in")
if (file.exists(prune_path)) {
  pruned <- fread(prune_path, header = FALSE, col.names = "SNP", showProgress = FALSE)
  pruned_set <- dat[pruned, on = "SNP", nomatch = 0L]
  if (nrow(pruned_set) < 30L) {
    coord <- tstrsplit(pruned$SNP, ":", fixed = TRUE)
    if (length(coord) >= 2L) {
      prune_chr <- suppressWarnings(as.integer(coord[[1]]))
      prune_bp <- suppressWarnings(as.integer(coord[[2]]))
      keep <- !is.na(prune_chr) & !is.na(prune_bp)
      pruned_set <- dat[.(prune_chr[keep], prune_bp[keep]), on = .(CHR, BP), nomatch = 0L]
    }
  }
  if (nrow(pruned_set) >= 30L) {
    zmat_s <- as.matrix(pruned_set[, .(Z_1, Z_2)])
    pmat_s <- as.matrix(pruned_set[, .(P_1, P_2)])
    VarZ_s <- var.placo(zmat_s, pmat_s, p.threshold = 1e-4)
    CorZ_s <- cor.pearson(zmat_s, pmat_s, p.threshold = 1e-4, returnMatrix = FALSE)
  } else {
    VarZ_s <- c(NA_real_, NA_real_)
    CorZ_s <- NA_real_
  }
} else {
  pruned_set <- dat[0]
  VarZ_s <- c(NA_real_, NA_real_)
  CorZ_s <- NA_real_
}

extreme <- dat[abs(Z_1)^2 > 80 | abs(Z_2)^2 > 80]
primary <- dat[!(abs(Z_1)^2 > 80 | abs(Z_2)^2 > 80)]

calc_one <- function(i, frame = primary) {
  ans <- tryCatch(placo.plus(Z = c(frame$Z_1[[i]], frame$Z_2[[i]]), VarZ = VarZ, CorZ = CorZ), error = function(e) NULL)
  if (is.null(ans)) return(NA_real_)
  as.numeric(ans$p.placo.plus)
}

idx <- seq_len(nrow(primary))
if (length(idx)) {
  chunks <- split(idx, cut(seq_along(idx), breaks = max(1L, min(threads * 4L, length(idx))), labels = FALSE))
  if (threads > 1L && .Platform$OS.type != "windows") {
    pvals <- unlist(parallel::mclapply(chunks, function(chunk) vapply(chunk, calc_one, numeric(1)), mc.cores = threads, mc.preschedule = FALSE), use.names = FALSE)
  } else {
    pvals <- vapply(idx, calc_one, numeric(1))
  }
} else {
  pvals <- numeric(0)
}
primary[, P_PLACO_PLUS := pvals]
primary[, Z_PRODUCT := Z_1 * Z_2]
primary[, direction_product := fifelse(BETA_1 * BETA_2 > 0, "POSITIVE_EFFECT_PRODUCT", fifelse(BETA_1 * BETA_2 < 0, "NEGATIVE_EFFECT_PRODUCT", "ZERO_EFFECT_PRODUCT"))]
primary[, pair := pair]
extreme[, Z_PRODUCT := Z_1 * Z_2]
extreme[, direction_product := fifelse(BETA_1 * BETA_2 > 0, "POSITIVE_EFFECT_PRODUCT", fifelse(BETA_1 * BETA_2 < 0, "NEGATIVE_EFFECT_PRODUCT", "ZERO_EFFECT_PRODUCT"))]
extreme[, dual_gwas_significant := P_1 < 5e-8 & P_2 < 5e-8]
extreme[, pair := pair]

keep_primary <- c("pair", "SNP", "CHR", "BP", "A1", "A2", "BETA_1", "BETA_2", "SE_1", "SE_2", "Z_1", "Z_2", "P_1", "P_2", "EAF_1", "EAF_2", "N_1", "N_2", "Z_PRODUCT", "P_PLACO_PLUS", "direction_product")
keep_extreme <- c("pair", "SNP", "CHR", "BP", "A1", "A2", "BETA_1", "BETA_2", "Z_1", "Z_2", "P_1", "P_2", "Z_PRODUCT", "direction_product", "dual_gwas_significant")
keep_primary <- intersect(keep_primary, names(primary))
keep_extreme <- intersect(keep_extreme, names(extreme))
write_tsv_gz(primary[, ..keep_primary], file.path(results_dir, "pair_runs", paste0(pair, ".placo.tsv.gz")))
write_tsv_gz(extreme[, ..keep_extreme], file.path(results_dir, "pair_runs", paste0(pair, ".extreme.tsv.gz")))

null_mask <- pmat[, 1] >= 1e-4 | pmat[, 2] >= 1e-4
sign_row <- data.table(
  pair = pair,
  eligible_variants = nrow(dat),
  primary_variants = nrow(primary),
  extreme_z_variants = nrow(extreme),
  z1_positive = sum(dat$Z_1 > 0, na.rm = TRUE), z1_negative = sum(dat$Z_1 < 0, na.rm = TRUE), z1_zero = sum(dat$Z_1 == 0, na.rm = TRUE),
  z2_positive = sum(dat$Z_2 > 0, na.rm = TRUE), z2_negative = sum(dat$Z_2 < 0, na.rm = TRUE), z2_zero = sum(dat$Z_2 == 0, na.rm = TRUE),
  z1_nonfinite = sum(!is.finite(dat$Z_1)), z2_nonfinite = sum(!is.finite(dat$Z_2)),
  placo_nonfinite = sum(!is.finite(primary$P_PLACO_PLUS)),
  null_variants = sum(null_mask, na.rm = TRUE),
  ld_pruned_null_variants = nrow(pruned_set),
  VarZ1 = VarZ[[1]], VarZ2 = VarZ[[2]], CorZ = CorZ,
  VarZ1_ld_pruned = VarZ_s[[1]], VarZ2_ld_pruned = VarZ_s[[2]], CorZ_ld_pruned = CorZ_s
)
fwrite(sign_row, file.path(results_dir, "pair_runs", paste0(pair, ".sign_qc.tsv")), sep = "\t", quote = FALSE)

param_row <- data.table(pair = pair, null_variants = sum(null_mask, na.rm = TRUE), VarZ1 = VarZ[[1]], VarZ2 = VarZ[[2]], CorZ = CorZ, ld_pruned_null_variants = nrow(pruned_set), VarZ1_ld_pruned = VarZ_s[[1]], VarZ2_ld_pruned = VarZ_s[[2]], CorZ_ld_pruned = CorZ_s)
fwrite(param_row, file.path(results_dir, "pair_runs", paste0(pair, ".params.tsv")), sep = "\t", quote = FALSE)
cat(pair, "complete", nrow(primary), "primary variants\n")
