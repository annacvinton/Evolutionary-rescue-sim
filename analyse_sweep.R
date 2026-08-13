#!/usr/bin/env Rscript
# =============================================================================
# Evolutionary rescue: phase decomposition analysis
#
#   Rscript analyse_sweep.R [run_summary.csv]
#
# Reads the run-level summary produced by reduce_sweep.py and writes:
#   models_table.csv     joint p-values, every parameter x every outcome
#   effects_table.csv    means by treatment level for each outcome
#   fig1_phase_signatures.pdf
#   fig2_comparison_matrix.pdf
#   fig3_dispersal.pdf
#   fig4_slope_severity.pdf
#   fig5_draw_robustness.pdf
#
# Base R only. No packages required.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args)) args[1] else "run_summary.csv"
d <- read.csv(infile, stringsAsFactors = FALSE)
cat(sprintf("read %s: %d runs\n", infile, nrow(d)))

# ---- derive phase metrics ---------------------------------------------------
d$cell <- interaction(d$slope, d$disp, d$mutsd, d$pert_treat,
  d$patch_sd, d$ac, d$draw,
  drop = TRUE
)
# note: runs are identified by rep AND batch; cell deliberately excludes both
d$resist <- d$trough_n / d$base_n # depth: higher = shallower crash
d$recover <- (d$peak_n - d$trough_n) / d$base_n
d$rec_rel <- d$peak_n / d$base_n
d$rec_time <- d$peak_t - d$trough_t
d <- d[is.finite(d$resist) & d$base_n > 0, ]
rownames(d) <- NULL # row names must be 1..n for the cluster index below

s <- d[d$survived == 1 & d$rec_time > 0, ]
rownames(s) <- NULL
s$y_recover <- log1p(pmax(s$recover, -0.999))
s$y_rec_rel <- log1p(pmax(s$rec_rel, -0.999))
s$y_time <- log(s$rec_time)
cat(sprintf("  %d runs, %d survivors with a rebound\n", nrow(d), nrow(s)))

for (v in c("slope", "disp", "mutsd", "pert_treat", "patch_sd", "ac", "draw")) {
  d[[paste0("f_", v)]] <- factor(d[[v]])
  s[[paste0("f_", v)]] <- factor(s[[v]])
}

# ---- cluster-robust covariance (no packages) --------------------------------
vcovCL <- function(fit, cluster) {
  X <- model.matrix(fit)
  u <- residuals(fit, type = "response")
  if (inherits(fit, "glm")) {
    u <- u * fit$family$mu.eta(fit$linear.predictors) /
      fit$family$variance(fitted(fit))
  }
  idx <- as.integer(rownames(X))
  stopifnot(max(idx) <= length(cluster)) # guard against the subsetting trap
  cl <- factor(cluster[idx])
  bread <- summary(fit)$cov.unscaled
  if (inherits(fit, "glm")) bread <- vcov(fit)
  meat <- matrix(0, ncol(X), ncol(X))
  for (g in levels(cl)) {
    i <- which(cl == g)
    if (!length(i)) next
    sg <- crossprod(X[i, , drop = FALSE], u[i])
    meat <- meat + tcrossprod(sg)
  }
  G <- nlevels(cl)
  n <- nrow(X)
  k <- ncol(X)
  adj <- (G / (G - 1)) * ((n - 1) / (n - k))
  if (inherits(fit, "glm")) {
    return(adj * bread %*% meat %*% bread)
  }
  sig <- bread # (X'X)^-1
  adj * sig %*% meat %*% sig
}

# joint Wald test. `main` = the factor's own terms only (no colon);
# `inter` = terms containing both factors. Matching is on term structure,
# not a loose regex, because R orders interaction names unpredictably.
selMain <- function(nm, f) grep(paste0("^", f), nm)[!grepl(":", nm[grep(paste0("^", f), nm)])]
selInter <- function(nm, f1, f2) which(grepl(":", nm) & grepl(f1, nm) & grepl(f2, nm))
jointP <- function(fit, V, idx) {
  b <- coef(fit)
  if (!length(idx)) {
    return(NA_real_)
  }
  bb <- b[idx]
  VV <- V[idx, idx, drop = FALSE]
  ok <- is.finite(bb) & is.finite(diag(VV))
  if (!any(ok)) {
    return(NA_real_)
  }
  bb <- bb[ok]
  VV <- VV[ok, ok, drop = FALSE]
  W <- tryCatch(as.numeric(t(bb) %*% solve(VV) %*% bb), error = function(e) NA_real_)
  if (!is.finite(W)) {
    return(NA_real_)
  }
  pchisq(W, df = length(bb), lower.tail = FALSE)
}

# ---- model set --------------------------------------------------------------
FORM <- paste(
  "~ f_slope * f_pert_treat + f_ac * f_pert_treat +",
  "f_patch_sd * f_pert_treat + f_disp + f_mutsd"
)
getIdx <- function(nm, what) {
  switch(what,
    slope = selMain(nm, "f_slope"),
    autocorrelation = selMain(nm, "f_ac"),
    patch_sd = selMain(nm, "f_patch_sd"),
    dispersal = selMain(nm, "f_disp"),
    mutation = selMain(nm, "f_mutsd"),
    `ac x severity` = selInter(nm, "f_ac", "f_pert_treat"),
    `slope x severity` = selInter(nm, "f_slope", "f_pert_treat")
  )
}
term_names <- c(
  "slope", "autocorrelation", "patch_sd", "dispersal", "mutation",
  "ac x severity", "slope x severity"
)

outcomes <- list(
  resistance = list(y = "resist", dat = "s"),
  `recovery magnitude` = list(y = "y_recover", dat = "s"),
  `recovery vs baseline` = list(y = "y_rec_rel", dat = "s"),
  `recovery time` = list(y = "y_time", dat = "s"),
  persistence = list(y = "survived", dat = "d")
)

res <- matrix(NA_real_, length(outcomes), length(term_names),
  dimnames = list(names(outcomes), term_names)
)
fits <- list()
# Autocorrelation is undefined when patch SD is 0 (a flat field has no spatial
# structure), so landscape terms are fitted on patchy landscapes only. Slope,
# dispersal and mutation are fitted on everything.
landscape_terms <- c("autocorrelation", "patch_sd", "ac x severity")
for (o in names(outcomes)) {
  spec <- outcomes[[o]]
  dat <- get(spec$dat)
  f <- as.formula(paste(spec$y, FORM))
  fit <- if (o == "persistence") glm(f, binomial, dat) else lm(f, dat)
  V <- vcovCL(fit, dat$cell)
  nm <- names(coef(fit))
  pat <- dat[dat$patch_sd > 0, ]
  rownames(pat) <- NULL
  fitP <- if (o == "persistence") glm(f, binomial, pat) else lm(f, pat)
  VP <- vcovCL(fitP, pat$cell)
  nmP <- names(coef(fitP))
  for (tn in term_names) {
    res[o, tn] <- if (tn %in% landscape_terms) {
      jointP(fitP, VP, getIdx(nmP, tn))
    } else {
      jointP(fit, V, getIdx(nm, tn))
    }
  }
  fits[[o]] <- fit
  cat(sprintf("  fitted %-22s n=%d (patchy subset n=%d)\n", o, nobs(fit), nobs(fitP)))
}
write.csv(round(res, 8), "models_table.csv")

cat("\n================ JOINT p-VALUES ================\n")
pr <- res
show <- apply(pr, c(1, 2), function(p) {
  if (is.na(p)) "  na  " else if (p > .05) "  ---  " else if (p < 1e-16) "<1e-16" else formatC(p, format = "e", digits = 1)
})
print(noquote(show))
cat("\n  '---' = p > 0.05\n")

# ---- effect tables ----------------------------------------------------------
eff <- function(var, dat_s, dat_d, restrict = NULL) {
  a <- if (is.null(restrict)) dat_s else dat_s[restrict(dat_s), ]
  b <- if (is.null(restrict)) dat_d else dat_d[restrict(dat_d), ]
  m <- aggregate(cbind(resist, recover, rec_rel, rec_time) ~ a[[var]], a, mean)
  names(m)[1] <- var
  p <- aggregate(survived ~ b[[var]], b, mean)
  names(p) <- c(var, "persistence")
  merge(m, p, by = var)
}
patchy <- function(z) z$patch_sd > 0
tabs <- list(
  autocorrelation = eff("ac", s, d, patchy),
  patch_sd        = eff("patch_sd", s, d),
  slope           = eff("slope", s, d),
  dispersal       = eff("disp", s, d),
  mutation        = eff("mutsd", s, d)
)
out <- do.call(rbind, lapply(names(tabs), function(n) {
  cbind(parameter = n, setNames(tabs[[n]], c("level", names(tabs[[n]])[-1])))
}))
write.csv(out, "effects_table.csv", row.names = FALSE)
cat("\n================ EFFECT SIZES ================\n")
print(out, digits = 3)


cat("\nwrote models_table.csv, effects_table.csv")
