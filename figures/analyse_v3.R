#!/usr/bin/env Rscript
# ============================================================================
# Analysis for the v3 dataset. Reproduces every number in
# methods_and_results_summary_v3: viability classification, core-cell filter,
# phase profiles, models with cell-clustered SEs and joint Wald tests,
# variance explained, and the abundance-mediation check. Base R only.
#
#   Rscript analyse_v3.R          (expects run_summary_v3.csv, controls_v3.csv)
#
# Outputs: viability_final.csv, phase_profiles_v3.csv, models_v3.txt
# ============================================================================

## ---- 1. Viability classification from the unperturbed controls ------------
ctrl <- read.csv("controls_v3.csv")
K5 <- c("slope","disp","mutsd","patch_sd","ac")
# one row per run at each checkpoint -> wide by t
w <- reshape(ctrl, idvar = c(K5,"draw","rep"), timevar = "t", direction = "wide")
w$alive1500 <- !is.na(w$n.1500) & w$n.1500 > 10
via <- aggregate(cbind(surv = w$alive1500, n800 = w$n.800, n1500 = w$n.1500),
                 w[K5], function(x) if (is.logical(x)) mean(x) else median(x, na.rm = TRUE))
via$surv  <- aggregate(w$alive1500, w[K5], mean)$x
via$n800  <- aggregate(w$n.800,  w[K5], median, na.rm = TRUE)$x
via$n1500 <- aggregate(w$n.1500, w[K5], function(x) if (all(is.na(x))) NA else median(x, na.rm = TRUE))$x
via$trend <- (via$n1500 - via$n800) / via$n800
via$klass <- ifelse(via$surv < 0.8 | (!is.na(via$n1500) & via$n1500 < 50), "non-viable",
             ifelse(!is.na(via$trend) & via$trend < -0.10, "declining", "viable"))
write.csv(via[, c(K5,"surv","n800","n1500","trend","klass")], "viability_final.csv", row.names = FALSE)
cat("viability (slopes 0.6-1.0):",
    paste(names(table(via$klass[via$slope %in% c(.6,.8,1)])),
          table(via$klass[via$slope %in% c(.6,.8,1)]), collapse="  "), "\n")

## ---- 2. Core dataset: balanced factorial, persisting cells ----------------
d <- read.csv("run_summary_v3.csv")
d <- d[d$slope > 0, ]
d <- merge(d, via[, c("slope","disp","mutsd","patch_sd","ac","surv","klass")],
           by = c("slope","disp","mutsd","patch_sd","ac"))
core <- d[d$slope < 1.2 & d$surv >= 0.8, ]
core$resist  <- core$trough_n / core$base_n
core$recover <- (core$peak_n - core$trough_n) / core$base_n
core$cell <- interaction(core$slope, core$disp, core$mutsd, core$pert_treat,
                         core$patch_sd, core$ac, core$draw, drop = TRUE)
s <- core[core$survived == 1 & core$base_n > 0 & core$peak_t > core$trough_t, ]
s$y <- log1p(s$recover)
h  <- s[s$patch_sd > 0, ];  hd <- core[core$patch_sd > 0, ]
cat("core:", nrow(core), "runs;", nlevels(droplevels(factor(
    interaction(core$slope,core$disp,core$mutsd,core$patch_sd,core$ac)))), "cells;",
    nrow(s), "survivors\n\n")

## ---- 3. Phase profiles (cluster means +/- SE) -----------------------------
prof <- function(dat, fac, val) {
  cm <- aggregate(dat[[val]], list(level = dat[[fac]], cell = dat$cell), mean)
  ag <- aggregate(cm$x, list(level = cm$level),
                  function(z) c(m = mean(z), se = sd(z)/sqrt(length(z))))
  data.frame(level = ag$level, mean = ag$x[, "m"], se = ag$x[, "se"])
}
rows <- list()
flab <- c(disp = "dispersal", slope = "slope", mutsd = "mutation", ac = "autocorrelation")
olab <- c(resist = "resistance", recover = "recovery")
for (fc in c("disp","slope","mutsd","ac")) {
  for (v in c("resist","recover")) {
    dat <- if (fc == "ac") h[h$patch_sd == 2, ] else s
    p <- prof(dat, fc, v); p$factor <- flab[fc]; p$outcome <- olab[v]; rows[[length(rows)+1]] <- p
  }
  dat <- if (fc == "ac") hd[hd$patch_sd == 2, ] else core
  p <- prof(dat, fc, "survived"); p$factor <- flab[fc]; p$outcome <- "persistence"
  rows[[length(rows)+1]] <- p
}
profiles <- do.call(rbind, rows)
write.csv(profiles, "phase_profiles_v3.csv", row.names = FALSE)

## ---- 4. Models: clustered SEs, joint Wald tests, variance explained -------
# cluster-robust covariance (CR1) for lm/glm
clustered_vcov <- function(m, cl) {
  X <- model.matrix(m); u <- residuals(m, type = "response")
  if (inherits(m, "glm")) u <- residuals(m, type = "pearson") * sqrt(summary(m)$dispersion)
  cl <- droplevels(factor(cl)); M <- nlevels(cl); N <- nrow(X); K <- ncol(X)
  meat <- Reduce(`+`, lapply(split(seq_len(N), cl), function(ix) {
    g <- colSums(X[ix, , drop = FALSE] * u[ix]); tcrossprod(g) }))
  bread <- if (inherits(m, "glm")) chol2inv(chol(crossprod(X * sqrt(m$weights)))) else
           chol2inv(chol(crossprod(X)))
  adj <- M/(M-1) * (N-1)/(N-K)
  bread %*% meat %*% bread * adj
}
joint_wald <- function(m, V, pattern) {
  b <- coef(m); ix <- grep(pattern, names(b)); ix <- ix[!grepl(":", names(b)[ix])]
  if (!length(ix)) return(NA)
  W <- t(b[ix]) %*% solve(V[ix, ix]) %*% b[ix]
  pchisq(as.numeric(W), df = length(ix), lower.tail = FALSE)
}
F0 <- y ~ factor(slope)*factor(pert_treat) + factor(ac) + factor(patch_sd) +
          factor(disp) + factor(mutsd) + factor(slope):factor(disp)
sink("models_v3.txt")
cat("=== v3 models: balanced core, patchy landscapes, cell-clustered SEs ===\n\n")
for (spec in list(list(y="resist", lab="resistance"), list(y="y", lab="log1p recovery"))) {
  f <- as.formula(sub("^y", spec$y, paste(deparse(F0), collapse = " ")))
  m <- lm(f, h); V <- clustered_vcov(m, h$cell)
  cat(spec$lab, ": joint p --",
      " slope", format.pval(joint_wald(m,V,"factor\\(slope\\)"), digits=2),
      " disp",  format.pval(joint_wald(m,V,"factor\\(disp\\)"),  digits=2),
      " mut",   format.pval(joint_wald(m,V,"factor\\(mutsd\\)"), digits=2),
      " AC",    format.pval(joint_wald(m,V,"factor\\(ac\\)"),    digits=2), "\n")
}
mp <- glm(update(F0, survived ~ .), binomial, hd); Vp <- clustered_vcov(mp, hd$cell)
cat("persistence   : joint p --",
    " disp", format.pval(joint_wald(mp,Vp,"factor\\(disp\\)"), digits=2),
    " AC",   format.pval(joint_wald(mp,Vp,"factor\\(ac\\)"),   digits=2), "\n\n")
# variance explained (recovery), dropping factor + its interactions
full <- lm(F0, h)
for (t in c("slope","disp","mutsd","ac")) {
  tt <- terms(F0); keep <- !grepl(paste0("factor\\(", t, "\\)"), attr(tt, "term.labels"))
  red <- lm(reformulate(attr(tt, "term.labels")[keep], "y"), h)
  ss <- sum(residuals(red)^2) - sum(residuals(full)^2)
  cat("eta2 recovery", t, ":", round(100*ss/(ss+sum(residuals(full)^2)),1), "%\n")
}
# mediation: AC given log(baseline)
mb <- lm(update(F0, . ~ . + log(base_n)), h); Vb <- clustered_vcov(mb, h$cell)
cat("\nmediation: AC joint p with log(baseline) in model:",
    format.pval(joint_wald(mb,Vb,"factor\\(ac\\)"), digits=2),
    "; log(baseline) coef", round(coef(mb)["log(base_n)"],3), "\n")
sink()
cat(readLines("models_v3.txt"), sep = "\n")

## ---- 5. Interaction screen: drop-one from the full two-way model -----------
# Same convention for all three outcomes: fit all 15 pairwise interactions, drop
# one, record the loss (partial eta2 for the phases; delta McFadden pseudo-R2 for
# persistence). Matches the table in the main-messages document.
F <- c("slope","disp","mutsd","pert_treat","patch_sd","ac")
pairs <- combn(F, 2)
ints  <- apply(pairs, 2, function(p) paste0("factor(", p[1], "):factor(", p[2], ")"))
mainsF <- paste0("factor(", F, ")")
fullf <- function(y, drop = NULL) reformulate(c(mainsF, setdiff(ints, drop)), y)
sink("models_v3.txt", append = TRUE)
cat("\n=== Interaction screen (drop-one from full two-way model) ===\n")
scr <- data.frame(interaction = apply(pairs, 2, paste, collapse = " x "))
for (yv in c("resist","y")) {
  fl <- lm(fullf(yv), h); ssf <- sum(residuals(fl)^2)
  scr[[yv]] <- sapply(ints, function(t) { r <- lm(fullf(yv, t), h)
    ss <- sum(residuals(r)^2) - ssf; round(100*ss/(ss+ssf), 1) })
}
pfull <- glm(fullf("survived"), binomial, hd)
pr2 <- function(m) 1 - m$deviance/m$null.deviance
scr$persist_dR2 <- sapply(ints, function(t) round(100*(pr2(pfull) - pr2(glm(fullf("survived", t), binomial, hd))), 2))
names(scr)[2:3] <- c("resist_eta2","recover_eta2")
print(scr[order(-scr$recover_eta2), ], row.names = FALSE)
sink()
cat(readLines("models_v3.txt"), sep = "\n")

## ---- 6. Figure inputs: every CSV the plot_*.R scripts read -----------------
# interactions_final.csv / maineffects_eta2.csv
write.csv(data.frame(row.names = scr$interaction, resistance = scr$resist_eta2,
                     recovery = scr$recover_eta2, persistence_dR2 = scr$persist_dR2),
          "interactions_final.csv")
me <- data.frame(row.names = F)
for (yv in c("resist","y")) {
  fl <- lm(reformulate(mainsF, yv), h); ssf <- sum(residuals(fl)^2)
  me[[yv]] <- sapply(F, function(f) { r <- lm(reformulate(setdiff(mainsF, paste0("factor(",f,")")), yv), h)
    ss <- sum(residuals(r)^2) - ssf; round(100*ss/(ss+ssf), 1) })
}
names(me) <- c("resistance","recovery"); write.csv(me, "maineffects_eta2.csv")

# cells_persist_base_v3.csv (cell x draw medians/means) and fig3_panels.csv
cb <- aggregate(cbind(base = core$base_n, persist = core$survived),
                core[c(K5,"draw")], function(x) if (all(x %in% 0:1)) mean(x) else median(x))
cb$base    <- aggregate(core$base_n,   core[c(K5,"draw")], median)$x
cb$persist <- aggregate(core$survived, core[c(K5,"draw")], mean)$x
write.csv(cb, "cells_persist_base_v3.csv", row.names = FALSE)
q <- cb[cb$patch_sd > 0, ]
f3 <- aggregate(q$base, q[c("patch_sd","ac")], function(x) c(base = median(x),
        base_lo = quantile(x,.25), base_hi = quantile(x,.75)))
f3 <- data.frame(f3[1:2], f3$x); names(f3)[3:5] <- c("base","base_lo","base_hi")
pp <- aggregate(q$persist, q[c("patch_sd","ac")], function(x) c(m = mean(x), se = sd(x)/sqrt(length(x))))
f3$persist <- pp$x[, "m"]; f3$pse <- pp$x[, "se"]
write.csv(f3, "fig3_panels.csv", row.names = FALSE)

# fig5_data.csv (homogeneous vs gradient)
raw <- read.csv("run_summary_v3.csv"); s0 <- raw[raw$slope == 0, ]
f5 <- rbind(data.frame(slope = "0 (homogeneous)", baseline = median(s0$base_n), persistence = mean(s0$survived)),
            do.call(rbind, lapply(c(0.6,0.8,1.0), function(sl) { q <- core[core$slope == sl, ]
              data.frame(slope = as.character(sl), baseline = median(q$base_n), persistence = mean(q$survived)) })))
write.csv(f5, "fig5_data.csv", row.names = FALSE)

# fig6_draws.csv (recovery x AC per landscape draw, SD 2)
s2 <- s[s$patch_sd == 2, ]
write.csv(aggregate(s2$recover, s2[c("draw","ac")], mean) |> setNames(c("draw","ac","recover")),
          "fig6_draws.csv", row.names = FALSE)

# fig6_kernel.csv from the extended kernel check (if present)
if (!file.exists("fig6_kernel.csv")) {
  k <- read.csv(file.path(getwd(), "Supplement", "kernel_check2.csv")); k <- k[k$patch_sd > 0, ]
  k$run <- interaction(k$kernel, k$ac, k$slope, k$disp, k$mutsd, k$pert, k$draw, k$rep, drop = TRUE)
  per <- do.call(rbind, lapply(split(k, k$run), function(g) {
    g <- g[order(g$t), ]; pre <- g[g$t <= 250, ]; post <- g[g$t > 250, ]
    if (!nrow(pre) || !nrow(post)) return(NULL)
    tr <- post[which.min(post$n), ]; aft <- post[post$t >= tr$t, ]
    data.frame(kernel = g$kernel[1], cond = paste0("AC", g$ac[1]),
               recover = (max(aft$n) - tr$n)/median(pre$n), surv = as.integer(tail(aft$n,1) > 10)) }))
  out <- aggregate(cbind(recover, surv) ~ kernel + cond, per, function(x) mean(x))
  out$recover <- sapply(seq_len(nrow(out)), function(i) { z <- per[per$kernel==out$kernel[i] & per$cond==out$cond[i],]
    mean(z$recover[z$surv == 1]) })
  names(out)[4] <- "persist"; write.csv(out, "fig6_kernel.csv", row.names = FALSE)
}
vl <- via; vl$label <- ifelse(vl$patch_sd == 0, "flat", paste0("SD", vl$patch_sd, "/AC", vl$ac))
write.csv(vl, "viability_labeled.csv", row.names = FALSE)
cat("figure inputs written\n")
