#!/usr/bin/env Rscript
# =============================================================================
# Trajectory figure: what a rescue curve looks like, and how spatial
# autocorrelation changes it.
#
#   Rscript plot_trajectories.R [traj_subset.csv]
#
# Expects a slice of all_results_combined.csv holding one treatment
# combination across landscape conditions. Base R only.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(if (length(args)) args[1] else "traj_subset.csv")

TPERT <- min(d$t[d$pert_value > 0])          # when the shift fires

# patch SD 0 is flat, so its autocorrelation label is meaningless: one condition
d$cond <- ifelse(d$patch_sd == 0, "flat", paste0("ac", d$ac))
d$run  <- interaction(d$draw, d$batch, d$rep, drop = TRUE)

conds  <- c("flat", "ac0", "ac2", "ac4")
titles <- c("flat landscape", "rough, uncorrelated", "rough, intermediate",
            "rough, correlated")
pal    <- c("#666666", "#D95F02", "#B8A02E", "#1B9E9E")

pdf("fig6_trajectories.pdf", width = 10, height = 5.4)
par(mfcol = c(2, 4), mar = c(3.6, 4.2, 2.4, 0.8), mgp = c(2.4, .6, 0),
    tcl = -.3, las = 1, bty = "n", oma = c(0, 0, 2, 0))

for (k in seq_along(conds)) {
  s <- d[d$cond == conds[k], ]

  # ---- top row: abundance, log scale ----
  plot(NA, xlim = c(0, 810), ylim = c(1, 2000), log = "y",
       xlab = "", ylab = if (k == 1) "population size" else "",
       main = titles[k], cex.main = 1, col.main = "#333333")
  grid(NA, NULL, col = "#eeeeee", lty = 1)
  abline(v = TPERT, col = "#cc3333", lwd = 1.2, lty = 2)
  for (r in unique(s$run)) {
    q <- s[s$run == r, ]; q <- q[order(q$t), ]
    lines(q$t, pmax(q$n, 1), col = adjustcolor(pal[k], .10), lwd = .7)
  }
  med <- tapply(s$n, s$t, median)
  lines(as.numeric(names(med)), pmax(med, 1), col = pal[k], lwd = 2.4)
  if (k == 1) {
    text(TPERT + 20, 1600, "perturbation", col = "#cc3333", cex = .72, adj = 0)
    legend("bottomright", c("individual runs", "median"), col = c(adjustcolor(pal[k], .5), pal[k]),
           lwd = c(.8, 2.4), bty = "n", cex = .7)
  }

  # ---- bottom row: trait variation ----
  plot(NA, xlim = c(0, 810), ylim = c(0, 42),
       xlab = "time", ylab = if (k == 1) "trait SD" else "")
  grid(NA, NULL, col = "#eeeeee", lty = 1)
  abline(v = TPERT, col = "#cc3333", lwd = 1.2, lty = 2)
  for (r in unique(s$run)) {
    q <- s[s$run == r, ]; q <- q[order(q$t), ]
    lines(q$t, q$u_sd, col = adjustcolor(pal[k], .10), lwd = .7)
  }
  medu <- tapply(s$u_sd, s$t, median, na.rm = TRUE)
  lines(as.numeric(names(medu)), medu, col = pal[k], lwd = 2.4)

  # survival annotation
  last <- do.call(rbind, lapply(split(s, s$run), function(q) q[which.max(q$t), ]))
  mtext(sprintf("%.0f%% persist", 100 * mean(last$n > 10)), side = 3,
        line = -1.3, cex = .65, col = pal[k])
}
mtext("Population size (top) and trait variation (bottom) through a perturbation",
      outer = TRUE, cex = .95, col = "#333333")
dev.off()
cat("wrote fig6_trajectories.pdf\n")
