#!/usr/bin/env Rscript
# =============================================================================
# Supplementary figure: landscape viability across heterogeneity amplitude
# and spatial autocorrelation, with no perturbation.
#
#   Rscript plot_boundary_map.R [boundary_map.csv]
#
# Left panel: median equilibrium abundance (t = 800). Right: fraction of runs
# surviving. 12 unperturbed runs per cell (3 landscape draws x 4 replicates),
# slope 1.0, dispersal 3, mutation on. Base R only.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(if (length(args)) args[1] else "boundary_map.csv")
d$alive <- as.numeric(d$final_n > 10)

sds <- sort(unique(d$sd)); acs <- sort(unique(d$ac))
M_n <- tapply(d$final_n, list(d$sd, d$ac), median)
M_a <- tapply(d$alive,   list(d$sd, d$ac), mean)

pdf("figS1_viability_map.pdf", width = 9.2, height = 4.3)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.8, 5.6), mgp = c(2.5, .7, 0),
    tcl = -.3, las = 1)

draw_map <- function(M, main, pal, legend_lab, breaks) {
  image(x = seq_along(sds), y = seq_along(acs), z = M,
        col = pal, breaks = breaks, axes = FALSE,
        xlab = "landscape heterogeneity (patch SD)",
        ylab = "spatial autocorrelation", main = main, cex.main = 1)
  axis(1, seq_along(sds), sds); axis(2, seq_along(acs), acs)
  for (i in seq_along(sds)) for (j in seq_along(acs))
    text(i, j, if (M[i, j] >= 10) round(M[i, j]) else sprintf("%.2f", M[i, j]),
         cex = .62, col = ifelse(M[i, j] > max(breaks) * .55, "white", "grey15"))
  # legend strip
  usr <- par("usr")
  xl <- usr[2] + .35; xr <- usr[2] + .65
  yy <- seq(usr[3] + .3, usr[4] - .3, length.out = length(pal) + 1)
  for (k in seq_along(pal))
    rect(xl, yy[k], xr, yy[k + 1], col = pal[k], border = NA, xpd = TRUE)
  text(xr + .12, c(yy[1], yy[length(yy)]),
       round(range(breaks)), cex = .6, adj = 0, xpd = TRUE)
  text((xl + xr)/2, usr[4] - .02, legend_lab, cex = .62, adj = c(.5, -1), xpd = TRUE)
}

pal_n <- hcl.colors(24, "YlGnBu", rev = TRUE)
draw_map(M_n, "equilibrium abundance (median, t = 800)",
         pal_n, "n", seq(0, max(M_n, na.rm = TRUE), length.out = 25))

pal_a <- hcl.colors(24, "YlGnBu", rev = TRUE)
M_a_txt <- M_a
image(x = seq_along(sds), y = seq_along(acs), z = M_a,
      col = pal_a, breaks = seq(0, 1, length.out = 25), axes = FALSE,
      xlab = "landscape heterogeneity (patch SD)",
      ylab = "spatial autocorrelation",
      main = "fraction of populations persisting", cex.main = 1)
axis(1, seq_along(sds), sds); axis(2, seq_along(acs), acs)
for (i in seq_along(sds)) for (j in seq_along(acs))
  text(i, j, sprintf("%.2f", M_a[i, j]), cex = .62,
       col = ifelse(M_a[i, j] > .55, "white", "grey15"))
# threshold contour: where survival crosses 0.5
contour(x = seq_along(sds), y = seq_along(acs), z = M_a, levels = .5,
        add = TRUE, drawlabels = FALSE, lwd = 2.2, col = "#cc3333", lty = 1)
usr <- par("usr")
xl <- usr[2] + .35; xr <- usr[2] + .65
yy <- seq(usr[3] + .3, usr[4] - .3, length.out = length(pal_a) + 1)
for (k in seq_along(pal_a))
  rect(xl, yy[k], xr, yy[k + 1], col = pal_a[k], border = NA, xpd = TRUE)
text(xr + .12, c(yy[1], yy[length(yy)]), c("0", "1"), cex = .6, adj = 0, xpd = TRUE)

dev.off()
cat("wrote figS1_viability_map.pdf\n")
