# =============================================================================
# TODO
# =============================================================================
#' 1. Plot landscapes to show effects of simulation parameters and batches
#' 2. Plot a few abundance curves to show how the simulation works over time
#' 3. Re-do plots below with unified colour scheme

library(ggplot2)
library(viridis)

# =============================================================================
# DATA
# =============================================================================
source("analyse_sweep.R")
simsteps_df <- read.csv("all_results_combined.csv")

head(simsteps_df)
simsteps_df$ID <- with(simsteps_df, paste(slope, disp, mutsd, pert_treat, patch_sd, ac, draw, batch, rep, sep = "_"))

## find ID for which we have the most data points
TargetID <- names(which.max(table(simsteps_df$ID)))

plot_df <- simsteps_df[simsteps_df$ID == TargetID, ]
plot(plot_df$t, plot_df$n, type = "l", xlab = "time", ylab = "abundance")



# =============================================================================
# PLOTS
# =============================================================================

## Landscapes ------------------------------
LS_fs <- list.files(file.path(getwd(), "simulation", "landscapes"), pattern = ".txt")

lapply(LS_fs, function(i) {
    landscape <- read.table(file.path(getwd(), "simulation", "landscapes", i), header = FALSE)
    landscape$x <- rep(1:100, each = 100)
    landscape$y <- rep(1:100, times = 100)
    # print(summary(landscape))
    # # # convert 10000 entry long vector into 100x100 matrix
    # landscape_matrix <- matrix(landscape$V1, nrow = 100, ncol = 100, byrow = TRUE)
    # # image(landscape_matrix, axes = FALSE, col = hcl.colors(24, "Blues", rev = TRUE), main = title)
    # file name contains information about the simulation parameters like so "L_ac0_sd0_r0.txt" for AC, SD, and R values of 0. We split these out to use as the title of the plot and name of the produced png file
    title <- tools::file_path_sans_ext(gsub("L_", "", i))
    ## plotting with ggplot2
    p <- ggplot(landscape, aes(x = x, y = y, fill = V1)) +
        geom_tile() +
        scale_fill_viridis() +
        labs(title = title, x = NULL, y = NULL) +
        theme_bw() +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    save_path <- file.path(getwd(), "simulation", "landscapes", paste0(tools::file_path_sans_ext(i), ".png"))
    ggsave(save_path, p, width = 5, height = 5, dpi = 300)
})





pal <- c("#1B9E9E", "#E08214", "#8DA0CB")
setpar <- function(...) {
    par(
        mar = c(4.2, 4.4, 2.6, 1), mgp = c(2.5, .7, 0),
        tcl = -.3, las = 1, bty = "n", ...
    )
}
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
bars <- function(x, m, e, col) {
    arrows(x, m - e, x, m + e,
        angle = 90,
        code = 3, length = .03, col = col
    )
}

# --- Fig 1: the phase signatures -------------------------------------------
pdf("fig1_phase_signatures.pdf", 8.5, 4.6)
setpar(mfrow = c(1, 2))
hs <- s[s$patch_sd > 0, ]
for (j in 1:2) {
    v <- c("ac", "patch_sd")[j]
    dat <- if (j == 1) hs else s
    lv <- sort(unique(dat[[v]]))
    r_m <- sapply(lv, function(l) mean(dat$resist[dat[[v]] == l], na.rm = TRUE))
    r_e <- sapply(lv, function(l) sem(dat$resist[dat[[v]] == l]))
    c_m <- sapply(lv, function(l) mean(dat$recover[dat[[v]] == l], na.rm = TRUE))
    c_e <- sapply(lv, function(l) sem(dat$recover[dat[[v]] == l]))
    rs <- r_m / mean(r_m)
    cs <- c_m / mean(c_m)
    plot(seq_along(lv), rs,
        type = "o", pch = 16, col = pal[1], lwd = 2,
        ylim = range(c(
            rs - r_e / mean(r_m), cs - c_e / mean(c_m),
            rs + r_e / mean(r_m), cs + c_e / mean(c_m)
        )),
        xaxt = "n", xlab = c("spatial autocorrelation", "patch SD")[j],
        ylab = "effect, scaled to its own mean",
        main = c(
            "autocorrelation acts on recovery",
            "patch SD acts on resistance"
        )[j], cex.main = 1
    )
    axis(1, seq_along(lv), lv)
    bars(seq_along(lv), rs, r_e / mean(r_m), pal[1])
    lines(seq_along(lv), cs, type = "o", pch = 17, col = pal[2], lwd = 2)
    bars(seq_along(lv), cs, c_e / mean(c_m), pal[2])
    abline(h = 1, lty = 3, col = "grey60")
    if (j == 1) {
        legend("topleft", c("resistance", "recovery"),
            col = pal[1:2],
            pch = c(16, 17), lwd = 2, bty = "n", cex = .85
        )
    }
}
dev.off()

# --- Fig 2: comparison matrix ----------------------------------------------
pdf("fig2_comparison_matrix.pdf", 8.2, 4.4)
setpar(mar = c(7, 10, 3, 5))
M <- -log10(pmax(res, 1e-17))
M[is.na(M)] <- 0
image(t(M[nrow(M):1, ]),
    col = hcl.colors(24, "Blues", rev = TRUE),
    axes = FALSE, main = "evidence by outcome and parameter"
)
axis(1, seq(0, 1, length.out = ncol(M)), colnames(M), las = 2, cex.axis = .8)
axis(2, seq(0, 1, length.out = nrow(M)), rev(rownames(M)), las = 1, cex.axis = .8)
for (i in seq_len(nrow(M))) {
    for (j in seq_len(ncol(M))) {
        lab <- if (res[i, j] > .05 || is.na(res[i, j])) "ns" else if (res[i, j] < 1e-16) "***" else if (res[i, j] < 1e-4) "**" else "*"
        text((j - 1) / (ncol(M) - 1), (nrow(M) - i) / (nrow(M) - 1), lab,
            cex = .8, col = ifelse(M[i, j] > 8, "white", "grey20")
        )
    }
}
mtext("-log10 p;  *** p<1e-16   ** p<1e-4   * p<0.05   ns not significant",
    side = 1, line = 5.4, cex = .72
)
dev.off()

# --- Fig 3: dispersal -------------------------------------------------------
pdf("fig3_dispersal.pdf", 9.5, 3.6)
setpar(mfrow = c(1, 3))
lv <- sort(unique(s$disp))
for (k in 1:3) {
    v <- c("resist", "recover", NA)[k]
    if (k < 3) {
        m <- sapply(lv, function(l) mean(s[[v]][s$disp == l], na.rm = TRUE))
        e <- sapply(lv, function(l) sem(s[[v]][s$disp == l]))
        ttl <- c("resistance", "recovery")[k]
    } else {
        m <- sapply(lv, function(l) mean(d$survived[d$disp == l]))
        e <- sapply(lv, function(l) sem(d$survived[d$disp == l]))
        ttl <- "persistence"
    }
    plot(seq_along(lv), m,
        type = "o", pch = 16, col = pal[1], lwd = 2,
        ylim = range(c(m - e, m + e)), xaxt = "n",
        xlab = "dispersal SD", ylab = ttl, main = ttl, cex.main = 1
    )
    axis(1, seq_along(lv), lv)
    bars(seq_along(lv), m, e, pal[1])
}
dev.off()

# --- Fig 4: slope x severity ------------------------------------------------
pdf("fig4_slope_severity.pdf", 8.5, 4.6)
setpar(mfrow = c(1, 2))
sl <- sort(unique(s$slope))
pt <- sort(unique(s$pert_treat))
for (k in 1:2) {
    v <- c("resist", "recover")[k]
    plot(NA,
        xlim = range(pt), ylim = range(sapply(sl, function(L) {
            sapply(pt, function(P) {
                mean(s[[v]][s$slope == L & s$pert_treat == P],
                    na.rm = TRUE
                )
            })
        }), na.rm = TRUE),
        xlab = "perturbation magnitude", ylab = c("resistance", "recovery")[k],
        main = c("resistance", "recovery")[k], cex.main = 1
    )
    for (i in seq_along(sl)) {
        m <- sapply(pt, function(P) {
            mean(s[[v]][s$slope == sl[i] & s$pert_treat == P],
                na.rm = TRUE
            )
        })
        lines(pt, m, type = "o", pch = 15 + i, col = pal[i], lwd = 2)
    }
    if (k == 1) {
        legend("topright", paste("slope", sl),
            col = pal, pch = 16:18,
            lwd = 2, bty = "n", cex = .85
        )
    }
}
dev.off()

# --- Fig 5: draw robustness -------------------------------------------------
pdf("fig5_draw_robustness.pdf", 8.5, 4.4)
setpar(mfrow = c(1, 2))
acl <- sort(unique(hs$ac))
dws <- sort(unique(hs$draw))
hd <- d[d$patch_sd > 0, ]
for (k in 1:2) {
    ylab <- c("recovery", "persistence")[k]
    plot(NA,
        xlim = c(1, length(acl)),
        ylim = if (k == 1) {
            range(sapply(dws, function(w) {
                sapply(acl, function(a) {
                    mean(hs$recover[hs$ac == a & hs$draw == w], na.rm = TRUE)
                })
            }))
        } else {
            range(sapply(dws, function(w) {
                sapply(acl, function(a) {
                    mean(hd$survived[hd$ac == a & hd$draw == w])
                })
            }))
        },
        xaxt = "n", xlab = "spatial autocorrelation", ylab = ylab,
        main = paste(ylab, "in each landscape draw"), cex.main = 1
    )
    axis(1, seq_along(acl), acl)
    for (i in seq_along(dws)) {
        m <- if (k == 1) {
            sapply(acl, function(a) {
                mean(hs$recover[hs$ac == a & hs$draw == dws[i]], na.rm = TRUE)
            })
        } else {
            sapply(acl, function(a) {
                mean(hd$survived[hd$ac == a & hd$draw == dws[i]])
            })
        }
        lines(seq_along(acl), m, type = "o", pch = 15 + i, col = pal[i], lwd = 2)
    }
    if (k == 1) {
        legend("bottomright", paste("draw", dws),
            col = pal, pch = 16:18,
            lwd = 2, bty = "n", cex = .85
        )
    }
}
dev.off()
