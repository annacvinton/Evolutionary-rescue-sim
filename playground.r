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

### Fig 5: draw robustness --------
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

# ---- ggplot version of Fig 5 (kept separate from the base-R plot above, mirrors
#      the "Fig 1: the phase signatures" and "Fig 4: slope x severity" plot_grid
#      patterns from viz.r) ----
# library(ggplot2)
# library(dplyr)
# library(cowplot)
#
acl <- sort(unique(hs$ac))
dws <- sort(unique(hs$draw))
hd <- d[d$patch_sd > 0, ]

make_draw_df <- function(dat, v, use_survived = FALSE) {
    out <- bind_rows(lapply(dws, function(w) {
        data.frame(
            ac = acl,
            value = if (use_survived) {
                sapply(acl, function(a) mean(dat$survived[dat$ac == a & dat$draw == w]))
            } else {
                sapply(acl, function(a) mean(dat[[v]][dat$ac == a & dat$draw == w], na.rm = TRUE))
            },
            draw = w
        )
    })) %>%
        mutate(draw = factor(draw, levels = dws))
    out
}

p_recover <- ggplot(
    make_draw_df(hs, "recover"),
    aes(x = factor(ac, levels = acl), y = value, colour = draw, shape = draw, group = draw)
) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_colour_manual(values = pal[seq_along(dws)], name = "draw") +
    scale_shape_manual(values = 15 + seq_along(dws), name = "draw") +
    labs(
        x = "spatial autocorrelation",
        y = "recovery",
        title = "recovery in each landscape draw"
    ) +
    mytheme +
    theme(
        legend.position = "top",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.text = element_text(size = 11),
        legend.key.size = unit(2, "lines"),
        legend.box = "horizontal"
    ) +
    guides(
        colour = guide_legend(title = "draw", nrow = 1),
        shape = guide_legend(title = "draw", nrow = 1)
    )

p_persist <- ggplot(
    make_draw_df(hd, "survived", use_survived = TRUE),
    aes(x = factor(ac, levels = acl), y = value, colour = draw, shape = draw, group = draw)
) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_colour_manual(values = pal[seq_along(dws)], name = "draw") +
    scale_shape_manual(values = 15 + seq_along(dws), name = "draw") +
    labs(
        x = "spatial autocorrelation",
        y = "persistence",
        title = "persistence in each landscape draw"
    ) +
    mytheme

legend <- cowplot::get_legend(p_recover)

p5 <- plot_grid(
    plot_grid(
        p_recover + theme(legend.position = "none"),
        p_persist + theme(legend.position = "none"),
        ncol = 2
    ),
    legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
)
p5

# ggsave(p5, filename = "fig5_draw_robustness_ggplot.pdf", width = 8.5, height = 4.4)
