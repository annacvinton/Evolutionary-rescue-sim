#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes]
#' CONTENTS:
#'  - Visualisation
#'  DEPENDENCIES:
#'  - analyse_sweep.R
#' AUTHOR: [Erik Kusch, Anna Vinton]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list = ls())

## Directories ------------------------------------------------------------
Dir.Base <- getwd()
Dir.simulation <- file.path(Dir.Base, "simulation")
Dir.writeup <- file.path(Dir.Base, "writeup")
Dir.Supplement <- file.path(Dir.Base, "Supplement")
Dir.Landscapes <- file.path(Dir.simulation, "landscapes_v2")
sapply(c(Dir.simulation, Dir.writeup, Dir.Supplement, Dir.Landscapes), function(x) if (!dir.exists(x)) dir.create(x))

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, repos = "http://cran.us.r-project.org")
    }
    require(x, character.only = TRUE)
}
package_vec <- c(
    "ggplot2", # for plotting
    "viridis", # colour palettes
    "cowplot", # for arranging plots
    "patchwork", # for arranging plots
    "dplyr", # for data wrangling
    "tidyr" # for data wrangling
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
# none needed for now

## Plot and theme definition ----------------------------------------------
### colour palettes -----------
## take three colours from viridis palette for the three outcomes (resistance, recovery, persistence) from the first 75% of the palette
pal <- viridis::viridis(3, begin = 0, end = 0.75)
# pal <- c("#1B9E9E", "#E08214", "#8DA0CB")

### ggplot theme -----------
theme_set(theme_bw())
mytheme <- theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        # legend.position = c(0.08, 0.85),
        legend.background = element_blank(),
        panel.grid.minor = element_blank()
    )

# DATA ====================================================================
message("Reading in data...")
## Individual Runs --------------------------------------------------------
simsteps_df <- read.csv("all_results_v2.csv")

# PLOTTING ================================================================
message("Plotting...")

## Trajectories -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(if (length(args)) args[1] else "traj_subset_v2.csv")

TPERT <- min(d$t[d$pert_value > 0]) # when the shift fires

# patch SD 0 is flat, so its autocorrelation label is meaningless: one condition
d$cond <- ifelse(d$patch_sd == 0, "flat", paste0("ac", d$ac))
d$run <- interaction(d$draw, d$batch, d$rep, drop = TRUE)

conds <- c("flat", "ac0", "ac2", "ac4")
titles <- c(
    "flat landscape", "rough, uncorrelated", "rough, intermediate",
    "rough, correlated"
)
pal_traj <- c("#666666", pal)

## gplot counterpart
d$cond_f <- factor(d$cond, levels = conds, labels = titles)
d$n_plot <- pmax(d$n, 1)
pal_named <- setNames(pal_traj, titles)

med_n <- aggregate(n_plot ~ cond_f + t, data = d, FUN = median)
med_u <- aggregate(u_sd ~ cond_f + t, data = d, FUN = median, na.rm = TRUE)

last <- do.call(rbind, lapply(
    split(d, interaction(d$run, d$cond_f, drop = TRUE)),
    function(q) q[which.max(q$t), ]
))
persist <- aggregate(cbind(persisted = n > 10) ~ cond_f, data = last, FUN = mean)
persist$label <- sprintf("%.0f%% persist", 100 * persist$persisted)

pert_label <- data.frame(
    cond_f = factor(titles[1], levels = titles),
    t = TPERT + 20, y = 1600, lab = "perturbation"
)

p_top <- ggplot(d, aes(t, n_plot, color = cond_f)) +
    geom_vline(xintercept = TPERT, color = "#cc3333", linewidth = 0.6, linetype = 2) +
    geom_line(aes(group = run, alpha = "individual runs"), linewidth = 0.35, show.legend = TRUE) +
    geom_line(data = med_n, aes(t, n_plot, alpha = "median"), linewidth = 1) +
    geom_text(
        data = pert_label, aes(t, y, label = lab), color = "#cc3333",
        hjust = 0, size = 2.6, inherit.aes = FALSE
    ) +
    facet_wrap(~cond_f, nrow = 1) +
    scale_y_log10(limits = c(1, 2000)) +
    scale_x_continuous(limits = c(0, 810)) +
    scale_color_manual(values = pal_named, guide = "none") +
    scale_alpha_manual(
        name = NULL, values = c("individual runs" = .15, "median" = 1),
        guide = guide_legend(override.aes = list(linewidth = c(.6, 1.2), color = "grey30"))
    ) +
    labs(x = NULL, y = "population size") +
    theme_bw(base_size = 11) +
    theme(
        strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank(),
        legend.position = "inside", legend.position.inside = c(.99, .02),
        legend.justification = c(1, 0), legend.background = element_blank()
    )

p_bottom <- ggplot(d, aes(t, u_sd, color = cond_f)) +
    geom_vline(xintercept = TPERT, color = "#cc3333", linewidth = 0.6, linetype = 2) +
    geom_line(aes(group = run), alpha = .15, linewidth = 0.35) +
    geom_line(data = med_u, aes(t, u_sd), linewidth = 1) +
    geom_label(data = persist, aes(x = 750, y = 39, label = label), size = 2.6, hjust = 1) +
    facet_wrap(~cond_f, nrow = 1) +
    scale_y_continuous(limits = c(0, 42)) +
    scale_x_continuous(limits = c(0, 810)) +
    scale_color_manual(values = pal_named, guide = "none") +
    labs(x = "time", y = "trait SD") +
    theme_bw(base_size = 11) +
    theme(strip.text = element_blank(), panel.grid.minor = element_blank())

fig <- p_top / p_bottom +
    plot_annotation(
        title = "Population size (top) and trait variation (bottom) through a perturbation",
        theme = theme(plot.title = element_text(hjust = .5, size = 11, color = "#333333"))
    )
# fig
ggsave(fig, filename = file.path(Dir.writeup, "fig6_trajectories.png"), width = 9.5, height = 6, dpi = 300)

## Landscapes -------------------------------------------------------------
LS_fs <- list.files(Dir.Landscapes, pattern = ".txt")
## make dataframe that stores file names and simulation settings
LS_df <- data.frame(
    file = LS_fs,
    ac = as.numeric(gsub("ac", "", sapply(strsplit(gsub("L_|.txt", "", LS_fs), "_"), function(x) x[1]))),
    patch_sd = as.numeric(gsub("sd", "", sapply(strsplit(gsub("L_|.txt", "", LS_fs), "_"), function(x) x[2]))),
    replicate = as.numeric(gsub("r", "", sapply(strsplit(gsub("L_|.txt", "", LS_fs), "_"), function(x) x[3])))
)
#' ac is spatial autocorrelation,
#' patch_sd is the landscape heterogeneity
#' replicate is the replicate number of the landscape
LS_df$ID <- with(LS_df, paste(ac, patch_sd, sep = "_")) # add an ID for treatment to simplify later grouping and plotting

## plotting triptychs of landscapes for each combination of ac and patch_sd acrross replicates
lapply(unique(LS_df$ID), FUN = function(i) {
    sub_df <- LS_df[LS_df$ID == i, ]
    plots_ls <- lapply(sub_df$file, function(j) {
        landscape <- read.table(file.path(Dir.Landscapes, j), header = FALSE)
        landscape$x <- rep(1:100, each = 100)
        landscape$y <- rep(1:100, times = 100)
        # title <- tools::file_path_sans_ext(gsub("L_", "", j))
        p <- ggplot(landscape, aes(x = x, y = y, fill = V1)) +
            geom_tile() +
            scale_fill_viridis() +
            # labs(title = title, x = NULL, y = NULL) +
            mytheme +
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    })
    p <- plot_grid(
        ggplot() +
            theme_void() +
            labs(title = paste0("ac = ", sub_df$ac[1], ", patch SD = ", sub_df$patch_sd[1])) +
            theme(
                plot.title = element_text(
                    hjust = 0.5,
                    face = "bold",
                    size = 14
                )
            ),
        plot_grid(
            plotlist = plots_ls,
            nrow = 1,
            ncol = length(plots_ls),
            align = "h",
            labels = paste0("replicate ", sub_df$replicate),
            label_size = 12,
            label_y = 1, # move labels upward
            vjust = -0.5, # tweak vertical placement
            hjust = -2.5, # center label horizontally
            rel_widths = c(1, rep(1, length(plots_ls)))
        ),
        ncol = 1,
        rel_heights = c(0.25, 2)
    )
    ggsave(p, filename = file.path(Dir.Supplement, paste0("figS_landscape_triptych_", i, ".png")), width = 13, height = 5, dpi = 300)
})

## Model Outcomes ---------------------------------------------------------
invisible(source("analyse_sweep.R"))
#' registers:
#' - d: run-level summary of all simulations
#' - s: run-level summary of simulations where populations persisted
#' - res: model results
#' - out: effects table

### Fig 1: the phase signatures --------
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
hs <- s[s$patch_sd > 0, ]

make_phase_df <- function(dat, xvar) {
    lv <- sort(unique(dat[[xvar]]))
    out <- dat %>%
        group_by(.data[[xvar]]) %>%
        summarise(
            resist_mean = mean(resist, na.rm = TRUE),
            recover_mean = mean(recover, na.rm = TRUE),
            resist_se = sem(resist),
            recover_se = sem(recover),
            .groups = "drop"
        ) %>%
        mutate(
            resist_scaled = resist_mean / mean(resist_mean),
            recover_scaled = recover_mean / mean(recover_mean),
            resist_se_scaled = resist_se / mean(resist_mean),
            recover_se_scaled = recover_se / mean(recover_mean)
        ) %>%
        rename(level = all_of(xvar)) %>%
        pivot_longer(
            cols = c(resist_scaled, recover_scaled),
            names_to = "outcome",
            values_to = "scaled"
        ) %>%
        mutate(
            outcome = recode(outcome, resist_scaled = "resistance", recover_scaled = "recovery"),
            se = case_when(
                outcome == "resistance" ~ resist_se_scaled,
                outcome == "recovery" ~ recover_se_scaled,
                TRUE ~ NA_real_
            ),
            level = factor(level, levels = lv)
        ) %>%
        select(level, outcome, scaled, se)
    out
}

p_ac <- ggplot(
    make_phase_df(hs, "ac"),
    aes(x = level, y = scaled, colour = outcome, group = outcome)
) +
    geom_line(linewidth = 0.8, aes(linetype = outcome)) +
    geom_point(aes(shape = outcome), size = 2) +
    geom_errorbar(aes(ymin = scaled - se, ymax = scaled + se), width = 0.08) +
    geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
    scale_colour_manual(values = pal[1:2], name = "") +
    scale_shape_manual(values = c(16, 17), name = "") +
    scale_linetype_manual(values = c(1, 1), name = "") +
    labs(
        x = "Spatial Autocorrelation",
        y = "Effect (scaled to its own mean)",
        title = "by Spatial Autocorrelation"
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
        colour = guide_legend(title = "", nrow = 1),
        shape = guide_legend(title = "", nrow = 1),
        linetype = guide_legend(title = "", nrow = 1)
    )

p_sd <- ggplot(
    make_phase_df(s, "patch_sd"),
    aes(x = level, y = scaled, colour = outcome, group = outcome, pch = outcome)
) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = scaled - se, ymax = scaled + se), width = 0.08) +
    geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
    scale_colour_manual(values = pal[1:2]) +
    labs(
        x = "Patch SD",
        y = "Effect (scaled to its own mean)",
        colour = "",
        title = "by Patch SD"
    ) +
    mytheme

legend <- cowplot::get_legend(p_ac)
p1 <- plot_grid(
    plot_grid(
        p_ac + theme(legend.position = "none"),
        p_sd + theme(legend.position = "none"),
        ncol = 2
    ),
    ggplot() +
        theme_void(),
    legend,
    ncol = 1,
    rel_heights = c(1, 0.05, 0.1)
)
# p1
ggsave(p1, filename = file.path(Dir.writeup, "fig1_phase_signatures.png"), width = 8.5, height = 4.6, dpi = 300)

### Fig 2: comparison matrix --------
M <- -log10(pmax(res, 1e-17))
M[is.na(M)] <- 0

mat_df <- expand.grid(
    row = rownames(M),
    col = colnames(M),
    stringsAsFactors = FALSE
) %>%
    dplyr::mutate(
        score = as.vector(M),
        pval = as.vector(res),
        sig = dplyr::case_when(
            is.na(pval) ~ "ns",
            pval > 0.05 ~ "ns",
            pval < 1e-16 ~ "***",
            pval < 1e-4 ~ "**",
            TRUE ~ "*"
        ),
        row = factor(row, levels = rev(rownames(M))),
        col = factor(col, levels = colnames(M))
    )

sig_breaks <- c(-log10(0.05), -log10(1e-4), -log10(1e-16))
sig_labels <- c("* (p<0.05)", "** (p<1e-4)", "*** (p<1e-16)")

p2 <- ggplot(mat_df, aes(x = col, y = row, fill = score)) +
    geom_tile(width = 0.98, height = 0.98) +
    geom_label(aes(label = sig), size = 4, show.legend = FALSE, fill = "white", color = "black") +
    scale_fill_viridis_c(
        option = "viridis", begin = 0, end = 0.7, limits = c(0, max(M, na.rm = TRUE)),
        breaks = sig_breaks, labels = sig_labels,
        guide = guide_colourbar(
            barheight = unit(10, "lines"), barwidth = unit(1.2, "lines"),
            title.position = "top"
        )
    ) +
    scale_colour_manual(values = c("grey20", "white")) +
    labs(
        title = "Evidence by Outcome and Parameter",
        x = NULL,
        y = NULL,
        fill = "-log10 p"
    ) +
    mytheme +
    coord_fixed(ratio = 1) +
    theme(
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11, hjust = 0.5)
    )

ggsave(p2, filename = file.path(Dir.writeup, "fig2_comparison_matrix.png"), width = 12, height = 6, dpi = 300)


### Fig 3: dispersal --------
lv <- sort(unique(s$disp))

disp_df <- bind_rows(
    lapply(c("resist", "recover"), function(v) {
        data.frame(
            disp  = lv,
            m     = sapply(lv, function(l) mean(s[[v]][s$disp == l], na.rm = TRUE)),
            e     = sapply(lv, function(l) sem(s[[v]][s$disp == l])),
            panel = if (v == "resist") "resistance" else "recovery"
        )
    })
) %>%
    bind_rows(
        data.frame(
            disp  = lv,
            m     = sapply(lv, function(l) mean(d$survived[d$disp == l])),
            e     = sapply(lv, function(l) sem(d$survived[d$disp == l])),
            panel = "persistence"
        )
    ) %>%
    mutate(panel = factor(panel, levels = c("resistance", "recovery", "persistence")))

p3 <- ggplot(disp_df, aes(x = factor(disp), y = m, group = panel, colour = panel)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2, shape = 16) +
    geom_errorbar(aes(ymin = m - e, ymax = m + e), width = 0.1) +
    scale_colour_manual(values = pal[c(2, 1, 3)], guide = "none") +
    facet_wrap(~panel, scales = "free_y", nrow = 1, strip.position = "top") +
    labs(x = "dispersal SD", y = NULL) +
    mytheme

ggsave(p3, filename = file.path(Dir.writeup, "fig3_dispersal.png"), width = 9.5, height = 3.6, dpi = 300)

### Fig 4: slope x severity --------
sl <- sort(unique(s$slope))
pt <- sort(unique(s$pert_treat))

make_slope_df <- function(dat, v) {
    out <- bind_rows(lapply(sl, function(L) {
        data.frame(
            pert_treat = pt,
            value = sapply(pt, function(P) mean(dat[[v]][dat$slope == L & dat$pert_treat == P], na.rm = TRUE)),
            slope = L
        )
    })) %>%
        mutate(slope = factor(slope, levels = sl))
    out
}

p_resist <- ggplot(
    make_slope_df(s, "resist"),
    aes(x = pert_treat, y = value, colour = slope, shape = slope, group = slope)
) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_colour_manual(values = pal[seq_along(sl)], name = "slope") +
    scale_shape_manual(values = 15 + seq_along(sl), name = "slope") +
    labs(
        x = "perturbation magnitude",
        y = "resistance",
        title = "resistance"
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
        colour = guide_legend(title = "slope", nrow = 1),
        shape = guide_legend(title = "slope", nrow = 1)
    )

p_recover <- ggplot(
    make_slope_df(s, "recover"),
    aes(x = pert_treat, y = value, colour = slope, shape = slope, group = slope)
) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_colour_manual(values = pal[seq_along(sl)], name = "slope") +
    scale_shape_manual(values = 15 + seq_along(sl), name = "slope") +
    labs(
        x = "perturbation magnitude",
        y = "recovery",
        title = "recovery"
    ) +
    mytheme

legend <- cowplot::get_legend(p_resist)

p4 <- plot_grid(
    plot_grid(
        p_resist + theme(legend.position = "none"),
        p_recover + theme(legend.position = "none"),
        ncol = 2
    ),
    legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
)
# p4

ggsave(p4, filename = file.path(Dir.writeup, "fig4_slope_severity.png"), width = 8.5, height = 4.6, dpi = 300)

### Fig 5: draw robustness --------
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

ggsave(p5, filename = file.path(Dir.writeup, "fig5_draw_robustness.png"), width = 8.5, height = 4.6, dpi = 300)


## Boundary Map -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(file.path(Dir.Supplement, if (length(args)) args[1] else "boundary_map.csv"))
d$alive <- as.numeric(d$final_n > 10)

sds <- sort(unique(d$sd))
acs <- sort(unique(d$ac))
M_n <- tapply(d$final_n, list(d$sd, d$ac), median)
M_a <- tapply(d$alive, list(d$sd, d$ac), mean)

agg <- aggregate(alive ~ sd + ac, data = d, FUN = mean)
agg$med_n <- aggregate(final_n ~ sd + ac, data = d, FUN = median)$final_n
agg$n_lab <- ifelse(agg$med_n >= 10, sprintf("%.0f", agg$med_n), sprintf("%.2f", agg$med_n))
agg$a_lab <- sprintf("%.2f", agg$alive)

tile_w <- min(diff(sds))
tile_h <- min(diff(acs))

p_n <- ggplot(agg, aes(sd, ac, fill = med_n)) +
    geom_tile(width = tile_w, height = tile_h, color = "white") +
    geom_label(aes(label = n_lab, color = med_n > max(agg$med_n) * .55),
        size = 3, show.legend = FALSE, colour = "black", fill = "white"
    ) +
    # scale_color_manual(values = c(`FALSE` = "white", `TRUE` = "grey15"), guide = "none") +
    scale_fill_viridis_c(name = "n", end = 0.75) +
    scale_x_continuous(breaks = sds) +
    scale_y_continuous(breaks = acs) +
    labs(
        x = "Landscape Heterogeneity (patch SD)", y = "Spatial Autocorrelation",
        title = "Equilibrium Abundance (median, t = 800)"
    ) +
    mytheme +
    theme(
        legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1.5, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal"
    )
# p_n

p_a <- ggplot(agg, aes(sd, ac, fill = alive)) +
    geom_tile(width = tile_w, height = tile_h, color = "white") +
    geom_label(aes(label = a_lab, color = alive > .55), size = 3, show.legend = FALSE, fill = "white", colour = "black") +
    # scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
    # geom_contour(aes(z = alive), breaks = .5, color = "#cc3333", linewidth = .9) +
    scale_fill_gradientn(colors = pal, limits = c(0, 1), name = "fraction\npersisting   ") +
    scale_x_continuous(breaks = sds) +
    scale_y_continuous(breaks = acs) +
    labs(
        x = "landscape heterogeneity (patch SD)", y = "spatial autocorrelation",
        title = "fraction of populations persisting"
    ) +
    mytheme +
    theme(
        legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1.5, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal"
    )
# p_a

fig <- p_n + p_a
ggsave(fig, filename = file.path(Dir.Supplement, "figS_boundary_map.png"), width = 9.5, height = 4.6, dpi = 300)
