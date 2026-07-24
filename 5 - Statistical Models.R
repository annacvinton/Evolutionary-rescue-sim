#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes]
#' CONTENTS:
#'  - Statistical Models of Evolutionary Rescue Moments
#'  DEPENDENCIES:
#'  - BEAST_ChangePoints.RData produced by "5 - BEAST-TimeSeries.r"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list = ls())

## Directories ------------------------------------------------------------
### Define dicrectories in relation to project directory
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "Data")
Dir.Exports <- file.path(Dir.Base, "Exports")
### Create directories which aren't present yet
Dirs <- c(Dir.Data, Dir.Exports)
CreateDir <- sapply(Dirs, function(x) if (!dir.exists(x)) dir.create(x))

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, repos = "http://cran.us.r-project.org")
    }
    require(x, character.only = TRUE)
}
package_vec <- c(
    "dplyr", "sandwich", "lmtest",
    "tidyr", "ggplot2"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
load(file.path(Dir.Exports, "MODELMetrics.RData"))
colnames(MODEL_Metrics) <- gsub("pert.name", "pert", colnames(MODEL_Metrics))

#' first filtering of data to
#' 1. remove runs that did not survive the perturbation
#' 2. keep only the last 5 perturbation levels (9-13)
#' 3. remove runs with missing change-point estimates
RunLevel <- MODEL_Metrics %>%
    # parameter columns may be stored as character/factor
    mutate(across(
        c(pert, rep, AC, DI, MU, SL, VA),
        ~ as.numeric(as.character(.))
    )) %>%
    filter(pert >= 9, pert <= 13, survival == TRUE, !is.na(CP1_n), !is.na(CP2_n))

# should be ~14,762 rows
nrow(RunLevel)
# CP1 lands at t = 470 in ~99.96% of runs which si the harcorded time-step at which the perturbation is introduced, so Res_Speed is pinned at the sampling interval.
mean(RunLevel$CP1_t == 470, na.rm = TRUE)


#' centre predictors, transform ratio outcomes, define clusters following https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2010.00012.x
#' Centring makes main effects interpretable at the middle of the design rather than at zero. Ratio outcomes are right-skewed, so log-transform.
Dat <- RunLevel %>%
    mutate(
        SL_c   = SL - 1,
        pert_c = pert - 11,
        AC_c   = AC - 2,
        VA_s   = log(VA) - log(50),
        DI_c   = DI - 1.75,
        y_rec  = log1p(Rec_Magnitude),
        y_spd  = log(pmax(Rec_Speed, 10)),
        y_pre  = log1p(pmax(Rec_MagnitudePre, -0.999)),
        cell   = interaction(AC, DI, MU, SL, VA, pert, drop = TRUE) # essential the same as the ID column, but removing the "rep" part, which is not relevant for clustering standard errors
    )


# MODELS ==================================================================
#' one interaction model, not a loop
# The original fitted lm(Outcome ~ 0 + Driver) separately at each
# perturbation level. That gives group means but cannot test whether an
# effect *changes* with severity, which is the actual question. Fit one
# model with severity as a predictor interacted with the driver.
#
# Standard errors are clustered by parameter combination, since runs
# sharing a combination are not independent.

fit_clustered <- function(formula, data) {
    m <- lm(formula, data = data)
    ct <- coeftest(m, vcov. = vcovCL(m, cluster = data$cell))
    list(model = m, coefs = ct)
}

# --- main models: does the slope effect depend on severity and mutation? ---
m_res <- fit_clustered(Res_Magnitude ~ SL_c * pert_c * MU + AC_c + VA_s + DI_c, Dat)
m_rec <- fit_clustered(y_rec ~ SL_c * pert_c * MU + AC_c + VA_s + DI_c, Dat)
m_spd <- fit_clustered(y_spd ~ SL_c * pert_c * MU + AC_c + VA_s + DI_c, Dat)
m_pre <- fit_clustered(y_pre ~ SL_c * pert_c * MU + AC_c + VA_s + DI_c, Dat)

m_rec$coefs # SL_c:pert_c:MU is the three-way term of interest

# --- same models fitted separately within each mutation regime ---
for (mu in c(1, 0)) {
    d <- filter(Dat, MU == mu)
    cat("\n---- MU =", mu, "----\n")
    print(fit_clustered(y_rec ~ SL_c * pert_c + AC_c + VA_s + DI_c, d)$coefs)
}


# ========================================================================
# ADDITION 1 — full interaction sweep for autocorrelation and variance
# ------------------------------------------------------------------------
# Tests AC and VA against each other, against slope, and against severity.
# The AC x severity term is the one that comes out non-zero, and only
# where new mutations are available.
# ========================================================================

sweep_form <- ~ AC_c * VA_s + AC_c * SL_c + AC_c * pert_c +
    VA_s * SL_c + VA_s * pert_c + SL_c * pert_c + DI_c

for (mu in c(1, 0)) {
    d <- filter(Dat, MU == mu)
    cat("\n---- AC/VA sweep, MU =", mu, "----\n")
    print(fit_clustered(update(sweep_form, y_rec ~ .), d)$coefs)
}


# ========================================================================
# ADDITION 2 — shape of the AC effect (monotonic vs humped vs bimodal)
# ------------------------------------------------------------------------
# With 5 levels a linear term cannot see a hump, so compare it against a
# free-shape (factor) fit. Same logic for VA, which has 3 levels.
# ========================================================================

S <- filter(Dat, MU == 1)

ac_lin <- lm(y_rec ~ AC_c * pert_c + SL_c * pert_c + VA_s + DI_c, data = S)
ac_fac <- lm(y_rec ~ factor(AC) * pert_c + SL_c * pert_c + VA_s + DI_c, data = S)
anova(ac_lin, ac_fac) # non-significant => linear is sufficient

va_lin <- lm(y_rec ~ VA_s * pert_c + SL_c * pert_c + AC_c + DI_c, data = S)
va_fac <- lm(y_rec ~ factor(VA) * pert_c + SL_c * pert_c + AC_c + DI_c, data = S)
anova(va_lin, va_fac) # non-significant => no humped variance effect


# ========================================================================
# ADDITION 3 — is the resistance/recovery trade-off a measurement artefact?
# ------------------------------------------------------------------------
# Rec_Magnitude is a ratio off the crash low, so a population that falls
# less has less room to climb. If the trade-off were only that, controlling
# for crash depth would remove it. Rec_MagnitudePre carries no such link.
# ========================================================================

a <- fit_clustered(y_rec ~ SL_c * pert_c + AC_c + VA_s + DI_c + MU, Dat)
b <- fit_clustered(y_rec ~ SL_c * pert_c + Res_Magnitude + AC_c + VA_s + DI_c + MU, Dat)
c(uncontrolled = coef(a$model)["SL_c"], controlled = coef(b$model)["SL_c"])


# # ========================================================================
# # ADDITION 4 (optional) — propagate change-point uncertainty properly; cannot do this with the new setup where posterior draws are already combined into a single dataset, but could be done if the posterior draws were kept separate.
# # ------------------------------------------------------------------------
# # Rather than discarding the posterior draws, fit the model once per draw
# # and combine with Rubin's rules. Total variance T = W + (1 + 1/M) * B,
# # so the intervals get *wider*, which is the point.
# # ========================================================================

# MODEL_Metrics <- MODEL_Metrics %>%
#     mutate(across(c(pert, rep, AC, DI, MU, SL, VA), ~ as.numeric(as.character(.)))) %>%
#     group_by(pert, rep, AC, DI, MU, SL, VA) %>%
#     mutate(draw = row_number()) %>%
#     ungroup()

# M <- 100
# sel <- sample(seq_len(1000), M)
# est_lst <- vector("list", M)
# se_lst <- vector("list", M)

# for (i in seq_along(sel)) {
#     d <- MODEL_Metrics %>%
#         filter(draw == sel[i], pert >= 9, pert <= 13) %>%
#         mutate(
#             SL_c = SL - 1, pert_c = pert - 11, AC_c = AC - 2,
#             VA_s = log(VA) - log(50), DI_c = DI - 1.75,
#             y_rec = log1p(pmax(Rec_Magnitude, -0.999)),
#             cell = interaction(AC, DI, MU, SL, VA, pert, drop = TRUE)
#         )
#     m <- lm(y_rec ~ SL_c * pert_c * MU + AC_c + VA_s + DI_c, data = d)
#     ct <- coeftest(m, vcov. = vcovCL(m, cluster = d$cell))
#     est_lst[[i]] <- ct[, "Estimate"]
#     se_lst[[i]] <- ct[, "Std. Error"]
# }

# est_mat <- do.call(rbind, est_lst)
# se_mat <- do.call(rbind, se_lst)

# theta <- colMeans(est_mat) # pooled estimate
# W <- colMeans(se_mat^2) # within-draw variance
# B <- apply(est_mat, 2, var) # between-draw variance
# Tvar <- W + (1 + 1 / M) * B # Rubin total variance
# se <- sqrt(Tvar)

# data.frame(
#     estimate = theta,
#     se = se,
#     z = theta / se,
#     p = 2 * (1 - pnorm(abs(theta / se)))
# )

# VISUALISATIONS ==========================================================
.fig_results <- character(0)
CHK <- function(label, observed, expected, tol = 0.01) {
    obs <- as.numeric(observed)
    exp <- as.numeric(expected)
    ok <- !is.na(obs) && abs(obs - exp) <= tol * max(1, abs(exp))
    cat(sprintf(
        "    %-42s obs %11.4f   exp %11.4f   %s\n",
        label, obs, exp, if (ok) "PASS" else "** FAIL **"
    ))
    .fig_results <<- c(.fig_results, if (ok) "PASS" else label)
}
HDR <- function(x) cat("\n", strrep("=", 74), "\n ", x, "\n", strrep("=", 74), "\n", sep = "")

Dir.Figs <- if (exists("Dir.Exports")) file.path(Dir.Exports, "Figures") else "Figures"
if (!dir.exists(Dir.Figs)) dir.create(Dir.Figs, recursive = TRUE)

# shared styling ----------------------------------------------------------
sev_cols <- setNames(hcl.colors(5, "viridis"), as.character(9:13))
mu_cols <- c("1" = "#1b7837", "0" = "#c0392b")
mu_labs <- c("1" = "With new mutations", "0" = "Existing variation only")
theme_paper <- theme_minimal(base_size = 9) +
    theme(
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        strip.text.y = element_text(face = "bold", size = 8),
        legend.position = "bottom", legend.key.size = unit(0.4, "cm")
    )

# ========================================================================
HDR("STEP 0 — check the coefficient names marginal_effect() relies on")
# ========================================================================
.chk <- lm(y_rec ~ AC_c * pert_c + SL_c * pert_c + VA_s + DI_c, data = filter(Dat, MU == 1))
cat("  interaction terms present:\n    ",
    paste(grep(":", names(coef(.chk)), value = TRUE), collapse = " | "), "\n",
    sep = ""
)
cat("  marginal_effect() looks up 'AC_c:pert_c' and 'SL_c:pert_c'.\n")
cat("  If they appear reversed above (e.g. 'pert_c:AC_c'), change the\n",
    "  `inter =` arguments in the Figure 4/5 calls to match.\n",
    sep = ""
)

# panel-data helper for Figures 1 and 3 -----------------------------------
fig_panel_data <- function(data, xvar) {
    data %>%
        select(MU, pert, all_of(xvar),
            `Crash depth` = Res_Magnitude,
            `Recovery magnitude` = Rec_Magnitude,
            `Recovery time` = Rec_Speed
        ) %>%
        pivot_longer(c(`Crash depth`, `Recovery magnitude`, `Recovery time`),
            names_to = "measure", values_to = "value"
        ) %>%
        group_by(MU, pert, measure, x = .data[[xvar]]) %>%
        summarise(med = median(value, na.rm = TRUE), n = n(), .groups = "drop") %>%
        filter(n >= 5) %>%
        mutate(
            measure = factor(measure,
                levels = c("Crash depth", "Recovery magnitude", "Recovery time")
            ),
            MU_label = factor(mu_labs[as.character(MU)], levels = unname(mu_labs)),
            pert_f = factor(pert, levels = 9:13)
        )
}

# ========================================================================
HDR("FIGURE 1 — outcomes across gradient slope")
# ========================================================================
d1 <- fig_panel_data(Dat, "SL")
cat("  rows feeding plot:", nrow(d1), "\n")
gv <- function(df, mu, meas, xx) df$med[df$MU == mu & df$measure == meas & df$pert == 9 & abs(df$x - xx) < 1e-6][1]
CHK("MU=1 SL0.8 sev9 crash depth", gv(d1, 1, "Crash depth", 0.8), -0.892)
CHK("MU=0 SL0.8 sev9 crash depth", gv(d1, 0, "Crash depth", 0.8), -0.921)
CHK("MU=1 SL1.2 sev9 recovery magnitude", gv(d1, 1, "Recovery magnitude", 1.2), 3.364)

p1 <- ggplot(d1, aes(x, med, colour = pert_f, group = pert_f)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.2) +
    # facet_grid(MU_label ~ measure, scales = "free_y", switch = "y") +
    facet_wrap(~ MU_label + measure, scales = "free_y") +
    scale_colour_manual(values = sev_cols, name = "Severity") +
    scale_x_continuous(breaks = c(0.8, 1.0, 1.2)) +
    labs(x = "Gradient slope", y = NULL) +
    theme_paper
ggsave(file.path(Dir.Figs, "Figure1_slope.png"), p1, width = 8.2, height = 4.8, dpi = 300)
cat("  saved Figure1_slope.png\n")


# ========================================================================
HDR("FIGURE 2 — trait variation, population size, completeness of recovery")
# ========================================================================
d2 <- Dat %>%
    select(SL,
        `Trait SD before shift` = pre_u_sd,
        `Population size before shift` = pre_n,
        `Recovery relative to original` = Rec_MagnitudePre
    ) %>%
    pivot_longer(-SL, names_to = "measure", values_to = "value") %>%
    group_by(measure, SL) %>%
    summarise(
        med = median(value, na.rm = TRUE),
        q1 = quantile(value, .25, na.rm = TRUE),
        q3 = quantile(value, .75, na.rm = TRUE), .groups = "drop"
    ) %>%
    mutate(measure = factor(measure, levels = c(
        "Trait SD before shift", "Population size before shift",
        "Recovery relative to original"
    )))
gm <- function(meas, s) d2$med[d2$measure == meas & abs(d2$SL - s) < 1e-6]
CHK("trait SD at SL0.8", gm("Trait SD before shift", 0.8), 23.7, tol = .02)
CHK("trait SD at SL1.2", gm("Trait SD before shift", 1.2), 36.9, tol = .02)
CHK("population size at SL0.8", gm("Population size before shift", 0.8), 1287, tol = .02)
CHK("recovery vs original at SL1.2", gm("Recovery relative to original", 1.2), -0.677, tol = .03)

p2 <- ggplot(d2, aes(SL, med)) +
    geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.03, linewidth = 0.4, colour = "#31688e") +
    geom_line(linewidth = 0.6, colour = "#31688e") +
    geom_point(size = 2, colour = "#31688e") +
    facet_wrap(~measure, scales = "free_y") +
    scale_x_continuous(breaks = c(0.8, 1.0, 1.2)) +
    labs(x = "Gradient slope", y = NULL) +
    theme_paper +
    theme(legend.position = "none")
ggsave(file.path(Dir.Figs, "Figure2_mechanism.png"), p2, width = 8.4, height = 2.9, dpi = 300)
cat("  saved Figure2_mechanism.png\n")


# ========================================================================
HDR("FIGURE 3 — outcomes across spatial clustering")
# ========================================================================
d3 <- fig_panel_data(Dat, "AC")
cat("  rows feeding plot:", nrow(d3), "\n")
g3 <- function(mu, meas, sev, xx) d3$med[d3$MU == mu & d3$measure == meas & d3$pert == sev & abs(d3$x - xx) < 1e-6][1]
CHK("MU=1 AC0 sev9 recovery magnitude", g3(1, "Recovery magnitude", 9, 0), 5.556)
CHK("MU=1 AC4 sev13 recovery magnitude", g3(1, "Recovery magnitude", 13, 4), 3.353)

p3 <- ggplot(d3, aes(x, med, colour = pert_f, group = pert_f)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.2) +
    # facet_grid(MU_label ~ measure, scales = "free_y", switch = "y") +
    facet_wrap(~ MU_label + measure, scales = "free_y") +
    scale_colour_manual(values = sev_cols, name = "Severity") +
    scale_x_continuous(breaks = 0:4) +
    labs(x = "Spatial clustering (autocorrelation)", y = NULL) +
    theme_paper
ggsave(file.path(Dir.Figs, "Figure3_clustering.png"), p3, width = 8.2, height = 4.8, dpi = 300)
cat("  saved Figure3_clustering.png\n")


# marginal-effect helper for Figures 4 and 5 ------------------------------
marginal_effect <- function(model, vc, main, inter, mult = 1, sev = seq(9, 13, by = 0.1)) {
    b <- coef(model)
    d <- sev - 11
    est <- (b[[main]] + b[[inter]] * d) * mult
    v <- (vc[main, main] + d^2 * vc[inter, inter] + 2 * d * vc[main, inter]) * mult^2
    data.frame(pert = sev, est = est, lo = est - 1.96 * sqrt(v), hi = est + 1.96 * sqrt(v))
}
effect_by_regime <- function(formula, main, inter, mult) {
    bind_rows(lapply(c(1, 0), function(mu) {
        d <- filter(Dat, MU == mu)
        m <- lm(formula, data = d)
        marginal_effect(m, vcovCL(m, cluster = d$cell), main, inter, mult) %>%
            mutate(MU = as.character(mu))
    }))
}
plot_effect <- function(df, ylab) {
    ggplot(df, aes(pert, est, colour = MU, fill = MU)) +
        geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
        geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
        geom_line(linewidth = 0.8) +
        scale_colour_manual(values = mu_cols, labels = mu_labs, name = NULL) +
        scale_fill_manual(values = mu_cols, labels = mu_labs, name = NULL) +
        scale_x_continuous(breaks = 9:13) +
        labs(x = "Severity of environmental shift", y = ylab) +
        theme_paper
}

# ========================================================================
HDR("FIGURE 4 — clustering effect across severity (AC 0 -> AC 4)")
# ========================================================================
d4 <- effect_by_regime(y_rec ~ AC_c * pert_c + SL_c * pert_c + VA_s + DI_c,
    main = "AC_c", inter = "AC_c:pert_c", mult = 4
)
e4 <- function(mu, sev) d4$est[d4$MU == mu & abs(d4$pert - sev) < 1e-6][1]
CHK("MU=1 clustering effect at sev9", e4(1, 9), -0.0645)
CHK("MU=1 clustering effect at sev11", e4(1, 11), 0.0624)
CHK("MU=1 clustering effect at sev13", e4(1, 13), 0.1893)
p4 <- plot_effect(d4, "Effect of fully clustered vs.\nunclustered landscape on recovery")
ggsave(file.path(Dir.Figs, "Figure4_clustering_effect.png"), p4, width = 4.4, height = 3.2, dpi = 300)
cat("  saved Figure4_clustering_effect.png\n")


# ========================================================================
HDR("FIGURE 5 — slope effect across severity (+0.2 slope)")
# ========================================================================
d5 <- effect_by_regime(y_rec ~ SL_c * pert_c + AC_c + VA_s + DI_c,
    main = "SL_c", inter = "SL_c:pert_c", mult = 0.2
)
e5 <- function(mu, sev) d5$est[d5$MU == mu & abs(d5$pert - sev) < 1e-6][1]
CHK("MU=1 slope effect at sev9", e5(1, 9), -0.3522)
CHK("MU=1 slope effect at sev13", e5(1, 13), -0.4532)
CHK("MU=0 slope effect at sev9", e5(0, 9), -0.1257)
CHK("MU=0 slope effect at sev13", e5(0, 13), 0.2413)
p5 <- plot_effect(d5, "Effect of +0.2 gradient slope on\nlog(1 + recovery magnitude)")
ggsave(file.path(Dir.Figs, "Figure5_slope_effect.png"), p5, width = 4.4, height = 3.2, dpi = 300)
cat("  saved Figure5_slope_effect.png\n")


# ========================================================================
HDR("SUMMARY")
# ========================================================================
fails <- setdiff(.fig_results, "PASS")
cat(sprintf(
    "\n  %d figure checks run, %d passed, %d failed\n",
    length(.fig_results), sum(.fig_results == "PASS"), length(fails)
))
if (length(fails)) {
    cat("\n  FAILED:\n")
    for (f in fails) cat("   -", f, "\n")
} else {
    cat("\n  All figure checks passed, and 5 PNGs written to ", Dir.Figs, "\n", sep = "")
}
