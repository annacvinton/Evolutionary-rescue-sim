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
    "brms",
    "ggplot2",
    "tidybayes",
    "plyr",
    "cowplot",
    "pbapply",
    "betareg",
    "MuMIn",
    "ggplot2",
    "statmod",
    "numDeriv",
    "formula.tools",
    "performance",
    "dplyr",
    "tidyr",
    "ggdist",
    "MASS",
    "stringr",
    "gridExtra"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================

## Data Preparation -------------------------------------------------------
if (file.exists(file.path(Dir.Exports, "MODEL_Metrics.RData"))) {
    load(file.path(Dir.Exports, "MODEL_Metrics.RData"))
} else {
    ## Simulation Outcomes ----------------------------------------------------
    EVORES_Metrics <- readRDS(file.path(Dir.Exports, "ModelData.rds"))
    # EVORES_Metrics <- EVORES_Metrics[EVORES_Metrics$pert.name >= 9, ]
    # EVORES_Metrics$EvoRes[!EVORES_Metrics$SuffDip] <- "Insufficient Population Crash"
    # EVORES_Metrics$survival <- as.numeric(EVORES_Metrics$survival)

    ## Population-Level Summaries ---------------------------------------------
    load(file.path(Dir.Exports, "POPULATION_TimeStep.RData"))

    ### Break Points --------------
    load(file.path(Dir.Exports, "BEAST_ChangePoints.RData"))
    # head(BEAST_ChangePoints_df)

    ## Modelling Data Frame----------------------------------------------------
    EVORES_Metrics$MatchID <- with(EVORES_Metrics, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_"))
    Data_df$MatchID <- with(Data_df, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_")) # population-level metrics
    BEAST_ChangePoints_df$MatchID <- with(BEAST_ChangePoints_df, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_")) # change points
    # DISTRIBUTION_Metrics$MatchID <- with(DISTRIBUTION_Metrics, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_"))

    ## make model data frame for each breakpoint
    MatchID_vec <- unique(intersect(unique(BEAST_ChangePoints_df$MatchID), unique(EVORES_Metrics$MatchID)))
    MatchID_vec <- unique(intersect(MatchID_vec, unique(Data_df$MatchID)))

    Combination_ls <- pblapply(MatchID_vec,
        cl = parallel::detectCores(),
        FUN = function(x) {
            # x <- MatchID_vec[976]
            # print(x)
            Iter_df <- BEAST_ChangePoints_df[BEAST_ChangePoints_df$MatchID == x, ]

            if (sum(table(Iter_df$changePoint)) != 2000) {
                return(NULL)
                stop("Not all iterations have 2000 posterior samples for both breakpoints, happening for MatchID:", x)
            }

            #' 1. For each posterior sample in Iter_df, need to combine breakpoint 1 and breakpoint 2 metrics
            Cp1 <- Iter_df[Iter_df$changePoint == 1, ]
            Cp2 <- Iter_df[Iter_df$changePoint == 2, ]
            Cps <- data.frame(
                CP1_t = Cp1$posteriorSample,
                CP1_n = Cp1$n,
                CP2_t = Cp2$posteriorSample,
                CP2_n = Cp2$n,
                AC = Cp1$AC,
                DI = Cp2$DI,
                MU = Cp1$MU,
                SL = Cp1$SL,
                VA = Cp1$VA,
                pert = Cp1$pert.name,
                rep = Cp1$rep
            )

            #' 2. Adding EvoRes metrics to Cps
            Cps$pre_t <- 460
            Cps$pre_n <- EVORES_Metrics[EVORES_Metrics$MatchID == x, c("n_pre")]
            Cps$total_t <- EVORES_Metrics[EVORES_Metrics$MatchID == x, c("t")]

            #' 3. Adding population-level metrics to Cps
            Pop_iter <- Data_df[Data_df$MatchID == x, ]
            Cps_ls <- lapply(1:nrow(Cps), function(i) {
                # i = 1
                cols <- c("u_mean", "u_sd", "u_skewness", "u_kurtosis", "MalAdap_mean", "MalAdap_sd", "MalAdap_skewness", "MalAdap_kurtosis", "NNDist_mean", "NNDist_sd", "NNDist_skewness", "NNDist_kurtosis")
                Pre_df <- Pop_iter[which(Pop_iter$t == 460), cols]
                colnames(Pre_df) <- paste("pre", cols, sep = "_")
                CP1_df <- Pop_iter[which(Pop_iter$t == Cps[i, "CP1_t"]), cols]
                colnames(CP1_df) <- paste("CP1", cols, sep = "_")
                CP2_df <- Pop_iter[which(Pop_iter$t == Cps[i, "CP2_t"]), cols]
                colnames(CP2_df) <- paste("CP2", cols, sep = "_")
                cbind(Cps[i, ], Pre_df, CP1_df, CP2_df)
            })

            #' 4. return back
            return(do.call(rbind, Cps_ls))
        }
    )

    #' Combining all iterations into one data frame
    # remove NULL elements from list
    KeepElems <- sapply(Combination_ls, function(x) if (!is.null(x)) TRUE else FALSE)
    MODEL_Metrics <- do.call(rbind, Combination_ls)

    #' 6. Relative Metrics
    MODEL_Metrics$Res_Magnitude <- (MODEL_Metrics$CP1_n - MODEL_Metrics$pre_n) / MODEL_Metrics$pre_n
    MODEL_Metrics$Res_Speed <- MODEL_Metrics$CP1_t - MODEL_Metrics$pre_t
    MODEL_Metrics$Rec_Magnitude <- (MODEL_Metrics$CP2_n - MODEL_Metrics$CP1_n) / MODEL_Metrics$CP1_n
    MODEL_Metrics$Rec_Speed <- MODEL_Metrics$CP2_t - MODEL_Metrics$CP1_t
    MODEL_Metrics$Rec_MagnitudePre <- (MODEL_Metrics$CP2_n - MODEL_Metrics$pre_n) / MODEL_Metrics$pre_n
    MODEL_Metrics$Rec_SpeedPre <- MODEL_Metrics$CP2_t - MODEL_Metrics$pre_t

    #' 7. Saving final file
    save(MODEL_Metrics, file = file.path(Dir.Exports, "MODEL_Metrics.RData"))
}

#' ####################################################################### #
#' Changes to "6 - Statistical Models.R" used for the resistance/recovery
#' reanalysis. Picks up after MODEL_Metrics has been built/loaded.
#'
#' NOTE: written to match the existing script's conventions but not executed
#' in R, so please run it once and sanity-check the output before relying on it.
#' ####################################################################### #

library(dplyr)
library(sandwich) # clustered standard errors
library(lmtest) # coeftest


# ========================================================================
# CHANGE 1 (the important one) — collapse posterior draws to one row per run
# ------------------------------------------------------------------------
# MODEL_Metrics holds one row per BEAST posterior draw, so each simulation
# appears ~1000 times. Fitting lm() to all rows treats draws from the same
# run as independent observations and deflates the SEs by roughly sqrt(1000).
# Reduce each run to a single value first.
# ========================================================================

RunLevel <- MODEL_Metrics %>%
    # parameter columns may be stored as character/factor
    mutate(across(
        c(pert, rep, AC, DI, MU, SL, VA),
        ~ as.numeric(as.character(.))
    )) %>%
    group_by(pert, rep, AC, DI, MU, SL, VA) %>%
    summarise(
        n_draws = n(),
        Res_Magnitude = median(Res_Magnitude, na.rm = TRUE),
        Rec_Magnitude = median(Rec_Magnitude, na.rm = TRUE),
        Rec_Speed = median(Rec_Speed, na.rm = TRUE),
        Rec_MagnitudePre = median(Rec_MagnitudePre, na.rm = TRUE),
        CP1_t = median(CP1_t, na.rm = TRUE),
        pre_n = median(pre_n, na.rm = TRUE),
        pre_u_sd = median(pre_u_sd, na.rm = TRUE),
        pre_MalAdap_mean = median(pre_MalAdap_mean, na.rm = TRUE),
        pre_NNDist_mean = median(pre_NNDist_mean, na.rm = TRUE),
        .groups = "drop"
    )

# should be ~14,983 rows, with n_draws == 1000 throughout
nrow(RunLevel)
table(RunLevel$n_draws)


# ========================================================================
# CHANGE 2 — keep both mutation treatments
# ------------------------------------------------------------------------
# The original filtered to MU == 0. MU enters the model as a factor instead,
# since the two regimes differ qualitatively rather than only in magnitude.
#
# CHANGE 3 — drop Res_Speed
# ------------------------------------------------------------------------
# CP1 lands at t = 470 in ~99.98% of runs against the fixed t = 460
# reference, so Res_Speed is pinned at the sampling interval.
# ========================================================================

# quick check of the Res_Speed problem
mean(RunLevel$CP1_t == 470)

Dat <- RunLevel %>% filter(pert >= 9, pert <= 13)


# ========================================================================
# CHANGE 4 — centre predictors, transform ratio outcomes, define clusters
# ------------------------------------------------------------------------
# Centring makes main effects interpretable at the middle of the design
# rather than at zero. Ratio outcomes are right-skewed, so log-transform.
# ========================================================================

Dat <- Dat %>%
    mutate(
        SL_c   = SL - 1,
        pert_c = pert - 11,
        AC_c   = AC - 2,
        VA_s   = log(VA) - log(50),
        DI_c   = DI - 1.75,
        y_rec  = log1p(Rec_Magnitude),
        y_spd  = log(pmax(Rec_Speed, 10)),
        y_pre  = log1p(pmax(Rec_MagnitudePre, -0.999)),
        cell   = interaction(AC, DI, MU, SL, VA, pert, drop = TRUE)
    )


# ========================================================================
# CHANGE 5 (the other important one) — one interaction model, not a loop
# ------------------------------------------------------------------------
# The original fitted lm(Outcome ~ 0 + Driver) separately at each
# perturbation level. That gives group means but cannot test whether an
# effect *changes* with severity, which is the actual question. Fit one
# model with severity as a predictor interacted with the driver.
#
# Standard errors are clustered by parameter combination, since runs
# sharing a combination are not independent.
# ========================================================================

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


# ========================================================================
# ADDITION 4 (optional) — propagate change-point uncertainty properly
# ------------------------------------------------------------------------
# Rather than discarding the posterior draws, fit the model once per draw
# and combine with Rubin's rules. Total variance T = W + (1 + 1/M) * B,
# so the intervals get *wider*, which is the point.
# ========================================================================

MODEL_Metrics <- MODEL_Metrics %>%
    mutate(across(c(pert, rep, AC, DI, MU, SL, VA), ~ as.numeric(as.character(.)))) %>%
    group_by(pert, rep, AC, DI, MU, SL, VA) %>%
    mutate(draw = row_number()) %>%
    ungroup()

M <- 100
sel <- sample(seq_len(1000), M)
est_lst <- vector("list", M)
se_lst <- vector("list", M)

for (i in seq_along(sel)) {
    d <- MODEL_Metrics %>%
        filter(draw == sel[i], pert >= 9, pert <= 13) %>%
        mutate(
            SL_c = SL - 1, pert_c = pert - 11, AC_c = AC - 2,
            VA_s = log(VA) - log(50), DI_c = DI - 1.75,
            y_rec = log1p(pmax(Rec_Magnitude, -0.999)),
            cell = interaction(AC, DI, MU, SL, VA, pert, drop = TRUE)
        )
    m <- lm(y_rec ~ SL_c * pert_c * MU + AC_c + VA_s + DI_c, data = d)
    ct <- coeftest(m, vcov. = vcovCL(m, cluster = d$cell))
    est_lst[[i]] <- ct[, "Estimate"]
    se_lst[[i]] <- ct[, "Std. Error"]
}

est_mat <- do.call(rbind, est_lst)
se_mat <- do.call(rbind, se_lst)

theta <- colMeans(est_mat) # pooled estimate
W <- colMeans(se_mat^2) # within-draw variance
B <- apply(est_mat, 2, var) # between-draw variance
Tvar <- W + (1 + 1 / M) * B # Rubin total variance
se <- sqrt(Tvar)

data.frame(
    estimate = theta,
    se = se,
    z = theta / se,
    p = 2 * (1 - pnorm(abs(theta / se)))
)
