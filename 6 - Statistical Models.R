#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes]
#' CONTENTS:
#'  - Statistical Models of
#'    - Survival Likelihood
#'    - Evolutionary Rescue Likelihood
#'    - Evolutionary Rescue Success
#'  DEPENDENCIES:
#'  - EVORES_Metrics.RData produced by "4 - Evolutionary Rescue Success Metrics.R"
#'  - DISTRIBUTIONS_Spatial.RData produced by "5 - DistributionComparison.R"
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
    "MASS"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
if (file.exists(file.path(Dir.Exports, "MODEL_Metrics.RData"))) {
    load(file.path(Dir.Exports, "MODEL_Metrics.RData"))
} else {
    ## Simulation Outcomes ----------------------------------------------------
    EVORES_Metrics <- readRDS(file.path(Dir.Exports, "ModelData.rds"))
    EVORES_Metrics <- EVORES_Metrics[EVORES_Metrics$pert.name >= 9, ]
    # EVORES_Metrics$EvoRes[!EVORES_Metrics$SuffDip] <- "Insufficient Population Crash"
    # EVORES_Metrics$survival <- as.numeric(EVORES_Metrics$survival)

    ## Population-Level Summaries ---------------------------------------------
    load(file.path(Dir.Exports, "POPULATION_TimeStep.RData"))

    ## Distribution Comparisons -----------------------------------------------
    ### Non-Spatially Explicit ----
    load(file.path(Dir.Exports, "DISTRIBUTIONS_NonSpatial.RData"))
    DISTRIBUTIONS_NonSpatial <- DISTRIBUTIONS_NonSpatial[DISTRIBUTIONS_NonSpatial$Pert >= 9, ]

    ### Spatially Explicit ----
    load(file.path(Dir.Exports, "DISTRIBUTIONS_Spatial.RData"))
    DISTRIBUTIONS_Spatial <- DISTRIBUTIONS_Spatial[["Summary_df"]]
    rownames(DISTRIBUTIONS_Spatial) <- c()
    DISTRIBUTIONS_Spatial <- DISTRIBUTIONS_Spatial[DISTRIBUTIONS_Spatial$Pert >= 9, ]

    ### Fusing ----
    DISTRIBUTION_Metrics <- base::merge(DISTRIBUTIONS_Spatial, DISTRIBUTIONS_NonSpatial)
    DISTRIBUTION_Metrics$Comparisons <- paste(
        paste("[Environment]", DISTRIBUTION_Metrics$Envir),
        "vs.",
        paste("[Individuals]", DISTRIBUTION_Metrics$Indivs)
    )
    DISTRIBUTION_Metrics$Comparisons <-
        factor(DISTRIBUTION_Metrics$Comparisons,
            levels = c(
                "[Environment] Pre vs. [Individuals] Pre",
                "[Environment] Post vs. [Individuals] PostMin",
                "[Environment] Post vs. [Individuals] PostMax",
                "[Environment] Pre vs. [Individuals] PostMin",
                "[Environment] Pre vs. [Individuals] PostMax"
            )
        )
    colnames(DISTRIBUTION_Metrics)[4] <- "pert.name"

    ## Modelling Data Frame----------------------------------------------------
    EVORES_Metrics$MatchID <- with(EVORES_Metrics, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_"))
    Data_df$MatchID <- with(Data_df, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_")) # population-level metrics
    DISTRIBUTION_Metrics$MatchID <- with(DISTRIBUTION_Metrics, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_"))

    ## attach population and distribution metrics to EVORES_Metrics
    Combination_ls <- pblapply(1:nrow(EVORES_Metrics),
        # cl = parallel::detectCores(),
        FUN = function(x) {
            # x <- 15011
            # print(x)
            Iter_df <- EVORES_Metrics[x, c("t_minpost", "t_maxpost", "MatchID")]

            ## adding population metrics
            Data_iter <- Data_df[Data_df$MatchID == Iter_df$MatchID, ]
            t_vec <- c(pre = 460, postmin = Iter_df$t_minpost, postmax = Iter_df$t_maxpost)
            Data_iter <- Data_iter[match(t_vec, Data_iter$t), ]
            fuse_df <- lapply(1:nrow(Data_iter), FUN = function(y) {
                ret_df <- Data_iter[y, colnames(Data_df) %nin% c("t", "pert.name", "rep", "AC", "DI", "MU", "SL", "VA", "ID", "MatchID")]
                colnames(ret_df) <- paste(names(t_vec)[y], colnames(ret_df), sep = "_")
                ret_df
            })
            Pop_add <- cbind(
                do.call(cbind, fuse_df)
                # , "MatchID" = Data_iter$MatchID[1]
            )

            Dist_iter <- DISTRIBUTION_Metrics[DISTRIBUTION_Metrics$MatchID == Iter_df$MatchID, ]
            if (nrow(Dist_iter) == 6) {
                print(x)
                Dist_iter <- Dist_iter[c(1, 3, 5), ]
            } # happens for three runs that we have two dsitribution overlap comparisons
            if (nrow(Dist_iter) == 0) {
                message(x)
                Dist_iter[1:5, ] <- NA
                Dist_iter$Comparisons <- rev(levels(DISTRIBUTION_Metrics$Comparisons))
            }
            fuse_df <- lapply(1:nrow(Dist_iter), FUN = function(y) {
                ret_df <- Dist_iter[y, colnames(Dist_iter) %nin% c("t", "pert.name", "rep", "AC", "DI", "MU", "SL", "VA", "ID", "MatchID", "Comparisons")]
                colnames(ret_df) <- paste(as.character(Dist_iter$Comparisons)[y], colnames(ret_df), sep = "_")
                ret_df
            })
            Dist_add <- cbind(do.call(cbind, fuse_df), "MatchID" = Data_iter$MatchID[1])

            cbind(Pop_add, Dist_add)
        }
    )
    Combination_df <- do.call(rbind, Combination_ls)
    MODEL_Metrics <- join(EVORES_Metrics, Combination_df)
    save(MODEL_Metrics, file = file.path(Dir.Exports, "MODEL_Metrics.RData"))
}
# ## express times between abundances as percentages
# MODEL_Metrics$perc_t_maxpost <- (MODEL_Metrics$t_maxpost - MODEL_Metrics$t_minpost) / (MODEL_Metrics$t - MODEL_Metrics$t_minpost)
# MODEL_Metrics$perc_t_minpost <- (MODEL_Metrics$t_minpost - 450) / 450
# ## Make recovery and resistance outcomes
# Bayes_df <- MODEL_Metrics[MODEL_Metrics$survival == 1, ] # select only survived runs, others cannot experience evores
# Bayes_df$Recovery <- (Bayes_df$perc_maxpostpre + Bayes_df$perc_t_maxpost) / 2 # recovery is mean of percentage change in abundance and time it took to get there
# Bayes_df$Resistance <- (Bayes_df$perc_minpost * Bayes_df$perc_t_minpost) / 2 # same as above

Bayes_df <- MODEL_Metrics[MODEL_Metrics$DI == 2 & MODEL_Metrics$pert.name >= 9 & MODEL_Metrics$SL == "1.0", ]
colnames(Bayes_df) <- gsub(pattern = " ", replacement = "_", colnames(Bayes_df))

## subset for MU
Bayes_df_mu0 <- Bayes_df[Bayes_df$MU == 0, ]
Bayes_df_mu1 <- Bayes_df[Bayes_df$MU == 1, ]
Data_ls <- list(Mu0 = Bayes_df_mu0, Mu1 = Bayes_df_mu1)

# MODELS ==================================================================

## Initial Models of Simulation Settings on Primary Outcomes --------------
Outcomes_Primary <- c(
    "Res_Mag_Immediate",
    "Res_Mag_Delayed",
    "Res_Speed_Delayed",
    "Rec_Mag_PreBase",
    "Rec_Speed_PreBase",
    "Rec_Mag_PostBase",
    "Rec_Speed_PostBase"
)

InitialModels_ls <- lapply(names(Data_ls), FUN = function(Mu_iter) {
    # Mu_iter <- "Mu0"
    Data_iter <- Data_ls[[Mu_iter]]
    Outcome_draws <- lapply(Outcomes_Primary, FUN = function(Outcome) {
        # Outcome <- "Res_Mag_Immediate"
        colnames(Data_iter)[colnames(Data_iter) == Outcome] <- "Outcome"
        print(Outcome)
        # print(summary(Data_iter$Outcome))
        # Pert_draws <- lapply(sort(unique(Data_iter$pert.name)), FUN = function(pert) {
        #     print(pert)
        # pert <- 14
        # Data_iter <- Data_iter[Data_iter$pert.name == pert, ]
        if (grepl(Outcome, pattern = "Speed") | grepl(Outcome, pattern = "Res")) {
            Data_iter <- Data_iter[Data_iter$Outcome <= 1, ] # these are definitive outliers
            fit <- betareg(Outcome ~ pert.name + AC * VA, data = Data_iter)
            R2 <- summary(fit)$pseudo.r.squared
        } else {
            fit <- glm(Outcome ~ pert.name + AC * VA, data = Data_iter, family = Gamma(link = "log"))
            ll_full <- logLik(fit)
            ll_null <- logLik(glm(Outcome ~ 1, data = Data_iter, family = Gamma(link = "log")))
            R2 <- 1 - as.numeric(ll_full / ll_null)
        }
        # Extract coefficients and covariance matrix
        coef_est <- coef(fit)
        vcov_mat <- vcov(fit)
        # Simulate 1000 draws from multivariate normal
        set.seed(42)
        sim_coefs <- MASS::mvrnorm(5000, mu = coef_est, Sigma = vcov_mat) %>%
            as.data.frame()
        # Convert to long format for ggdist
        sim_long <- sim_coefs %>%
            pivot_longer(everything(), names_to = "term", values_to = "value")

        if (class(fit)[1] != "betareg") { # adjustment needed for log link
            sim_long <- sim_long %>%
                mutate(value = exp(value))
        }
        list(coeffs = sim_long, R2 = R2)
        # })
        # names(Pert_draws) <- sort(unique(Data_iter$pert.name))
        # Pert_draws
    })
    names(Outcome_draws) <- Outcomes_Primary
    Outcome_draws
})
names(InitialModels_ls)

lapply(InitialModels_ls[[1]], "[[", "R2")
lapply(InitialModels_ls[[2]], "[[", "R2")


# ## Logistic Models of Evolutionary Rescue -----------------------
# Logit_df <- Bayes_df[Bayes_df$SuffDip != "Insufficient Population Crash", ]
# Logit_df$EvoRes <- Logit_df$EvoRes == "TRUE"
# logit_model <- glm(EvoRes ~ AC + DI + MU + SL + VA, data = Logit_df, family = binomial)
# summary(logit_model)

# ## Continuous Models of Resistance and Recovery -----------------------
# ### Model Building Functions ------------
# generate_valid_subsets <- function(vec) {
#     base_terms <- vec[!grepl(":", vec)] # Get individual terms
#     all_subsets <- list()

#     for (term in vec) {
#         new_vec <- setdiff(vec, term) # Remove one term

#         # Identify removed base terms
#         removed_base <- setdiff(base_terms, new_vec)

#         # Remove any interaction terms that contain removed base terms
#         valid_vec <- new_vec[!sapply(new_vec, function(x) any(strsplit(x, ":")[[1]] %in% removed_base))]

#         all_subsets <- append(all_subsets, list(valid_vec))
#     }

#     return(all_subsets)
# }

# do_model_selection <- function(model, terms, data = Bayes_df, mode = c("FORWARD", "BACKWARD")) {
#     BestBIC <- CurBIC <- BIC(model)

#     while (TRUE) {
#         CurBIC <- BIC(model)

#         if (mode == "FORWARD") {
#             trialterms <- terms[terms %nin% attr(terms(model), "term.labels")]
#             newformulae <- as.list(
#                 c(
#                     paste(as.character(formula(model)), trialterms, sep = "+"),
#                     paste(as.character(formula(model)), trialterms, sep = "*")
#                 )
#             )
#         }
#         if (mode == "BACKWARD") {
#             Curterms <- attr(terms(model), "term.labels")
#             Curterms <- generate_valid_subsets(Curterms)
#             newformulae <- lapply(Curterms, function(termsiter) {
#                 paste(all.vars(formula(model))[1], paste(termsiter, collapse = "+"), sep = "~")
#             })
#         }

#         newmodels <- lapply(newformulae, FUN = function(formulatry) {
#             updated_model <- update(model, as.formula(formulatry))
#         })

#         NewModPos <- which.min(unlist(lapply(newmodels, BIC)))

#         BestBIC <- BIC(newmodels[[NewModPos]])
#         if (BestBIC > CurBIC) {
#             print("Best model reached")
#             break()
#         }
#         model <- newmodels[[NewModPos]]
#         print(as.character(formula(model)))
#     }
#     return(model)
# }

# ### Resistance ------------
# base_formula <- "pert.name + AC + SL + DI + VA +
#                 pert.name:AC + pert.name:SL + AC:SL + pert.name:DI +
#                 AC:DI + SL:DI + pert.name:VA + AC:VA + SL:VA + DI:VA +
#                 pert.name:AC:SL + pert.name:AC:DI + pert.name:SL:DI +
#                 AC:SL:DI + pert.name:AC:VA + pert.name:SL:VA + AC:SL:VA +
#                 pert.name:DI:VA + AC:DI:VA + SL:DI:VA + pert.name:AC:SL:DI +
#                 pert.name:AC:SL:VA + pert.name:AC:DI:VA + pert.name:SL:DI:VA +
#                 AC:SL:DI:VA + pert.name:AC:SL:DI:VA"


# #### MU0 ###########
# message("Resistance MU 1")
# resist_fwd_mu0 <- do_model_selection(
#     data = Bayes_df_mu0,
#     terms = all.vars(formula(paste0("test ~", base_formula)))[-1],
#     model = betareg(formula(paste("Resistance ~", base_formula)),
#         data = Bayes_df_mu0, na.action = "na.fail"
#     ),
#     mode = "BACKWARD"
# )

# resist_bkw_mu0 <- do_model_selection(
#     data = Bayes_df_mu0,
#     terms = all.vars(formula(paste0("test ~", base_formula)))[-1],
#     model = betareg(formula(paste("Resistance ~", 1)),
#         data = Bayes_df_mu0, na.action = "na.fail"
#     ),
#     mode = "FORWARD"
# )

# resist_mu0 <- list(FWD = resist_fwd_mu0, BKW = resist_bkw_mu0)[[
#     which.min(unlist(lapply(list(FWD = resist_fwd_mu0, BKW = resist_bkw_mu0), BIC)))
# ]]
# print(r2(resist_mu0))

# #### MU1 ###########
# message("Resistance MU 1")
# resist_fwd_mu1 <- do_model_selection(
#     data = Bayes_df_mu1,
#     terms = all.vars(formula(paste0("test ~", base_formula)))[-1],
#     model = betareg(formula(paste("Resistance ~", base_formula)),
#         data = Bayes_df_mu1, na.action = "na.fail"
#     ),
#     mode = "BACKWARD"
# )

# resist_bkw_mu1 <- do_model_selection(
#     data = Bayes_df_mu1,
#     terms = all.vars(formula(paste0("test ~", base_formula)))[-1],
#     model = betareg(formula(paste("Resistance ~", 1)),
#         data = Bayes_df_mu1, na.action = "na.fail"
#     ),
#     mode = "FORWARD"
# )

# resist_mu1 <- list(FWD = resist_fwd_mu1, BKW = resist_bkw_mu1)[[
#     which.min(unlist(lapply(list(FWD = resist_fwd_mu1, BKW = resist_bkw_mu1), BIC)))
# ]]
# print(r2(resist_mu1))

# ### Recovery ------------
# formula_obj <- as.formula(Test ~ pert.name * AC * SL * DI * VA *
#     pre_u_sd * Resistance * c(pre_u_sd - postmin_u_sd))

# # [Environment] Post vs. [Individuals] PostMin_Closest_dist_mean
# # [Environment] Pre vs. [Individuals] Pre_Closest_dist_mean

# terms_obj <- terms(formula_obj)
# # Get the terms as a character vector
# all_terms <- attr(terms_obj, "term.labels")

# base_formula <- paste(all_terms, collapse = "+")

# #### MU0 ###########
# message("Recovery MU 0")
# recov_fwd_mu0 <- do_model_selection(
#     data = Bayes_df_mu0,
#     terms = all.vars(formula(paste0("test ~", base_formula)))[-1],
#     glm(formula(paste("Recovery ~", base_formula)), data = Bayes_df_mu0, family = Gamma(link = "log")),
#     mode = "BACKWARD"
# )

# recov_bkw_mu0 <- do_model_selection(
#     data = Bayes_df_mu0,
#     terms = all.vars(formula(paste0("test ~", base_formula)))[-1],
#     glm(formula(paste("Recovery ~", 1)), data = Bayes_df_mu0, family = Gamma(link = "log")),
#     mode = "FORWARD"
# )

# recov_mu0 <- list(FWD = recov_fwd_mu0, BKW = recov_bkw_mu0)[[
#     which.min(unlist(lapply(list(FWD = recov_fwd_mu0, BKW = recov_bkw_mu0), BIC)))
# ]]
# print(r2(recov_mu0))

# #### MU1 ###########
# message("Recovery MU 1")
# recov_fwd_mu1 <- do_model_selection(
#     data = Bayes_df_mu1,
#     terms = all.vars(formula(paste0("test ~", base_formula)))[-1],
#     glm(formula(paste("Recovery ~", base_formula)), data = Bayes_df_mu1, family = Gamma(link = "log")),
#     mode = "BACKWARD"
# )

# recov_bkw_mu1 <- do_model_selection(
#     data = Bayes_df_mu1,
#     terms = all.vars(formula(paste0("test ~", base_formula)))[-1],
#     glm(formula(paste("Recovery ~", 1)), data = Bayes_df_mu1, family = Gamma(link = "log")),
#     mode = "FORWARD"
# )

# recov_mu1 <- list(FWD = recov_fwd_mu1, BKW = recov_bkw_mu1)[[
#     which.min(unlist(lapply(list(FWD = recov_fwd_mu1, BKW = recov_bkw_mu1), BIC)))
# ]]
# print(r2(recov_mu1))


# # VISUALISATION =============================
# stop("INSEPCT")

# ## plots of coefficients
# library(tidyr)

# ModelCoeff_ls <- list(
#     Mu0 = summary(recov_mu0[[1]])$coefficients,
#     Mu1 = summary(recov_mu1[[1]])$coefficients
# )

# PlotCoeffs_ls <- lapply(1:length(ModelCoeff_ls), FUN = function(mat) {
#     Name <- names(ModelCoeff_ls)[mat]
#     mat <- ModelCoeff_ls[[mat]]
#     ret <- as.data.frame(mat) %>%
#         tibble::rownames_to_column("row_name") %>% # Add row names as a column
#         pivot_longer(
#             cols = -row_name, # All columns except the row_name column
#             names_to = "column_name", # Name for new column storing column names
#             values_to = "value" # Name for new column storing values
#         )
#     ret$model <- Name
#     ret
# })
# PlotCoeffs_df <- do.call(rbind, PlotCoeffs_ls)

# # Filter rows for "Estimate" and "Std. Error", then reshape the data
# library(dplyr)
# library(tidyr)
# plot_data <- PlotCoeffs_df %>%
#     filter(column_name %in% c("Estimate", "Std. Error")) %>%
#     pivot_wider(names_from = column_name, values_from = value) %>%
#     mutate(
#         lower_bound = Estimate - 0.5 * `Std. Error`, # Lower range
#         upper_bound = Estimate + 0.5 * `Std. Error` # Upper range
#     )

# # Plot the estimates with error bars
# ggplot(plot_data, aes(x = Estimate, y = model)) +
#     geom_point(size = 5) + # Point for the estimate
#     geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound), width = 0.2) + # Error bars
#     # labs(
#     #   title = "Estimate with Range (± Half Std. Error)",
#     #   y = "Parameter",
#     #   x = "Estimate"
#     # ) +
#     geom_vline(xintercept = 0) +
#     facet_wrap(~row_name, ncol = 1, scales = "free_x") +
#     theme_minimal()
# ggsave(file.path(getwd(), "Exports/GLMCoeffs_Recovery.jpg"), width = 9, height = 90, units = "cm")
