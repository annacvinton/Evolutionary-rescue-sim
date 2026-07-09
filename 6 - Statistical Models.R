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
MODEL_Metrics$pert <- as.numeric(MODEL_Metrics$pert)
MODEL_Metrics <- MODEL_Metrics[MODEL_Metrics$Rec_Magnitude >= 0, ] # limit to only those runs where second change-point indicators abundance recovery

## Data Limiting ----------------------------------------------------------
MODEL_Metrics <- MODEL_Metrics[MODEL_Metrics$pert >= 9, ]
MODEL_Metrics <- MODEL_Metrics[MODEL_Metrics$MU == 0, ]

# MODEL_Metrics <- MODEL_Metrics[MODEL_Metrics$AC == 1, ]
# MODEL_Metrics <- MODEL_Metrics[MODEL_Metrics$DI == 2, ]

# MODELS ==================================================================
Outcomes <- c("Res_Magnitude", "Res_Speed", "Rec_Magnitude", "Rec_Speed")
Parameters <- c("AC", "DI", "SL", "VA")
Combinations <- c(
    "DI = 1.5, SL = 1.2, VA = 50",
    "DI = 1.5, SL = 1.2, VA = 100",
    "DI = 2, SL = 1, VA = 100",
    "DI = 2, AC = 1|2, VA = 25|50|100",
    "AC = 0, DI = 2, SL = 1"
)

subset_by_combination <- function(df, combination) {
    parts <- strsplit(combination, split = ",\\s*")[[1]]
    keep <- rep(TRUE, nrow(df))

    for (part in parts) {
        kv <- strsplit(part, split = "\\s*=\\s*")[[1]]
        parameter <- trimws(kv[1])
        values <- trimws(strsplit(kv[2], split = "\\|")[[1]])
        values <- suppressWarnings(as.numeric(values))
        values <- values[!is.na(values)]

        keep <- keep & df[[parameter]] %in% values
    }

    df[keep, , drop = FALSE]
}

## Iterate over combinations and identify outcome as the parameter not listed
IterData_ls <- lapply(Combinations, FUN = function(Combination) {
    # Combination <- Combinations[4]
    print(Combination)

    ## identify driver parameter for model
    DriverParam <- Parameters[!(sapply(Parameters, FUN = function(Parameter) {
        # Parameter <- Parameters[1]
        grepl(Parameter, Combination)
    }))]

    IterData <- subset_by_combination(MODEL_Metrics, Combination)

    ## Iterate over outcome variables and produce plots facetted by perturbation level
    OutcomePlots_ls <- pblapply(Outcomes, FUN = function(Outcome) {
        # Outcome <- "Rec_Magnitude"
        # print(Outcome)
        Coeffs_df <- lapply(sort(unique(MODEL_Metrics$pert)), FUN = function(pert) {
            # pert <- 15
            # print(pert)
            IterData <- IterData[IterData$pert == pert, ]
            n_runs <- length(unique(IterData$rep))
            LinMod <- tryCatch(
                {
                    LinMod <- lm(
                        as.formula(paste0(Outcome, " ~ 0 + ", DriverParam)),
                        data = IterData
                    )
                },
                error = function(e) {
                    # message("Error in lm for pert: ", pert, " - ", e$message)
                    return(NULL)
                }
            )

            if (is.null(LinMod)) {
                return(NULL)
            }

            # extract coefficients and p-values
            LinMod_summary <- summary(LinMod)
            CoeffNames <- rownames(LinMod_summary$coefficients) # Coefficient names
            estimates <- LinMod_summary$coefficients[, 1] # Coefficients
            StdError <- LinMod_summary$coefficients[, 2] # Standard errors
            pvals <- LinMod_summary$coefficients[, 4] # p-values

            plot_df <- data.frame(
                Coefficients = CoeffNames, Estimates = estimates, P_values = pvals, StError = StdError,
                pert = paste("Perturbation: ", str_pad(pert, 2, "left", "0"), "; # of Simulations: ", n_runs, sep = "")
            )

            plot_df
        })
        Coeffs_df <- do.call(rbind, Coeffs_df[!unlist(lapply(Coeffs_df, is.null))])

        ggplot(Coeffs_df, aes(x = Coefficients, y = Estimates)) +
            geom_bar(stat = "identity", aes(fill = P_values < 0.05)) +
            labs(x = "Driver", y = "Estimate", title = paste("Outcome:", Outcome, "; Combination:", Combination)) +
            geom_errorbar(aes(ymin = Estimates - StError, ymax = Estimates + StError), width = 0.5) +
            scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "darkred"), name = "Significant") +
            # lims(y = c(0, NA)) +
            facet_wrap(~pert, ncol = 1) +
            theme_minimal() +
            theme(legend.position = "bottom", text = element_text(size = 20))
    })

    ## save plots into one PDF file
    ggsave(
        filename = file.path(Dir.Exports, paste0("OutcomePlots_", gsub("[^A-Za-z0-9]", "_", Combination), ".pdf")),
        plot = marrangeGrob(grobs = OutcomePlots_ls, nrow = 2, ncol = 2),
        width = 42, height = 42
    )
})
