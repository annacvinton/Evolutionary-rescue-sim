#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes]
#' CONTENTS:
#'  - Bayesian Estimation of Changepoints in abundance Time Series
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
Dir.Beast <- file.path(Dir.Exports, "BEAST")
### Create directories which aren't present yet
Dirs <- c(Dir.Data, Dir.Exports, Dir.Beast)
CreateDir <- sapply(Dirs, function(x) if (!dir.exists(x)) dir.create(x))

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, repos = "http://cran.us.r-project.org")
    }
    require(x, character.only = TRUE)
}
package_vec <- c(
    "Rbeast",
    "ggplot2",
    "pbapply",
    "parallel"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
load(file.path(Dir.Exports, "POPULATION_TimeStep.RData")) # loads Data_df
ID_vec <- unique(Data_df$ID)

## Conditions, time and survival of runs -----------------------------------
conditions <- Data_df[, c("pert.name", "rep", "AC", "DI", "SL", "VA", "MU")]
conditions <- conditions[which(!duplicated(conditions)), ]
# runtimes <- aggregate(Data_df, t ~ pert.name + rep + AC + DI + MU + SL + VA, FUN = max)
runtimes <- aggregate(Data_df, t ~ ID, FUN = max)
runtimes <- cbind(runtimes, do.call(
    rbind,
    pblapply(runtimes$ID, function(id) {
        # id <- runtimes$ID[1]
        x <- Data_df[which(Data_df$ID == id), c("pert.name", "rep", "AC", "DI", "SL", "VA", "MU", "n", "t", "u_sd")]
        if (460 %in% x$t) {
            x$pre_n <- x$n[which(x$t == 460)]
            x$pre_u_sd <- x$u_sd[which(x$t == 460)]
        } else {
            x$pre_n <- NA
            x$pre_u_sd <- NA
        }
        x <- x[1, c("pert.name", "rep", "AC", "DI", "SL", "VA", "MU", "pre_n", "pre_u_sd")]
        x$survival <- ifelse(runtimes$t[runtimes$ID == id] < 1110, FALSE, TRUE)
        return(x)
    })
))

# MODELS ==================================================================
cl <- parallel::makeCluster(ifelse(parallel::detectCores() > 30, 30, parallel::detectCores()))
clusterExport(cl = cl, varlist = c("Dir.Data", "Dir.Exports", "Dir.Beast", "install.load.package", "package_vec", "Data_df", "ID_vec"))
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

Model_ls <- pblapply(
    ID_vec, # [-1:-47832],
    cl = cl,
    function(id) {
        # id <- ID_vec[47832]
        print(paste(which(ID_vec == id), "/", length(ID_vec)))
        # print(paste("Fitting BEAST model for ID", id))
        x <- Data_df[Data_df$ID == id, ]
        x_ts <- x[order(x$t), ]
        if (nrow(x_ts) == 1) {
            return(NULL)
        }
        n_ts <- ts(
            data = as.numeric(x_ts$n),
            start = min(x_ts$t),
            frequency = 1 / median(diff(x_ts$t))
        )

        if (max(time(n_ts)) < 460 | min(time(n_ts)) > 460 | length(n_ts) == 2) {
            return(NULL)
        }

        FNAME <- file.path(Dir.Beast, paste0("BEAST_", id, ".RData"))
        if (file.exists(FNAME)) {
            load(FNAME)
        } else {
            beast_mod <- Rbeast::beast(
                n_ts,
                tcp.minmax = c(0, 3),
                period = "none",
                quiet = TRUE,
                print.progress = FALSE
            )
            save(beast_mod, file = FNAME)
            Sys.sleep(1)
            # plot(beast_mod)
        }

        ## extract relevant information about change-points
        cp <- beast_mod$trend$cp # extract change-points
        cpPr <- beast_mod$trend$cpPr # extract change-point probabilities
        cpCI <- beast_mod$trend$cpCI # extract change-point credible intervals
        p <- beast_mod$trend$cpOccPr # extract posterior distribution of change-point occurrence
        time <- as.numeric(time(n_ts))

        ## remove NA values
        valid <- !is.na(cp)
        cp <- cp[valid]
        cpPr <- cpPr[valid]
        if (is.null(nrow(cpCI))) {
            cpCI <- data.frame(cpCI[1], cpCI[2])
        } else {
            cpCI <- cpCI[valid, ]
        }


        ## subset to change-points only starting when the perturbation is introduced at t = 460
        after460 <- cp > 460 # a changepoint at 460 denotes basically immediate extinction
        if (sum(after460) == 0) {
            return(NULL)
        }
        cp <- cp[after460]
        cpPr <- cpPr[after460]
        if (is.null(nrow(cpCI))) {
            cpCI <- data.frame(cpCI[1], cpCI[2])
        } else {
            cpCI <- cpCI[after460, ]
        }

        ## order change-points by their occurrence in time
        ord <- order(cp)
        cp <- cp[ord]
        cpPr <- cpPr[ord]
        if (is.null(nrow(cpCI))) {
            cpCI <- data.frame(cpCI[1], cpCI[2])
        } else {
            cpCI <- cpCI[ord, ]
        }

        ## extract the intervals of the two earliest change-points after 460 - these should be the resistance and recovery points and build their posterior distributions
        if (nrow(cpCI) >= 2) {
            numCI <- 1:2
        } else {
            numCI <- 1
        }
        changePoints_ls <- lapply(numCI, function(i) {
            # i = 1
            idx <- which(time >= cpCI[i, 1] & time <= cpCI[i, 2])
            post_cp <- p[idx]
            post_cp <- post_cp / sum(post_cp)
            time_cp <- time[idx]
            if (length(time_cp) == 1) {
                samples_cp <- rep(time_cp, 1000)
            } else {
                samples_cp <- sample(time_cp, 1000, replace = TRUE, prob = post_cp)
            }
            data.frame(changePoint = i, posteriorSample = samples_cp)
        })
        changePoints_df <- do.call(rbind, changePoints_ls)
        changePoints_df$n <- x_ts$n[match(changePoints_df$posteriorSample, x_ts$t)]
        changePoints_df$AC <- strsplit(id, split = "-")[[1]][1]
        changePoints_df$DI <- strsplit(id, split = "-")[[1]][2]
        changePoints_df$MU <- strsplit(id, split = "-")[[1]][3]
        changePoints_df$SL <- strsplit(id, split = "-")[[1]][4]
        changePoints_df$VA <- strsplit(id, split = "-")[[1]][5]
        changePoints_df$pert.name <- strsplit(id, split = "-")[[1]][6]
        changePoints_df$rep <- strsplit(id, split = "-")[[1]][7]

        return(changePoints_df)
    }
)
BEAST_ChangePoints_df <- do.call(rbind, Model_ls[sapply(Model_ls, function(x) if (!is.null(x)) TRUE else FALSE)])
save(BEAST_ChangePoints_df, file = file.path(Dir.Exports, "BEAST_ChangePoints.RData"))
# load(file.path(Dir.Exports, "BEAST_ChangePoints.RData"))

BEAST_ChangePoints_df$ID <- with(BEAST_ChangePoints_df, paste(AC, DI, MU, SL, VA, pert.name, rep, sep = "-"))

med_n_chps <- aggregate(data = BEAST_ChangePoints_df, n ~ changePoint + ID, FUN = median)
med_t_chps <- aggregate(data = BEAST_ChangePoints_df, posteriorSample ~ changePoint + ID, FUN = median)
sd_n_chps <- aggregate(data = BEAST_ChangePoints_df, n ~ changePoint + ID, FUN = sd)
sd_t_chps <- aggregate(data = BEAST_ChangePoints_df, posteriorSample ~ changePoint + ID, FUN = sd)

Chps_summary_df <- pblapply(runtimes$ID, function(id) {
    # id <- runtimes$ID[2]
    if (!id %in% med_n_chps$ID) {
        return(
            data.frame(
                CP1_n = NA,
                CP1_n_sd = NA,
                CP1_t = NA,
                CP1_t_sd = NA,
                CP2_n = NA,
                CP2_n_sd = NA,
                CP2_t = NA,
                CP2_t_sd = NA,
                ID = id
            )
        )
    }
    CP1_n <- med_n_chps[which(med_n_chps$ID == id & med_n_chps$changePoint == 1), "n"]
    CP2_n <- med_n_chps[which(med_n_chps$ID == id & med_n_chps$changePoint == 2), "n"]
    CP1_n_sd <- sd_n_chps[which(sd_n_chps$ID == id & sd_n_chps$changePoint == 1), "n"]
    CP2_n_sd <- sd_n_chps[which(sd_n_chps$ID == id & sd_n_chps$changePoint == 2), "n"]
    CP1_t <- med_t_chps[which(med_t_chps$ID == id & med_t_chps$changePoint == 1), "posteriorSample"]
    CP2_t <- med_t_chps[which(med_t_chps$ID == id & med_t_chps$changePoint == 2), "posteriorSample"]
    CP1_t_sd <- sd_t_chps[which(sd_t_chps$ID == id & sd_t_chps$changePoint == 1), "posteriorSample"]
    CP2_t_sd <- sd_t_chps[which(sd_t_chps$ID == id & sd_t_chps$changePoint == 2), "posteriorSample"]
    data.frame(
        CP1_n = ifelse(length(CP1_n) == 0, NA, CP1_n),
        CP1_n_sd = ifelse(length(CP1_n_sd) == 0, NA, CP1_n_sd),
        CP1_t = ifelse(length(CP1_t) == 0, NA, CP1_t),
        CP1_t_sd = ifelse(length(CP1_t_sd) == 0, NA, CP1_t_sd),
        CP2_n = ifelse(length(CP2_n) == 0, NA, CP2_n),
        CP2_n_sd = ifelse(length(CP2_n_sd) == 0, NA, CP2_n_sd),
        CP2_t = ifelse(length(CP2_t) == 0, NA, CP2_t),
        CP2_t_sd = ifelse(length(CP2_t_sd) == 0, NA, CP2_t_sd),
        ID = id
    )
})
Chps_summary_df <- do.call(rbind, Chps_summary_df)

MODEL_Metrics <- merge(
    x = runtimes,
    y = Chps_summary_df,
    by = "ID",
    all.x = TRUE
)
colnames(MODEL_Metrics)[which(colnames(MODEL_Metrics) == "t")] <- "t_max"

MODEL_Metrics$Res_Magnitude <- (MODEL_Metrics$CP1_n - MODEL_Metrics$pre_n) / MODEL_Metrics$pre_n
MODEL_Metrics$Res_Speed <- MODEL_Metrics$CP1_t - 460
MODEL_Metrics$Rec_Magnitude <- (MODEL_Metrics$CP2_n - MODEL_Metrics$CP1_n) / MODEL_Metrics$CP1_n
MODEL_Metrics$Rec_Speed <- MODEL_Metrics$CP2_t - MODEL_Metrics$CP1_t
MODEL_Metrics$Rec_MagnitudePre <- (MODEL_Metrics$CP2_n - MODEL_Metrics$pre_n) / MODEL_Metrics$pre_n
MODEL_Metrics$Rec_SpeedPre <- MODEL_Metrics$CP2_t - 460

save(MODEL_Metrics, file = file.path(Dir.Exports, "ModelMetrics.RData"))
