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
    "ggplot2"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
load(file.path(Dir.Exports, "POPULATION_TimeStep.RData")) # loads Data_df
ID_vec <- unique(Data_df$ID)

# MODELS ==================================================================
Model_ls <- lapply(ID_vec, function(id) {
    # id <- ID_vec[20]
    print(paste(which(ID_vec == id), "/", length(ID_vec)))
    # print(paste("Fitting BEAST model for ID", id))
    x <- Data_df[Data_df$ID == id, ]
    x_ts <- x[order(x$t), ]
    n_ts <- ts(
        data = as.numeric(x_ts$n),
        start = min(x_ts$t),
        frequency = 1 / median(diff(x_ts$t))
    )

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
    after450 <- cp > 450 # a changepoint at 460 denotes basically immediate extinction
    if (sum(after450) == 0) {
        return(NULL)
    }
    cp <- cp[after450]
    cpPr <- cpPr[after450]
    if (is.null(nrow(cpCI))) {
        cpCI <- data.frame(cpCI[1], cpCI[2])
    } else {
        cpCI <- cpCI[after450, ]
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
    Sys.sleep(60)
})
