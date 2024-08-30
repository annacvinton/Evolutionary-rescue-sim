#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Calculation of time-series moments (highest/lowest population size, time to lowest population size, etc.)
#'  - Classification of simulation runs as either (1) extinct, (2) not rescued, (3) evolutionary rescued
#'  DEPENDENCIES:
#'  - 2 - SummaryMetricsExtract.R must have been run and produced POPULATION_TimeStep.RData
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list=ls())

## Directories ------------------------------------------------------------
### Define dicrectories in relation to project directory
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "Data")
Dir.Exports <- file.path(Dir.Base, "Exports")
### Create directories which aren't present yet
Dirs <- c(Dir.Data, Dir.Exports)
CreateDir <- sapply(Dirs, function(x) if(!dir.exists(x)) dir.create(x))

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}
package_vec <- c(
  "pbapply",
  "parallel",
  "ggplot2",
  "cowplot",
  "dplyr"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
## Loading Data -----------------------------------------------------------
load(file.path(Dir.Exports, "POPULATION_TimeStep.RData"))

## Conditions for runs and time of runs -----------------------------------
conditions <- Data_df[,c("pert.name", "rep", "AC", "DI", "SL", "VA", "MU")]
conditions <- conditions[which(!duplicated(conditions)), ]
runtimes <- aggregate(Data_df, t ~ pert.name + rep + AC + DI + MU + SL + VA, FUN = max)

## Survival Classification ------------------------------------------------
runtimes$survival <- TRUE
runtimes$survival[runtimes$t < 1110] <- FALSE
runtimes$ID <- 1:nrow(runtimes)

## Non-Survival Classification --------------------------------------------
noext <- runtimes[runtimes$survival, ]

cl <- makeCluster(detectCores())
clusterExport(cl, c("conditions", "Data_df", "%nin%", "noext"))

evoressuc_ls <- pblapply(1:nrow(noext), 
                         cl = cl,
                         FUN = function(l){
                           ### Data subsetting ----
                           # message(l)
                           iter2_df <- Data_df[
                             Data_df$pert.name == noext[l, "pert.name"] &
                               Data_df$rep == noext[l, "rep"] &
                               Data_df$AC == noext[l, "AC"] &
                               Data_df$DI == noext[l, "DI"] &
                               Data_df$MU == noext[l, "MU"] &
                               Data_df$SL == noext[l, "SL"] &
                               Data_df$VA == noext[l, "VA"]
                             ,   
                           ]
                           ### not enough data or misaligned time-step of perturbation switch ----
                           iter_df <- iter2_df[iter2_df$t>460, ]
                           if(nrow(iter_df) < 65 | 460 %nin% iter2_df$t){
                             bind_df <- data.frame(
                               ID = noext$ID[l], 
                               n_pre = NA,
                               n_minpost = NA,
                               t_minpost = NA,
                               n_maxpost = NA,
                               t_maxpost = NA,
                               perc_minpost = NA,
                               perc_maxpostmin = NA,
                               perc_maxpostpre = NA,
                               survival = NA,
                               SuffDip = NA,
                               SuffReb = NA,
                               EvoRes = NA
                             )
                           }else{
                             ### Success Metrics ----
                             pre_n <- Data_df[ ## pre-perturbation
                               Data_df$pert.name == noext[l, "pert.name"] &
                                 Data_df$rep == noext[l, "rep"] &
                                 Data_df$AC == noext[l, "AC"] &
                                 Data_df$DI == noext[l, "DI"] &
                                 Data_df$MU == noext[l, "MU"] &
                                 Data_df$SL == noext[l, "SL"] &
                                 Data_df$VA == noext[l, "VA"] &
                                 Data_df$t == 460
                               , 
                               "n"
                             ]
                             ## minimum population size post perturbation
                             minpost <- which.min(iter_df$n)
                             minpost_df <- iter_df[minpost, c("n", "t")]
                             minpost_df$ChangePre <- minpost_df$n/pre_n
                             
                             ## if minimum is not final step, calculate maximum post perturbation & post-minimum
                             if(minpost != nrow(iter_df)){
                               maxpost <- which.max(iter_df$n[(minpost+1):nrow(iter_df)])+minpost
                               maxpost_df <- iter_df[maxpost, c("n", "t")]
                               maxpost_df$ChangePre <- maxpost_df$n/pre_n
                               maxpost_df$ChangePost <- maxpost_df$n/minpost_df$n
                             }else{
                               maxpost_df <- data.frame(
                                 n = NA,
                                 t = NA,
                                 ChangePre = NA,
                                 ChangePost = NA
                               )
                             }

                             ## combine into data frame
                             bind_df <- data.frame(
                               n_pre = pre_n,
                               ID = noext$ID[l], 
                               n_minpost = minpost_df$n,
                               t_minpost = minpost_df$t,
                               n_maxpost = maxpost_df$n,
                               t_maxpost = maxpost_df$t,
                               perc_minpost = minpost_df$ChangePre,
                               perc_maxpostmin = maxpost_df$ChangePost,
                               perc_maxpostpre = maxpost_df$ChangePre,
                               survival = TRUE,
                               SuffDip = FALSE,
                               SuffReb = FALSE,
                               EvoRes = FALSE
                             )
                             
                             ## classification
                             DipCut <- 0.1 # go down to below 10% of pre-perturbation abundance
                             RebCut <- 0.5 # bounce back to at least 50% of pre-perturbation abundance
                             bind_df$SuffDip <- bind_df$perc_minpost < DipCut
                             bind_df$SuffReb <- bind_df$perc_maxpostpre >= RebCut
                             bind_df$EvoRes <- (bind_df$SuffDip+bind_df$SuffReb == 2)
                           }
                           bind_df
                         })
EVORES_Metrics <- do.call(rbind, evoressuc_ls)
EVORES_Metrics <- na.omit(EVORES_Metrics)
## merging survival runs
save_df <- base::merge(x = runtimes[], 
                       y = EVORES_Metrics[, colnames(EVORES_Metrics) != "survival"], 
                       by = "ID")
### bring back in non-survival runs
EVORES_Metrics <- bind_rows(save_df, runtimes[!runtimes$survival, ])
### saving data
save(EVORES_Metrics, file = file.path(Dir.Exports, "EVORES_Metrics.RData"))
