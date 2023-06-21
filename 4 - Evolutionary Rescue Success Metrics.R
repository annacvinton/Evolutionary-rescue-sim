#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Evolutionary Success Reading
#'  DEPENDENCIES:
#'  - 2 - SummaryMetricsExtract.R must have been run
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
  "parallel"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
load(file.path(Dir.Exports, "NN_Metrics.RData"))
conditions <- Data_df[,c("pert.name", "rep", "AC", "DI", "SL", "VA", "MU")]
conditions <- conditions[which(!duplicated(conditions)), ]
conditions

cl <- makeCluster(detectCores())
clusterExport(cl, c("conditions", "Data_df"))

EvoResSuc_ls <- pblapply(1:nrow(conditions), 
                       cl = cl,
                       FUN = function(x){

  iter_df <- Data_df[
    which(conditions[x,]$pert.name == Data_df$pert.name & 
            conditions[x,]$rep == Data_df$rep &
            conditions[x,]$AC == Data_df$AC &
            conditions[x,]$DI == Data_df$DI &
            conditions[x,]$SL == Data_df$SL &
            conditions[x,]$VA == Data_df$VA &
            conditions[x,]$MU == Data_df$MU)
    ,
  ]
  
  if(max(iter_df$t) < 470 | min(iter_df$t)>450){
    NULL
  }else{
    data.frame(Min = min(iter_df$n), # minimum population size
               t_Min = iter_df$t[which.min(iter_df$n)], # time of minimum population size
               n_pre = iter_df$n[iter_df$t == max(iter_df$t[iter_df$t<470])], # last time step before perturbation
               n_470 = ifelse(length(iter_df$n[iter_df$t == 470]) == 0, NA, iter_df$n[iter_df$t == 470]), # directly after perturbation
               n_600 = ifelse(length(iter_df$n[iter_df$t == 600]) == 0, NA, iter_df$n[iter_df$t == 600]), # quite some time after perturbation
               max_post = max(iter_df$n[iter_df$t > 460]), # maximum pop size under perturbation
               t_max_post = iter_df$t[iter_df$t > 460][which.max(iter_df$n[iter_df$t > 460])], # time of maximum pop size under perturbation
               max_t = max(iter_df$t) # time to which run was recorded
    )
  }
  
})
NonSkip <- which(!unlist(lapply(EvoResSuc_ls, is.null)))
EvoResSuc_df <- cbind(conditions[NonSkip,], do.call(rbind, EvoResSuc_ls[NonSkip]))
write.table(EvoResSuc_df, file = file.path(Dir.Exports, "EvoResSuccessMetrics.txt"))
