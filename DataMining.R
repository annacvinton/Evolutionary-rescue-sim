#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Data Mining
#'  DEPENDENCIES:
#'  - DataCleaning.R must have been run
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
  "ggplot2", # for plotting
  "viridis", # colour palettes
  "pbapply", # for apply functionality with progress bars
  "data.table", # for rbindlist
  "tidybayes" # for some nice distribution visualisations with the ggplot engine
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
cl <- parallel::makeCluster(parallel::detectCores())
parallel::clusterExport(cl,
                        varlist = c("Dir.Data", "install.load.package", "package_vec"),
                        envir = environment()
)
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

## Loading & Extraction ---------------------------------------------------
Sims_fs <- list.files(Dir.Data, ".rds")

DataMine_ls <- pblapply(Sims_fs, cl = cl, function(x){
  # message(x)
  data_df <- readRDS(file.path(Dir.Data, x))
  Ident_vec <- unlist(strsplit(tools::file_path_sans_ext(x), "_")) # vector of identifiers
  Ident_vec <- as.numeric(gsub("[^0-9.-]", "", Ident_vec)) # retain only numbers
  
  ## total number of samples
  samples <- data.frame(samples = nrow(data_df),
                        Autocorrelation = Ident_vec[1],
                        Dispersal = Ident_vec[2],
                        Mutation = Ident_vec[3],
                        Slope = Ident_vec[4],
                        Variance = Ident_vec[5])
  
  
  ## how many samples per treatment
  SamplesPerTreatment <- as.data.frame(table(data_df$pert.name))
  colnames(SamplesPerTreatment) <- c("pert.name", "samples")
  SamplesPerTreatment <- cbind(SamplesPerTreatment, 
                               data.frame(Autocorrelation = rep(Ident_vec[1], nrow(SamplesPerTreatment)),
                                          Dispersal = rep(Ident_vec[2], nrow(SamplesPerTreatment)),
                                          Mutation = rep(Ident_vec[3], nrow(SamplesPerTreatment)),
                                          Slope = rep(Ident_vec[4], nrow(SamplesPerTreatment)),
                                          Variance = rep(Ident_vec[5], nrow(SamplesPerTreatment))))
  
  ## timesteps per treatment (replicates count separately)
  timepertreatandrep <- data.frame(pert = data_df$pert.name, 
                                   timerep = paste(data_df$rep, data_df$t, sep = "_"))
  TimestepsPerPert <- aggregate(timerep ~ pert, data = timepertreatandrep, FUN = function(x) length(unique(x)))
  colnames(TimestepsPerPert) <- c("pert.name", "timesteps")
  TimestepsPerPert <- cbind(TimestepsPerPert, 
                            data.frame(Autocorrelation = rep(Ident_vec[1], nrow(TimestepsPerPert)),
                                       Dispersal = rep(Ident_vec[2], nrow(TimestepsPerPert)),
                                       Mutation = rep(Ident_vec[3], nrow(TimestepsPerPert)),
                                       Slope = rep(Ident_vec[4], nrow(TimestepsPerPert)),
                                       Variance = rep(Ident_vec[5], nrow(TimestepsPerPert))))
  
  ## timesteps following perturbation
  timeafterpert <- data_df[c(
    which(!duplicated(paste(data_df$pert.name, data_df$rep, data_df$pert.value, sep = "_")))[-1],
    nrow(data_df)
  ), ]
  NonPertRuns <- which(c(head(timeafterpert$pert.value,-1) == tail(timeafterpert$pert.value,-1)))
  timeafterpert <- timeafterpert[-NonPertRuns, ]
  if(nrow(timeafterpert) %% 2 != 0){
    timeafterpert <- timeafterpert[-nrow(timeafterpert),]
  }
  timeafterpert <- lapply(seq(from = 1, to = nrow(timeafterpert), by = 2), function(counter){
    data.frame(pert = timeafterpert[counter,"pert.value"],
               tsteps = length(unique(data_df[rownames(timeafterpert)[counter]:rownames(timeafterpert)[counter+1],"t"]))
    )
  })
  TimestepsAfterPert <- do.call(rbind, timeafterpert)
  colnames(TimestepsAfterPert) <- c("pert.value", "timesteps")
  TimestepsAfterPert <- cbind(TimestepsAfterPert, 
                              data.frame(Autocorrelation = rep(Ident_vec[1], nrow(TimestepsAfterPert)),
                                         Dispersal = rep(Ident_vec[2], nrow(TimestepsAfterPert)),
                                         Mutation = rep(Ident_vec[3], nrow(TimestepsAfterPert)),
                                         Slope = rep(Ident_vec[4], nrow(TimestepsAfterPert)),
                                         Variance = rep(Ident_vec[5], nrow(TimestepsAfterPert))))
  ## export mined data
  list(samples = samples,
       SamplesPerTreatment = SamplesPerTreatment,
       TimestepsPerPert = TimestepsPerPert,
       TimestepsAfterPert = TimestepsAfterPert)  
})
save(DataMine_ls, file = file.path(Dir.Exports, "DataMine_ls.RData"))

parallel::stopCluster(cl)
closeAllConnections()

## data extraction and reformatting
samples_df <- do.call(rbind, lapply(DataMine_ls, "[[", "samples"))
samples_df <- aggregate(samples ~., data = samples_df, FUN = sum)
write.csv(samples_df, file = file.path(Dir.Exports, "Samples.csv"))

SamplesPerTreatment_df <- do.call(rbind, lapply(DataMine_ls, "[[", "SamplesPerTreatment"))




ggplot(SamplesPerTreatment_df, 
       aes(x = samples, 
           y = factor(pert.name, levels = c(0, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)), 
           col = factor(Mutation)
           )
       ) + 
  geom_boxplot(position = "dodge") + 
  facet_wrap(~Slope+Variance, scales = "free") + 
  theme_bw()


stop("Autocorrelation and Dispersal")

SamplesPerTreatment_df <- aggregate(samples ~., data = SamplesPerTreatment_df, FUN = sum)

TimestepsPerPert_df <- do.call(rbind, lapply(DataMine_ls, "[[", "TimestepsPerPert"))
TimestepsPerPert_df <- aggregate(timesteps ~., data = TimestepsPerPert_df, FUN = sum)

TimestepsAfterPert_df <- do.call(rbind, lapply(DataMine_ls, "[[", 4))
TimestepsAfterPert_df <- aggregate(timesteps ~., data = TimestepsAfterPert_df, FUN = sum)














































