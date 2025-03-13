#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Extraction of population-level metrics at each time-step for each simulation run
#'  DEPENDENCIES:
#'  - 1 - DataCleaning.R must have been run
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
  "parallel",
  "e1071", # for skewness and kurtosis calculation
  "sp",
  # "rgeos",
  "nabor" # for knn nearest neighbours
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
## Parallelisation --------------------------------------------------------
cl <- parallel::makeCluster(parallel::detectCores())
parallel::clusterExport(cl,
                        varlist = c("Dir.Data", "Dir.Exports",
                                    "install.load.package", "package_vec"),
                        envir = environment()
)
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

## Loading & Extraction ---------------------------------------------------
Sims_fs <- list.files(Dir.Data, ".rds") ## data files
Sims_fs <- Sims_fs[order(file.size(file.path(Dir.Data, Sims_fs)), decreasing = FALSE)]
## loop over data files
Data_ls <- pblapply(Sims_fs, 
                    cl = cl,
                    function(x){
                      ### Load data ----
                      # x <- Sims_fs[1]
                      message(x)
                      
                      if(file.exists(file.path(Dir.Exports, paste0("TEMP_TStepPop_", x, ".RData")))){
                        load(file.path(Dir.Exports, paste0("TEMP_TStepPop_", x, ".RData")))
                      }else{
                        data_df <- readRDS(file.path(Dir.Data, x))
                        keeprows <- (data_df$t <= 1110)
                        data_df <- data_df[keeprows, ]
                        duplicates_check <- with(data_df, paste(
                          pert.name, rep, t, sep = "_"
                        ))
                        
                        ## Simulation Settings ----
                        Ident_vec <- unlist(strsplit(tools::file_path_sans_ext(x), "_")) # vector of identifiers
                        Ident_vec <- as.numeric(gsub("[^0-9.-]", "", Ident_vec)) # retain only numbers
                        AC <- Ident_vec[1]
                        DI <- Ident_vec[2]
                        MU <- Ident_vec[3]
                        SL <- Ident_vec[4]
                        VA <- Ident_vec[5]

                        ## Population size ----
                        N <- aggregate(n ~ t+pert.name+rep, data = data_df, FUN = mean)
                        N$AC <- AC
                        N$DI <- DI
                        N$MU <- MU
                        N$SL <- SL
                        N$VA <- VA

                        ## Trait measures ----
                        u_mean <- aggregate(u ~ t+pert.name+rep, data = data_df,
                                            FUN = mean)$u
                        u_sd <- aggregate(u ~ t+pert.name+rep, data = data_df,
                                          FUN = sd)$u
                        u_skewness <- aggregate(u ~ t+pert.name+rep, data = data_df,
                                                FUN = e1071::skewness)$u
                        u_kurtosis <- aggregate(u ~ t+pert.name+rep, data = data_df,
                                                FUN = e1071::kurtosis)$u
                        ## Adaptedness measures ----
                        AdaptFunc <- function(SL, x ,u, patch){
                          abs(SL*x+patch - u)
                        }
                        data_df$MalAdaptedness <- do.call(function(SL, x ,u, patch) AdaptFunc(SL, x, u, patch),
                                                          args = list(data_df$x,
                                                                      data_df$u,
                                                                      data_df$patch,
                                                                      SL)
                        )
                        MalAdap_mean <- aggregate(MalAdaptedness ~ t+pert.name+rep, data = data_df,
                                                  FUN = mean)$MalAdaptedness
                        MalAdap_sd <- aggregate(MalAdaptedness ~ t+pert.name+rep, data = data_df,
                                                FUN = sd)$MalAdaptedness
                        MalAdap_skewness <- aggregate(MalAdaptedness ~ t+pert.name+rep, data = data_df,
                                                      FUN = e1071::skewness)$MalAdaptedness
                        MalAdap_kurtosis <- aggregate(MalAdaptedness ~ t+pert.name+rep, data = data_df,
                                                      FUN = e1071::kurtosis)$MalAdaptedness

                        ## Spatial measures ----
                        min.d <- c()
                        SimSteps <- unique(duplicates_check) # loop over all individual combinations of timesteps respective to simulation runs
                        pb <- txtProgressBar(max = length(SimSteps), style = 3)
                        k <- 1
                        for(i in SimSteps){
                          sp.mydata <- data_df[duplicates_check == i,]
                          coordinates(sp.mydata) <- ~y+ x
                          if(nrow(sp.mydata) < 2){
                            nearestdist <- rep(NA, nrow(sp.mydata))
                          }else{
                            nearestdist <- knn(coordinates(sp.mydata), coordinates(sp.mydata), k = 2)$nn.dist[,2] # nearest neighbour distance
                          }
                          min.d <- c(min.d, nearestdist)
                          setTxtProgressBar(pb, k)
                          k <- k+1
                        }
                        data_df$NN_Distance <- min.d

                        NNDist_mean <- aggregate(NN_Distance ~ t+pert.name+rep, data = data_df,
                                                 FUN = mean, na.action = NULL)$NN_Distance
                        NNDist_sd <- aggregate(NN_Distance ~ t+pert.name+rep, data = data_df,
                                               FUN = sd, na.action = NULL)$NN_Distance
                        NNDist_skewness <- aggregate(NN_Distance ~ t+pert.name+rep, data = data_df,
                                                     FUN = e1071::skewness, na.action = NULL)$NN_Distance
                        NNDist_kurtosis <- aggregate(NN_Distance ~ t+pert.name+rep, data = data_df,
                                                     FUN = e1071::kurtosis, na.action = NULL)$NN_Distance

                        ## Final data frame ----
                        return_df <- cbind(N,
                                           data.frame(
                                             u_mean = u_mean,
                                             u_sd = u_sd,
                                             u_skewness = u_skewness,
                                             u_kurtosis = u_kurtosis,
                                             MalAdap_mean = MalAdap_mean,
                                             MalAdap_sd = MalAdap_sd,
                                             MalAdap_skewness = MalAdap_skewness,
                                             MalAdap_kurtosis = MalAdap_kurtosis,
                                             NNDist_mean = NNDist_mean,
                                             NNDist_sd = NNDist_sd,
                                             NNDist_skewness = NNDist_skewness,
                                             NNDist_kurtosis = NNDist_kurtosis
                                           )
                        )
                        save(return_df, file = file.path(Dir.Exports, paste0("TEMP_TStepPop_", x, ".RData")))
                      }
                      return_df
                    })
## Saving & Quitting ------------------------------------------------------
Data_df <- do.call(rbind, Data_ls)
Data_df$ID <- with(Data_df, paste(AC, DI, MU, SL, VA, pert.name, rep, sep = "-"))
save(Data_df, file = file.path(Dir.Exports, "POPULATION_TimeStep.RData"))
unlink(list.files(Dir.Exports, pattern = "TEMP_TStepPop_", full.names = TRUE))

stopCluster(cl)
