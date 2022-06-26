#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Resistance Analyses
#'  DEPENDENCIES:
#'  - None
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
  "data.table" # for rbindlist
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
## Loading & Extraction ---------------------------------------------------
message("Currently only analysing first 50 simulation runs")
Sims_fs <- list.files(Dir.Data, ".txt")[1:50]

message("duplicate reps need to be accounted for/removed. Is there a replicate? Yes - does it go to maximum time (that we set) or does population size go to 0? Can we just pick the longer of duplicate replicates?")

message("some files have characters in the pert.value column - e.g: Sims_fs[50] in row 931860")

## reading data and extracting relevant information
Pert_ls <- pblapply(Sims_fs, function(x){
  ## identifiers for simulation run
  Ident_vec <- unlist(strsplit(tools::file_path_sans_ext(x), "_")) # vector of identifiers
  Ident_vec <- as.numeric(gsub("[^0-9.-]", "", Ident_vec)) # retain only numbers
  
  ## data within simulation data
  data <- read.csv(file = file.path(Dir.Data, x)) # reading data
  names(data) <- c("pert.value", "pert.name", "rep", "t","n","x","y","u","id","patch") # assigning variable names
  PertPos <- which(c(FALSE, tail(data$pert.value,-1) > head(data$pert.value,-1))) # identify where perturbations increase
  PercLoss <- data$n[PertPos]/data$n[PertPos-1] # loss of population size in % of time-step before perturbation hit
  PertInc <- as.numeric(data$pert.value[PertPos])-
    as.numeric(data$pert.value[PertPos-1]) # increase in perturbation
  
  ## exporting of extracted data
  data.frame(Loss = PercLoss,
             Perturbation = PertInc,
             Autocorrelation = rep(Ident_vec[1], length(PercLoss)),
             Dispersal = rep(Ident_vec[2], length(PercLoss)),
             Mutation = rep(Ident_vec[3], length(PercLoss)),
             Slope = rep(Ident_vec[4], length(PercLoss)),
             Variance = rep(Ident_vec[5], length(PercLoss))
  )
})
Pert_df <- rbindlist(Pert_ls)
save(Pert_df, file = file.path(Dir.Data, "ResistanceDFFirst50.RData"))

load(file.path(Dir.Data, "ResistanceDFFirst50.RData"))

Pert_df <- Pert_df[Pert_df$Perturbation < 17, ]
Pert_df <- Pert_df[Pert_df$Loss < 1, ]

ggplot(Pert_df, aes(x = Perturbation, y = Loss, col = factor(Mutation))) + 
  geom_point() + 
  geom_smooth() + 
  theme_bw()

























