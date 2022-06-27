#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Data Cleaning
#'  DEPENDENCIES:
#'  - RenamingScript.R (located in Data folder executed)
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
  "parallel" # for parallel computation
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================

## Part 2 files -----------------------------------------------------------
part2_fs <- list.files(Dir.Data, pattern = "part2.txt")
if(length(part2_fs)>0){
  for(name in part2_fs){
    
    name2 <- list.files(Dir.Data, pattern = strsplit(name, split = ".part2")[[1]][1])[2]
    print(paste("Fusing", name2, "and", name))
    data <- read.csv(file = file.path(Dir.Data, name2))
    data_2 <- read.csv(file = file.path(Dir.Data, name))
    names(data) <- names(data_2) <- c("pert.value", "pert.name", "rep", "t","n","x","y","u","id","patch") # assigning variable names
    
    if(min(unique(data_2$rep)) < 100){data_2$rep <- data_2$rep+100} 
    data <- rbind(data, data_2)
    
    write.table(data, file = file.path(Dir.Data, name2))
    unlink(file.path(Dir.Data, name))
  }
  
}

## Duplicate Replicates ---------------------------------------------------
cl <- parallel::makeCluster(parallel::detectCores())
parallel::clusterExport(cl,
                        varlist = c("Dir.Data", "install.load.package", "package_vec"),
                        envir = environment()
)
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

Sims_fs <- list.files(Dir.Data, ".txt")

pblapply(Sims_fs, 
         cl = cl,
         function(x){
  # x <- Sims_fs[1]
           # message(x)
  
  data_df <- na.omit(read.table(file = file.path(Dir.Data, x)))
  if(ncol(data_df) != 10){
    data_df <- na.omit(read.csv(file = file.path(Dir.Data, x), header = FALSE))
    names(data_df) <- c("pert.value", "pert.name", "rep", "t","n","x","y","u","id","patch") # assigning variable names
  }
  
  ## Row Removal ----
  RowPertNA <- which(is.na(as.numeric(data_df$pert.value))) # remove rows where pert.name is a character value
  RowPertTooBig <- which(as.numeric(data_df$pert.value) > 16) # remove rows where pert.value is too big (these are usually) shifted columns
  if(length(c(RowPertNA, RowPertTooBig)) != 0){
    data_df <- data_df[-c(RowPertNA, RowPertTooBig), ] # remove error rows
  }

  ## Duplicate Removal ----
  RowsDup <- c()
  PosNewSim <- which(c(TRUE, tail(data_df$t,-1) < head(data_df$t,-1))) # identify where time decreases
  NewSim_df <- data_df[PosNewSim, ]
  Identifiersdup <- paste(NewSim_df$pert.name, NewSim_df$rep, sep="_")
  duplicates <- which(duplicated(Identifiersdup))
  if(length(duplicates) != 0){
    duplicates <- unique(Identifiersdup[duplicates])
    for(i in duplicates){
      ## identify instances of duplicate
      PosNewSim <- which(Identifiersdup == i)
      ## identify run time of each duplicate
      finRow <- as.numeric(rownames(NewSim_df)[(PosNewSim+1)])-1 # final row of each duplicates run
      if(is.na(finRow[length(finRow)])){finRow[length(finRow)] <- nrow(data_df)}
      ts <- data_df$t[finRow] # maximum run time of duplicate replicates
      ## identify shorter time
      longestrep <- which.max(ts) # if these are the same, the first will be selected
      ## identify row in newsim data frame that corresponds to replicates that ought to be deleted
      begin <- as.numeric(rownames(NewSim_df)[PosNewSim[-longestrep]])
      end <- as.numeric(rownames(NewSim_df)[(PosNewSim[-longestrep]+1)])-1
      if(is.na(end[length(end)])){end[length(end)] <- nrow(data_df)} # this catches when we have to delete the last sim in the file
      RowsDupApp <- unlist(sapply(1:length(begin), function(l){
        begin[l]:end[l]
      }))
      ## export rows that need to be removed
      RowsDup <- c(RowsDup, RowsDupApp)
    }
    data_df <- data_df[-RowsDup,] # removing duplicate rows
  }
  saveRDS(data_df, file = file.path(Dir.Data, 
                                    paste0(tools::file_path_sans_ext(x), ".rds"))
          )
  unlink(file.path(Dir.Data, x))
  })

parallel::stopCluster(cl)
closeAllConnections()
