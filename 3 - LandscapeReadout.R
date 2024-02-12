#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Landscape metric reading; simulations where initialised with specific parameters for landscape characteristics, but randomness in simulation modulated them; analyses should be run on the actual values
#'  DEPENDENCIES:
#'  - *.txt files containing environment specifications in Data folder
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
  "raster",
  "pbapply",
  "readr"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
Enviro_fs <- list.files(Dir.Data, pattern = ".txt") # .txt files containing values of cells in Environment
Slopes <- c(0.8, 1, 1.2) # slope values that were set for modulation of environments and alter cell values stored in .txt files

## loop over .txt files and make a data frame out of them
Enviro_table <- do.call(rbind, 
        lapply(X = Enviro_fs, FUN = function(i){
          ## reading data in ----
          dat_df <- read.table(file.path(Dir.Data, i))
          
          ## make data into a matrix and fix differences in indexing between C and R ----
          dat_mat <- matrix(NA, ncol = 100, nrow = 100, byrow = FALSE)
          for(asd in 0:99){
            dat_mat[nrow(dat_mat)-asd, ] <- dat_df$V1[(asd*100+1):(asd*100+100)]
          }
          dat_rast <- raster(dat_mat, xmn = -0.5, xmx = 99.5, ymn = -0.5, ymx = 99.5)
          
          ## write "fixed" raster ----
          writeRaster(dat_rast, filename = file.path(Dir.Data, tools::file_path_sans_ext(i)), format = "CDF", overwrite = TRUE)
          
          ## loop over slopes and extract landscape metrics and make a data frame ----
          do.call(rbind, 
                  lapply(X = Slopes, FUN = function(SL){
                    ### slopes affect cell values
                    sl_vec <- (1:100)*SL
                    sl_mat <- matrix(rep(sl_vec, 100), nrow = 100, ncol = 100, byrow = TRUE)
                    sl_rast <- raster(sl_mat, xmn=-0.5, xmx=99.5, ymn=-0.5, ymx=99.5)
                    enviro_rast <- dat_rast+sl_rast
                    
                    ## get metrics
                    (autocor_val <- Moran(enviro_rast))
                    (variance_val <- var(as.vector(enviro_rast))) 
                    
                    ## save temporary files with raster values
                    saveRDS(values(enviro_rast), file = file.path(Dir.Exports, paste(as.numeric(gsub(".*?([0-9]+).*", "\\1", i)),
                                                                                     as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                                                                                     gsub(pattern = ".*var", replacement = "", x = i))
                                                                                     ),
                                                                                     SL, ".rds",
                                                                                     sep = "_")
                    ))
                    
                    ## return landscape metrics as data frame
                    data.frame(name = i,
                               SL = as.numeric(SL),
                               AC_input = as.numeric(gsub(".*?([0-9]+).*", "\\1", i)),
                               AC_data = as.numeric(autocor_val),
                               VA_input = as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                                          gsub(pattern = ".*var", replacement = "", x = i))
                                                     ),
                               VA_data = as.numeric(variance_val))
                    
                  })
                  )
        })
      )
## save environment characteristics table
save(Enviro_table, file = file.path(Dir.Exports, "ENVIRONMENT_Parameters.RData"))

## load raster values; save as list; then delete temporary files
Histo_fs <- list.files(Dir.Exports, pattern = "_.rds")
Enviro_Cells_ls <- lapply(Histo_fs, FUN = function(x){
  readRDS(file.path(Dir.Exports,x))
})
names(Enviro_Cells_ls) <- Histo_fs # names are AC_VA_SL_rds
save(Enviro_Cells_ls, file = file.path(Dir.Exports, "ENVIRONMENT_Cells.RData"))
unlink(file.path(Dir.Exports, Histo_fs), recursive = TRUE)
