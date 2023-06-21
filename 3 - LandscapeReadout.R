#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Landscape metric reading
#'  DEPENDENCIES:
#'  - 
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
Enviro_fs <- list.files(Dir.Data, pattern = ".txt")
Slopes <- c(0.8, 1, 1.2)

Enviro_table <- do.call(rbind, 
        lapply(X = Enviro_fs, FUN = function(i){
          
          dat_df <- read.table(file.path(Dir.Data, i))
          dat_mat <- matrix(NA, ncol = 100, nrow = 100, byrow = FALSE)
          for(asd in 0:99){
            dat_mat[nrow(dat_mat)-asd, ] <- dat_df$V1[(asd*100+1):(asd*100+100)]
          }
          dat_rast <- raster(dat_mat, xmn = -0.5, xmx = 99.5, ymn = -0.5, ymx = 99.5)
          
          # Enviro_ras <- raster(nrows=100, ncols=100, 
          #                      xmn=-0.5, xmx=99.5, ymn=-0.5, ymx=99.5, 
          #                      vals = base::rev(dat_df$V1)
          # )
          # Enviro_ras <- flip(Enviro_ras, direction = "x")
          writeRaster(dat_rast, filename = file.path(Dir.Data, tools::file_path_sans_ext(i)), format = "CDF", overwrite = TRUE)
          
          do.call(rbind, 
                  lapply(X = Slopes, FUN = function(SL){
                    
                    sl_vec <- (1:100)*SL
                    sl_mat <- matrix(rep(sl_vec, 100), nrow = 100, ncol = 100, byrow = TRUE)
                    # sl_rast <- raster(sl_mat, xmn = 0, xmx = 99, ymn = 0, ymx = 99)
                    sl_rast <- raster(sl_mat, xmn=-0.5, xmx=99.5, ymn=-0.5, ymx=99.5)
                    # enviro_rast <- dat_rast+sl_rast
                    enviro_rast <- dat_rast+sl_rast
                    (autocor_val <- Moran(enviro_rast))
                    (variance_val <- var(as.vector(enviro_rast))) 
                    
                    saveRDS(values(enviro_rast), file = file.path(Dir.Exports, paste(as.numeric(gsub(".*?([0-9]+).*", "\\1", i)),
                                                                                     as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                                                                                     gsub(pattern = ".*var", replacement = "", x = i))
                                                                                     ),
                                                                                     SL, ".rds",
                                                                                     sep = "_")
                    ))
                    
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
Enviro_table

boxplot(AC_data ~ AC_input, data = Enviro_table)
boxplot(VA_data ~ VA_input, data = Enviro_table)

Histo_fs <- list.files(Dir.Exports, pattern = ".rds")
EnvirDist_ls <- lapply(Histo_fs, FUN = function(x){
  readRDS(file.path(Dir.Exports,x))
})
names(EnvirDist_ls) <- Histo_fs # names are AC_VA_SL_rds
save(EnvirDist_ls, file = file.path(Dir.Exports, "EnviroDist_ls.RData"))
unlink(file.path(Dir.Exports, Histo_fs), recursive = TRUE)

## Exemplary quantification of patch value from coordinates
# dat_df <- read.table(file.path(Dir.Data, "ac0var25.txt"))
# ID_df <- readRDS(file.path(Dir.Data, "AC0_DI2_MU1_SL1_VA25.rds"))
# ID_df <- ID_df[ID_df$t==310,]
# summary(ID_df$patch)
# 
# cell10_df <- ID_df[ID_df$x< 10 & ID_df$x>9.5 & ID_df$y< 10 & ID_df$y>9.5, ]
# hist(cell10_df$patch)
# 
# firstx <- newy <- 10 
# newx <- firstx+newy*100
# if(newx == 1e4){newx <- 9999}
# (Text <- dat_df$V1[newx+1]) # plus 1 because c++ index starts at 0
# unique(cell10_df$patch)

