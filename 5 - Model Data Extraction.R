#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes]
#' CONTENTS:
#'  - Comparison of landscape parameter distributions and individual trait distributions
#'  DEPENDENCIES:
#'  - EVORES_Metrics.RData produced by "4 - Evolutionary Rescue Success Metrics.R"
#'  - ENVIRONMENT_Cells.RData produced by "3 - Landscape Readout.R"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list = ls())

## Directories ------------------------------------------------------------
### Define dicrectories in relation to project directory
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "Data")
Dir.Exports <- file.path(Dir.Base, "Exports")
### Create directories which aren't present yet
Dirs <- c(Dir.Data, Dir.Exports)
CreateDir <- sapply(Dirs, function(x) if (!dir.exists(x)) dir.create(x))

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
  #   "ggplot2",
  "overlapping",
  "lattice",
  "DescTools",
  "sp",
  "raster",
  "emdist",
  "e1071"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
Base_rast <- raster(list.files(Dir.Data, pattern = ".nc", full.names = TRUE)[1])
load(file.path(Dir.Exports, "EVORES_Metrics.RData")) # load EVORES_Metrics
EVORES_Metrics <- EVORES_Metrics[, -20:-22]
EvoResSuc_df <- EVORES_Metrics[EVORES_Metrics$survival, ]
load(file.path(Dir.Exports, "ENVIRONMENT_Cells.RData")) # load Enviro_Cells_ls

# EvoResSuc_df <- EvoResSuc_df[1:200, ] # just for testing purposes

## Helper: load population file for a row of EvoResSuc_df -----------------
# Accepts either a numeric index (i) or a single-row data.frame/list matching EvoResSuc_df columns
get_pop_df_from_row <- function(row) {
  # resolve numeric index to a single-row data.frame
  if (is.numeric(row)) {
    if (row < 1 || row > nrow(EvoResSuc_df)) stop("row index out of range in get_pop_df_from_row")
    row <- EvoResSuc_df[row, , drop = FALSE]
  }

  # construct file selection using the same logic as the original inline block
  Pop_f <- list.files(Dir.Data, pattern = paste0("AC", row$AC, "_"))
  Pop_f <- Pop_f[grep(paste0("VA", row$VA, ".rds"), Pop_f)]
  SLCheck <- ifelse(row$SL == 0.8, ".8", row$SL)
  SLCheck <- ifelse(row$SL == 1, "1", SLCheck)
  Pop_f <- Pop_f[grep(paste0("SL", SLCheck, "_"), Pop_f)]
  Pop_f <- Pop_f[grep(paste0("DI", row$DI, "_"), Pop_f)]
  Pop_f <- Pop_f[grep(paste0("MU", row$MU, "_"), Pop_f)]

  if (length(Pop_f) != 1) {
    stop(paste(
      "file selection went wrong for",
      if (!is.null(row$ID)) row$ID else paste0("AC=", row$AC, " VA=", row$VA, " SL=", row$SL, " DI=", row$DI, " MU=", row$MU)
    ))
  }
  readRDS(file.path(Dir.Data, Pop_f))
}

# Directory to store per-row temporary cache files
Dir.Temp <- file.path(Dir.Exports, "TEMP_CACHE")
if (!dir.exists(Dir.Temp)) dir.create(Dir.Temp, recursive = TRUE)

# Multi-Cores
nC <- 30 # floor(parallel::detectCores() * 0.5)
# if(Sys.info()['sysname'] == "Linux"){cl = nC}
cl <- makeCluster(nC)
clusterExport(cl = cl, varlist = c("Dir.Data", "install.load.package", "package_vec", "get_pop_df_from_row", "Base_rast", "Dir.Temp"))
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

# ANALYSES ================================================================
## Ouctome Metrics --------------------------------------------------------
### Resistance
#### Immediate
EvoResSuc_df$Res_Mag_Immediate <- EvoResSuc_df$n_post/EvoResSuc_df$n_pre
#### Delayed
EvoResSuc_df$Res_Mag_Delayed <- EvoResSuc_df$n_minpost/EvoResSuc_df$n_pre
EvoResSuc_df$Res_Speed_Delayed <- abs(EvoResSuc_df$t_minpost-460)/EvoResSuc_df$t

### Recovery
#### Pre-Baseline
EvoResSuc_df$Rec_Mag_PreBase <- EvoResSuc_df$n_maxpost/EvoResSuc_df$n_pre
EvoResSuc_df$Rec_Speed_PreBase <- abs(EvoResSuc_df$t_maxpost-460)/(EvoResSuc_df$t-460)
#### Post-Baseline
EvoResSuc_df$Rec_Mag_PostBase <- EvoResSuc_df$n_maxpost/EvoResSuc_df$n_minpost
EvoResSuc_df$Rec_Speed_PostBase <- abs(EvoResSuc_df$t_maxpost-EvoResSuc_df$t_minpost)/(EvoResSuc_df$t-EvoResSuc_df$t_minpost)

EvoResSuc_df <- EvoResSuc_df[,c("ID", "pert.name", "rep", "AC", "DI", "MU", "SL", "VA", "t", "survival", "n_pre", "n_post", "n_minpost", "t_minpost", "n_maxpost", "t_maxpost", "Res_Mag_Immediate", "Res_Mag_Delayed", "Res_Speed_Delayed", "Rec_Mag_PreBase", "Rec_Speed_PreBase", "Rec_Mag_PostBase", "Rec_Speed_PostBase")]

clusterExport(cl = cl, varlist = c("Enviro_Cells_ls", "EvoResSuc_df"))

## Genetic Diversity ------------------------------------------------------
### Traits of Individuals alive
message("Traits of individuals alive")
Traits_ls <- pblapply(seq_len(nrow(EvoResSuc_df)), cl = cl, FUN = function(i) {
  # i <- 1
  cache_file <- file.path(Dir.Temp, paste0("Traits_row_", i, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  ## file selection (delegated to helper)
  Pop_df <- get_pop_df_from_row(i)

  ## selecting data in file
  Popb_df <- Pop_df[Pop_df$rep == EvoResSuc_df$rep[i] & Pop_df$pert.name == EvoResSuc_df$pert.name[i], ]

  ## extracting trait information
  time_points <- c(pre = 460, minpost = EvoResSuc_df[i, "t_minpost"], maxpost = EvoResSuc_df[i, "t_maxpost"])
  ret_ls <- lapply(seq_along(time_points), FUN = function(tp) {
    # tp <- 1
    t_iter <- time_points[tp]
    Indivs_t <- Popb_df[Popb_df$t == t_iter, ]
    c(SD = sd(Indivs_t$u), Unique = length(unique(Indivs_t$u)) / nrow(Indivs_t))
  })
  names(ret_ls) <- names(time_points)
  ret_vec <- unlist(ret_ls)
  names(ret_vec) <- paste0("Trait.", names(ret_vec))

  saveRDS(ret_vec, cache_file)
  ret_vec
})
EvoResSuc_df <- cbind(EvoResSuc_df, do.call(rbind, Traits_ls))

### Trait landscape (mean traits of individuals per cell in environment); Variance, Autocorrelation, Slope
message("Traits landscapes")
TraitLandscape_ls <- pblapply(seq_len(nrow(EvoResSuc_df)), cl = cl, FUN = function(i) {
  # i <- 1
  print(i)
  cache_file <- file.path(Dir.Temp, paste0("TraitLandscape_row_", i, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  ## file selection (delegated to helper)
  Pop_df <- get_pop_df_from_row(i)

  ## selecting data in file
  Popb_df <- Pop_df[Pop_df$rep == EvoResSuc_df$rep[i] & Pop_df$pert.name == EvoResSuc_df$pert.name[i], ]

  ## extracting trait information
  time_points <- c(pre = 460, minpost = EvoResSuc_df[i, "t_minpost"], maxpost = EvoResSuc_df[i, "t_maxpost"])
  ret_ls <- lapply(seq_along(time_points), FUN = function(tp) {
    # tp <- 1
    ## subsetting
    t_iter <- time_points[tp]
    Indivs_t <- Popb_df[Popb_df$t == t_iter, ]

    ## rasterising
    coordinates(Indivs_t) <- ~ x + y
    Indivs_ras <- rasterize(x = Indivs_t[, "u"], y = Base_rast, fun = mean)$u

    ## calculating metrics
    autocor_val <- Moran(Indivs_ras)
    variance_val <- var(as.vector(Indivs_ras), na.rm = TRUE)
    slope_val <- coef(lm(u ~ x + y, data = Indivs_t))["x"]
    c(AC = autocor_val, VA = variance_val, SL = slope_val)
  })
  names(ret_ls) <- names(time_points)
  ret_vec <- unlist(ret_ls)
  names(ret_vec) <- paste0("TraitScape.", names(ret_vec))
  saveRDS(ret_vec, cache_file)
  ret_vec
})
EvoResSuc_df <- cbind(EvoResSuc_df, do.call(rbind, TraitLandscape_ls))

## Genetic Diversity in Environment ---------------------------------------
### Average dif between ind trait and env value
message("Traits vs Environments")
Adaptedness_ls <- pblapply(seq_len(nrow(EvoResSuc_df)), # cl = cl,
  FUN = function(i) {
    # i <- 1
    # print(i)
    cache_file <- file.path(Dir.Temp, paste0("Adaptedness_row_", i, ".rds"))
    if (file.exists(cache_file)) {
      return(readRDS(cache_file))
    }
    ## file selection (delegated to helper)
    Pop_df <- get_pop_df_from_row(i)

    ## selecting data in file
    Popb_df <- Pop_df[Pop_df$rep == EvoResSuc_df$rep[i] & Pop_df$pert.name == EvoResSuc_df$pert.name[i], ]

    ## selecting environment raster
    dat_rast <- raster(file.path(Dir.Data, paste0("ac", EvoResSuc_df$AC[i], "var", EvoResSuc_df$VA[i], ".nc")))
    sl_vec <- (1:100) * EvoResSuc_df$SL[i]
    sl_mat <- matrix(rep(sl_vec, 100), nrow = 100, ncol = 100, byrow = TRUE)
    sl_rast <- raster(sl_mat, xmn = -0.5, xmx = 99.5, ymn = -0.5, ymx = 99.5)
    Enviro_ras <- dat_rast + sl_rast

    ## extracting trait information
    time_points <- c(pre = 460, minpost = EvoResSuc_df[i, "t_minpost"], maxpost = EvoResSuc_df[i, "t_maxpost"])
    ret_ls <- lapply(seq_along(time_points), FUN = function(tp) {
      # tp <- 1
      ## subsetting
      t_iter <- time_points[tp]
      Indivs_t <- Popb_df[Popb_df$t == t_iter, ]

      ## rasterising
      coordinates(Indivs_t) <- ~ x + y

      if (grepl(pattern = "post", names(time_points))[tp]) {
        Enviro_ras2 <- Enviro_ras + EvoResSuc_df$pert.name[i]
      } else {
        Enviro_ras2 <- Enviro_ras
      }

      extracted_env_values <- extract(Enviro_ras2, Indivs_t)
      mean(abs(Indivs_t$u - extracted_env_values))
    })
    names(ret_ls) <- names(time_points)
    ret_vec <- unlist(ret_ls)
    names(ret_vec) <- paste0("Maladap.", names(ret_vec))
    saveRDS(ret_vec, cache_file)
    ret_vec
  }
)
EvoResSuc_df <- cbind(EvoResSuc_df, do.call(rbind, Adaptedness_ls))

### Earth mover distance (pre pert genetic values and post pert env values)
message("Earth Mover Distance")
EMD_ls <- pblapply(seq_len(nrow(EvoResSuc_df)), cl = cl, FUN = function(i) {
  # i <- 1
  # print(i)
  cache_file <- file.path(Dir.Temp, paste0("EMD_row_", i, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  ## file selection (delegated to helper)
  Pop_df <- get_pop_df_from_row(i)

  ## selecting data in file
  Popb_df <- Pop_df[Pop_df$rep == EvoResSuc_df$rep[i] & Pop_df$pert.name == EvoResSuc_df$pert.name[i], ]

  ## selecting environment raster
  dat_rast <- raster(file.path(Dir.Data, paste0("ac", EvoResSuc_df$AC[i], "var", EvoResSuc_df$VA[i], ".nc")))
  sl_vec <- (1:100) * EvoResSuc_df$SL[i]
  sl_mat <- matrix(rep(sl_vec, 100), nrow = 100, ncol = 100, byrow = TRUE)
  sl_rast <- raster(sl_mat, xmn = -0.5, xmx = 99.5, ymn = -0.5, ymx = 99.5)
  Enviro_ras <- dat_rast + sl_rast + EvoResSuc_df$pert.name[i]

  ## extracting trait information
  Indivs_t <- Popb_df[Popb_df$t == 460, ]

  ## rasterising
  coordinates(Indivs_t) <- ~ x + y
  Indivs_ras <- rasterize(x = Indivs_t[, "u"], y = Base_rast, fun = mean)$u

  # Patch_ras <- rasterize(x = Indivs_t[, "patch"], y = Base_rast, fun = mean)[[2]]
  # Patch_ras[is.na(Patch_ras)] <- 0

  # EMD
  set.seed(42)
  EMD2D <- emd2d(matrix(Indivs_ras, ncol = 100), matrix(Enviro_ras #+ Patch_ras
    ,
    ncol = 100
  ), max.iter = 5e2) # may want to set max.iter = 1 for testing

  saveRDS(EMD2D, cache_file)
  EMD2D
})
EvoResSuc_df$EMD <- unlist(EMD_ls)

## Genetic Diversity over Time --------------------------------------------
### Recovery evolution
message("Recovery evolution")
RecovEvo_ls <- pblapply(seq_len(nrow(EvoResSuc_df)), 
cl = cl, 
FUN = function(i) {
  # i <- 1
  # print(i)
  cache_file <- file.path(Dir.Temp, paste0("RecoveryEvo_row_", i, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  ## file selection (delegated to helper)
  Pop_df <- get_pop_df_from_row(i)

  ## selecting data in file
  Popb_df <- Pop_df[Pop_df$rep == EvoResSuc_df$rep[i] & Pop_df$pert.name == EvoResSuc_df$pert.name[i], ]

  time_points <- data.frame(
    Comparison = c("GenDivChange.preVSmaxpost", "GenDivChange.minpostVSmaxpost"),
    baseline = c(460, EvoResSuc_df[i, "t_minpost"]),
    compare = rep(EvoResSuc_df[i, "t_maxpost"], 2)
  )
  ret_ls <- lapply(1:nrow(time_points), FUN = function(tp) {
    # tp <- 1
    ## subsetting
    t_iter <- time_points[tp, ]
    # message(t_iter$compare)

    Indivs_baseline <- Popb_df[Popb_df$t == t_iter$baseline, "u"]
    Indivs_compare <- Popb_df[Popb_df$t == t_iter$compare, "u"]

    c(
      SD = abs(sd(Indivs_compare) - sd(Indivs_baseline)) / (t_iter$compare - t_iter$baseline),
      mean = abs(mean(Indivs_compare) - mean(Indivs_baseline)) / (t_iter$compare - t_iter$baseline)
    )
  })
  names(ret_ls) <- time_points$Comparison
  ret_vec <- unlist(ret_ls)
  saveRDS(ret_vec, cache_file)
  ret_vec
})
EvoResSuc_df <- cbind(EvoResSuc_df, do.call(rbind, RecovEvo_ls))

### Resistance evolution
message("Resistance evolution")
ResisEvo_ls <- pblapply(seq_len(nrow(EvoResSuc_df)), cl = cl, FUN = function(i) {
  # i <- 1
  # print(i)
  cache_file <- file.path(Dir.Temp, paste0("ResistanceEvo_row_", i, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  ## file selection (delegated to helper)
  Pop_df <- get_pop_df_from_row(i)

  ## selecting data in file
  Popb_df <- Pop_df[Pop_df$rep == EvoResSuc_df$rep[i] & Pop_df$pert.name == EvoResSuc_df$pert.name[i], ]

  time_points <- data.frame(
    Comparison = c("GenDivChange.preVSpost", "GenDivChange.preVSminpost"),
    baseline = rep(460, 2),
    compare = c(470, EvoResSuc_df[i, "t_minpost"])
  )
  ret_ls <- lapply(1:nrow(time_points), FUN = function(tp) {
    # tp <- 1
    ## subsetting
    t_iter <- time_points[tp, ]

    Indivs_baseline <- Popb_df[Popb_df$t == t_iter$baseline, "u"]
    Indivs_compare <- Popb_df[Popb_df$t == t_iter$compare, "u"]

    c(
      SD = abs(sd(Indivs_compare) - sd(Indivs_baseline)) / (t_iter$compare - t_iter$baseline),
      mean = abs(mean(Indivs_compare) - mean(Indivs_baseline)) / (t_iter$compare - t_iter$baseline)
    )
  })
  names(ret_ls) <- time_points$Comparison
  ret_vec <- unlist(ret_ls)
  saveRDS(ret_vec, cache_file)
  ret_vec
})
EvoResSuc_df <- cbind(EvoResSuc_df, do.call(rbind, ResisEvo_ls))


## Spatial Population Measures --------------------------------------------
message("Spatial Population Measures")
SpatPop_ls <- pblapply(
  seq_len(nrow(EvoResSuc_df)), 
  cl = cl,
  FUN = function(i) {
  # i <- 1
  # print(i)
  cache_file <- file.path(Dir.Temp, paste0("SpatialPopEvo_row_", i, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  ## file selection (delegated to helper)
  Pop_df <- get_pop_df_from_row(i)

  ## selecting data in file
  Popb_df <- Pop_df[Pop_df$rep == EvoResSuc_df$rep[i] & Pop_df$pert.name == EvoResSuc_df$pert.name[i], ]

  ## base environment raster
  dat_rast <- raster(file.path(Dir.Data, paste0("ac", EvoResSuc_df$AC[i], "var", EvoResSuc_df$VA[i], ".nc")))
  sl_vec <- (1:100) * EvoResSuc_df$SL[i]
  sl_mat <- matrix(rep(sl_vec, 100), nrow = 100, ncol = 100, byrow = TRUE)
  sl_rast <- raster(sl_mat, xmn = -0.5, xmx = 99.5, ymn = -0.5, ymx = 99.5)
  Enviro_ras <- dat_rast + sl_rast

  ## distances
  time_points <- c(pre = 460, minpost = EvoResSuc_df[i, "t_minpost"], maxpost = EvoResSuc_df[i, "t_maxpost"])
  ret_ls <- lapply(seq_along(time_points), FUN = function(tp) {
    # tp <- 1
    ## subsetting
    t_iter <- time_points[tp]
    Indivs_t <- Popb_df[Popb_df$t == t_iter, ]
    Indivs_sp <- Indivs_t

    ## rasterising
    coordinates(Indivs_sp) <- ~ x + y
    Indivs_ras <- rasterize(x = Indivs_sp[, "u"], y = Base_rast, fun = mean)$u

    if (grepl(pattern = "post", names(time_points))[tp]) {
      Enviro_ras2 <- Enviro_ras + EvoResSuc_df$pert.name[i]
    } else {
      Enviro_ras2 <- Enviro_ras
    }

    ## distance to best environment per individual
    Dist.indivVSenvir2 <- sapply(1:nrow(Indivs_t), FUN = function(ind) {
      # ind <- 1
      indiv_coord <- Indivs_t[ind, c("x", "y")]
      indiv_trait <- Indivs_t[ind, "u"]

      # cells with matching environment
      match_cell <- which.min(abs(values(Enviro_ras2) - indiv_trait))[1]

      match_coords <- xyFromCell(Enviro_ras2, match_cell)
      dists <- spDistsN1(match_coords, as.numeric(indiv_coord), longlat = FALSE)
      dists
    })
    Dist.indivVSenvir <- c(
      Mean = mean(Dist.indivVSenvir2),
      SD = sd(Dist.indivVSenvir2),
      Skew = skewness(Dist.indivVSenvir2)
    )

    ## distance to nearest neighbour per individual
    # Compute pairwise distances
    D <- as.matrix(dist(Indivs_t[, c("x", "y")]))

    # Replace diagonal (self-distance) with NA or Inf
    diag(D) <- NA

    # For each individual, get nearest neighbour and distance
    nearest_distance <- apply(D, 1, min, na.rm = TRUE)

    Dist.indivVSindiv <- c(
      Mean = mean(nearest_distance),
      SD = sd(nearest_distance),
      Skew = skewness(nearest_distance)
    )

    ## cell filling
    populated <- sum(!is.na(values(Indivs_ras))) / ncell(Indivs_ras)

    vals <- getValues(Indivs_ras)
    cells <- which(!is.na(vals))

    # Get xy coordinates of those cells
    xy <- xyFromCell(Indivs_ras, cells)

    centroid <- c(mean(xy[, 1]), mean(xy[, 2]), sd(xy[, 1]))

    SpatArrang <- c(
      Populated = populated,
      Centroid.x.mean = centroid[1],
      Centroid.y.mean = centroid[2],
      Centroid.x.sd = centroid[3],
      Centroid.y.sd = centroid[4]
    )


    ## return
    list(
      Dist.indivVSenvir = Dist.indivVSenvir,
      Dist.indivVSindiv = Dist.indivVSindiv,
      SpatArrang = SpatArrang
    )

  })

  names(ret_ls) <- names(time_points)
  ret_vec <- unlist(ret_ls)
  saveRDS(ret_vec, cache_file)
  ret_vec
})
EvoResSuc_df <- cbind(EvoResSuc_df, do.call(rbind, SpatPop_ls))

saveRDS(EvoResSuc_df, file = file.path(Dir.Exports, "ModelData.rds"))
head(EvoResSuc_df)
