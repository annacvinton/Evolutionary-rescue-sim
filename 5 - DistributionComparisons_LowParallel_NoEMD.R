#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Comparison of landscape parameter distributions and individual trait distributions
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
# rm(list=ls())

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
  "overlapping",
  "lattice",
  "DescTools",
  "sp",
  "raster"
  # ,
  # "emdist"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
EvoResSuc_df <- read.csv(file = file.path(Dir.Exports, "EvoResSuccessMetrics.csv"))
load(file = file.path(Dir.Exports, "EnviroDist_ls.RData")) # load EnvirDist_ls

## Non-Spatially Explicit Overlaps ----------------------------------------
if(file.exists(file.path(Dir.Exports, "Distrib_NonSpatial_df.csv"))){
  Overlap_df <- read.csv(file.path(Dir.Exports, "Distrib_NonSpatial_df.csv"))
}else{
  nC <- parallel::detectCores()
  cl <- makeCluster(nC)
  clusterExport(cl = cl, varlist = c("Dir.Data", "install.load.package", "package_vec"))
  clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))
  clusterExport(cl = cl, varlist = c("EnvirDist_ls", "EvoResSuc_df"))
  Overlap_ls <- pblapply(X = names(EnvirDist_ls),
                         cl = cl,
                         FUN = function(x){
                           # x <- names(EnvirDist_ls)[1]
                           message(paste("EnvirDistName =", x))
                           
                           # Environment Data
                           AC <- as.numeric(unlist(strsplit(x = x, split = "_"))[1])
                           VA <- as.numeric(unlist(strsplit(x = x, split = "_"))[2])
                           SL <- as.numeric(unlist(strsplit(x = x, split = "_"))[3])
                           
                           # Individual Data
                           Pop_fs <- list.files(Dir.Data, pattern = ".rds")
                           Iter_fs <- Pop_fs[grep(paste0("AC", AC, "_"), Pop_fs)]
                           Iter_fs <- Iter_fs[grep(paste0("VA", VA, ".rds"), Iter_fs)]
                           Iter_fs <- Iter_fs[grep(paste0("SL", ifelse(SL == 0.8, ".8", SL), "_"), Iter_fs)]
                           
                           Iter_ls <- lapply(Iter_fs, FUN = function(y){
                             # y <- Iter_fs[1] # simulation run loop here
                             print(paste("File =", y))
                             Pop_df <- readRDS(file.path(Dir.Data, y))
                             Characteristics <- as.numeric(unlist(regmatches(unlist(strsplit(x = y, split = "_")),
                                                                             gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                                                      unlist(strsplit(x = y, split = "_"))))
                             ))
                             DI <- Characteristics[2]
                             MU <- Characteristics[3]
                             
                             Pert_ls <- lapply(unique(Pop_df$pert.name), FUN = function(z){
                               # z <- unique(Pop_df$pert.name)[1] # perturbation magnitude loop here
                               print(paste("Pert =", z))
                               
                               EvoResSuc_iter <- EvoResSuc_df[EvoResSuc_df$pert.name == z & 
                                                                EvoResSuc_df$AC == AC & 
                                                                EvoResSuc_df$SL == SL & 
                                                                EvoResSuc_df$VA == VA & 
                                                                EvoResSuc_df$DI == DI &
                                                                EvoResSuc_df$MU == MU,]
                               Pop_df <- Pop_df[Pop_df$pert.name == z,] 
                               
                               rep_ls <- lapply(unique(EvoResSuc_iter$rep), FUN = function(b){
                                 # b <- EvoResSuc_iter$rep[1]
                                 print(paste("Rep =", b))
                                 
                                 Popb_df <- Pop_df[Pop_df$rep == b,]
                                 
                                 t_min <- EvoResSuc_iter$t_minpost[which(EvoResSuc_iter$rep == b)]
                                 t_max <- EvoResSuc_iter$t_maxpost[which(EvoResSuc_iter$rep == b)]
                                 comparison <- data.frame(Envir = c("Pre", "Post", "Post"),
                                                          Traits = c("Pre", "PostMin", "PostMax"))
                                 
                                 plot_df <- data.frame(Value = c(
                                   EnvirDist_ls[[x]],
                                   EnvirDist_ls[[x]]+z,
                                   Popb_df$u[Popb_df$t == 460],
                                   Popb_df$u[Popb_df$t == t_min],
                                   Popb_df$u[Popb_df$t == t_max]
                                 ),
                                 Origin = c(
                                   rep("Environment", (length(EnvirDist_ls[[x]]))*2),
                                   rep("Individuals", length(Popb_df$u[Popb_df$t == 460])),
                                   rep("Individuals", length(Popb_df$u[Popb_df$t == t_min])),
                                   rep("Individuals", length(Popb_df$u[Popb_df$t == t_max]))
                                 ),
                                 Time = c(
                                   rep("Pre", length(EnvirDist_ls[[x]])),
                                   rep("Post", length(EnvirDist_ls[[x]])),
                                   rep("Pre", length(Popb_df$u[Popb_df$t == 460])),
                                   rep("PostMin", length(Popb_df$u[Popb_df$t == t_min])),
                                   rep("PostMax", length(Popb_df$u[Popb_df$t == t_max]))
                                 )
                                 )
                                 
                                 OVERLAP <- c()
                                 for(i in 1:nrow(comparison)){
                                   if(length(plot_df$Value[plot_df$Origin == "Individuals" & plot_df$Time == comparison$Traits[i]])<2 # These come from cases where the population went down to just one individual
                                      | 
                                      comparison$Traits[i] == "PostMax" & is.na(t_max) # No postmaximu identified
                                   ){
                                     OV <- NA 
                                   }else{
                                     OV <- overlap(list(X2 = plot_df$Value[plot_df$Origin == "Environment" & plot_df$Time == comparison$Envir[i]],
                                                        X1 = plot_df$Value[plot_df$Origin == "Individuals" & plot_df$Time == comparison$Traits[i]]), 
                                                   plot = FALSE)$OV 
                                   }
                                   OVERLAP <- c(OVERLAP, 
                                                OV
                                   )
                                 }
                                 
                                 # 2 if trait value is within optimal range given the environment, 1 if not
                                 OV_perc <- c()
                                 for(i in 1:nrow(comparison)){
                                   CheckOverlap <- (plot_df$Value[plot_df$Origin == "Individuals" & plot_df$Time == comparison$Traits[i]] <= 
                                                      max(plot_df$Value[plot_df$Origin == "Environment" & plot_df$Time == comparison$Envir[i]])) + 
                                     (plot_df$Value[plot_df$Origin == "Individuals" & plot_df$Time == comparison$Traits[i]] >= 
                                        min(plot_df$Value[plot_df$Origin == "Environment" & plot_df$Time == comparison$Envir[i]]))
                                   
                                   OV_perc <- c(OV_perc, 
                                                sum(CheckOverlap == 2)/length(CheckOverlap)        
                                   )
                                 }
                                 
                                 data.frame(
                                   AC = rep(AC, nrow(comparison)),
                                   VA = rep(VA, nrow(comparison)),
                                   SL = rep(SL, nrow(comparison)),
                                   Pert = rep(z, nrow(comparison)),
                                   rep = rep(b, nrow(comparison)),
                                   OV_Dens = OVERLAP,
                                   OV_Perc = OV_perc,
                                   Envir = comparison$Envir,
                                   Indivs = comparison$Traits,
                                   ID = EvoResSuc_iter$ID[which(EvoResSuc_iter$rep == b)]
                                 )
                               })
                               do.call(rbind, rep_ls)
                             })
                             do.call(rbind, Pert_ls)
                           })
                           do.call(rbind, Iter_ls)
                         })
  Overlap_df <- do.call(rbind, Overlap_ls)
  write.csv(Overlap_df, file = file.path(Dir.Exports, "Distrib_NonSpatial_df.csv"))
}




## Spatially Explicit Overlaps --------------------------------------------
# start_n = 1
# direction_n = -1
order <- names(EnvirDist_ls)#[start_n:length(EnvirDist_ls)]
# if(direction_n == -1){order <- rev(order)}
Overlap_ls <- lapply(X = order,
                     # cl = cl,
                     FUN = function(x){
                       # x <- order[2]
                       # Environment Data
                       AC <- as.numeric(unlist(strsplit(x = x, split = "_"))[1])
                       VA <- as.numeric(unlist(strsplit(x = x, split = "_"))[2])
                       SL <- as.numeric(unlist(strsplit(x = x, split = "_"))[3])
                       
                       # Establishing environmental raster
                       dat_rast <- raster(file.path(Dir.Data, paste0("ac", AC, "var", VA, ".nc")))
                       sl_vec <- (1:100)*SL
                       sl_mat <- matrix(rep(sl_vec, 100), nrow = 100, ncol = 100, byrow = TRUE)
                       sl_rast <- raster(sl_mat, xmn=-0.5, xmx=99.5, ymn=-0.5, ymx=99.5)
                       Enviro_ras <- dat_rast+sl_rast
                       
                       # Individual Data
                       Pop_fs <- list.files(Dir.Data, pattern = ".rds")
                       Iter_fs <- Pop_fs[grep(paste0("AC", AC), Pop_fs)]
                       Iter_fs <- Iter_fs[grep(paste0("VA", VA), Iter_fs)]
                       Iter_fs <- Iter_fs[grep(paste0("SL", ifelse(SL == 0.8, ".8", SL)), Iter_fs)]
                       
                       Iter_ls <- lapply(Iter_fs, FUN = function(y){
                         # y <- Iter_fs[1]
                         print(paste("File =", y))
                         Pop_df <- readRDS(file.path(Dir.Data, y))
                         Characteristics <- as.numeric(unlist(regmatches(unlist(strsplit(x = y, split = "_")),
                                                                         gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                                                  unlist(strsplit(x = y, split = "_"))))
                         ))
                         DI <- Characteristics[2]
                         MU <- Characteristics[3]
                         
                         #### Perturbation Loop ----
                         message("## Perturbation-Loop")
                         Pert_ls <- lapply(unique(Pop_df$pert.name),
                                           FUN = function(z){ #[1:2]
                                             # z <- unique(Pop_df$pert.name)[1] # perturbation magnitude loop here
                                             print(paste("Perturbation:", z))
                                             Pop_df <- Pop_df[Pop_df$pert.name == z,] 
                                             EvoResSuc_iter <- EvoResSuc_df[EvoResSuc_df$pert.name == z & 
                                                                              EvoResSuc_df$AC == AC & 
                                                                              EvoResSuc_df$SL == SL & 
                                                                              EvoResSuc_df$VA == VA & 
                                                                              EvoResSuc_df$DI == DI &
                                                                              EvoResSuc_df$MU == MU,]
                                             #### Repetition Loop ----
                                             # if(file.exists(file.path(Dir.Exports, paste0("TEMP_", y, "_",  z, ".RData")))){
                                             #   load(file.path(Dir.Exports, paste0("TEMP_", y, "_",  z, ".RData")))
                                             # }else{
                                             
                                             rep_ls <- lapply(unique(EvoResSuc_iter$rep),
                                                              # cl = cl,
                                                              FUN = function(b){
                                                                # b <- EvoResSuc_iter$rep[1]
                                                                print(paste("Perturbation:", z, "| Rep =", b))
                                                                
                                                                Popb_df <- Pop_df[Pop_df$rep == b,]
                                                                
                                                                t_min <- EvoResSuc_iter$t_minpost[which(EvoResSuc_iter$rep == b)]
                                                                t_max <- EvoResSuc_iter$t_maxpost[which(EvoResSuc_iter$rep == b)]
                                                                comparison <- data.frame(Envir = c("Pre", "Post", "Post"),
                                                                                         Traits = c("Pre", "PostMin", "PostMax"))
                                                                Enviro_ls <- list(Pre = Enviro_ras,
                                                                                  Post = Enviro_ras+z)
                                                                
                                                                #### Comparison Loop ----
                                                                if(file.exists(file.path(Dir.Exports, paste0("TEMP_", y, "_",  z, "_", b, "_NoEMD.RData")))){
                                                                  load(file.path(Dir.Exports, paste0("TEMP_", y, "_",  z, "_", b, "_NoEMD.RData")))
                                                                }else{
                                                                  nC <- nrow(comparison)
                                                                  cl <- makeCluster(nC)
                                                                  clusterExport(cl = cl, varlist = c("Dir.Data", "install.load.package", "package_vec",
                                                                                                     "EnvirDist_ls", "EvoResSuc_df", "z", "Pop_df", "EvoResSuc_iter", "Enviro_ras",
                                                                                                     "DI", "MU", "AC", "SL", "VA", "y", "x", "EvoResSuc_df", "%nin%",
                                                                                                     "t_min", "t_max", "comparison", "Enviro_ls", "Popb_df", "b"),
                                                                                envir=environment())
                                                                  clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))
                                                                  comp_ls <- pblapply(1:nrow(comparison), 
                                                                                      cl = cl,
                                                                                      FUN = function(comp_iter){
                                                                                        # comp_iter <- 2
                                                                                        t_iter <- c(460, t_min, t_max)[comp_iter]
                                                                                        comparison_iter <- comparison[comp_iter,]
                                                                                        print(paste("Comparison:", paste(comparison_iter, collapse = " vs. ")))
                                                                                        
                                                                                        ## select correct environment layer
                                                                                        Enviro_iter <- Enviro_ls[[which(names(Enviro_ls) == comparison_iter$Envir)]]
                                                                                        
                                                                                        ## make spatialpoints object and raster of individuals
                                                                                        Indivs_sp <- Pop_df[Pop_df$t == t_iter, ] ## right before perturbation
                                                                                        if(nrow(Indivs_sp) == 0){
                                                                                          Indivs_stack <- NULL
                                                                                          Ras_abs <- NULL
                                                                                          Best_df <- NULL
                                                                                          Closest_df <- NULL
                                                                                          # EMD2D <- NULL
                                                                                        }else{
                                                                                          coordinates(Indivs_sp) <- ~x+y
                                                                                          Indivs_ras <- rasterize(x = Indivs_sp[,"u"], y = Enviro_ras, fun = mean)$u
                                                                                          IndivsSD_ras <- rasterize(x = Indivs_sp[,"u"], y = Enviro_ras, fun = sd)$u
                                                                                          
                                                                                          # ## Contrast spatially explicit distributions
                                                                                          # plot(Indivs_ras, colNA = "grey")
                                                                                          # plot(Enviro_ras)
                                                                                          # 
                                                                                          # plot(abs(Indivs_ras-Enviro_ras), colNA = "grey")
                                                                                          # hist(abs(values(Indivs_ras-Enviro_ras)))
                                                                                          
                                                                                          # Raster differences
                                                                                          Ras_abs <- abs(Indivs_ras-Enviro_iter)
                                                                                          # plot(Ras_abs, colNA = "grey")
                                                                                          
                                                                                          # # EMD
                                                                                          # set.seed(42)
                                                                                          # EMD2D <- emd2d(matrix(Indivs_ras, ncol = 100), matrix(Enviro_iter, ncol = 100), max.iter = 1e3)
                                                                                          
                                                                                          # Best-suited environment
                                                                                          Best_df <- apply(Pop_df, MARGIN = 1, FUN = function(p){
                                                                                            Indiv_iter <- p
                                                                                            ## Cell ID
                                                                                            Best_cell <- which.min(abs(values(Enviro_iter)-Indiv_iter["u"]))
                                                                                            ## Difference
                                                                                            Best_diff <- abs(values(Enviro_iter)[Best_cell]-Indiv_iter["u"])
                                                                                            ## Distance to
                                                                                            Best_dist <- pointDistance(p1 = Indiv_iter[c("x", "y")], # individual location, 
                                                                                                                       p2 = raster::xyFromCell(
                                                                                                                         Enviro_iter, 
                                                                                                                         Best_cell
                                                                                                                       ), # centroid of cell with best conditions
                                                                                                                       lonlat = FALSE, allpairs = FALSE)
                                                                                            ## reporting back
                                                                                            data.frame(
                                                                                              Best_cell = Best_cell,
                                                                                              Best_diff = Best_diff,
                                                                                              Best_dist = Best_dist,
                                                                                              ID = Indiv_iter["id"])
                                                                                          })
                                                                                          Best_df <- do.call(rbind, Best_df)
                                                                                          
                                                                                          # Closest Distance to Matching Environment
                                                                                          Closest_ls <- apply(Pop_df, MARGIN = 1, FUN = function(p){
                                                                                            Indiv_iter <- p
                                                                                            ## Cell IDs
                                                                                            # message("Matching environment is abs(trait - envir) < 3")
                                                                                            Closest_cells <- which(abs(values(Enviro_iter)-Indiv_iter["u"]) < 3)
                                                                                            ## Differences
                                                                                            Closest_diffs <- abs(values(Enviro_iter)[Closest_cells]-Indiv_iter["u"])
                                                                                            ## Distance to
                                                                                            Closest_dists <- pointDistance(p1 = Indiv_iter[c("x", "y")], # individual location, 
                                                                                                                           p2 = raster::xyFromCell(
                                                                                                                             Enviro_iter, 
                                                                                                                             Closest_cells
                                                                                                                           ), # centroid of cell with best conditions
                                                                                                                           lonlat = FALSE, allpairs = FALSE)
                                                                                            ## reporting back
                                                                                            list(
                                                                                              Closest_cells = Closest_cells,
                                                                                              Closest_diffs = Closest_diffs,
                                                                                              Closest_dists = Closest_dists,
                                                                                              ID = Indiv_iter["id"]
                                                                                            )
                                                                                          })
                                                                                          Closest_df <- data.frame(Closest_diffs = unlist(lapply(Closest_ls, "[[", "Closest_diffs")), 
                                                                                                                   Closest_dists = unlist(lapply(Closest_ls, "[[", "Closest_dists"))) 
                                                                                        }
                                                                                        # Saving summary metrics
                                                                                        list(
                                                                                          Summary_df = data.frame(
                                                                                            AC = AC,
                                                                                            VA = VA,
                                                                                            SL = SL,
                                                                                            Pert = z,
                                                                                            rep = b, 
                                                                                            Best_dist_med = ifelse(is.null(median(Best_df$Best_dist)), NA, median(Best_df$Best_dist)),
                                                                                            Best_dist_sd = ifelse(is.null(sd(Best_df$Best_dist)), NA, sd(Best_df$Best_dist)),
                                                                                            Best_diff_med = ifelse(is.null(median(Best_df$Best_diff)), NA, median(Best_df$Best_diff)),
                                                                                            Best_diff_sd = ifelse(is.null(sd(Best_df$Best_diff)), NA, sd(Best_df$Best_diff)),
                                                                                            Closest_dist_med = ifelse(is.null(median(Closest_df$Closest_dists)), NA, median(Closest_df$Closest_dists)),
                                                                                            Closest_dist_sd = ifelse(is.null(sd(Closest_df$Closest_dists)), NA, sd(Closest_df$Closest_dists)),
                                                                                            Closest_diff_med = ifelse(is.null(median(Closest_df$Closest_diffs)), NA, median(Closest_df$Closest_diffs)),
                                                                                            Closest_diff_sd = ifelse(is.null(sd(Closest_df$Closest_diffs)), NA, sd(Closest_df$Closest_diffs)),
                                                                                            # EMD2D = EMD2D,
                                                                                            Envir = comparison_iter$Envir,
                                                                                            Indivs = comparison_iter$Traits,
                                                                                            ID = EvoResSuc_iter$ID[which(EvoResSuc_iter$rep == b)]
                                                                                          )
                                                                                          ,
                                                                                          RasterDiff = Ras_abs,
                                                                                          Best_df = Best_df,
                                                                                          Closest_df = Closest_df
                                                                                        )
                                                                                      })
                                                                  names(comp_ls) <- paste(comparison$Envir, comparison$Traits, sep = "-")
                                                                  stopCluster(cl)
                                                                  Up_ls <- list(Summary_df = do.call(rbind, lapply(comp_ls, "[[", "Summary_df"))
                                                                                # ,
                                                                                # Comparisons_ls = lapply(comp_ls, "[", -1)
                                                                  )
                                                                  save(Up_ls, file = file.path(Dir.Exports, paste0("TEMP_", y, "_",  z, "_", b, "_NoEMD.RData")))
                                                                }
                                                                Up_ls
                                                              })
                                             names(rep_ls) <- unique(EvoResSuc_iter$rep)
                                             Up_ls <- list(Summary_df = do.call(rbind, lapply(rep_ls, "[[", "Summary_df"))
                                                           # ,
                                                           # Comparisons_ls = lapply(lapply(rep_ls, "[", -1), "[[", "Comparisons_ls")
                                             )
                                             # names(Up_ls) <- c("Summary_df", z)
                                             # stopCluster(cl)
                                             # save(Up_ls, file = file.path(Dir.Exports, paste0("TEMP_", y, "_",  z, ".RData")))
                                             # }
                                             Up_ls
                                           }) # Pert_ls loop
                         Up_ls <- 
                           # list(Summary_df = 
                           do.call(rbind, lapply(Pert_ls, "[[", "Summary_df"))
                         # ,
                         # Perturbations_ls = lapply(lapply(Pert_ls, "[", -1), "[[", 1)
                         # )
                         # names(Up_ls) <- c("Summary_df", z)
                         # save(Up_ls, file = file.path(Dir.Exports, paste0("TEMP_", Iter_fs, "_",  z, ".RData")))
                         Up_ls
                       }) # Iter_ls loop
                       Up_ls <- do.call(rbind, Iter_ls)
                       Up_ls
                     }) # Overlap_ls loop
SpatialOver_df <- do.call(rbind, Overlap_ls)
write.csv(SpatialOver_df, file = file.path(Dir.Exports, "Distrib_Spatial_df_NoEMD.csv"))