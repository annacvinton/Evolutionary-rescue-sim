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
  "overlapping",
  "lattice",
  "DescTools",
  "sp",
  "raster",
  "emdist"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
load(file.path(Dir.Exports, "EVORES_Metrics.RData")) # load EVORES_Metrics
EvoResSuc_df <- EVORES_Metrics[EVORES_Metrics$survival, ]
load(file.path(Dir.Exports, "ENVIRONMENT_Cells.RData")) # load Enviro_Cells_ls

## Non-Spatially Explicit Overlaps ----------------------------------------
if(file.exists(file.path(Dir.Exports, "DISTRIBUTIONS_NonSpatial.RData"))){
  load(file.path(Dir.Exports, "DISTRIBUTIONS_NonSpatial.RData")) # loads DISTRIBUTIONS_NonSpatial
}else{
  nC <- parallel::detectCores()
  # if(Sys.info()['sysname'] == "Linux"){cl = nC}
  cl <- makeCluster(nC)
  clusterExport(cl = cl, varlist = c("Dir.Data", "install.load.package", "package_vec"))
  clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))
  clusterExport(cl = cl, varlist = c("Enviro_Cells_ls", "EvoResSuc_df"))
  Overlap_ls <- pblapply(X = names(Enviro_Cells_ls),
                         # cl = cl,
                         FUN = function(x){
                           # x <- names(Enviro_Cells_ls)[1]

                           # Environment Data
                           AC <- as.numeric(unlist(strsplit(x = x, split = "_"))[1])
                           VA <- as.numeric(unlist(strsplit(x = x, split = "_"))[2])
                           SL <- as.numeric(unlist(strsplit(x = x, split = "_"))[3])

                           message(paste("AC =", AC, "; VA =", VA, "; SL =", SL))
                           
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
                                   Enviro_Cells_ls[[x]],
                                   Enviro_Cells_ls[[x]]+z,
                                   Popb_df$u[Popb_df$t == 460],
                                   Popb_df$u[Popb_df$t == t_min],
                                   Popb_df$u[Popb_df$t == t_max]
                                 ),
                                 Origin = c(
                                   rep("Environment", (length(Enviro_Cells_ls[[x]]))*2),
                                   rep("Individuals", length(Popb_df$u[Popb_df$t == 460])),
                                   rep("Individuals", length(Popb_df$u[Popb_df$t == t_min])),
                                   rep("Individuals", length(Popb_df$u[Popb_df$t == t_max]))
                                 ),
                                 Time = c(
                                   rep("Pre", length(Enviro_Cells_ls[[x]])),
                                   rep("Post", length(Enviro_Cells_ls[[x]])),
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
                                   OV_Dens = OVERLAP, # shared area between both density curves divided by their combined total area 
                                   OV_Perc = OV_perc, # percentage of trait values who fall within the range of the environment values
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
  DISTRIBUTIONS_NonSpatial <- do.call(rbind, Overlap_ls)
  save(DISTRIBUTIONS_NonSpatial, file = file.path(Dir.Exports, "DISTRIBUTIONS_NonSpatial.RData"))
  closeAllConnections()
}

## Spatially Explicit Overlaps --------------------------------------------
if(file.exists(file.path(Dir.Exports, "DISTRIBUTIONS_Spatial.RData"))){
  load(file.path(Dir.Exports, "DISTRIBUTIONS_Spatial.RData"))
}else{
  order <- names(Enviro_Cells_ls)#[start_n:length(Enviro_Cells_ls)]
  
  nC <- ifelse(length(order)>parallel::detectCores(), parallel::detectCores(), length(order))
  cl <- makeCluster(nC)
  clusterExport(cl = cl, varlist = c("Dir.Data", "Dir.Exports", "install.load.package", "package_vec", "%nin%",
                                     "EvoResSuc_df"),
                envir=environment())
  clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))
  
  Overlap_ls <- pblapply(X = order,
                         cl = cl,
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
                           SLCheck <- ifelse(SL == 0.8, ".8", SL)
                           SLCheck <- ifelse(SL == 1, "1_", SL)
                           Iter_fs <- Iter_fs[grep(paste0("SL", SLCheck), Iter_fs)]
                           
                           Iter_ls <- lapply(Iter_fs, FUN = function(y){
                             # y <- Iter_fs[1]
                             print(paste("File =", y, "|", which(Iter_fs == y), "/", length(Iter_fs)))
                             Pop_df <- readRDS(file.path(Dir.Data, y))
                             Characteristics <- as.numeric(unlist(regmatches(unlist(strsplit(x = y, split = "_")),
                                                                             gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                                                      unlist(strsplit(x = y, split = "_"))))
                             ))
                             DI <- Characteristics[2]
                             MU <- Characteristics[3]
                             
                             #### Perturbation Loop ----
                             # message("## Perturbation-Loop")
                             perts_vec <- unique(Pop_df$pert.name)[unique(Pop_df$pert.name)>= 9]
                             Pert_ls <- lapply(perts_vec,
                                               # cl = cl,
                                               FUN = function(z){ #[1:2]
                                                 # z <- unique(Pop_df$pert.name)[1] # perturbation magnitude loop here
                                                 # print(paste("Perturbation:", z))
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
                                                                    message(paste("Perturbation:", z, 
                                                                                  "| Rep =", which(unique(EvoResSuc_iter$rep) == b), 
                                                                                  "/", length(unique(EvoResSuc_iter$rep))))
                                                                    
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
                                                                      comp_ls <- lapply(1:nrow(comparison), 
                                                                                        # cl = cl,
                                                                                        FUN = function(comp_iter){
                                                                                          # comp_iter <- 1
                                                                                          t_iter <- c(460, t_min, t_max)[comp_iter]
                                                                                          comparison_iter <- comparison[comp_iter,]
                                                                                          # message("Comparison: ", paste(comparison_iter, collapse = " vs. "))
                                                                                          
                                                                                          ## select correct environment layer
                                                                                          Enviro_iter <- Enviro_ls[[which(names(Enviro_ls) == comparison_iter$Envir)]]
                                                                                          
                                                                                          ## make spatialpoints object and raster of individuals
                                                                                          Indivs_sp <- Popb_df[Popb_df$t == t_iter, ] ## right before perturbation
                                                                                          if(nrow(Indivs_sp) == 0){ # this should not be triggered anymore
                                                                                            Indivs_stack <- NULL
                                                                                            Ras_abs <- NULL
                                                                                            Best_df <- NULL
                                                                                            Closest_df <- NULL
                                                                                            Weighted_maladap <- NULL
                                                                                            # EMD2D <- NULL
                                                                                          }else{
                                                                                            coordinates(Indivs_sp) <- ~x+y
                                                                                            Indivs_ras <- rasterize(x = Indivs_sp[,"u"], y = Enviro_ras, fun = mean)$u
                                                                                            IndivsSD_ras <- rasterize(x = Indivs_sp[,"u"], y = Enviro_ras, fun = sd)$u
                                                                                            
                                                                                            # Raster differences
                                                                                            Ras_abs <- abs(Indivs_ras-Enviro_iter)
                                                                                            
                                                                                            # Best-suited environment
                                                                                            # print("Best cell distance and difference")
                                                                                            Best_df <- apply(data.frame(Indivs_sp), MARGIN = 1, FUN = function(p){
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
                                                                                            # print("Closest suitable distance and difference")
                                                                                            Closest_ls <- apply(data.frame(Indivs_sp), MARGIN = 1, FUN = function(p){
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
                                                                                            
                                                                                            # Weighted distance by adaptdeness
                                                                                            # print("Distance-weighted adaptedness")
                                                                                            Weighted_ls <- apply(data.frame(Indivs_sp), MARGIN = 1, FUN = function(p){
                                                                                              Indiv_iter <- p
                                                                                              ## Cell IDs
                                                                                              # message("Matching environment is abs(trait - envir) < 3")
                                                                                              ## Differences
                                                                                              Weighted_maladap <- abs(values(Enviro_iter)-as.numeric(Indiv_iter["u"]))
                                                                                              ## Distances
                                                                                              Weighted_dists <- pointDistance(p1 = Indiv_iter[c("x", "y")], # individual location, 
                                                                                                                              p2 = coordinates(Enviro_iter), # centroids of all cells
                                                                                                                              lonlat = FALSE, allpairs = FALSE)
                                                                                              
                                                                                              
                                                                                              ## Distance to
                                                                                              
                                                                                              ## reporting back
                                                                                              list(
                                                                                                Weighted_maladap = Weighted_maladap/Weighted_dists,
                                                                                                ID = Indiv_iter["id"]
                                                                                              )
                                                                                            })
                                                                                            Weighted_maladap <- as.numeric(unlist(lapply(Weighted_ls, "[[", "Weighted_maladap")))
                                                                                          }
                                                                                          # Saving summary metrics
                                                                                          list(
                                                                                            Summary_df = data.frame(
                                                                                              AC = AC,
                                                                                              VA = VA,
                                                                                              SL = SL,
                                                                                              Pert = z,
                                                                                              rep = b, 
                                                                                              ## best environmental cell
                                                                                              Best_dist_mean = ifelse(is.null(mean(Best_df$Best_dist)), NA, mean(Best_df$Best_dist)),
                                                                                              Best_dist_sd = ifelse(is.null(sd(Best_df$Best_dist)), NA, sd(Best_df$Best_dist)),
                                                                                              Best_diff_mean = ifelse(is.null(mean(Best_df$Best_diff)), NA, mean(Best_df$Best_diff)),
                                                                                              Best_diff_sd = ifelse(is.null(sd(Best_df$Best_diff)), NA, sd(Best_df$Best_diff)),
                                                                                              ## closest cell whose value is less than 3 units away from individual traits
                                                                                              Closest_dist_mean = ifelse(is.null(Closest_df), NA, mean(Closest_df$Closest_dists)),
                                                                                              Closest_dist_sd = ifelse(is.null(Closest_df), NA, sd(Closest_df$Closest_dists)),
                                                                                              Closest_diff_mean = ifelse(is.null(Closest_df), NA, mean(Closest_df$Closest_diffs)),
                                                                                              Closest_diff_sd = ifelse(is.null(Closest_df), NA, sd(Closest_df$Closest_diffs)),
                                                                                              ## Weighted maladaptedness in full environment (weighted by distance)
                                                                                              Weighted_maladap_mean = ifelse(is.null(Weighted_maladap), NA, mean(Weighted_maladap)),
                                                                                              Weighted_maladap_sd = ifelse(is.null(Weighted_maladap), NA, sd(Weighted_maladap)),
                                                                                              # EMD2D = EMD2D,
                                                                                              Envir = comparison_iter$Envir,
                                                                                              Indivs = comparison_iter$Traits,
                                                                                              ID = EvoResSuc_iter$ID[which(EvoResSuc_iter$rep == b)]
                                                                                            ),
                                                                                            RasterDiff = values(Ras_abs)
                                                                                          )
                                                                                        })
                                                                      names(comp_ls) <- paste(comparison$Envir, comparison$Traits, sep = "-")
                                                                      Up_ls <- list(Summary_df = do.call(rbind, lapply(comp_ls, "[[", "Summary_df")),
                                                                                    RasterDiff = lapply(comp_ls, "[", "RasterDiff")
                                                                      )
                                                                      save(Up_ls, file = file.path(Dir.Exports, paste0("TEMP_", y, "_",  z, "_", b, "_NoEMD.RData")))
                                                                    }
                                                                    Up_ls
                                                                  })
                                                 names(rep_ls) <- unique(EvoResSuc_iter$rep)
                                                 Up_ls <- list(Summary_df = do.call(rbind, lapply(rep_ls, "[[", "Summary_df")),
                                                               RasterDiff = lapply(rep_ls, "[[", "RasterDiff"))
                                                 Up_ls
                                               }) # Pert_ls loop
                             names(Pert_ls) <- perts_vec
                             Pert_ls <- Pert_ls[which(!unlist(lapply(lapply(lapply(Pert_ls, "[[", "Summary_df"), nrow), is.null)))] # select only perturbations where survival occured
                             Up_ls <- list(Summary_df = do.call(rbind, lapply(Pert_ls, "[[", "Summary_df")),
                                           RasterDiff = lapply(Pert_ls, "[[", "RasterDiff")
                             )
                             Up_ls
                           }) # Iter_ls loop
                           names(Iter_ls) <- Iter_fs
                           Up_ls <- list(Summary_df = do.call(rbind, lapply(Iter_ls, "[[", "Summary_df")),
                                         RasterDiff = lapply(Iter_ls, "[[", "RasterDiff")
                           )
                           Up_ls
                         }) # Overlap_ls loop
  DISTRIBUTIONS_Spatial <- list(Summary_df = do.call(rbind, lapply(Overlap_ls, "[[", "Summary_df")),
                                RasterDiff = lapply(Overlap_ls, "[[", "RasterDiff"))
  save(DISTRIBUTIONS_Spatial, file = file.path(Dir.Exports, "DISTRIBUTIONS_Spatial.RData"))
  unlink(list.files(Dir.Exports, pattern = "TEMP_", full.names = TRUE))
  stopCluster(cl)
}

## Earth-Mover Distance ---------------------------------------------------
order <- names(Enviro_Cells_ls)#[start_n:length(Enviro_Cells_ls)]
nC <- ifelse(length(order)>parallel::detectCores(), parallel::detectCores(), length(order))
cl <- makeCluster(nC)
clusterExport(cl = cl, varlist = c("Dir.Data", "Dir.Exports", "install.load.package", "package_vec", "%nin%",
                                   "EvoResSuc_df"),
              envir=environment())
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

Overlap_ls <- pblapply(X = order,
                       cl = cl,
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
                         SLCheck <- ifelse(SL == 0.8, ".8", SL)
                         SLCheck <- ifelse(SL == 1, "1_", SL)
                         Iter_fs <- Iter_fs[grep(paste0("SL", SLCheck), Iter_fs)]
                         
                         Iter_ls <- lapply(Iter_fs, FUN = function(y){
                           # y <- Iter_fs[1]
                           print(paste("File =", y, "|", which(Iter_fs == y), "/", length(Iter_fs)))
                           Pop_df <- readRDS(file.path(Dir.Data, y))
                           Characteristics <- as.numeric(unlist(regmatches(unlist(strsplit(x = y, split = "_")),
                                                                           gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                                                    unlist(strsplit(x = y, split = "_"))))
                           ))
                           DI <- Characteristics[2]
                           MU <- Characteristics[3]
                           
                           #### Perturbation Loop ----
                           # message("## Perturbation-Loop")
                           perts_vec <- unique(Pop_df$pert.name)[unique(Pop_df$pert.name)>= 9]
                           Pert_ls <- lapply(perts_vec,
                                             # cl = cl,
                                             FUN = function(z){ #[1:2]
                                               # z <- unique(Pop_df$pert.name)[1] # perturbation magnitude loop here
                                               Pop_df <- Pop_df[Pop_df$pert.name == z,] 
                                               EvoResSuc_iter <- EvoResSuc_df[EvoResSuc_df$pert.name == z & 
                                                                                EvoResSuc_df$AC == AC & 
                                                                                EvoResSuc_df$SL == SL & 
                                                                                EvoResSuc_df$VA == VA & 
                                                                                EvoResSuc_df$DI == DI &
                                                                                EvoResSuc_df$MU == MU,]
                                               #### Repetition Loop ----
                                               rep_ls <- lapply(unique(EvoResSuc_iter$rep),
                                                                # cl = cl,
                                                                FUN = function(b){
                                                                  # b <- EvoResSuc_iter$rep[1]
                                                                  message(paste("Perturbation:", z, 
                                                                                "| Rep =", which(unique(EvoResSuc_iter$rep) == b), 
                                                                                "/", length(unique(EvoResSuc_iter$rep))))
                                                                  
                                                                  Popb_df <- Pop_df[Pop_df$rep == b,]
                                                                  
                                                                  t_min <- EvoResSuc_iter$t_minpost[which(EvoResSuc_iter$rep == b)]
                                                                  t_max <- EvoResSuc_iter$t_maxpost[which(EvoResSuc_iter$rep == b)]
                                                                  comparison <- data.frame(Envir = c("Pre", "Post", "Post"),
                                                                                           Traits = c("Pre", "PostMin", "PostMax"))
                                                                  Enviro_ls <- list(Pre = Enviro_ras,
                                                                                    Post = Enviro_ras+z)
                                                                  
                                                                  #### Comparison Loop ----
                                                                  if(file.exists(file.path(Dir.Exports, paste0("TEMP_", y, "_",  z, "_", b, "_EMD.RData")))){
                                                                    load(file.path(Dir.Exports, paste0("TEMP_", y, "_",  z, "_", b, "_EMD.RData")))
                                                                  }else{
                                                                    comp_ls <- lapply(1:nrow(comparison), 
                                                                                      # cl = cl,
                                                                                      FUN = function(comp_iter){
                                                                                        # comp_iter <- 1
                                                                                        t_iter <- c(460, t_min, t_max)[comp_iter]
                                                                                        comparison_iter <- comparison[comp_iter,]
                                                                                        
                                                                                        ## select correct environment layer
                                                                                        Enviro_iter <- Enviro_ls[[which(names(Enviro_ls) == comparison_iter$Envir)]]
                                                                                        
                                                                                        ## make spatialpoints object and raster of individuals
                                                                                        Indivs_sp <- Popb_df[Popb_df$t == t_iter, ] ## right before perturbation
                                                                                        if(nrow(Indivs_sp) == 0){ # this should not be triggered anymore
                                                                                          EMD2D <- NULL
                                                                                        }else{
                                                                                          coordinates(Indivs_sp) <- ~x+y
                                                                                          Indivs_ras <- rasterize(x = Indivs_sp[,"u"], y = Enviro_ras, fun = mean)$u
                                                                                          # EMD
                                                                                          set.seed(42)
                                                                                          EMD2D <- emd2d(matrix(Indivs_ras, ncol = 100), matrix(Enviro_iter, ncol = 100), max.iter = 5e3) # may want to set max.iter = 1 for testing purposes
                                                                                        }
                                                                                        data.frame(AC = AC,
                                                                                                   VA = VA,
                                                                                                   SL = SL,
                                                                                                   Pert = z,
                                                                                                   rep = b, 
                                                                                                   EMD2D = EMD2D,
                                                                                                   Envir = comparison_iter$Envir,
                                                                                                   Indivs = comparison_iter$Traits,
                                                                                                   ID = EvoResSuc_iter$ID[which(EvoResSuc_iter$rep == b)])
                                                                                        
                                                                                      })
                                                                    Up_ls <- do.call(rbind, comp_ls)
                                                                    save(Up_ls, file = file.path(Dir.Exports, paste0("TEMP_", y, "_",  z, "_", b, "_EMD.RData")))
                                                                  }
                                                                  Up_ls
                                                                })
                                               do.call(rbind, rep_ls)
                                             }) # Pert_ls loop
                           do.call(rbind, Pert_ls)
                         }) # Iter_ls loop
                         do.call(rbind, Iter_ls)
                       }) # Overlap_ls loop
stop("Merge with DISTRIBUTION_Spatial by ID to add EMD2D to previous data frame DISTRIBUTIONS_Spatial$Summary_df")
# DISTRIBUTIONS_Spatial <- list(Summary_df = do.call(rbind, lapply(Overlap_ls, "[[", "Summary_df")),
#                               RasterDiff = lapply(Overlap_ls, "[[", "RasterDiff"))
# save(DISTRIBUTIONS_Spatial, file = file.path(Dir.Exports, "DISTRIBUTIONS_Spatial.RData"))
# unlink(list.files(Dir.Exports, pattern = "TEMP_", full.names = TRUE))
# stopCluster(cl)