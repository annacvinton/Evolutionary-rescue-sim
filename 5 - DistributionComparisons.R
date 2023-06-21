#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Comparison of landscape parameter distributions and individual trait distributions
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

## Cores ------------------------------------------------------------------
nC <- parallel::detectCores()
cl <- makeCluster(nC)
clusterExport(cl = cl, varlist = c("Dir.Data", "install.load.package", "package_vec"))
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

# DATA ====================================================================
load(file = file.path(Dir.Exports, "EnviroDist_ls.RData")) # load EnvirDist_ls

## Non-Spatially Explicit Overlaps ----------------------------------------
clusterExport(cl = cl, varlist = c("EnvirDist_ls"))
Overlap_ls <- pblapply(X = names(EnvirDist_ls),
         cl = cl,
         FUN = function(x){
           # x <- names(EnvirDist_ls)[1]
           # print(x)
           
           # Environment Data
           AC <- as.numeric(unlist(strsplit(x = x, split = "_"))[1])
           VA <- as.numeric(unlist(strsplit(x = x, split = "_"))[2])
           SL <- as.numeric(unlist(strsplit(x = x, split = "_"))[3])
           
           # Individual Data
           Pop_fs <- list.files(Dir.Data, pattern = ".rds")
           Iter_fs <- Pop_fs[grep(paste0("AC", AC), Pop_fs)]
           Iter_fs <- Iter_fs[grep(paste0("VA", VA), Iter_fs)]
           Iter_fs <- Iter_fs[grep(paste0("SL", ifelse(SL == 0.8, ".8", SL)), Iter_fs)]
    
           Iter_ls <- lapply(Iter_fs, FUN = function(y){
             # y <- Iter_fs[1] # simulation run loop here
             # print(y)
             Pop_df <- readRDS(file.path(Dir.Data, y))
             
             Pert_ls <- lapply(unique(Pop_df$pert.name), FUN = function(z){
               # z <- unique(Pop_df$pert.name)[1] # perturbation magnitude loop here
               Pop_df <- Pop_df[Pop_df$pert.name == z,] 
               plot_df <- data.frame(Value = c(
                 EnvirDist_ls[[x]]+z,
                 Pop_df$u[Pop_df$t == 460]
               ),
               Origin = c(
                 rep("Environment", length(EnvirDist_ls[[x]])),
                 rep("Individuals", length(Pop_df$u[Pop_df$t == 460]))
               )
               )

               # ggplot(data = plot_df, aes(x = Value, fill = Origin)) +
               #   # geom_histogram(alpha = 0.3) +
               #   geom_density(alpha = 0.4) +
               #   theme_bw() + labs(title = paste("Perturbation", z, ";",
               #                                   "AC", AC, ";",
               #                                   "VA", VA, ";",
               #                                   "SL", SL
               #                                   ),
               #                     )
               # 
               # ggplot(data = plot_df, aes(x = Value, fill = Origin, colour = Origin)) +
               #   geom_histogram(aes(y= after_stat(density)), alpha=0.4, 
               #                  position = "identity", lwd = 0.2) +
               #   ggtitle("Normalized")

               # overlap of density areas
               if(!(460 %in% unique(Pop_df$t))){
                 OVERLAP <- NA
                 OV_perc <- NA
               }else{
                 OVERLAP <- overlap(list(X2 = EnvirDist_ls[[x]], 
                                         X1 = Pop_df$u[Pop_df$t == 460]), 
                                    plot = FALSE)$OV
                 
                 # Dens_Enviro <- density(plot_df$Value[plot_df$Origin == "Environment"])
                 # Dens_Indivs <- density(plot_df$Value[plot_df$Origin == "Individuals"])
                 # convolve(Dens_Enviro, Dens_Indivs, type = "open")
                 # sum(Dens_Enviro$x * Dens_Indivs$x)
                 # 
                 # emd(A = matrix(c(Dens_Enviro$y, Dens_Enviro$x), nrow = length(Dens_Enviro$x)), 
                 #     B = matrix(c(Dens_Indivs$y, Dens_Indivs$x), nrow = length(Dens_Indivs$x)), 
                 #     dist="euclidean")
                 
                 # 2 if trait value is within optimal range given the environment, 1 if not
                 CheckOverlap <- (Pop_df$u[Pop_df$t == 460] <= max(EnvirDist_ls[[x]])) + 
                   (Pop_df$u[Pop_df$t == 460] >= min(EnvirDist_ls[[x]]))
                 OV_perc <- sum(CheckOverlap == 2)/length(CheckOverlap)
               }
               
               data.frame(
                 AC = AC,
                 VA = VA,
                 SL = SL,
                 Pert = z,
                 OV_Dens = OVERLAP,
                 OV_Perc = OV_perc
               )
             })
             do.call(rbind, Pert_ls)
           })
           do.call(rbind, Iter_ls)
         })
Overlap_df <- do.call(rbind, Overlap_ls)
save(Overlap_df, file = file.path(Dir.Exports, "Overlap_df.RData"))
write.csv(Overlap_df, file = file.path(Dir.Exports, "Overlap_df.csv"))

ggplot(data = Overlap_df, aes(x = factor(VA), y = OV_Dens)) + 
  geom_boxplot() + 
  theme_bw()

## Spatially Explicit Overlaps ----------------------------------------
x <- names(EnvirDist_ls)[[4]]

# Environment Data
AC <- as.numeric(unlist(strsplit(x = x, split = "_"))[1])
VA <- as.numeric(unlist(strsplit(x = x, split = "_"))[2])
SL <- as.numeric(unlist(strsplit(x = x, split = "_"))[3])

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

y <- Iter_fs[1]
# Iter_ls <- lapply(Iter_fs, FUN = function(y){
  # print(y)
  # y <- Iter_fs[1] # simulation run loop here
  Pop_df <- readRDS(file.path(Dir.Data, y))
  
  #### Perturbation Loop ----
  message("## Perturbation-Loop")
  Pert_ls <- lapply(unique(Pop_df$pert.name), FUN = function(z){ #[1:2]
    # z <- unique(Pop_df$pert.name)[1] # perturbation magnitude loop here
    print(paste("Perturbation:", z))
    Pop_df <- Pop_df[Pop_df$pert.name == z,] 
    #### Time Loop ----
    message("## Time-Loop")
    Time_ls <- lapply(c(460, 540), FUN = function(t_iter){
      # t_iter <- 460
      print(paste("Time:", t_iter))
      
      ## make spatialpoints object and raster of individuals
      Indivs_sp <- Pop_df[Pop_df$t == t_iter, ] ## right before perturbation
      if(nrow(Indivs_sp) == 0){
        Ras_abs <- NULL
        Best_df <- NULL
        Closest_df <- NULL
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
        Ras_abs <- abs(Indivs_ras-Enviro_ras)
        # plot(Ras_abs, colNA = "grey")
        
        # EMD
        set.seed(42)
        EMD2D <- emd2d(matrix(Indivs_ras, ncol = 100), matrix(Enviro_ras, ncol = 100), max.iter = 1e2)
        
        # Best-suited environment
        Best_df <- pbapply(Pop_df, MARGIN = 1, FUN = function(p){
          Indiv_iter <- p
          ## Cell ID
          Best_cell <- which.min(abs(values(Enviro_ras)-Indiv_iter["u"]))
          ## Difference
          Best_diff <- abs(values(Enviro_ras)[Best_cell]-Indiv_iter["u"])
          ## Distance to
          Best_dist <- pointDistance(p1 = Indiv_iter[c("x", "y")], # individual location, 
                                     p2 = raster::xyFromCell(
                                       Enviro_ras, 
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
        Closest_ls <- pbapply(Pop_df, MARGIN = 1, FUN = function(p){
          Indiv_iter <- p
          ## Cell IDs
          # message("Matching environment is abs(trait - envir) < 3")
          Closest_cells <- which(abs(values(Enviro_ras)-Indiv_iter["u"]) < 3)
          ## Differences
          Closest_diffs <- abs(values(Enviro_ras)[Closest_cells]-Indiv_iter["u"])
          ## Distance to
          Closest_dists <- pointDistance(p1 = Indiv_iter[c("x", "y")], # individual location, 
                                         p2 = raster::xyFromCell(
                                           Enviro_ras, 
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
          t_indivs = t_iter,
          # t_enviro = 9999,
          Best_dist_med = median(Best_df$Best_dist),
          Best_dist_sd = sd(Best_df$Best_dist),
          Best_diff_med = median(Best_df$Best_diff),
          Best_diff_sd = sd(Best_df$Best_diff),
          Closest_dist_med = median(Closest_df$Closest_dists),
          Closest_dist_sd = sd(Closest_df$Closest_dists),
          Closest_diff_med = median(Closest_df$Closest_diffs),
          Closest_diff_sd = sd(Closest_df$Closest_diffs)
        ),
        RasterDiff = Ras_abs,
        EMD2D = EMD2D, 
        Best_df = Best_df,
        Closest_df = Closest_df
      )
    })
    Time_ls
  })
stop("summarise Pert_ls")    
    
    
    
    
    
    
    
  # })
# })







