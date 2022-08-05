#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Metric summaries
#'  DEPENDENCIES:
#'  - MetricExtrac.R must have been run
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
  "e1071" # for skewness and kurtosis calculation
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
load(file.path(Dir.Exports, "NN_Metrics.RData"))
# View(head(Data_df))
Data_df <- Data_df[Data_df$t >= 310, ]

Combos_df <- with(Data_df, expand.grid(unique(SL), unique(VA), unique(MU)))
colnames(Combos_df) <- c("SL", "VA", "MU")

PlotCols <- c("n", colnames(Data_df)[11:22])
for(PlotVar in PlotCols){
  message(PlotVar)
  Plot_df <- Data_df[, -match(PlotCols[which(PlotCols %nin% PlotVar)], colnames(Data_df))]
  colnames(Plot_df)[which(colnames(Plot_df) == PlotVar)] <- "PlotVar"
  # Plot_df$Facettes <- paste("Dispersal =", Plot_df$DI,
  #                           " | ",  
  #                           "Perturbation =", Plot_df$pert.name)
  pb <- txtProgressBar(max = nrow(Combos_df), style = 3)
  for(i in 1:nrow(Combos_df)){
    
    Iter_df <- Plot_df[
      Plot_df$MU == Combos_df[i, "MU"] &
        Plot_df$VA == Combos_df[i, "VA"] &
        Plot_df$SL == Combos_df[i, "SL"], ]
    ggplot(
      data = rbind(Iter_df, Iter_df), 
      aes(x = t, y = PlotVar, col = factor(AC), 
          fill = factor(AC))) + 
      stat_summary(fun.data = mean_se, na.rm = TRUE,
                   geom = "ribbon", alpha = 0.6) +
      stat_summary(fun = mean, geom = "line", na.rm = TRUE,
                   size = 1.2, alpha = 0.9) +
      theme_bw() + labs(fill = "Autocorrelation",
                        col = "Autocorrelation",
                        title = paste(
                          "Slope =", Combos_df[i, "SL"], " | ", 
                          "Variance =", Combos_df[i, "VA"], " | ", 
                          "Evolution =", Combos_df[i, "MU"] 
                        )) + 
      scale_color_viridis_d(direction = -1) + scale_fill_viridis_d(direction = -1) + 
      facet_wrap(~ pert.name + factor(DI), ncol = 4) + ylab(PlotVar)
    ggsave(filename = file.path(Dir.Exports, 
                                paste0(PlotVar,
                                       "_",
                                       paste(Combos_df[i,], collapse="-"), 
                                       ".png")),
           width = 32, height = 32, unit = "cm"
    )
    setTxtProgressBar(pb, i)
  }
}


      




