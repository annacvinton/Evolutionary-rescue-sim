#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Evolutionary Success Reading
#'  DEPENDENCIES:
#'  - 2 - SummaryMetricsExtract.R must have been run
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
  "cowplot"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

# DATA ====================================================================
load(file.path(Dir.Exports, "NN_Metrics.RData"))
conditions <- Data_df[,c("pert.name", "rep", "AC", "DI", "SL", "VA", "MU")]
conditions <- conditions[which(!duplicated(conditions)), ]

runtimes <- aggregate(Data_df, t ~ pert.name + rep + AC + DI + MU + SL + VA, FUN = max)
runtimes$survival <- TRUE
runtimes$survival[runtimes$t < 1000] <- FALSE
runtimes$ID <- 1:nrow(runtimes)
noext <- runtimes[runtimes$survival, ]

cl <- makeCluster(detectCores())
clusterExport(cl, c("conditions", "Data_df", "%nin%", "noext"))

evoressuc_ls <- pblapply(1:nrow(noext), 
                         cl = cl,
                         FUN = function(l){
                           # message(l)
                           iter2_df <- Data_df[
                             Data_df$pert.name == noext[l, "pert.name"] &
                               Data_df$rep == noext[l, "rep"] &
                               Data_df$AC == noext[l, "AC"] &
                               Data_df$DI == noext[l, "DI"] &
                               Data_df$MU == noext[l, "MU"] &
                               Data_df$SL == noext[l, "SL"] &
                               Data_df$VA == noext[l, "VA"]
                             ,   
                           ]
                           iter_df <- iter2_df[iter2_df$t>460, ]
                           if(nrow(iter_df) < 65 | 460 %nin% iter2_df$t){
                             bind_df <- data.frame(
                               ID = noext$ID[l], 
                               n_minpost = NA,
                               t_minpost = NA,
                               n_maxpost = NA,
                               t_maxpost = NA,
                               perc_minpost = NA,
                               perc_maxpostmin = NA,
                               perc_maxpostpre = NA
                             )
                           }else{
                             pre_n <- Data_df[
                               Data_df$pert.name == noext[l, "pert.name"] &
                                 Data_df$rep == noext[l, "rep"] &
                                 Data_df$AC == noext[l, "AC"] &
                                 Data_df$DI == noext[l, "DI"] &
                                 Data_df$MU == noext[l, "MU"] &
                                 Data_df$SL == noext[l, "SL"] &
                                 Data_df$VA == noext[l, "VA"] &
                                 Data_df$t == 460
                               , 
                               "n"
                             ]
                             
                             minpost <- which.min(iter_df$n)
                             minpost_df <- iter_df[minpost, c("n", "t")]
                             minpost_df$ChangePre <- minpost_df$n/pre_n
                             
                             if(minpost != nrow(iter_df)){
                               maxpost <- which.max(iter_df$n[(minpost+1):nrow(iter_df)])+minpost
                               maxpost_df <- iter_df[maxpost, c("n", "t")]
                               maxpost_df$ChangePre <- maxpost_df$n/pre_n
                               maxpost_df$ChangePost <- maxpost_df$n/minpost_df$n
                             }else{
                               maxpost_df <- data.frame(
                                 n = NA,
                                 t = NA,
                                 ChangePre = NA,
                                 ChangePost = NA
                               )
                             }
                             
                             bind_df <- data.frame(
                               ID = noext$ID[l], 
                               n_minpost = minpost_df$n,
                               t_minpost = minpost_df$t,
                               n_maxpost = maxpost_df$n,
                               t_maxpost = maxpost_df$t,
                               perc_minpost = minpost_df$ChangePre,
                               perc_maxpostmin = maxpost_df$ChangePost,
                               perc_maxpostpre = maxpost_df$ChangePre
                             )
                           }
                           bind_df
                         })

evoressuc_df <- do.call(rbind, evoressuc_ls)
# save(evoressuc_df, file = "evoressuc_df.RData")

print(paste(sum(is.na(evoressuc_df$n_minpost)), "of the", nrow(evoressuc_df), "potential evolutionary resuce runs had to be discarded due to data issues (mostly, runs missing writing events)"))

DipCut <- 0.1 # go down to below 10% of pre-perturbation abundance
RebCut <- 0.5 # bounce back to at least 50% of pre-perturbation abundance

PotEvoRes1_gg <- ggplot(data = evoressuc_df, aes(x = perc_minpost)) + 
  geom_histogram(bins = 1e3) + 
  geom_vline(aes(xintercept = DipCut), col = "darkred") + 
  annotate("text", x = DipCut+0.4, y = 200, label = paste0("Everything to the left of this cutoff (", DipCut, ") is considered as \n a population crash which may be subject to evolutionary rescue")) + 
  theme_bw() + 
  labs(title = "Abundance at minimum population size post-perturbation", 
       x = "% of pre-perturbation abundance", y = "Count")

maxtest_df <- evoressuc_df[which(evoressuc_df$perc_minpost < DipCut), ]
PotEvoRes2_gg <- ggplot(data = maxtest_df, aes(x = perc_maxpostpre)) + 
  geom_histogram(bins = 1e3) + 
  geom_vline(aes(xintercept = RebCut), col = "forestgreen") + 
  annotate("text", x = RebCut+0.7, y = 75, label = paste0("Everything to the right of this cutoff (", RebCut, ") is considered as \n a population rebound indicative of evolutionary rescue")) + 
  theme_bw() + 
  labs(title = "Abundance at maximum population size following the post-perturbation minimum", 
       x = "% of pre-perturbation abundance", y = "Count")

## run number summaries
n_totalruns <- nrow(runtimes)
n_extinctruns <- sum(!runtimes$survival)
n_potevoresruns <- length(na.omit(evoressuc_df$perc_minpost))
n_dipruns <- nrow(maxtest_df)
n_rebruns <- sum(na.omit(maxtest_df$perc_maxpostmin) > RebCut)

n_plot <- ggplot(data.frame(x = 1:10), aes(x = x, y = x/2)) + 
  geom_point(col = "white") + 
  annotate("text", x = 0, y = 5, label = paste("Total number of executed simulations =", n_totalruns), hjust = 0) + 
  annotate("text", x = 0, y = 4, label = paste("Number of executed simulations ending in extinction =", n_extinctruns), hjust = 0) + 
  annotate("text", x = 0, y = 3, label = paste("Number of executed simulations for which evolutionary resuce metrics can be computed =", n_potevoresruns), hjust = 0) + 
  annotate("text", x = 0, y = 2, label = paste("Number of executed simulations whose populations crash hard enough =", n_dipruns), hjust = 0) + 
  annotate("text", x = 0, y = 1, label = paste("Number of executed simulations who crash hard enough and rebound sufficiently =", n_rebruns), hjust = 0) + 
  theme_void()

plot_save <- plot_grid(PotEvoRes1_gg, PotEvoRes2_gg)
ggsave(plot_save, file = file.path(Dir.Exports, "EvoResCutOffs.jpg"), width = 32, height = 9, units = "cm")
ggsave(n_plot, file = file.path(Dir.Exports, "EvoResCutOffsRuns.jpg"), width = 22, height = 9, units = "cm")


## Making one big data frame for modelling
### appending columns for evores success metrics
runtimes$EvoRes <- FALSE
save_df <- base::merge(x = runtimes, y = evoressuc_df, by = "ID")
save_df$EvoRes[which(save_df$perc_minpost < DipCut & save_df$perc_maxpostpre > RebCut)] <- TRUE
### removing runs which have abnormalities as identified in evores extraction
save_df <- save_df[save_df$ID %nin% evoressuc_df$ID[is.na(evoressuc_df$n_minpost)], ]
write.csv(save_df, file = file.path(Dir.Exports, "EvoResSuccessMetrics.csv"))


## INDIVIDUAL run plotting according to evolutionary rescue or not
Plot_df <- merge(Data_df, save_df[,1:8])
Plot_df <- Plot_df[Plot_df$ID %in% save_df$ID[save_df$EvoRes], ]
counts <- data.frame(pert.name = names(table(save_df$pert.name[save_df$EvoRes])),
                     label = paste("n =", table(save_df$pert.name[save_df$EvoRes]))
)
counts$pert.name <- factor(counts$pert.name, levels = counts$pert.name)

ribbons_df <- aggregate(data = Plot_df, n ~ pert.name + t, FUN = mean)
ribbons_df$SD <- aggregate(data = Plot_df, n ~ pert.name + t, FUN = sd)$n

TS_gg <- ggplot(ribbons_df, aes(x = t, y = n)) + 
  geom_line() + 
  geom_ribbon(aes(y = n, ymin = n - SD, ymax = n + SD), alpha = .2) +
  # geom_smooth(col = "black") + 
  theme_bw() + 
  labs(x = "Time [t]", y = "Abundance [n]") + 
  facet_wrap(~factor(pert.name, levels = counts$pert.name), ncol = 5) + 
  expand_limits(y = 0) + 
  geom_text(
    data = counts,
    mapping = aes(x = max(ribbons_df$t), y = max(ribbons_df$n)+200, label = label),
    hjust   = 1,
    vjust   = 1
  )
TS_gg
ggsave(TS_gg, file = file.path(Dir.Exports, "EvoResRuns_TS.jpg"), width = 22, height = 9, units = "cm")
