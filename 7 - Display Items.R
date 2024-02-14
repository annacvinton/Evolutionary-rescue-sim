#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Display Item Creation
#'  DEPENDENCIES:
#'  - All previous scripts have been run
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
### CRAN ----
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}
package_vec <- c(
  "ggplot2",
  "remotes",
  "stringr",
  "dplyr",
  "cowplot",
  "ggpubr"
)
sapply(package_vec, install.load.package)

### NON-CRAN ----
if("ggsankey" %in% rownames(installed.packages()) == FALSE){
  remotes::install_github("davidsjoberg/ggsankey")
}
library(ggsankey)

## Functionality ----------------------------------------------------------
"%nin%" <- Negate("%in%") # a function for negation of %in% function

# DISPLAY ITMES ===========================================================

## Sankey Plot of Simulation Run Classifications --------------------------
# data loading
load(file.path(Dir.Exports, "EVORES_Metrics.RData"))
EvoResSuc_df <- EVORES_Metrics
EvoResSuc_df$EvoRes[!EvoResSuc_df$survival] <- "Not Possible"

# data extraction
sankey_df <- EvoResSuc_df[, c("pert.name", "survival", "EvoRes"
                              # , "AC", "DI", "MU", "SL", "VA"
                              )]
sankey_df$Run <- "Total Simulations"
colnames(sankey_df) <- c("Perturbation Magnitude", "Survival", 
                         "Evolutionary Rescue", "Simulation"
                         # , 
                         # "Spatial Autocorrelation", "Dispersal", "Mutation", 
                         # "Spatial Slope", "Spatial Variation"
                         )

# data reformatting for plotting
sankey_df <- sankey_df[sankey_df$`Perturbation Magnitude` >= 9, ]
sankey_df$`Perturbation Magnitude` <- str_pad(sankey_df$`Perturbation Magnitude`, width = 2, side = "left", pad = "0")
sankey_df$Survival <- str_replace_all(sankey_df$Survival, c("TRUE" = "Survival", "FALSE" = "Extinction"))

sankey_df$`Evolutionary Rescue` <- str_replace_all(sankey_df$`Evolutionary Rescue`, c("TRUE" = "Evolutionary Rescue", "FALSE" = "No Evolutionary Rescue"))

# data orientation for plotting
sankey_df <- make_long(sankey_df, 
                       "Simulation", 
                       # "Spatial Autocorrelation", "Spatial Slope", "Spatial Variation", "Dispersal", "Mutation",
                       "Perturbation Magnitude", 
                       "Survival", "Evolutionary Rescue"
                       )
sankey_count <- sankey_df%>%
  dplyr::group_by(node)%>%
  tally()
sankey_df <- merge(sankey_df, sankey_count, by.x = 'node', by.y = 'node', all.x = TRUE)

sankey_df$node <- factor(sankey_df$node, levels = c("09", "10", "11", "12", "13", "14", "15", "16", "Evolutionary Rescue", "Survival", "Extinction", "No Evolutionary Rescue", "Not Possible", "Total Simulations"))

sankey_df$next_node <- factor(sankey_df$next_node, levels = c("09", "10", "11", "12", "13", "14", "15", "16", "Evolutionary Rescue", "Survival", "Extinction", "No Evolutionary Rescue", "Not Possible", "Total Simulations"))

sankey_df$node[which(sankey_df$node == "Not Possible")] <- NA

sankey_gg <- ggplot(sankey_df[!is.na(sankey_df$node), ], aes(x = x
               , next_x = next_x
               , node = node
               , next_node = next_node
               , fill = factor(node)
               
               , label = paste0(node,"; n=", n)
               )) + 
  geom_sankey(flow.alpha = 0.5,  color = "gray40", show.legend = TRUE) +
  geom_sankey_label(size = 4, color = "white", fill= "gray40") +
  theme_bw() + theme(legend.position = "none") +  
  theme(axis.title = element_blank()
                  , axis.text.y = element_blank()
                  , axis.ticks = element_blank()  
                  , panel.grid = element_blank()) + 
  scale_fill_viridis_d(option = "F", direction = 1, na.translate = FALSE)
sankey_gg

ggsave(sankey_gg, filename = file.path(Dir.Exports, "PLOT_Sankey.png"), width = 16*2, height = 9*2, units = "cm")

## Combinations of Simulation Settings ------------------------------------
# data loading
load(file.path(Dir.Exports, "EVORES_Metrics.RData"))
EvoResSuc_df <- EVORES_Metrics
EvoResSuc_df$EvoRes[!EvoResSuc_df$survival] <- "Not Possible"
EvoResSuc_df[EvoResSuc_df$pert.name >= 9, ]

# data extraction
CombsBase_df <- EvoResSuc_df[, c("pert.name", "survival", "EvoRes", 
                             "AC", "DI", "MU", "SL", "VA")]
colnames(CombsBase_df) <- c("Perturbation Magnitude", "Survival",
                        "Total Simulations",
                        "Spatial Autocorrelation", "Dispersal", "Mutation",
                        "Spatial Slope", "Spatial Variation")

# tabulating differences
Combs_df <- aggregate(`Total Simulations` ~ `Spatial Autocorrelation` + `Spatial Slope` + 
                        `Spatial Variation` + 
                        Dispersal + Mutation, 
                      data = CombsBase_df, FUN = length)

CombsBase_df$`Total Simulations` <- as.logical(CombsBase_df$`Total Simulations`)
Combs_df$`Evolutionary Rescue` <- aggregate(`Total Simulations` ~ `Spatial Autocorrelation` + `Spatial Slope` + 
                        `Spatial Variation` + 
                        Dispersal + Mutation, 
                      data = CombsBase_df, FUN = sum)$`Total Simulations`

CombsBase_df$Survival <- as.logical(CombsBase_df$Survival)
Combs_df$Survival <- aggregate(Survival ~ `Spatial Autocorrelation` + `Spatial Slope` + 
                        `Spatial Variation` + 
                        Dispersal + Mutation, 
                      data = CombsBase_df, FUN = sum)$Survival

CombPlots_ls <- lapply(c("Total Simulations", "Survival", "Evolutionary Rescue"), function(Outcome_i){
  Combs_ls <- lapply(unique(Combs_df$`Spatial Variation`), function(VA_i){
    ggplot(data = Combs_df[Combs_df$`Spatial Variation` == VA_i, ],
           aes(x = `Spatial Autocorrelation`, 
               y = `Spatial Slope`, 
               fill = Outcome)) + 
      geom_tile(colour = "black") + 
      geom_label(aes(label = Outcome), fill = "white") + 
      theme_bw() + 
      facet_grid(Dispersal ~ Mutation, labeller = label_both) + 
      scale_fill_viridis_c(option = "F", direction = 1, 
                           limits = c(min(Combs_df$Outcome), 
                                      max(Combs_df$Outcome)), 
                           name = Outcome_i
                           ) + 
      labs(title = paste("Spatial Variation =", VA_i)) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(legend.position = "bottom", legend.key.width = unit(3, "cm"))
  })
  leg <- get_legend(Combs_ls[[1]])
  
  Combs_gg <- plot_grid(
    plot_grid(
      plotlist = lapply(Combs_ls, function(Plot_i){Plot_i + theme(legend.position = "none")}), ncol = 3),
    as_ggplot(leg), ncol = 1, rel_heights = c(1, 0.2))
  Combs_gg
  
})

ggsave(Combs_gg, filename = file.path(Dir.Exports, "PLOT_SettingCombinations.png"), width = 16*3, height = 9*2, units = "cm")


## Simulation Abundance Time-Series; Conceptual ---------------------------
# Data Loading
load(file.path(Dir.Exports, "SummaryTimeStep.RData"))
Abund_df <- Data_df
rm(Data_df)

### Overview of Runs for Specific Perturbation ----
# Assigning traceable IDs to abundance data
AbundIdents <- with(Abund_df, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_"))
EvoResIdents <- with(EvoResSuc_df, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_"))
Abund_df$ID <- EvoResSuc_df$ID[match(AbundIdents, EvoResIdents)]

# match outcomes to to IDs
Abund_df$survival <- EvoResSuc_df$survival[match(Abund_df$ID, EvoResSuc_df$ID)]
Abund_df$EvoRes <- EvoResSuc_df$EvoRes[match(Abund_df$ID, EvoResSuc_df$ID)]
Abund_df$Outcome <- str_replace_all(Abund_df$EvoRes, c("TRUE" = "Evolutionary Rescue", 
                                                       "FALSE" = "No Evolutionary Rescue",
                                                       "Not Possible" = "Extinction"))
Abund_df <- Abund_df[!is.na(Abund_df$Outcome), ]

# Subsetting for perturbation 10 due to good split of outcomes
Abundplot_df <- Abund_df[Abund_df$pert.name == 10, ]
Abundplot_df <- Abundplot_df[Abundplot_df$t <= 1000, ]
Abundplot_df <- Abundplot_df[Abundplot_df$t >= 310, ]
# message("ATTENTION ATTENTION - RE-INVESTIGATE HOW EXTINCTION IS CLASSIFIED/ASSIGNED. MAY HAVE TO CLEAN MORE FOR CRAPPED OUT RUNS")
# Abundplot_df$n[Abundplot_df$Outcome == "Extinction" & Abundplot_df$t > 470 & Abundplot_df$n > 100] <- 0

# making means and standard deviations for plotting
Abundplot_df <- data.frame(aggregate(x = n ~ t+Outcome, data = Abundplot_df, FUN = mean),
                           sd = aggregate(x = n ~ t+Outcome, data = Abundplot_df, FUN = sd)[,3])
Abundplot_df$sd[is.na(Abundplot_df$sd)] <- 0

# plotting
ConceptTime_gg <- ggplot(Abundplot_df, aes(x = t, y = n, fill = Outcome, color = Outcome)) + 
  geom_point(alpha = 0.1) +
  geom_line() + 
  geom_ribbon(aes(y = n, ymin = ifelse((n - sd) < 0, 0, n - sd), ymax = n + sd, fill = Outcome), alpha = .2) + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  ylim(c(0, NA)) +
  theme_bw()
ConceptTime_gg

### Overview of all Runs with Cutoffs and resulting classifications ----
message("can I somehow add the dpcut and rebcut?")
# DipCut <- 0.1 # go down to below 10% of pre-perturbation abundance
# RebCut <- 0.5 # bounce back to at least 50% of pre-perturbation abundance
# 
# print(paste(sum(is.na(EVORES_Metrics$n_minpost)), "of the", nrow(EVORES_Metrics), "potential evolutionary resuce runs had to be discarded due to data issues (mostly, runs missing writing events)"))
# 
# 
# 
# PotEvoRes1_gg <- ggplot(data = EVORES_Metrics, aes(x = perc_minpost)) + 
#   geom_histogram(bins = 1e3) + 
#   geom_vline(aes(xintercept = DipCut), col = "darkred") + 
#   annotate("text", x = DipCut+0.4, y = 200, label = paste0("Everything to the left of this cutoff (", DipCut, ") is considered as \n a population crash which may be subject to evolutionary rescue")) + 
#   theme_bw() + 
#   labs(title = "Abundance at minimum population size post-perturbation", 
#        x = "% of pre-perturbation abundance", y = "Count")
# 
# maxtest_df <- EVORES_Metrics[which(EVORES_Metrics$perc_minpost < DipCut), ]
# PotEvoRes2_gg <- ggplot(data = maxtest_df, aes(x = perc_maxpostpre)) + 
#   geom_histogram(bins = 1e3) + 
#   geom_vline(aes(xintercept = RebCut), col = "forestgreen") + 
#   annotate("text", x = RebCut+0.7, y = 75, label = paste0("Everything to the right of this cutoff (", RebCut, ") is considered as \n a population rebound indicative of evolutionary rescue")) + 
#   theme_bw() + 
#   labs(title = "Abundance at maximum population size following the post-perturbation minimum", 
#        x = "% of pre-perturbation abundance", y = "Count")
# 
# # ## run number summaries
# # n_totalruns <- nrow(runtimes)
# # n_extinctruns <- sum(!runtimes$survival)
# # n_potevoresruns <- length(na.omit(EVORES_Metrics$perc_minpost))
# # n_dipruns <- nrow(maxtest_df)
# # n_rebruns <- sum(na.omit(maxtest_df$perc_maxpostmin) > RebCut)
# # 
# # n_plot <- ggplot(data.frame(x = 1:10), aes(x = x, y = x/2)) + 
# #   geom_point(col = "white") + 
# #   annotate("text", x = 0, y = 5, label = paste("Total number of executed simulations =", n_totalruns), hjust = 0) + 
# #   annotate("text", x = 0, y = 4, label = paste("Number of executed simulations ending in extinction =", n_extinctruns), hjust = 0) + 
# #   annotate("text", x = 0, y = 3, label = paste("Number of executed simulations for which evolutionary resuce metrics can be computed =", n_potevoresruns), hjust = 0) + 
# #   annotate("text", x = 0, y = 2, label = paste("Number of executed simulations whose populations crash hard enough =", n_dipruns), hjust = 0) + 
# #   annotate("text", x = 0, y = 1, label = paste("Number of executed simulations who crash hard enough and rebound sufficiently =", n_rebruns), hjust = 0) + 
# #   theme_void()
# # 
# # plot_save <- plot_grid(PotEvoRes1_gg, PotEvoRes2_gg)
# # ggsave(plot_save, file = file.path(Dir.Exports, "EvoResCutOffs.jpg"), width = 32, height = 9, units = "cm")
# # ggsave(n_plot, file = file.path(Dir.Exports, "EvoResCutOffsRuns.jpg"), width = 22, height = 9, units = "cm")
# 
# 
# ## INDIVIDUAL run plotting according to evolutionary rescue or not
# Plot_df <- merge(Data_df, save_df[,1:8])
# Plot_df <- Plot_df[Plot_df$ID %in% save_df$ID[save_df$EvoRes], ]
# counts <- data.frame(pert.name = names(table(save_df$pert.name[save_df$EvoRes])),
#                      label = paste("n =", table(save_df$pert.name[save_df$EvoRes]))
# )
# counts$pert.name <- factor(counts$pert.name, levels = counts$pert.name)
# 
# 
# ribbons_df <- aggregate(data = Plot_df, n ~ pert.name + t, FUN = mean)
# ribbons_df$SD <- aggregate(data = Plot_df, n ~ pert.name + t, FUN = sd)$n
# 
# TS_gg <- ggplot(ribbons_df, aes(x = t, y = n)) + 
#   geom_line() + 
#   geom_ribbon(aes(y = n, ymin = n - SD, ymax = n + SD), alpha = .2) +
#   # geom_smooth(col = "black") + 
#   theme_bw() + 
#   labs(x = "Time [t]", y = "Abundance [n]") + 
#   facet_wrap(~factor(pert.name, levels = counts$pert.name), ncol = 5) + 
#   expand_limits(y = 0) + 
#   geom_text(
#     data = counts,
#     mapping = aes(x = max(ribbons_df$t), y = max(ribbons_df$n)+200, label = label),
#     hjust   = 1,
#     vjust   = 1
#   )
# TS_gg
# ggsave(TS_gg, file = file.path(Dir.Exports, "EvoResRuns_TS.jpg"), width = 22, height = 9, units = "cm")

## Environmental Parameterisation and Actual Environmental Values ---------