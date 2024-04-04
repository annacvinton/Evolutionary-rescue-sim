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
  "ggpubr",
  "pbapply",
  "zoo",
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

## Data Availability and Simulation Settings ------------------------------
### Sankey Plot of Simulation Run Classifications ----
# data loading
load(file.path(Dir.Exports, "EVORES_Metrics.RData"))
EvoResSuc_df <- EVORES_Metrics
EvoResSuc_df <- EvoResSuc_df[EvoResSuc_df$pert.name >= 9, ]
EvoResSuc_df$EvoRes[!EvoResSuc_df$SuffDip] <- "Insufficient Population Crash"

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

sankey_df$node <- factor(sankey_df$node, levels = c("09", "10", "11", "12", "13", "14", "15", "16", "Evolutionary Rescue", "Survival", "Extinction", "No Evolutionary Rescue", "Insufficient Population Crash", "Total Simulations"))

sankey_df$next_node <- factor(sankey_df$next_node, levels = c("09", "10", "11", "12", "13", "14", "15", "16", "Evolutionary Rescue", "Survival", "Extinction", "No Evolutionary Rescue", "Insufficient Population Crash", "Total Simulations"))

# sankey_df$node[which(sankey_df$node == "Not Possible")] <- NA

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

### Combinations of Simulation Settings ----
# data loading
load(file.path(Dir.Exports, "EVORES_Metrics.RData"))
EvoResSuc_df <- EVORES_Metrics
EvoResSuc_df <- EvoResSuc_df[EvoResSuc_df$pert.name >= 9, ]
EvoResSuc_df <- EvoResSuc_df[EvoResSuc_df$SuffDip | is.na(EvoResSuc_df$SuffDip), ]
EvoResSuc_df$EvoRes[!EvoResSuc_df$survival] <- FALSE
# EvoResSuc_df$EvoRes[!EvoResSuc_df$SuffDip] <- "insufficient crash"

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

CombsBase_df$Survival <- as.logical(CombsBase_df$Survival)
Combs_df$Survival <- aggregate(Survival ~ `Spatial Autocorrelation` + `Spatial Slope` + 
                        `Spatial Variation` + 
                        Dispersal + Mutation, 
                      data = CombsBase_df, FUN = sum)$Survival

CombsBase_df$`Total Simulations` <- as.logical(CombsBase_df$`Total Simulations`)
Combs_df$`Evolutionary Rescue` <- aggregate(`Total Simulations` ~ `Spatial Autocorrelation` + `Spatial Slope` + 
                                              `Spatial Variation` + 
                                              Dispersal + Mutation, 
                                            data = CombsBase_df, FUN = sum)$`Total Simulations`

Combs_df$`Evolutionary Rescue` <- round(Combs_df$`Evolutionary Rescue`/Combs_df$`Survival`, 3)*100
Combs_df$Survival <- round(Combs_df$Survival/Combs_df$`Total Simulations`, 3)*100


CombPlots_ls <- lapply(c("Total Simulations", "Survival", "Evolutionary Rescue"), function(Outcome_i){
  print(Outcome_i)
  colnames(Combs_df)[which(colnames(Combs_df) == Outcome_i)] <- "Outcome"
  Combs_ls <- lapply(unique(Combs_df$`Spatial Variation`), function(VA_i){
    ggplot(data = Combs_df[Combs_df$`Spatial Variation` == VA_i, ],
           aes(x = `Spatial Autocorrelation`, 
               y = `Spatial Slope`, 
               fill = Outcome)) + 
      geom_tile(colour = "black") + 
      geom_label(aes(label = Outcome), fill = "white", size = 2) + 
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
Combs_gg <- plot_grid(plotlist = CombPlots_ls, ncol = 1, labels = "AUTO")
Combs_gg

ggsave(
  Combs_gg, 
  filename = file.path(Dir.Exports, "PLOT_SettingCombinations.png"), 
  width = 16*2, height = 9*4, units = "cm")


## Time-Series Plotting ---------------------------------------------------
### Abundance; Conceptual ----
# Data Loading
load(file.path(Dir.Exports, "POPULATION_TimeStep.RData"))
Abund_df <- Data_df
rm(Data_df)
load(file.path(Dir.Exports, "EVORES_Metrics.RData"))
EvoResSuc_df <- EVORES_Metrics
rm(EVORES_Metrics)

EvoResSuc_df$EvoRes[!EvoResSuc_df$survival] <- "Not Possible"
EvoResSuc_df$EvoRes[!EvoResSuc_df$SuffDip] <- "Insufficient Crash"
EvoResSuc_df <- EvoResSuc_df[EvoResSuc_df$pert.name >= 9, ]

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

### Overview of Runs for Specific Perturbation ----
# Subsetting for perturbation 10 due to good split of outcomes
Abundplot_df <- Abund_df[Abund_df$pert.name == 10, ]
Abundplot_df <- Abundplot_df[Abundplot_df$t <= 1000, ]
Abundplot_df <- Abundplot_df[Abundplot_df$t >= 310, ]

# making means and standard deviations for plotting
Abundplot_df <- data.frame(aggregate(x = n ~ t+Outcome, data = Abundplot_df, FUN = mean),
                           sd = aggregate(x = n ~ t+Outcome, data = Abundplot_df, FUN = sd)[,3])
Abundplot_df$sd[is.na(Abundplot_df$sd)] <- 0

# plotting
ConceptTime_gg <- ggplot(Abundplot_df, aes(x = t, y = n, fill = Outcome, color = Outcome)) + 
  geom_point(alpha = 0.1) +
  geom_line() + 
  geom_ribbon(aes(y = n, ymin = ifelse((n - sd) < 0, 0, n - sd), ymax = n + sd, fill = Outcome), alpha = .2) + 
  scale_fill_manual(values = c("forestgreen", "darkred", "darkblue", "orange")) + 
  scale_color_manual(values = c("forestgreen", "darkred", "darkblue", "orange")) + 
  ylim(c(0, NA)) +
  theme_bw() + labs(y = "Abundance", x = "Time")
ConceptTime_gg

ggsave(
  ConceptTime_gg, 
  filename = file.path(Dir.Exports, "PLOT_ConceptualTimeSeries.png"), 
  width = 16*2, height = 9*2, units = "cm")

### Visualisation of Classification Cutoffs ----
# select single run for each outcome
RunsID <- c(unique(Abund_df$ID[Abund_df$pert.name == 10 & Abund_df$Outcome == "Extinction"])[1],
            unique(Abund_df$ID[Abund_df$pert.name == 10 & Abund_df$Outcome == "No Evolutionary Rescue"])[1],
            unique(Abund_df$ID[Abund_df$pert.name == 10 & Abund_df$Outcome == "Evolutionary Rescue"])[1],
            unique(Abund_df$ID[Abund_df$pert.name == 10 & Abund_df$Outcome == "Insufficient Crash"])[2])

# prepare runs for plotting
Abundplot_df <- Abund_df[Abund_df$ID %in% RunsID, ]
Abundplot_df <- Abundplot_df[Abundplot_df$t <= 1000, ]
Abundplot_df <- Abundplot_df[Abundplot_df$t >= 310, ]

# identify zones per run
Zones <- data.frame(Lower = Abundplot_df$n[Abundplot_df$t == 460]*0.1,
                    Upper = Abundplot_df$n[Abundplot_df$t == 460]*0.5,
                    Outcome = unique(Abundplot_df$Outcome))
AnnText <- data.frame(Y = c(1500, 500, 75, 2120, 2120),
                      X = c(980, 980, 980, 450, 750),
                      Outcome = "Extinction",
                      Colours = c("forestgreen", "orange", "darkblue", "black", "black"),
                      Labels = c("Evolutionary Rescue", "No Evolutionary Rescue",
                                 "Required Population Crash", "Pre-Perturbation", "Post-Perturbation"))
# plotting
ConceptZones_gg <- ggplot(Abundplot_df, aes(x = t, y = n)) + 
  geom_point(alpha = 0.1) +
  geom_line() + 
  ## Coloured Zones
  geom_rect(
    data = Zones,
    mapping = aes(xmin = 310, xmax = Inf, ymin = Lower, ymax = Upper,
                  x = NULL, y = NULL), fill = "orange", alpha = 0.2) +
  geom_rect(
    data = Zones,
    mapping = aes(xmin = 310, xmax = Inf, ymin = Upper, ymax = Inf,
                  x = NULL, y = NULL), fill = "forestgreen", alpha = 0.2) +
  geom_rect(
    data = Zones,
    mapping = aes(xmin = 310, xmax = Inf, ymin = 0, ymax = Lower,
                  x = NULL, y = NULL), fill = "darkblue", alpha = 0.2) +
  ## Labelling
  geom_vline(xintercept = 460, linetype="dotted", linewidth = 0.3) +
  geom_label(data = AnnText, hjust = 1, col = AnnText$Colours,
             aes(x = X, y = Y, label = Labels)) +
  ## Plot Limits and Facetting
  # ylim(c(0, 2200)) +
  xlim(c(310, NA)) +
  facet_wrap(~factor(Outcome, levels = c("Extinction", "Insufficient Crash", "No Evolutionary Rescue", "Evolutionary Rescue")), ncol = 1, scales = "free_y") + 
  theme_bw() + labs(y = "Abundance", x = "Time")
ConceptZones_gg

ggsave(
  ConceptZones_gg, 
  filename = file.path(Dir.Exports, "PLOT_ConceptualZones.png"), 
  width = 16*2, height = 9*5, units = "cm")

### Abundance; Full Data ----
# harmonise length of runs for plotting
Abundplot_df <- Abund_df[Abund_df$pert.name >= 9, ]
Abundplot_df <- Abundplot_df[Abundplot_df$t <= 1000, ]
Abundplot_df <- Abundplot_df[Abundplot_df$t >= 310, ]

# establish band ranges
Abundplot_df <- data.frame(
  aggregate(x = n ~ t+Outcome+pert.name+AC+DI+MU+SL+VA, data = Abundplot_df, FUN = mean),
  sd = aggregate(x = n ~ t+Outcome+pert.name+AC+DI+MU+SL+VA, data = Abundplot_df, FUN = sd)[,3])
Abundplot_df$sd[is.na(Abundplot_df$sd)] <- 0
colnames(Abundplot_df)[3:8] <- c("Perturbation Magnitude",
                                   "Spatial Autocorrelation", "Dispersal", "Mutation",
                                   "Spatial Slope", "Spatial Variation")


# plotting
## characteristics for data subsetting
AllChars <- Abundplot_df[
  !duplicated(Abundplot_df[, c("Perturbation Magnitude", "Spatial Autocorrelation", "Spatial Slope", "Spatial Variation")]), 
  c("Perturbation Magnitude", "Spatial Autocorrelation", "Spatial Slope", "Spatial Variation")]
AllChars <- AllChars[
  with(AllChars, 
       order(`Perturbation Magnitude`, 
             `Spatial Autocorrelation`,
             `Spatial Slope`,
             `Spatial Variation`
             )),
]

FullTime_ls <- pblapply(1:nrow(AllChars), FUN = function(x){
  Characteristics <- AllChars[x,]
  Subset <- which(
    apply(Abundplot_df[, c("Perturbation Magnitude", "Spatial Autocorrelation", "Spatial Slope", "Spatial Variation")], 
          MARGIN = 1, FUN = function(y){
            sum(y == Characteristics)
          }) == 4)
  FullTime_gg <- ggplot(Abundplot_df[Subset,], aes(x = t, y = n, 
                                                   fill = factor(Outcome, levels = c("Extinction", "Insufficient Crash", "No Evolutionary Rescue", "Evolutionary Rescue")), 
                                                   color = factor(Outcome, levels = c("Extinction", "Insufficient Crash", "No Evolutionary Rescue", "Evolutionary Rescue")))) + 
    geom_point(alpha = 0.1) +
    geom_line() + 
    geom_ribbon(aes(y = n, ymin = ifelse((n - sd) < 0, 0, n - sd), ymax = n + sd, fill = Outcome), alpha = .2) + 
    scale_fill_manual(values = c("forestgreen", "darkred", "darkblue", "orange"), name = "Outcome") + 
    scale_color_manual(values = c("forestgreen", "darkred", "darkblue", "orange"), name = "Outcome") + 
    facet_grid(Dispersal ~ Mutation, labeller = label_both) + 
    ylim(c(0, NA)) +
    theme_bw() + labs(y = "Abundance", x = "Time", 
                      title = paste("Perturbation = ", Characteristics[1],
                                    " ; Spatial Autocorrelation = ", Characteristics[2],
                                    " ; Spatial Slope = ", Characteristics[3],
                                    " ; Spatial Variation = ", Characteristics[4])) + 
    theme(legend.position = "bottom")
  FullTime_gg
})
pdf(file.path(Dir.Exports, "PLOT_FullTimeSeries.pdf"), width = 11, height = 8)
for(i in 1:length(FullTime_ls)){
  print(FullTime_ls[[i]])
}
dev.off()

## Environment Parameters  ------------------------------------------------
### Environmental Parameterisation and Actual Environmental Values --------
load(file.path(Dir.Exports, "ENVIRONMENT_Parameters.RData"))
SL_labels <- c(`0.8` = "Slope = 0.8",
               `1` = "Slope = 1.0",
               `1.2` = "Slope = 1.2")

AC_comparisons <- rollapply(as.character(sort(unique(Enviro_table$AC_input))), 2, c)
AC_comparisons <- split(AC_comparisons, row(AC_comparisons))
AC_gg <- ggplot(Enviro_table, aes(x = factor(AC_input), y = AC_data)) + 
  geom_boxplot() + 
  facet_wrap(~factor(SL), labeller = as_labeller(SL_labels)) + 
  labs(x = "Autocorrelation Parameterisation", y = "Simulated Autocorrelation") + 
  # stat_compare_means(comparisons = AC_comparisons, 
  #                    label = 'p.signif') + 
  theme_bw()
AC_gg2 <- ggplot(Enviro_table, aes(x = AC_input, y = AC_data, 
                                   group = SL
                                   , fill = factor(SL), col = factor(SL)
)) +
  geom_smooth() + 
  labs(x = "Autocorrelation Parameterisation", y = "Simulated Autocorrelation",
       fill = "Slope", col = "Slope") + 
  theme_bw() + theme(legend.position = "top")

VA_comparisons <- rollapply(as.character(sort(unique(Enviro_table$VA_input))), 2, c)
VA_comparisons <- split(VA_comparisons, row(VA_comparisons))
VA_gg <- ggplot(Enviro_table, aes(x = factor(VA_input), y = VA_data)) + 
  geom_boxplot() + 
  # facet_wrap(~factor(SL), labeller = as_labeller(SL_labels)) + 
  labs(x = "Autocorrelation Variation", y = "Simulated Variation") + 
  stat_compare_means(comparisons = VA_comparisons, 
                     label = 'p.signif') + 
  theme_bw() 
VA_gg2 <- ggplot(Enviro_table, aes(x = VA_input, y = VA_data, 
                                   group = SL
                                   , fill = factor(SL), col = factor(SL)
)) +
  geom_smooth() + 
  labs(x = "Autocorrelation Variation", y = "Simulated Variation",
       fill = "Slope", col = "Slope") + 
  theme_bw() + theme(legend.position = "top")

SpatParam_gg <- cowplot::plot_grid(
  cowplot::plot_grid(AC_gg, AC_gg2, ncol = 2, rel_widths = c(2, 1)),
  cowplot::plot_grid(VA_gg, VA_gg2, ncol = 2, rel_widths = c(2, 1)),
  ncol = 1, labels = "auto"
  )
SpatParam_gg
ggsave(SpatParam_gg, filename = file.path(Dir.Exports, "PLOT_SpatialParameterisation.png"), width = 16*2, height = 16*2, units = "cm")

### Distribution Comparisons ----
#### Non-Spatially Explicit ----
load(file.path(Dir.Exports, "DISTRIBUTIONS_NonSpatial.RData"))
DISTRIBUTIONS_NonSpatial
DISTRIBUTIONS_NonSpatial$Comparisons <- paste(
  paste("[Environment]", DISTRIBUTIONS_NonSpatial$Envir), 
  "vs.", 
  paste("[Individuals]", DISTRIBUTIONS_NonSpatial$Indivs)
)
DISTRIBUTIONS_NonSpatial$Comparisons <- 
  factor(DISTRIBUTIONS_NonSpatial$Comparisons,
         levels = c(
           "[Environment] Pre vs. [Individuals] Pre",
           "[Environment] Post vs. [Individuals] PostMin",
           "[Environment] Post vs. [Individuals] PostMax"
         ))
DISTRIBUTIONS_NonSpatial <- DISTRIBUTIONS_NonSpatial[DISTRIBUTIONS_NonSpatial$Pert>=9, ]

OV_Dens_gg <- pblapply(sort(unique(DISTRIBUTIONS_NonSpatial$AC)), FUN = function(AC_iter){
  AC_df <- DISTRIBUTIONS_NonSpatial[DISTRIBUTIONS_NonSpatial$AC == AC_iter,]
  SL_ls <- lapply(sort(unique(DISTRIBUTIONS_NonSpatial$SL)), FUN = function(SL_iter){
    SL_df <- AC_df[AC_df$SL == SL_iter,]
    VA_ls <- lapply(sort(unique(DISTRIBUTIONS_NonSpatial$VA)), FUN = function(VA_iter){
      VA_df <- SL_df[SL_df$VA == VA_iter,]
      ggplot(VA_df, aes(y = OV_Dens, x = Pert,
                        fill = Comparisons,
                        col = Comparisons)) + 
        geom_point(size = 3) + 
        geom_smooth(method="loess") + 
        labs(y = "Overlap of Density Distributions", x = "Perturbation Magnitude",
             title = paste(
               "Autocorrelation =", AC_iter,
               "; Slope =", SL_iter,
               "; Variation =", VA_iter
             )) + 
        facet_grid(DI ~ MU, labeller = label_both) + 
        theme_bw() + theme(legend.position = "bottom") + ylim(c(0,1)) + xlim(c(9,16))
    })
    cowplot::plot_grid(
      cowplot::plot_grid(plotlist = lapply(VA_ls, FUN = function(x){
        x+theme(legend.position = "none")
      }), ncol = 3), 
      as_ggplot(get_legend(VA_ls[[1]])), ncol = 1, rel_heights = c(1, 0.1)
    )
  })
  cowplot::plot_grid(plotlist = SL_ls, ncol = 1)
})
pdf(file.path(Dir.Exports, "PLOT_NonSpatialOverlapDensity.pdf"), 
    width = 16, height = 16)
for(i in 1:length(OV_Dens_gg)){
  print(OV_Dens_gg[[i]])
}
dev.off()

#### Spatially Explicit ----
load(file.path(Dir.Exports, "DISTRIBUTIONS_Spatial.RData"))
DISTRIBUTIONS_Spatial <- DISTRIBUTIONS_Spatial[["Summary_df"]]
rownames(DISTRIBUTIONS_Spatial) <- c()
head(DISTRIBUTIONS_Spatial)
stop("Continue here")

## Analyses Inputs and Outcomes -------------------------------------------
### Survival or Extinction ----
### Survival or Evolutionary Rescue ----
### Characteristics of Evolutionary Rescue ----