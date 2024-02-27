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
  "pbapply"
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

### Combinations of Simulation Settings ----
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
Combs_df$`Evolutionary Rescue` <- round(Combs_df$`Evolutionary Rescue`/Combs_df$`Total Simulations`, 3)*100

CombsBase_df$Survival <- as.logical(CombsBase_df$Survival)
Combs_df$Survival <- aggregate(Survival ~ `Spatial Autocorrelation` + `Spatial Slope` + 
                        `Spatial Variation` + 
                        Dispersal + Mutation, 
                      data = CombsBase_df, FUN = sum)$Survival
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
Combs_gg <- plot_grid(plotlist = CombPlots_ls, ncol = 1, labels = "auto")
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
EvoResSuc_df[EvoResSuc_df$pert.name >= 9, ]

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
  scale_fill_manual(values = c("forestgreen", "darkred", "orange")) + 
  scale_color_manual(values = c("forestgreen", "darkred", "orange")) + 
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
            unique(Abund_df$ID[Abund_df$pert.name == 10 & Abund_df$Outcome == "Evolutionary Rescue"])[1])

# prepare runs for plotting
Abundplot_df <- Abund_df[Abund_df$ID %in% RunsID, ]
Abundplot_df <- Abundplot_df[Abundplot_df$t <= 1000, ]
Abundplot_df <- Abundplot_df[Abundplot_df$t >= 310, ]

# identify zones per run
Zones <- data.frame(Lower = Abundplot_df$n[Abundplot_df$t == 460]*0.1,
                    Upper = Abundplot_df$n[Abundplot_df$t == 460]*0.5,
                    Outcome = unique(Abundplot_df$Outcome))

# plotting
ConceptZones_gg <- ggplot(Abundplot_df, aes(x = t, y = n)) + 
  geom_point(alpha = 0.1) +
  geom_line() + 
  geom_rect(
    data = Zones,
    mapping = aes(xmin = 310, xmax = Inf, ymin = Lower, ymax = Upper,
                  x = NULL, y = NULL), fill = "orange", alpha = 0.2) +
  geom_rect(
    data = Zones,
    mapping = aes(xmin = 310, xmax = Inf, ymin = Upper, ymax = Inf,
                  x = NULL, y = NULL), fill = "forestgreen", alpha = 0.2) +
  ylim(c(0, NA )) +
  xlim(c(310, NA )) +
  facet_wrap(~factor(Outcome, levels = c("Extinction", "No Evolutionary Rescue", "Evolutionary Rescue")), ncol = 1) + 
  theme_bw() + labs(y = "Abundance", x = "Time")
ConceptZones_gg

ggsave(
  ConceptZones_gg, 
  filename = file.path(Dir.Exports, "PLOT_ConceptualZones.png"), 
  width = 16*2, height = 9*4, units = "cm")

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
                                                   fill = factor(Outcome, levels = c("Extinction", "No Evolutionary Rescue", "Evolutionary Rescue")), 
                                                   color = factor(Outcome, levels = c("Extinction", "No Evolutionary Rescue", "Evolutionary Rescue")))) + 
    geom_point(alpha = 0.1) +
    geom_line() + 
    geom_ribbon(aes(y = n, ymin = ifelse((n - sd) < 0, 0, n - sd), ymax = n + sd, fill = Outcome), alpha = .2) + 
    scale_fill_manual(values = c("darkred", "orange", "forestgreen"), name = "Outcome") + 
    scale_color_manual(values = c("darkred", "orange", "forestgreen"), name = "Outcome") + 
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


### Other Population Metrics; Full Data ----


## Environment Parameters  ------------------------------------------------

### Environmental Parameterisation and Actual Environmental Values --------

### Distribution Comparisons ----
#### Non-Spatially Explicit ----
#### Spatially Explicit ----

## Analyses Inputs and Outcomes -------------------------------------------
### Survival or Extinction ----
### Survival or Evolutionary Rescue ----
### Characteristics of Evolutionary Rescue ----