#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Display Item Creation
#'  DEPENDENCIES:
#'  - All previous scripts have been run
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
  "dplyr"
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
EvoResSuc_df <- read.csv(file.path(Dir.Exports, "EvoResSuccessMetrics.csv"))
EvoResSuc_df$EvoRes[is.na(EvoResSuc_df$EvoRes)] <- "Not Possible"

# data extraction
sankey_df <- EvoResSuc_df[, c("pert.name", "survival", "EvoRes")]
sankey_df$Run <- "Total Simulations"
colnames(sankey_df) <- c("Perturbation Magnitude", "Survival", "Evolutionary Rescue", "Simulation")

# data reformatting for plotting
sankey_df <- sankey_df[sankey_df$`Perturbation Magnitude` >= 9, ]
sankey_df$`Perturbation Magnitude` <- str_pad(sankey_df$`Perturbation Magnitude`, width = 2, side = "left", pad = "0")
sankey_df$Survival <- str_replace_all(sankey_df$Survival, c("TRUE" = "Survival", "FALSE" = "Extinction"))
sankey_df$`Evolutionary Rescue` <- str_replace_all(sankey_df$`Evolutionary Rescue`, c("TRUE" = "Evolutionary Rescue", "FALSE" = "No Evolutionary Rescue"))

# data orientation for plotting
sankey_df <- make_long(sankey_df, "Simulation", "Perturbation Magnitude", "Survival", "Evolutionary Rescue")
sankey_count <- sankey_df%>%
  dplyr::group_by(node)%>%
  tally()
sankey_df <- merge(sankey_df, sankey_count, by.x = 'node', by.y = 'node', all.x = TRUE)

sankey_gg <- ggplot(sankey_df, aes(x = x
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
  scale_fill_viridis_d(option = "G", na.translate = FALSE)
sankey_gg

## Simulation Abundance Time-Series; Conceptual ---------------------------
