#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes] 
#' CONTENTS: 
#'  - Statistical Models of
#'    - Survival Likelihood
#'    - Evolutionary Rescue Likelihood
#'    - Evolutionary Rescue Success
#'  DEPENDENCIES:
#'  - EVORES_Metrics.RData produced by "4 - Evolutionary Rescue Success Metrics.R"
#'  - DISTRIBUTIONS_Spatial.RData produced by "5 - DistributionComparison.R"
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
  "brms",
  "ggplot2",
  "tidybayes",
  "plyr",
  "cowplot",
  "pbapply"
)
sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

LogisticModelDiagnostics <- function(Model_brm, title = "NULL TITLE"){
  ### Title ----
  title <- ggdraw() + 
    draw_label(
      title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  ### Model Diagnostics ----
  # model diagnostics
  # summary(Model_brm) # all coefficients significant
  if(any(abs(rhat(Model_brm)-1) > 1.05)){# everything converged
    message("Some rhat values are worryingly large")
  }
  # plot(Model_brm) # everything converged
  # neff_ratio(Model_brm) # looks alright
  ppcheck_gg <- pp_check(Model_brm, type = "bars", ndraws = 42)
  
  ### Model Outputs ----
  # #### Predicted survival
  # EVORES_Metrics |>  
  #   add_epred_draws(Model_brm, ndraws = 50) |> 
  #   ggplot() +
  #   geom_boxplot(aes(x = factor(pert.name), y = .epred), width = 0.5) + 
  #   facet_grid(DI ~ MU, labeller = label_both) +
  #   labs(x = "Perturbation Magnitude", y = "Probability of Survival") + 
  #   theme_bw()
  # sum(!is.na(EVORES_Metrics$ID))
  
  #### Coefficient Plotting
  ## extract posterior draws
  plot_coeffs <- Model_brm |> 
    gather_draws(`^b_.*`, regex = TRUE) |>
    mutate(.panel = ifelse((grepl(.variable, pattern = "Intercept") | grepl(.variable, pattern = "pert.name")),
                           "Intercept", "Coefficient")) 
  ## significance of coefficients
  signif_df <- data.frame(coeff = rownames(summary(Model_brm)[["fixed"]]),
                          signif = abs(sign(summary(Model_brm)[["fixed"]]$`l-95% CI`) + 
                                         sign(summary(Model_brm)[["fixed"]]$`u-95% CI`)) == 2)
  plot_coeffs$signif <- signif_df$signif[match(gsub(plot_coeffs$.variable, pattern = "b_", replacement = ""), signif_df$coeff)]
  ## define plotting panel names
  plot_coeffs$.panel[
    which(plot_coeffs$.variable %in% unique(grep(paste(c("MU", "DI"),collapse = "|"), plot_coeffs$.variable, value = TRUE))) 
  ] = "Population"
  plot_coeffs$.panel[
    which(plot_coeffs$.variable %in% unique(grep(paste(c("AC", "SL", "VA"),collapse = "|"), plot_coeffs$.variable, value = TRUE)))
  ] = "Environment"
  plot_coeffs$.variable <- gsub(plot_coeffs$.variable, pattern = "b_", replacement = "")
  plot_coeffs$.variable[which(plot_coeffs$.variable == "pert.name")] <- "PM"
  ## define order of coefficients in plotting
  varfacs <- sort(unique(plot_coeffs$.variable))
  varfacs <- c("Intercept", "PM", varfacs[!(varfacs %in% c("Intercept", "PM"))])
  ## plotting
  Coeffplots <- plot_coeffs |> 
    ggplot(aes(x = .value, y = factor(.variable, levels = rev(varfacs)), fill = signif)) +
    geom_vline(xintercept = 0) + 
    stat_halfeye(normalize = "xy") + 
    facet_wrap(~ factor(.panel, levels = c("Intercept", "Population", "Environment")),
               scales = "free_y",
               ncol = 1) +
    theme_bw() + theme(legend.position = "none") + 
    labs(x = "Log Odds / Log Odds Change", y = "Coefficient")
  
  #### Combined Plots
  cowplot::plot_grid(
    title,
    cowplot::plot_grid(ppcheck_gg, Coeffplots),
    rel_heights = c(0.1, 1), ncol = 1
    )
}

MultiLevelModelDiagnostics <- function(Model_brm, type = c("cat", "cont", "cont")){
  ### Model Diagnostics ----
  # model diagnostics
  # summary(Model_brm) # all coefficients significant
  if(any(rhat(Model_brm) > 1.05)){# everything converged
    message("Some rhat values are worryingly large")
  }
  
  plot_coeffs <- Model_brm |> 
    gather_draws(`^b_.*`, regex = TRUE) |>
    mutate(.panel = ifelse((grepl(.variable, pattern = "Intercept") | grepl(.variable, pattern = "pert.name")),
                           "Intercept", "Coefficient")) 
  plot_coeffs$.variable <- gsub(plot_coeffs$.variable, pattern = "mu1_", replacement = "")
  plot_coeffs$.variable <- gsub(plot_coeffs$.variable, pattern = "mu2_", replacement = "")
  ## significance of coefficients
  plot_coeffs$.value
  
  Quantiles <- aggregate(.value~.variable, data = plot_coeffs, FUN = quantile, probs = c(0.05, 0.95))
  Quantiles$signif <- abs(sign(Quantiles$.value[,1]) + sign(Quantiles$.value[,2])) == 2
  plot_coeffs$signif <- Quantiles$signif[match(gsub(plot_coeffs$.variable, pattern = "b_", replacement = ""), Quantiles$.variable)]
  ## define plotting panel names
  plot_coeffs$.panel[
    which(plot_coeffs$.variable %in% unique(grep(paste(c("MU", "DI"),collapse = "|"), plot_coeffs$.variable, value = TRUE))) 
  ] = "Population"
  plot_coeffs$.panel[
    which(plot_coeffs$.variable %in% unique(grep(paste(c("AC", "SL", "VA"),collapse = "|"), plot_coeffs$.variable, value = TRUE)))
  ] = "Environment"
  plot_coeffs$.variable <- gsub(plot_coeffs$.variable, pattern = "b_", replacement = "")
  plot_coeffs$.variable[which(plot_coeffs$.variable == "pert.name")] <- "PM"
  ## define order of coefficients in plotting
  varfacs <- sort(unique(plot_coeffs$.variable))
  varfacs <- c("Intercept", "PM", varfacs[!(varfacs %in% c("Intercept", "PM"))])
  
  plots_ls <- lapply(1:length(Model_brm[["formula"]][["responses"]]),
                     FUN = function(x){
                       print(paste("Hierarchy Layer", x))
                       
                       if(type[x] == "cat"){
                         pp_gg <- pp_check(Model_brm, ndraws = 50,
                                           resp = Model_brm[["formula"]][["responses"]][x],
                                           type = "bars") 
                       }else{
                         pp_gg <- pp_check(Model_brm, ndraws = 50,
                                           resp = Model_brm[["formula"]][["responses"]][x])
                       }
                       
                       plot_iter <- plot_coeffs[startsWith(plot_coeffs$.variable, Model_brm[["formula"]][["responses"]][x]), ]
                       plot_iter$.variable <- gsub(pattern = paste0(Model_brm[["formula"]][["responses"]][x], "_"), 
                                                   replacement = "", plot_iter$.variable)
                       coeff_gg <- plot_iter |> 
                         ggplot(aes(x = .value, 
                                    y = factor(.variable, levels = unique(rev(gsub(pattern = paste0(Model_brm[["formula"]][["responses"]][x], "_"),
                                                                                   replacement = "", varfacs)))), 
                                    fill = signif)) +
                         geom_vline(xintercept = 0) + 
                         stat_halfeye(normalize = "xy") + 
                         facet_wrap(~ factor(.panel, levels = c("Intercept", "Population", "Environment", "Coefficient")),
                                    scales = "free",
                                    ncol = 1) +
                         theme_bw() + theme(legend.position = "none") + 
                         labs(x = "Log Odds / Log Odds Change", y = "Coefficient")
                       
                       gg_grid <- cowplot::plot_grid(
                         ggdraw() + 
                           draw_label(
                             Model_brm[["formula"]][["responses"]][x],
                             fontface = 'bold',
                             x = 0,
                             hjust = 0
                           ) +
                           theme(
                             # add margin on the left of the drawing canvas,
                             # so title is aligned with left edge of first plot
                             plot.margin = margin(0, 0, 0, 7)
                           ),
                         cowplot::plot_grid(pp_gg, coeff_gg, rel_widths = c(0.4, 0.6)),
                         rel_heights = c(0.1, 1), ncol = 1)
                       print(gg_grid)
                     }) 
}

# DATA ====================================================================
if(file.exists(file.path(Dir.Exports, "MODEL_Metrics.RData"))){
  load(file.path(Dir.Exports, "MODEL_Metrics.RData"))
}else{
  ## Simulation Outcomes ----------------------------------------------------
  load(file.path(Dir.Exports, "EVORES_Metrics.RData"))
  EVORES_Metrics <- EVORES_Metrics[EVORES_Metrics$pert.name >= 9, ]
  EVORES_Metrics$EvoRes[!EVORES_Metrics$SuffDip] <- "Insufficient Population Crash"
  EVORES_Metrics$survival <- as.numeric(EVORES_Metrics$survival)
  
  ## Population-Level Summaries ---------------------------------------------
  load(file.path(Dir.Exports, "POPULATION_TimeStep.RData"))
  
  ## Distribution Comparisons -----------------------------------------------
  ### Non-Spatially Explicit ----
  load(file.path(Dir.Exports, "DISTRIBUTIONS_NonSpatial.RData"))
  DISTRIBUTIONS_NonSpatial <- DISTRIBUTIONS_NonSpatial[DISTRIBUTIONS_NonSpatial$Pert>=9, ]
  
  ### Spatially Explicit ----
  load(file.path(Dir.Exports, "DISTRIBUTIONS_Spatial.RData"))
  DISTRIBUTIONS_Spatial <- DISTRIBUTIONS_Spatial[["Summary_df"]]
  rownames(DISTRIBUTIONS_Spatial) <- c()
  DISTRIBUTIONS_Spatial <- DISTRIBUTIONS_Spatial[DISTRIBUTIONS_Spatial$Pert>=9, ]
  
  ### Fusing ----
  DISTRIBUTION_Metrics <- base::merge(DISTRIBUTIONS_Spatial, DISTRIBUTIONS_NonSpatial)
  DISTRIBUTION_Metrics$Comparisons <- paste(
    paste("[Environment]", DISTRIBUTION_Metrics$Envir), 
    "vs.", 
    paste("[Individuals]", DISTRIBUTION_Metrics$Indivs)
  )
  DISTRIBUTION_Metrics$Comparisons <- 
    factor(DISTRIBUTION_Metrics$Comparisons,
           levels = c(
             "[Environment] Pre vs. [Individuals] Pre",
             "[Environment] Post vs. [Individuals] PostMin",
             "[Environment] Post vs. [Individuals] PostMax"
           ))
  colnames(DISTRIBUTION_Metrics)[4] <- "pert.name"
  
  ## Modelling Data Frame----------------------------------------------------
  EVORES_Metrics$MatchID <- with(EVORES_Metrics, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_"))
  Data_df$MatchID <- with(Data_df, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_")) # population-level metrics
  DISTRIBUTION_Metrics$MatchID <- with(DISTRIBUTION_Metrics, paste(pert.name, rep, AC, DI, MU, SL, VA, sep = "_"))
  
  ## attach population and distribution metrics to EVORES_Metrics
  Combination_ls <- pblapply(1:nrow(EVORES_Metrics), cl = parallel::detectCores(), FUN = function(x){
    Iter_df <- EVORES_Metrics[x,c("t_minpost", "t_maxpost", "MatchID")]
    
    ## adding population metrics
    Data_iter <- Data_df[Data_df$MatchID == Iter_df$MatchID, ]
    t_vec <- c(pre = 460, postmin = Iter_df$t_minpost, postmax = Iter_df$t_maxpost)
    Data_iter <- Data_iter[match(t_vec, Data_iter$t), ]
    fuse_df <- lapply(1:nrow(Data_iter), FUN = function(y){
      ret_df <- Data_iter[y, colnames(Data_df) %nin% c("t", "pert.name", "rep", "AC", "DI", "MU", "SL", "VA", "ID", "MatchID")]
      colnames(ret_df) <- paste(names(t_vec)[y], colnames(ret_df), sep = "_")
      ret_df
    })
    Pop_add <- cbind(do.call(cbind, fuse_df)
                     # , "MatchID" = Data_iter$MatchID[1]
    )
    
    Dist_iter <- DISTRIBUTION_Metrics[DISTRIBUTION_Metrics$MatchID == Iter_df$MatchID,]
    if(nrow(Dist_iter)==6){Dist_iter <- Dist_iter[c(1,3,5), ]} # happens for three runs that we have two dsitribution overlap comparisons
    if(nrow(Dist_iter) == 0){
      Dist_iter[1:3, ] <- NA
      Dist_iter$Comparisons <- rev(levels(DISTRIBUTION_Metrics$Comparisons))
    }
    fuse_df <- lapply(1:nrow(Dist_iter), FUN = function(y){
      ret_df <- Dist_iter[y, colnames(Dist_iter) %nin% c("t", "pert.name", "rep", "AC", "DI", "MU", "SL", "VA", "ID", "MatchID", "Comparisons")]
      colnames(ret_df) <- paste(as.character(Dist_iter$Comparisons)[y], colnames(ret_df), sep = "_")
      ret_df
    })
    Dist_add <- cbind(do.call(cbind, fuse_df), "MatchID" = Data_iter$MatchID[1])
    
    cbind(Pop_add, Dist_add)
  })
  Combination_df <- do.call(rbind, Combination_ls)
  MODEL_Metrics <- join(EVORES_Metrics, Combination_df)
  save(MODEL_Metrics, file = file.path(Dir.Exports, "MODEL_Metrics.RData")) 
}

# MODEL_Metrics$Recovery <- sqrt((MODEL_Metrics$n_pre - MODEL_Metrics$n_maxpost)^2 + (MODEL_Metrics$t_maxpost - MODEL_Metrics$t_minpost)^2) # Higher value = lower resistance, because either big dip in abundance or long time of decline until minimum
# MODEL_Metrics$Resistance <- sqrt((MODEL_Metrics$n_pre - MODEL_Metrics$n_minpost)^2 + (MODEL_Metrics$t_minpost - 450)^2)

MODEL_Metrics$perc_t_maxpost <- (MODEL_Metrics$t_maxpost - MODEL_Metrics$t_minpost)/(MODEL_Metrics$t-MODEL_Metrics$t_minpost)
MODEL_Metrics$perc_t_minpost <- (MODEL_Metrics$t_minpost - 450)/450

load(file.path(Dir.Exports, "POPULATION_TimeStep.RData"))
POPULATION_Timestep <- Data_df
rm(Data_df)

# MODELS ==================================================================
### Intermediate Optimum Slope --------------------
#' how does slope affect trait values?
ggplot(POPULATION_Timestep[
  POPULATION_Timestep$MU == 0 &
  POPULATION_Timestep$DI == 2 &
  POPULATION_Timestep$t == 450
  ,
],
       aes(y = u_sd, x = factor(SL))) + 
  geom_boxplot() + 
  facet_grid(AC ~ VA, labeller = label_both) + 
  theme_bw() + 
  labs(title = "Pre-perturbation trait variation - MU = 0, DI = 2")

ggplot(POPULATION_Timestep[
  POPULATION_Timestep$MU == 1 &
    POPULATION_Timestep$DI == 1.5 &
    POPULATION_Timestep$t == 450
  ,
],
aes(y = u_sd, x = factor(SL))) + 
  geom_boxplot() + 
  facet_grid(AC ~ VA, labeller = label_both) + 
  theme_bw() + 
  labs(title = "Pre-perturbation trait variation - MU = 1, DI = 1.5")

### Resilience Regression -------------------------
stop("SPLIT OUT Speediness and Change in abundance seperately before attempting aggregate metrics")


if(file.exists(file.path(Dir.Exports, "ContEvoRes_MULTILEVEL_brm.RData"))){
  load(file.path(Dir.Exports, "ContEvoRes_MULTILEVEL_brm.RData"))
}else{
  Bayes_df <- MODEL_Metrics[MODEL_Metrics$survival == 1, ]
  Bayes_df$Recovery <- Bayes_df$perc_maxpostpre * Bayes_df$perc_t_maxpost
  Bayes_df$Resistance <- Bayes_df$perc_minpost * Bayes_df$perc_t_minpost
  
  ContEvoRes_MULTILEVEL_brm <- brm(
    ## model formulae
    # Recovery
    bf(Recovery ~ perc_minpost + postmin_u_sd/postmax_u_sd) +  
    # Resistance
    bf(Resistance  ~ pert.name + abs(SL-1) + AC + pre_u_sd,
         family = hurdle_gamma()
      )
    , 
    ## data
    data = Bayes_df,
    ## settings for run
    chains = 4,
    cores = parallel::detectCores(),
    iter = 1e4,
    warmup = 3e3,
    seed = 42)
  save(ContEvoRes_MULTILEVEL_brm, file = file.path(Dir.Exports, "ContEvoRes_MULTILEVEL_brm.RData"))
}
MultiLevelModelDiagnostics(ContEvoRes_MULTILEVEL_brm, type = c("cont", "cont"))

if(file.exists(file.path(Dir.Exports, "ContEvoRes_MULTILEVEL2_brm.RData"))){
  load(file.path(Dir.Exports, "ContEvoRes_MULTILEVEL2_brm.RData"))
}else{
  Bayes_df <- MODEL_Metrics[MODEL_Metrics$survival == 1, ] #  & MODEL_Metrics$MU == 0
  Bayes_df$Recovery <- rowMeans(cbind(Bayes_df$perc_maxpostpre, Bayes_df$perc_t_maxpost))
  Bayes_df$Resistance <- rowMeans(cbind(Bayes_df$perc_minpost, Bayes_df$perc_t_minpost))
  ContEvoRes_MULTILEVEL2_brm <- brm(
    ## model formulae
    # Recovery
    bf(Recovery ~ perc_minpost * postmin_u_sd/postmax_u_sd * MU,
       family = hurdle_gamma()
       ) +  
    # Resistance
    bf(Resistance ~ pert.name + abs(SL-1) + AC + pre_u_sd,
         family = hurdle_gamma()
      )
    , 
    ## data
    data = Bayes_df,
    ## settings for run
    chains = 4,
    cores = parallel::detectCores(),
    iter = 1e4,
    warmup = 3e3,
    seed = 42)
  save(ContEvoRes_MULTILEVEL2_brm, file = file.path(Dir.Exports, "ContEvoRes_MULTILEVEL2_brm.RData"))
}
MultiLevelModelDiagnostics(ContEvoRes_MULTILEVEL2_brm, type = c("cont", "cont"))



### Population Trajectory Identification ----------
#' Trajectories:
#'  - Extinction
#'  - No Rescue
#'  - Rescue
#'    --> Mutation-facilitated
#'    --> Non-Mutation facilitated










