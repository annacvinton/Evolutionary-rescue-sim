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

LogisticModelDiagnostics <- function(Model_brm){
  ### Model Diagnostics ----
  # model diagnostics
  # summary(Model_brm) # all coefficients significant
  if(any(abs(rhat(Model_brm)-1) > 1.05)){# everything converged
    message("Some rhat values are worryingly large")
  }
  # plot(Model_brm) # everything converged
  # neff_ratio(Model_brm) # looks alright
  ppcheck_gg <- pp_check(Model_brm, type = "bars") # looks alright
  
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
  print(cowplot::plot_grid(ppcheck_gg, Coeffplots))
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

# MODELS ==================================================================
## Survival ---------------------------------------------------------------
Dir.LogiSurv <- file.path(Dir.Exports, "Survival_Logistic")
if(!dir.exists(Dir.LogiSurv)){dir.create(Dir.LogiSurv)}

### Model Fitting ----
## split by Mutation, Dispersal, and perturbation magnitude, because of the questions we ask
## can I average model coefficients
## expand.grid(unique(EVORES_Metrics$pert.name), unique(EVORES_Metrics$DI), c(0,1))
## t-tests between posterior samples against splitting factors
## separate model of all data just without MU and DI (split by these)

##### SINGLE MODEL
if(file.exists(file.path(Dir.LogiSurv, "BiSurv_brm.RData"))){
  load(file.path(Dir.LogiSurv, "BiSurv_brm.RData"))
}else{
  BiSurv_brm <- brm(formula = survival ~ pert.name + SL + DI * MU,
                            data = MODEL_Metrics,
                            family = bernoulli(link = "logit"),
                            warmup = 3e3,
                            iter = 1e4,
                            chains = 4,
                            cores = 4,
                            seed = 42)
  save(BiSurv_brm, file = file.path(Dir.LogiSurv, "BiSurv_brm.RData"))
}
BiSurvFinal_brm <- BiSurv_brm
LogisticModelDiagnostics(BiSurvFinal_brm)

##### SPLITTING
ModelCombs <- with(na.omit(MODEL_Metrics), expand.grid(unique(DI), unique(MU)))
colnames(ModelCombs) <- c("DI", "MU")

BiSurv_ls <- lapply(1:nrow(ModelCombs), FUN = function(iter){
  FNAME <- file.path(Dir.LogiSurv, paste0("BiSurv_brm_DI-", ModelCombs$DI[iter], "_MU-", ModelCombs$MU[iter],".RData"))
  if(file.exists(FNAME)){
    load(FNAME)
  }else{
    BiSurv_brm <- brm(formula = survival ~ pert.name + VA * AC * SL,
                        data = MODEL_Metrics[
                          (MODEL_Metrics$DI == ModelCombs$DI[iter] &
                             MODEL_Metrics$MU == ModelCombs$MU[iter]),
                        ],
                        family = bernoulli(link = "logit"),
                        warmup = 3e3,
                        iter = 1e4,
                        chains = 4,
                        cores = 4,
                        seed = 42)
    # Sys.sleep(60*30)
    save(BiSurv_brm, file = FNAME)
  }
  LogisticModelDiagnostics(BiSurv_brm)
  BiSurv_brm
})

## Evolutionary Rescue ----------------------------------------------------
Dir.LogiEvoRes <- file.path(Dir.Exports, "EvoRes_Logistic")
if(!dir.exists(Dir.LogiEvoRes)){dir.create(Dir.LogiEvoRes)}

### Logistic Regression ---------------------------
EVORES_Metrics <- MODEL_Metrics[MODEL_Metrics$EvoRes != "Insufficient Population Crash", ]
EVORES_Metrics$EvoRes <- as.numeric(as.logical(EVORES_Metrics$EvoRes))

#### Model Fitting ----
##### SINGLE MODEL
if(file.exists(file.path(Dir.LogiEvoRes, "BiEvoRes_brm.RData"))){
  load(file.path(Dir.LogiEvoRes, "BiEvoRes_brm.RData"))
}else{
  BiEvoRes_brm <- brm(formula = EvoRes ~ pert.name + VA * AC * abs(SL-1) + DI + factor(MU),
                      data = EVORES_Metrics,
                      family = bernoulli(link = "logit"),
                      warmup = 3e3,
                      iter = 1e4,
                      chains = 4,
                      cores = 4,
                      seed = 42)
  save(BiEvoRes_brm, file = file.path(Dir.LogiEvoRes, "BiEvoRes_brm.RData"))
}
BiEvoResFinal_brm <- BiEvoRes_brm
LogisticModelDiagnostics(BiEvoResFinal_brm)

##### SPLITTING
ModelCombs <- with(na.omit(EVORES_Metrics), expand.grid(unique(DI), unique(MU)))
colnames(ModelCombs) <- c("DI", "MU")

BiEvoRes_ls <- lapply(1:nrow(ModelCombs), FUN = function(iter){
  FNAME <- file.path(Dir.LogiEvoRes, paste0("BiEvoRes_brm_DI-", ModelCombs$DI[iter], "_MU-", ModelCombs$MU[iter],".RData"))
  if(file.exists(FNAME)){
    load(FNAME)
  }else{
    BiEvoRes_brm <- brm(formula = EvoRes ~ pert.name + VA * AC * abs(SL-1),
                        data = EVORES_Metrics[
                          (EVORES_Metrics$DI == ModelCombs$DI[iter] &
                            EVORES_Metrics$MU == ModelCombs$MU[iter]),
                        ],
                        family = bernoulli(link = "logit"),
                        warmup = 3e3,
                        iter = 1e4,
                        chains = 4,
                        cores = 4,
                        seed = 42)
    # Sys.sleep(60*30)
    save(BiEvoRes_brm, file = FNAME)
  }
  LogisticModelDiagnostics(BiEvoRes_brm)
  BiEvoRes_brm
})

### MULTILEVEL TRIAL -----
### to try:
#' Pre-perturbation population metrics (u_sd, etc.) are driven landscape parameters?
#' Pre-perturbation population metrics then predict post-perturbation population metrics together with perturbation magnitude and population parameters?
#' Evolutionary rescue is driven post-perturbation population metrics as well as population and landscape parameters
#' Might want to consider splitting the resulting model(s) by mutation and/or dispersal
#' Might want to run one full interaction model
#' Maybe run initially for MU = 0 and DI = 2 as well as MU = 1 and DI = 1.5?

# Combine the formulas into a single model
multi_level_model2 <- brm(
  # model formulae
  bf(EvoRes ~ PoMiTrSD, family = bernoulli()) + 
  bf(PoMiTrSD ~ 0 + factor(MU) * abs(SL-1) * DI, family = Gamma()),
  # data
  data = EvoResMod_df[EvoResMod_df$PoMiTrSD != 0,],
  # settings for run
  chains = 4,
  cores = parallel::detectCores(),
  iter = 1e4,
  warmup = 3e3
)
save(multi_level_model2, file = "multi_level_model2.RData")

multi_level_model <- brm(
  # model formulae
  bf(EvoRes ~ PoMiTrSD, family = bernoulli()) + 
  bf(PoMiTrSD ~ 0 + factor(MU) * abs(SL-1) * DI, family = hurdle_gamma()),
  # data
  data = EvoResMod_df,
  # settings for run
  chains = 4,
  cores = parallel::detectCores(),
  iter = 1e4,
  warmup = 3e3
)
save(multi_level_model, file = "multi_level_model.RData")

# Summarize the model
summary(multi_level_model)

# Plot the results
plot(multi_level_model)



Model_brm <- multi_level_model

# model diagnostics
# summary(Model_brm) # all coefficients significant
if(any(abs(rhat(Model_brm)-1) > 1.05)){# everything converged
  message("Some rhat values are worryingly large")
}
# plot(Model_brm) # everything converged
# neff_ratio(Model_brm) # looks alright
ppcheckEvoRes_gg <- pp_check(Model_brm, resp = "EvoRes", type = "bars") # looks alright
ppcheckTraits_gg <- pp_check(Model_brm, resp = "PoMiTrSD") # needs modification, outcome is strictly zero-positive

#### Coefficient Plotting
## extract posterior draws
plot_coeffs <- Model_brm |> 
  gather_draws(`^b_.*`, regex = TRUE) |>
  mutate(.panel = ifelse(grepl(.variable, pattern = "EvoRes"), 
                         "Evolutionary Rescue", "Post Minimum Trait Spread")) 

## significance of coefficients
signif_df <- data.frame(coeff = rownames(summary(Model_brm)[["fixed"]]),
                        signif = abs(sign(summary(Model_brm)[["fixed"]]$`l-95% CI`) + 
                                       sign(summary(Model_brm)[["fixed"]]$`u-95% CI`)) == 2)
plot_coeffs$signif <- signif_df$signif[match(gsub(plot_coeffs$.variable, pattern = "b_", replacement = ""), signif_df$coeff)]

# ## define plotting panel names
# plot_coeffs$.panel[
#   which(plot_coeffs$.variable %in% unique(grep(paste(c("MU", "DI"),collapse = "|"), plot_coeffs$.variable, value = TRUE))) 
# ] = "Population"
# plot_coeffs$.panel[
#   which(plot_coeffs$.variable %in% unique(grep(paste(c("AC", "SL", "VA"),collapse = "|"), plot_coeffs$.variable, value = TRUE)))
# ] = "Environment"
# plot_coeffs$.variable <- gsub(plot_coeffs$.variable, pattern = "b_", replacement = "")
# plot_coeffs$.variable[which(plot_coeffs$.variable == "pert.name")] <- "PM"

# ## define order of coefficients in plotting
# varfacs <- sort(unique(plot_coeffs$.variable))
# varfacs <- c("Intercept", "PM", varfacs[!(varfacs %in% c("Intercept", "PM"))])

## plotting
Coeffplots <- plot_coeffs |> 
  ggplot(aes(x = .value, y = factor(.variable), fill = signif)) +
  geom_vline(xintercept = 0) + 
  stat_halfeye(normalize = "xy") + 
  facet_wrap(~ factor(.panel),
             scales = "free",
             ncol = 2) +
  theme_bw() + theme(legend.position = "none") + 
  labs(x = "Coefficient ", y = "Coefficient Estimate")
Coeffplots

#### Combined Plots
cowplot::plot_grid(Coeffplots, cowplot::plot_grid(ppcheckEvoRes_gg, ppcheckTraits_gg, ncol = 2), ncol = 1)









































### Success Regression ----------------------------

















