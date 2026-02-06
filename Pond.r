#' ####################################################################### #
#' PROJECT: [Evolutionary Rescue in Complex Landscapes]
#' CONTENTS:
#'  - Display Item Creation
#'  DEPENDENCIES:
#'  - All previous scripts have been run
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list = ls())

## Directories ------------------------------------------------------------
### Define dicrectories in relation to project directory
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "Data")
Dir.Exports <- file.path(Dir.Base, "Exports")
### Create directories which aren't present yet
Dirs <- c(Dir.Data, Dir.Exports)
CreateDir <- sapply(Dirs, function(x) if (!dir.exists(x)) dir.create(x))

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

# ### NON-CRAN ----
# if ("ggsankey" %in% rownames(installed.packages()) == FALSE) {
#     remotes::install_github("davidsjoberg/ggsankey")
# }
# library(ggsankey)

## Functionality ----------------------------------------------------------
"%nin%" <- Negate("%in%") # a function for negation of %in% function

# DISPLAY ITMES ===========================================================
## Time-Series Plotting ---------------------------------------------------
### Abundance; Conceptual ----
# Data Loading
load(file.path(Dir.Exports, "POPULATION_TimeStep.RData"))
Abund_df <- Data_df
Abund_df <- Abund_df[Abund_df$pert.name > 7, ]
Abund_df <- Abund_df[Abund_df$pert.name < 14, ]
rm(Data_df)

sapply(c(0, 1), FUN = function(Mu) {
    Abundplot_df <- data.frame(
        aggregate(x = n ~ t + AC + DI + SL + VA + pert.name, data = Abund_df[Abund_df$MU == Mu, ], FUN = mean),
        sd = aggregate(x = n ~ t + AC + DI + SL + VA + pert.name, data = Abund_df[Abund_df$MU == Mu, ], FUN = sd)[, 7]
    )

    ConceptTime_gg <- ggplot(
        Abundplot_df,
        aes(x = t, y = n, colour = factor(pert.name), fill = factor(pert.name))
    ) +
        # geom_point(alpha = 0.1) +
        geom_line() +
        geom_ribbon(aes(y = n, ymin = ifelse((n - sd) < 0, 0, n - sd), ymax = n + sd), alpha = .2) +
        facet_grid(SL * VA ~ AC * DI, scales = "free", labeller = label_both) +
        lims(x = c(310, 1000), y = c(0, NA)) +
        theme_bw() +
        labs(y = "Abundance", x = "Time")
    ConceptTime_gg

    ggsave(ConceptTime_gg, file = paste0("Pond_MU", Mu, ".png"), width = 20, height = 15)
})

### Abundance; Individual ----
Mu <- 0
SimParams <- c("AC", "DI", "SL", "VA")
Abundplot_df <- data.frame(
    aggregate(x = n ~ t + AC + DI + SL + VA + pert.name, data = Abund_df[Abund_df$MU == Mu, ], FUN = mean),
    sd = aggregate(x = n ~ t + AC + DI + SL + VA + pert.name, data = Abund_df[Abund_df$MU == Mu, ], FUN = sd)[, 7]
)

for (Param in SimParams) {
    other_params <- setdiff(SimParams, Param)

    # split data by all combinations of the OTHER parameters
    df_list <- split(
        Abundplot_df,
        Abundplot_df[other_params],
        drop = TRUE
    )

    for (i in seq_along(df_list)) {
        df_sub <- df_list[[i]]

        combo_label <- paste(
            paste(other_params,
                sapply(df_sub[1, other_params], as.character),
                sep = " = "
            ),
            collapse = ", "
        )

        ## annotation data: number of simulations per facet
        annot_df <- aggregate(
            rep ~ get(Param) + pert.name,
            data = Abund_df[
                Abund_df$MU == Mu &
                    Reduce(`&`, Map(
                        `==`,
                        Abund_df[other_params],
                        df_sub[1, other_params]
                    )),
            ],
            FUN = function(x) length(unique(x))
        )
        colnames(annot_df) <- c(Param, "pert.name", "n_rep")
        annot_df$x <- 400
        annot_df$label <- paste("n =", annot_df$n_rep)
        base_y <- 500
        y_step <- 50
        annot_df$y <- base_y -
            (as.numeric(factor(annot_df$pert.name)) - 1) * y_step

        p <- ggplot(
            df_sub,
            aes(
                x = t, y = n,
                colour = factor(pert.name),
                fill = factor(pert.name)
            )
        ) +
            geom_line() +
            geom_ribbon(
                aes(
                    ymin = pmax(0, n - sd),
                    ymax = n + sd
                ),
                alpha = 0.2
            ) +
            facet_wrap(
                as.formula(paste("~", Param)),
                labeller = label_both
            ) +
            geom_text(
                data = annot_df,
                aes(x = x, y = y, label = label, colour = factor(pert.name)),
                inherit.aes = FALSE,
                size = 3,
                show.legend = FALSE
            ) +
            lims(x = c(310, 1000), y = c(0, NA)) +
            theme_bw() +
            labs(
                title = paste("Facetted by", Param),
                subtitle = combo_label,
                y = "Abundance",
                x = "Time"
            )

        print(p)

        pb <- ggplot_build(p)
        layout <- pb$layout$layout

        n_cols <- length(unique(layout$COL))
        n_rows <- length(unique(layout$ROW))

        panel_width <- 7
        panel_height <- 5

        plot_width <- n_cols * panel_width
        plot_height <- n_rows * panel_height

        ## save plot
        file_stub <- gsub("[^A-Za-z0-9_=]", "_", combo_label)
        ggsave(
            filename = paste0("Pond_", Param, "_", file_stub, ".png"),
            plot = p,
            width = plot_width,
            height = plot_height,
            dpi = 300,
            limitsize = FALSE
        )
    }
}
