# ---- Code to reproduce Figure S5c (bottom; DL-invasive) ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(
    glue,
    data.table,
    tidyverse,
    stringr,
    ggplot2,
    GaitiLabUtils,
    readxl,
    ggplot2,
    ggrastr,
    ggtext
)

# Required inputs
params <- list(
    plot_dir = "output/submission/figures",
    interactions_path = "misc/Table S2.xlsx", # can be downloaded online
    is_directed = TRUE,
    condition_varname = "Region",
    receiver = "Differentiated_like",
    sender = "Myeloid_Immunosuppressive",
    region_oi = "TC"
)

GaitiLabUtils::create_dir(params$plot_dir)

source_target_oi <- paste0(params$sender, "__", params$receiver)
source_target_oi_formatted <- str_replace_all(
    source_target_oi,
    c("/" = "_", " " = "_")
)

# Load detected interactions...")
detected_interactions <- readxl::read_excel(
    params$interactions_path,
    sheet = "All_predictions",
    skip = 1
) %>%
    # Only keep significant interactions based on p-adj, keep pair and region of interest
    filter(
        pval_adj < 0.05,
        source_target == c(paste0(params$sender, "__", params$receiver)),
        Region == params$region_oi
    ) %>%
    mutate(log10p = -log10(pval_adj)) %>%
    arrange(desc(log10p), .by_group = TRUE) %>%
    distinct(complex_interaction, .keep_all = TRUE) %>%
    separate(
        complex_interaction,
        into = c("ligand_complex", "receptor_complex"),
        sep = "__",
        remove = FALSE
    ) %>%
    mutate(
        # Add rank basd on on -log10(p-adj)
        order_id = row_number(),
        # Formatting strings for plot
        complex_interaction = str_replace_all(complex_interaction, "__", "-")
    )


interactions_oi <- c(
    "SPP1-CD44",
    "TGFB1-EGFR",
    "TGFB1-LPP"
)
p <- ggplot(
    data = detected_interactions,
    aes(x = order_id, y = log10p)
) +
    geom_point(stroke = NA) +
    GBMutils::GBM_theme() +
    labs(
        y = "-log10(p-adj)",
        x = "rank by -log10(p-adj)",
        title = str_replace_all(
            glue("{params$sender} - {params$receiver}"),
            "_",
            "-"
        ),
    ) +
    # Label interactions of interest
    ggrepel::geom_text_repel(
        data = detected_interactions %>%
            filter(complex_interaction %in% interactions_oi),
        aes(
            label = complex_interaction
        ),
        min.segment.length = unit(0.1, "lines"),
        show.legend = FALSE
    ) +
    theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.y.right = element_text(angle = 0, hjust = 0)
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    theme(aspect.ratio = 1) +
    coord_equal()

ggsave(
    plot = p,
    filename = glue(
        "FigS5c_rank_plot_{source_target_oi_formatted}_in_{params$region_oi}.pdf"
    ),
    height = 10,
    width = 10,
    path = params$plot_dir
)
