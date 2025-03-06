# ---- Code to reproduce Figure 2b ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(dplyr, ggplot2, ggrepel)

# Required inputs
params <- list(
    plot_dir = "output/submission/figures",
    interactions_path = "misc/Table S2.xlsx" # can be downloaded online
)

# Creating directories needed for outputs
GaitiLabUtils::create_dir(params$plot_dir)

# ---- Load data & data wrangling ---- #

interactions <- readxl::read_excel(
    params$interactions_path,
    sheet = "Glutamatergic and Invasive-high",
    skip = 1
)

interactions <- interactions %>%
    filter(Interaction_type != "NA") %>%
    group_by(Interaction_type) %>%
    summarise(n = n())

interactions$rank <- rank(-interactions$n, ties.method = "first")
interactions <- interactions %>% arrange(rank)
interactions$pathway_label <- ifelse(
    interactions$rank <= 4,
    interactions$Interaction_type,
    NA
)

# ---- Create figure ---- #
p <- ggplot(interactions, aes(x = rank, y = n)) +
    geom_point(color = "black", stroke = 1, shape = 21, size = 2) +
    geom_text_repel(
        aes(label = pathway_label),
        nudge_x = 5,
        nudge_y = 2,
        na.rm = TRUE
    ) +
    GBMutils::GBM_theme() +
    labs(
        y = "Number of Interactions Identified Between\n\nGlutamatergic Neuron - Invasive-high OPC/NPC1-like",
        x = "Rank"
    )
ggsave(plot = p, filename = "Fig2b.pdf", path = params$plot_dir)
