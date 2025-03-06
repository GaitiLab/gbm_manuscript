# ---- Code to reproduce Figure S4e ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(
    Seurat,
    data.table,
    tidyverse,
    stringr,
    ggplot2,
    ggpubr,
    argparse,
    ggExtra,
    patchwork,
    ggrepel,
    cowplot,
    Signac,
    varhandle,
    log4r,
    scales,
    dplyr,
    stats,
    RColorBrewer,
    GBMutils,
    circlize,
    ComplexHeatmap
)

# Required inputs
params <- list(
    # Path to Seurat object generated using this manuscript's data, raw data and final metadata can be downloaded online, see publication
    seurat_obj_path = "",
    degs_table_path = "misc/Table S2.xlsx", # can be downloaded online, see publication
    plot_dir = "output/figures"
)

GaitiLabUtils::create_dir(params$plot_dir)

# Load seurat object
curr_seurat_data <- readRDS(params$seurat_obj_path)
DefaultAssay(curr_seurat_data) <- "RNA"

# Invasive signature
degs <- readxl::read_excel(
    params$degs_table_path,
    sheet = "DEGs",
    skip = 1
) %>%
    data.frame()
inv_up_list <- degs %>%
    filter(Direction == "Upregulated in PT OPC/NPC1-like cells") %>%
    pull(gene) %>%
    list()
inv_down_list <- degs %>%
    filter(
        Direction == "Upregulated in tumor bulk (TE+TC) OPC/NPC1-like cells"
    ) %>%
    pull(gene) %>%
    list()


# Add invasive signature score to the Seurat object
curr_seurat_data <- AddModuleScore(
    curr_seurat_data,
    features = inv_up_list,
    name = "inv_up",
    nbin = 15
)
curr_seurat_data <- AddModuleScore(
    curr_seurat_data,
    features = inv_down_list,
    name = "inv_down",
    nbin = 15
)
curr_seurat_data$invasive_score <- curr_seurat_data$inv_up1 -
    curr_seurat_data$inv_down1

# Get metadata and calculate mean invasive score
metadata <- curr_seurat_data[[]]
metadata_malignant <- metadata %>%
    select(
        Sample,
        invasivity,
        invasive_score,
        CellClass_L1,
        Region,
        Patient
    ) %>%
    filter(CellClass_L1 == "Malignant") %>%
    mutate(
        Region = case_when(
            Region == "PT" ~ "Peri-tumoral",
            Region == "TE" ~ "Tumor edge",
            Region == "TC" ~ "Tumor core",
        )
    )

color_palette <- c(
    "Peri-tumoral" = "#0173b2",
    "Tumor edge" = "#de8f05",
    "Tumor core" = "#029e73"
)

metadata_malignant <- metadata_malignant %>%
    group_by(Sample, Region) %>%
    summarise(across(c(invasivity, invasive_score), mean, na.rm = TRUE))

p <- ggscatter(
    metadata_malignant,
    x = "invasivity",
    y = "invasive_score",
    xlab = "Mean invasivity signature module score",
    ylab = "Mean Neuronal signature module score",
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "spearman"
) +
    geom_point(aes(color = Region)) +
    scale_color_manual(values = color_palette)
ggsave(
    plot = p,
    filename = "FigS4e_invasivity_invasive_score_correlation_malignant.pdf",
    width = 12,
    height = 12,
    units = "cm",
    path = params$plot_dir
)
