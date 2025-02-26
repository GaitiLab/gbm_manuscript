# ---- Code to reproduce Figure S4e ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# Load packages
# TODO @Yiyan-YW add pacman::p_load() here
pacman::p_load(Seurat, data.table, tidyverse, stringr, ggplot2, ggpubr)

params <- list(
    seurat_obj_path = "",
    path_to_de_genes_table = "misc/SuppTables/Table S3.xlsx",
    plot_dir = "output/submission/figures"
)

GaitiLabUtils::create_dir(params$plot_dir)

# Load seurat object
curr_seurat_data <- readRDS(params$seurat_obj_path)
DefaultAssay(curr_seurat_data) <- "RNA"

# Invasive signature
# TODO @Yiyan-YW is below correct?
# inv_sig <- read.csv("misc/data/inv_sig.csv") %>% column_to_rownames("X")
# inv_up_list <- list(inv_sig$inv_up)
# inv_down_list <- list(inv_sig$inv_down)

degs <- readxl::read_excel(
    params$path_to_de_genes_table,
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
    p,
    filename = "FigS4e_invasivity_invasive_score_correlation_malignant.pdf",
    width = 12,
    height = 12,
    units = "cm",
    path = params$plot_dir
)
