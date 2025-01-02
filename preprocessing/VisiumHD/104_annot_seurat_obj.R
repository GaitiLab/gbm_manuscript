# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Annotate Seurat object",
    )
    parser$add_argument("--seurat_obj", type = "character", help = "Path to Seurat object")
    parser$add_argument("--look_up_table", type = "character", help = "Path to lookup table with cluster-celltype mappings, Excel file with at least one required column 'cluster'.")
    parser$add_argument("--cluster_name", type = "character", help = "Name of cluster")
    parser$add_argument("--filename", type = "character", help = "Filename of output Seurat object. (default='object')", default = "object")
    parser$add_argument("--sheet_name", default = "Sheet1", help = "Sheet to use in Excel file (default='Sheet1')")
    params <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    params <- list()
    params$log_level <- 5
    params$output_dir <- glue("{here::here()}/output/")
    params$seurat_obj <- "output/stardist_segmentation/cell/clustering/Gaiti_Yiyan__6425_cortex_1_A1.rds"
    params$look_up_table <- "misc/annot_cell_segm.xlsx"
    params$cluster_name <- "BANKSY_snn_res.2"
}

# Set up logging
logr <- init_logging(log_level = params$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(params$output_dir)

# Load additional libraries
pacman::p_load(Seurat, readxl)

log_info(("Load Seurat object..."))
seurat_obj <- readRDS(params$seurat_obj)

log_info("Load look-up table with cluster-to-cell type mapping...")
look_up_table <- read_excel(params$look_up_table, sheet = params$sheet_name) %>% mutate(cluster = as.character(cluster))

log_info("Map clusters to label...")
metadata <- seurat_obj@meta.data %>%
    # TODO temporary for updating annotation labels
    select(-!!sym(paste0(params$cluster_name, "_annot"))) %>%
    mutate(cluster = as.character(!!sym(params$cluster_name))) %>%
    left_join(look_up_table)

log_info("Add updated metadata dataframe to Seurat object...")
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)

log_info("Save Seurat object...")
saveRDS(seurat_obj, file = file.path(params$output_dir, glue("{params$filename}.rds")))


saveRDS(seurat_obj@meta.data,
    file = file.path(
        params$output_dir,
        glue("{params$filename}__meta.rds")
    )
)

log_info("Finished")
