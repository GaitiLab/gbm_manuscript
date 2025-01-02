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
        description = "DE analysis",
    )
    parser$add_argument("--sample_id", type = "character", default = NULL, help = "Sample ID (default = NULL)")
    parser$add_argument("--input_file", type = "character", help = "Path to seurat object")
    parser$add_argument("--cluster_name", type = "character", help = "Variable containing the cluster identities (default='seurat_clusters')", default = "seurat_clusters")

    params <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    params <- list()
    params$log_level <- 5
    params$output_dir <- glue("{here::here()}/output/")
    params$cluster_name <- "7__sub.cluster"
    params$input_file <- "output/nuclei_segmentation_approach/v1/v1.1/v1.1.2/501_subclustering/Gaiti_Yiyan__6425_cortex_3_D1.rds"
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
pacman::p_load(Seurat)

log_info("Load Seurat object...")
seurat_obj <- readRDS(params$input_file)
DefaultAssay(seurat_obj) <- "RNA"

Idents(seurat_obj) <- params$cluster_name

log_info("Find all markers...")
markers <- FindAllMarkers(
    subset(seurat_obj, subset = BANKSY_snn_res.2_annot != "Undetermined"),
    only.pos = TRUE, assay = "RNA", test.use = "MAST"
)

saveRDS(markers, file.path(params$output_dir, glue("{params$sample_id}_markers_by_celltype")))

log_info("Finished!")
