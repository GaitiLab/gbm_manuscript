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
        description = "Perform QC & Normalize",
    )
    parser$add_argument("--input_file", type = "character", default = NULL, help = "Path to Seurat object")
    parser$add_argument("--sample_id", type = "character", default = NULL, help = "Sample ID (default = NULL)")

    parser$add_argument("--assay", type = "character", default = "RNA", help = "Assay to use (default='RNA')")

    # QC filters
    parser$add_argument("--min_counts", type = "numeric", default = 1, help = "Min. number of counts per bin or nucleus (default=1)")
    parser$add_argument("--min_features", type = "numeric", default = 1, help = "Min. number of features per bin or nucleus (default=1)")
    parser$add_argument("--max_ratio_mt", type = "numeric", default = .25, help = "Min. number of features per bin or nucleus (default=.25)")
    parser$add_argument("--min_area", type = "numeric", default = 0, help = "Min. area (default = 0)")
    parser$add_argument("--max_area", type = "numeric", default = NULL, help = "Max area (default = NULL), if NULL then don't apply filter area < max_area")

    params <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    params <- list()
    params$log_level <- 5
    params$output_dir <- glue("{here::here()}/output")
    params$input_file <- "output/nuclei_segmentation_approach/005_clustering/with_elbow/Gaiti_Yiyan__6425_cortex_1_A1.rds"
    params$min_counts <- 1
    params$min_features <- 1
    params$max_ratio_mt <- .25
    # Assay: RNA or RNA.008um
    params$assay <- "RNA"
}

# Set up logging
logr <- init_logging(log_level = params$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

# Load additional libraries
library(Seurat)
library(dplyr)

# QC variables/metrics
basevarnames <- c("nCount", "nFeature")
varnames <- paste0(basevarnames, "_", params$assay)
names(varnames) <- basevarnames

log_info("Load Seurat object...")
object <- readRDS(params$input_file)

log_info("Apply QC filters using user-given thresholds...")
cells_to_keep <- object@meta.data %>%
    filter(
        !!sym(varnames[["nCount"]]) > params$min_counts,
        !!sym(varnames[["nFeature"]]) > params$min_features,
        mitoRatio < params$max_ratio_mt
    ) %>%
    rownames(.)

log_info(glue("Number of cells to keep: {length(cells_to_keep)}/{nrow(object@meta.data)}"))
object <- subset(object, cells = cells_to_keep)
DefaultAssay(object) <- params$assay

log_info("Normalize: 1. log-normalization & SCTransform")
object <- normalizedata_wrapper(seurat_object = object, apply_sct_transform = TRUE)
DefaultAssay(object) <- "RNA"

log_info("Create output directory...")
output_dir <- file.path(params$output_dir)
create_dir(output_dir)

log_info("Save Seurat object...")
saveRDS(object, file = file.path(output_dir, glue("{params$sample_id}.rds")))
log_info("Save metadata separately...")
saveRDS(object@meta.data, file = file.path(output_dir, glue("{params$sample_id}__meta.rds")))
log_info("Finished!")
