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
        description = "Create raw seurat object for Xenium outputs",
    )
    parser$add_argument("--input_dir", type = "character", help = "Directory with Xenium outputs")
    parser$add_argument("--sample_id", type = "character", help = "Id of sample used for saving object", default = "")
    args <- parser$parse_args()
} 

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Load additional libraries
log_info("Load additional libraries...")
pacman::p_load(Seurat)

log_info("Create a Xenium object...")
xenium_obj <- LoadXenium(args$input_dir, fov = "fov")

# Add additional metadata, i.e. area + coordinates
log_info("Load cell information...")
cells_info <- data.frame(fread(glue("{args$input_dir}/cells.csv.gz")))
rownames(cells_info) <- cells_info$cell_id
head(cells_info)
log_info("Add additional metadata...")
xenium_obj <- AddMetaData(xenium_obj, metadata = cells_info)

log_info("Determine mito. ratio...")
xenium_obj$mitoRatio <- PercentageFeatureSet(object = xenium_obj, pattern = "^MT-")
xenium_obj$mitoRatio <- xenium_obj@meta.data$mitoRatio / 100

# Compute cell-blank
xenium_obj$pct_blank <- xenium_obj$nCount_BlankCodeword / (xenium_obj$nCount_Xenium + xenium_obj$nCount_BlankCodeword) * 100

log_info("Save object as RDS...")
sample_id <- ifelse(args$sample_id == "", str_split(basename(args$input_dir), "__", simplify = TRUE)[, 3], args$sample_id)
saveRDS(xenium_obj, glue("{args$output_dir}/{sample_id}__raw.rds"))

log_info("Save object...")
saveRDS(xenium_obj@meta.data, glue("{args$output_dir}/{sample_id}__raw_metadata.rds"))

log_info("Finished!")
