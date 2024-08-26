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
        description = "Extract gene expression data",
    )
    parser$add_argument("--input_file", type = "character", help = "Input file (Seurat object)")
    parser$add_argument("--assay", type = "character", help = "Assay to use", default = "RNA")
    parser$add_argument("--layer", type = "character", help = "Layer to use", default = "counts")
    parser$add_argument("--sample_id", type = "character", help = "Sample ID", default = NULL)
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

pacman::p_load(Seurat)

log_info("Load input file...")
obj <- readRDS(args$input_file)
DefaultAssay(obj) <- args$assay
print(obj)

if (args$assay == "RNA") {
    obj <- NormalizeData(obj)
    obj <- ScaleData(obj)
}

log_info("Extract sample id from filename...")
sample_id <- ifelse(is.null(args$sample_id), str_split(basename(args$input_file), "__", simplify = TRUE)[1], args$sample_id)

log_info("Save counts/normalized data...")
saveRDS(
    LayerData(obj, assay = args$assay, layer = args$layer),
    file = glue("{args$output_dir}/{sample_id}__{args$assay}_{args$layer}.rds")
)

log_info("Finished!")
