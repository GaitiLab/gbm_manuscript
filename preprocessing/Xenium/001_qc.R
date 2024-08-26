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
        description = "Apply Quality Control",
    )
    parser$add_argument("--input_file", type = "character", default = "", help = "Path to seurat object (RDS)")
    parser$add_argument("--sample_id", type = "character", default = "", help = "Sample ID")
    parser$add_argument("--meta", type = "character", default = "", help = "Path to metadata object (RDS)")
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

sink(glue("{args$output_dir}/{args$sample_id}__qc.txt"))
log_info(glue("Sample ID: {args$sample_id}"))

log_info("Load additional libraries...")
pacman::p_load(Seurat)

# ---- FILTERS ---- #
min_counts <- 10

log_info("Set filters...")
log_info(glue("Min counts: {min_counts}"))

log_info("Load metadata...")
meta <- readRDS(args$meta)

log_info("Load Xenium Seurat object...")
seurat_obj <- readRDS(args$input_file)

log_info("Apply filters...")
log_info(glue("Total number of cells: {nrow(meta)}"))
# 1) Filtering counts
filter_counts <- meta %>% filter(
    nCount_Xenium >= min_counts
)
ncells <- filter_counts %>% nrow()
log_info(
    (glue("Total number of cells after filtering COUNTS: {ncells}"))
)

log_info("Subset Seurat object...")
seurat_obj <- subset(seurat_obj, subset = cell_id %in% rownames(filter_counts))

log_info("Save object and metadata...")

sample_id <- ifelse(args$sample_id == "", str_split(basename(args$input_dir), "__", simplify = TRUE)[, 3], args$sample_id)

saveRDS(seurat_obj, glue("{args$output_dir}/{args$sample_id}__qc.rds"))
saveRDS(filter_counts, glue("{args$output_dir}/{args$sample_id}__qc_metadata.rds"))
sink() # returns output to the console

log_info("Finished!")
