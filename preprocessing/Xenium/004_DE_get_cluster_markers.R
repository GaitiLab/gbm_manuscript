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
        description = "Perform DGE on Xenium",
    )
    parser$add_argument("--input_file", type = "character", help = "Path to xenium sample directory")
    parser$add_argument("--cluster_label", type = "character", help = "Cluster label in Seurat")
    parser$add_argument("--sample_id", type = "character", help = "Sample ID", default = NULL)
    parser$add_argument("--assay", type = "character", help = "Assay to use for DE testing", default = "SCT")

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
pacman::p_load(Seurat)

log_info("Load input file...")
obj <- readRDS(args$input_file)

log_info("Extract sample id from filename...")
if (is.null(args$sample_id)) {
    sample_id <- str_split(basename(args$input_file), "__", simplify = TRUE)[1]
} else {
    sample_id <- args$sample_id
}

DefaultAssay(obj) <- args$assay
Idents(object = obj) <- args$cluster_label

log_info("Do differential testing...")
de_markers <- FindAllMarkers(obj,
    test.use = "MAST", assay = args$assay, only.pos = TRUE
)

log_info("Save markers...")
saveRDS(de_markers, file = glue("{args$output_dir}/{sample_id}__DE_markers.rds"))

log_info("Finished!")
