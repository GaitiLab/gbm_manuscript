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

    params <- parser$parse_args()
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

log_info("Load input file...")
obj <- readRDS(params$input_file)

DefaultAssay(obj) <- params$assay
Idents(obj) <- params$cluster_label

log_info("Do differential testing...")
de_markers <- FindAllMarkers(
    subset(obj, subset = cell_type != "Undetermined"),
    assay = params$assay,
    test.use = "MAST", only.pos = TRUE
)

log_info("Save markers...")
saveRDS(de_markers, file = glue("{params$output_dir}/{params$sample_id}__DE_markers.rds"))

log_info("Finished!")
