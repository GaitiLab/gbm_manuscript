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
        description = "Pre-processing & Annotation",
    )
    parser$add_argument("--input_file", type = "character", help = "Path to xenium sample directory")
    parser$add_argument("--sample_id", type = "character", help = "Sample ID")
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
pacman::p_load(Seurat, ggtext)

log_info("Load Xenium Data...")
xenium_obj <- readRDS(args$input_file)

log_info("SCT normalization...")
xenium_obj <- SCTransform(xenium_obj, assay = "Xenium")

# Get PCs and make elbow plot
log_info("Perform PCA...")
xenium_obj <- RunPCA(xenium_obj, npcs = 50, features = rownames(xenium_obj))

log_info("Determine number of optimal PCs...")
npcs <- min(get_pcs(xenium_obj))
p <- ElbowPlot(xenium_obj, ndims = 50) + geom_vline(xintercept = npcs, lty = "dashed") +
    annotate("text", x = npcs + 3.5, y = 6, label = glue("Optimal PCs={npcs}"))
ggsave(plot = p, filename = glue("{args$output_dir}/{args$sample_id}_preproc_elbow_plt.pdf"), height = 5, width = 6)

log_info(glue("Number of optimal PCs: {npcs}"))

log_info("Find Neighbors...")
xenium_obj <- FindNeighbors(xenium_obj, reduction = "pca", dims = 1:npcs)

log_info("Run UMAP...")
xenium_obj <- RunUMAP(xenium_obj, dims = 1:npcs)

log_info("Clustering using Leiden...")

# Use Leiden (algorithm = 4)
xenium_obj <- FindClusters(xenium_obj,
    resolution = 0.5,
    algorithm = 4,
    cluster.name = glue("res_leiden_0_5"),
    method = "igraph"
)

log_info("Rename assay from Xenium to RNA...")
xenium_obj <- RenameAssays(xenium_obj, "Xenium" = "RNA")

log_info("Save RDS...")
saveRDS(xenium_obj, glue("{args$output_dir}/{args$sample_id}__preproc.rds"))
saveRDS(xenium_obj@meta.data, glue("{args$output_dir}/{args$sample_id}__preproc_metadata.rds"))

log_info("Finished!")
