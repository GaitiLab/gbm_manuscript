# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
require(GBMutils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Run BANKSY",
    )
    parser$add_argument("--lambda", type = "numeric", help = "Lambda parameter (0-1), lower = cell-typing, higher = spatial domains", default = 0.3)
    parser$add_argument("--input_file", type = "character", help = "Seurat object")
    parser$add_argument("--k_geom", type = "numeric", help = "Number of neighbors used in BANKSY default: 10", default = 10)
    parser$add_argument("--cluster_res", type = "numeric", help = "Resolution used for Leiden clustering", default = 0.5)
    parser$add_argument("--features_to_use", type = "character", help = "all or variable used in BANKSY", default = "all")
    parser$add_argument("--sample_id", type = "character", default = "", help = "Sample ID")
    parser$add_argument("--assay", type = "character", default = "RNA", help = "Sample ID")
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
pacman::p_load(Banksy, Seurat, SeuratWrappers)

log_info("Load seurat object (Xenium)...")
xenium_obj <- readRDS(params$input_file)
DefaultAssay(xenium_obj) <- params$assay
xenium_obj[["BANKSY"]] <- NULL


log_info("Run BANKSY...")
xenium_obj <- RunBanksy(xenium_obj,
    lambda = params$lambda, verbose = TRUE,
    assay = params$assay, slot = "data",
    features = params$features_to_use,
    k_geom = params$k_geom, dimx = "x_centroid", dimy = "y_centroid", spatial_mode = "kNN_r"
)

# Run PCA and UMAP
log_info("Determine number of optimal PCs...")
xenium_obj <- RunPCA(xenium_obj, assay = "BANKSY", features = rownames(xenium_obj), npcs = 50)

# TODO replace later with `elbow_method_wrapper`
npcs <- min(get_pcs(xenium_obj))
p <- ElbowPlot(xenium_obj, ndims = 50) + geom_vline(xintercept = npcs, lty = "dashed") +
    annotate("text", x = npcs + 3.5, y = 6, label = glue("Optimal PCs={npcs}"))
ggsave(plot = p, filename = glue("{params$sample_id}_banksy__elbow_plt.pdf"), height = 5, width = 6, path = params$output_dir)
log_info(glue("Number of optimal PCs: {npcs}"))
xenium_obj <- RunUMAP(xenium_obj, dims = 1:npcs)

# Clustering
log_info("Find Neighbors...")
xenium_obj <- FindNeighbors(xenium_obj, dims = 1:npcs)

log_info("Clustering using Leiden...")
xenium_obj <- FindClusters(xenium_obj,
    resolution = params$cluster_res, algorithm = 4,
    method = "igraph"
)

log_info("Save Seurat object...")
saveRDS(xenium_obj, file = file.path(params$output_dir, glue("{params$sample_id}.rds")))

log_info("Save metadata...")
saveRDS(xenium_obj@meta.data, file = file.path(params$output_dir, glue("{params$sample_id}__meta.rds")))

log_info("Finished!")
