# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Create Seurat object from gene-count matrix & log-normalization",
    )
    parser$add_argument("--sample_id", type = "character", default = NULL, help = "Sample ID (default = NULL)")
    parser$add_argument("--input_dir", type = "character", default = "", help = "Sample directory w/ '/outs/' folder")
    parser$add_argument("--is_spatial", type = "numeric", default = TRUE, help = "Load data with 'Load10X_Spatial' (1) or load Seurat object (0)")
    parser$add_argument("--assay", type = "character", default = "RNA", help = "Assay to use (default='RNA')")
    parser$add_argument("--bin_size", type = "numeric", default = 2, help = "Bin size to use for creating Seurat object when is_spatial = TRUE (using direct spaceranger outputs) (default=8)")
    parser$add_argument("--scale_factors_path", type = "character", help = "Path to scalefactors JSON file")

    params <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    params <- list()
    params$log_level <- 5
    params$output_dir <- glue("{here::here()}/output/")
    params$input_dir <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/gbm_visiumHD/output/nuclei_segmentation_approach/gene_count_mat"
    params$input_dir <- "/Users/joankant/Desktop/gaitigroup/Data/GBM_VisiumHD/20240708_LH00244_0121_B22MN7GLT3_Gaiti_Yiyan_Visium/Gaiti_Yiyan__6425_cortex_1_A1"
    params$sample_id <- "Gaiti_Yiyan__6425_cortex_1_A1"
    params$is_spatial <- 1
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

if (params$is_spatial) {
    log_info("Create Seurat object from Spaceranger outputs...")
    # Loading the data this way will result in an assay called 'RNA.008um' (depending on bin size)
    object <- Load10X_Spatial(data.dir = file.path(params$input_dir, "outs"), bin.size = params$bin_size, assay = "RNA")
    assay_name <- params$assay
    suffix <- str_split(assay_name, "\\.", simplify = TRUE)[2]

    # Extract centroids of bins from object
    coords <- GetTissueCoordinates(object)

    # Add spatial information to metadata
    meta <- object@meta.data %>%
        mutate(cell = rownames(.)) %>%
        left_join(coords, by = "cell") %>%
        rename(x_centroid = x, y_centroid = y)
} else {
    log_info("Load gene-count matrix...")
    counts <- Read10X_h5(file.path(params$input_dir, glue("{params$sample_id}__filtered_feature_bc_matrix_w_nuclei_seg.h5")))

    log_info("Create Seurat object from counts...")
    object <- CreateSeuratObject(counts = counts)

    log_info("Load spatial information (centroids)...")
    cell_nuclei_boundaries <- data.table::fread(file.path(params$input_dir, glue("{params$sample_id}__nucleus_boundaries.csv.gz"))) %>%
        data.frame() %>%
        select(id, x_centroid, y_centroid, area) %>%
        distinct()

    # Add spatial information to metadata
    meta <- object@meta.data %>%
        mutate(id = rownames(.)) %>%
        left_join(cell_nuclei_boundaries, by = "id")
}

log_info("Convert coordinates in px to microns...")
meta <- meta %>%
    rename(x_centroid_in_px = x_centroid, y_centroid_in_px = y_centroid) %>%
    mutate(
        x_centroid = x_centroid_in_px * scalefactors$microns_per_pixel,
        y_centroid = y_centroid_in_px * scalefactors$microns_per_pixel
    )

log_info("Add coordinates of centroids to metadata...")
object <- AddMetaData(object, meta)

log_info("Compute mito Ratio...")
object$mitoRatio <- PercentageFeatureSet(object = object, pattern = "^MT-")
object$mitoRatio <- object@meta.data$mitoRatio / 100

DefaultAssay(object) <- params$assay

log_info("Normalize data using log-normalization...")
object <- normalizedata_wrapper(seurat_object = object)

# Get PCs and make elbow plot
DefaultAssay(object) <- params$assay

log_info("Save Seurat object...")
saveRDS(object, file = file.path(params$output_dir, glue("{params$sample_id}.rds")))

log_info("Save metadata...")
saveRDS(object@meta.data, file = file.path(params$output_dir, glue("{params$sample_id}__meta.rds")))

log_info("Finished!")
