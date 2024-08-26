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
        description = "Determine neighbors for each cells based on nearest neighbors ('knn') or max distance ('radius')",
    )
    parser$add_argument("--annotated_cells_df", type = "character", help = "path to annotated_cells_df data file (dataframe with at least 'x_centroid', 'y_centroid' and a column with label/annotation)")
    parser$add_argument("--k_neighbors", type = "numeric", help = "Number of nearest neighbors to obtain", default = 30)
    parser$add_argument("--radius", type = "numeric", help = "Radius for searching for neighbors", default = 30)
    parser$add_argument("--approach", type = "character", help = "Approach: 'radius' or 'knn'", default = "knn")
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
library(dbscan)

# Ensure a valid approach is chosen
if (!(args$approach %in% c("knn", "radius"))) {
    stop("Choose valid approach 'knn' or 'radius'.")
}

log_info("Load annotated_cells_df data...")
annotated_cells_df <- readRDS(args$annotated_cells_df)

# Only select required columns for determining neighbors
annotated_cells_df <- annotated_cells_df %>%
    column_to_rownames("cell_id") %>%
    select(x_centroid, y_centroid)

if (args$approach == "knn") {
    log_info("Find k nearest neighbors for each cell")
    nn <- dbscan::kNN(
        x = annotated_cells_df,
        k = args$k_neighbors
    )
    log_info("Save dbscan output as rds...")
    saveRDS(nn, file = glue("{args$output_dir}/{args$sample_id}__k{args$k_neighbors}nn_dbscan.rds"))
} else if (args$approach == "radius") {
    log_info("Find fixed radius nearest neighbors for each cell")
    nn <- dbscan::frNN(
        x = annotated_cells_df,
        eps = args$radius
    )
    log_info("Save dbscan output as rds...")
    saveRDS(nn, file = glue("{args$output_dir}/{args$sample_id}__n_within_r{args$radius}_dbscan.rds"))
} else {
    stop("Choose valid approach 'knn' or 'radius'.")
}
log_info("Finished!")
