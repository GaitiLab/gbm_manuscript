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
        description = "Compute metric for n-nearest neighbors",
    )
    parser$add_argument("--meta", type = "character", help = "path to metadata file")
    parser$add_argument("--k_neighbors", type = "numeric", help = "Number of nearest neighbors to obtain", default = NULL)
    parser$add_argument("--radius", type = "numeric", help = "Radius for searching for neighbors", default = NULL)
    parser$add_argument("--approach", type = "character", help = "Approach: 'radius' or 'knn'", default = NULL)
    parser$add_argument("--group_varname", type = "character", help = "Variable to use with cell type labels")
    parser$add_argument("--sample_id", type = "character", help = "Sample ID")
    parser$add_argument("--n_cores", type = "numeric", help = "Number of cores")
    parser$add_argument("--dbscan", type = "character", help = "path to file with dbscan results", default = 1)
    parser$add_argument("--ix", type = "numeric", default = NULL, help = "Set to integer if you want to use shuffling")
    parser$add_argument("--n_iter", type = "numeric", default = 100, help = "Number of iterations, for shuffling (only used when 'ix' > 0) (default = 100)")

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

log_info("Load dbscan results...")
dbscan_res <- readRDS(args$dbscan)

log_info("Load metadata...")
# Only select cell id and labels
meta <- readRDS(args$meta) %>% select(cell_id, !!sym(args$group_varname))

# If missing set to 'Undetermined'
meta[is.na(meta[, args$group_varname]), args$group_varname] <- "Undetermined"

if (!is.null(args$ix)) {
    cell_oi_with_neighbors_df <- do.call(rbind, lapply(seq(args$n_iter),
        function(ix,
                 approach, dbscan_res, group_varname, df_with_labels, n_cores = 1, k_neighbors = NA) {
            message(glue::glue("Iteration: {ix}..."))
            # log_info("Shuffle labels...")
            labels <- meta %>% dplyr::pull(!!dplyr::sym(group_varname))
            # Add shuffled labels
            meta[, paste0("shuffled_, ", group_varname)] <- sample(labels, size = length(labels))
            neighbors_df <- compute_neighbor_metrics(
                approach = approach,
                dbscan_res = dbscan_res,
                # Use shuffled labels to check neighbors
                group_varname = paste0("shuffled_, ", group_varname),
                df_with_labels = meta,
                n_cores = n_cores,
                # k_neighbors is NULL when radius
                k_neighbors = k_neighbors
            ) %>%
                # Iteration
                dplyr::mutate(ix = ix)
        },
        approach = args$approach,
        dbscan_res = dbscan_res,
        group_varname = args$group_varname,
        df_with_labels = meta,
        n_cores = args$n_cores,
        # k_neighbors is NULL when radius
        k_neighbors = args$k_neighbors
    )) %>%
        dplyr::mutate(
            # Actual run
            run = args$ix,
            Sample = args$sample_id
        )
} else {
    cell_oi_with_neighbors_df <- compute_neighbor_metrics(
        approach = args$approach,
        dbscan_res = dbscan_res,
        group_varname = args$group_varname,
        df_with_labels = meta,
        n_cores = args$n_cores,
        # k_neighbors is NULL when radius
        k_neighbors = args$k_neighbors
    )
}

log_info("Save as rds...")
# TODO adapt for cell type pair x interaction
if (args$approach == "knn") {
    saveRDS(cell_oi_with_neighbors_df,
        file =
            ifelse(is.null(args$ix),
                glue("{args$output_dir}/{args$sample_id}__{args$group_varname}__k{args$k_neighbors}nn_mean.rds"),
                glue("{args$output_dir}/{args$sample_id}__{args$group_varname}__k{args$k_neighbors}nn_mean_ix{args$ix}.rds")
            )
    )
} else if (args$approach == "radius") {
    saveRDS(
        object = cell_oi_with_neighbors_df,
        file = ifelse(is.null(args$ix),
            glue("{args$output_dir}/{args$sample_id}__{args$group_varname}__n_within_r{args$radius}_mean.rds"),
            glue("{args$output_dir}/{args$sample_id}__{args$group_varname}__n_within_r{args$radius}_mean_ix{args$ix}.rds")
        )
    )
}
log_info("Finished!")
