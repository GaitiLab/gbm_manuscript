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
        description = "Combine shuffling results and compute p-values",
    )
    parser$add_argument("--radius", type = "numeric", help = "Radius for searching for neighbors")
    parser$add_argument("--k_neighbors", type = "numeric", help = "Number of nearest neighbors to obtain")
    parser$add_argument("--group_varname", type = "character", help = "Variable to use with cell type labels")
    parser$add_argument("--obs_df", type = "character", help = "Path to dataframe with observed results from dbscan")
    parser$add_argument("--perm_dir", type = "character", help = "Directory with results fro shuffling")
    parser$add_argument("--sample_id", type = "character", help = "Sample ID")
    parser$add_argument("--n_cores", type = "numeric", help = "Number of cores")
    parser$add_argument("--metric", type = "character", help = "Metric to use for testing", default = "mean_frac")
    parser$add_argument("--approach", type = "character", help = "Approach: 'radius' or 'knn'", default = NULL)
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

# TODO adapt for cell type pair x interaction
if (args$approach == "knn") {
    shuffle_file_pattern <- glue("{args$sample_id}__{args$group_varname}__k{args$k_neighbors}nn_mean_ix")
    shuffle_mean_outfile <- glue("{args$output_dir}/{args$sample_id}__{args$group_varname}__k{args$k_neighbors}nn_shuffled_mean.rds")
    # TODO adapt for gene-pairs
    outfile <- glue("{args$output_dir}/{args$sample_id}__{args$group_varname}__k{args$k_neighbors}nn_stats_res.rds")
} else if (args$approach == "radius") {
    shuffle_file_pattern <- glue("{args$sample_id}__{args$group_varname}__n_within_r{args$radius}_mean_ix")
    shuffle_mean_outfile <- glue("{args$output_dir}/{args$sample_id}__{args$group_varname}__n_within_r{args$radius}_shuffled_mean.rds")
    outfile <- glue("{args$output_dir}/{args$sample_id}__{args$group_varname}__n_within_r{args$radius}_stats_res.rds")
}

# Load additional libraries
log_info("Load observed results...")
neighbors_obs_df <- readRDS(args$obs_df)

# TODO uncomment later, just temporary
log_info("Get all files from shuffling...")
files <- list.files(args$perm_dir, full.names = TRUE)

missing_files <- setdiff(seq_len(args$n_iter), as.numeric(str_remove(str_split(files, "ix", simplify = TRUE)[, 2], ".rds")))

log_info("Missing files:")
print(missing_files)

files_mean <- files[str_detect(files, shuffle_file_pattern)]
log_info(glue("Number of files: {length(files_mean)}"))

log_info("Combine all shuffling results...")
neighbors_shuffled_df <- do.call(rbind, parallel::mclapply(files_mean, readRDS, mc.cores = args$n_cores)) %>%
    mutate(Sample = args$sample_id)

log_info("Save shuffled - computed means as rds...")
saveRDS(
    obj = neighbors_shuffled_df,
    file = shuffle_mean_outfile
)

# Get all cell-type pairs
all_celltype_pairs_shuffled <- neighbors_shuffled_df %>%
    pull(source_target) %>%
    unique()
all_celltype_pairs_obs <- neighbors_obs_df %>%
    pull(source_target) %>%
    unique()

all_celltype_pairs <- intersect(all_celltype_pairs_shuffled, all_celltype_pairs_obs)


log_info("Permutation testing...")
res_df <- do.call(
    rbind,
    parallel::mclapply(
        all_celltype_pairs,
        compute_perm_pval,
        shuffled_df = neighbors_shuffled_df,
        obs_df = neighbors_obs_df,
        metric = args$metric,
        mc.cores = args$n_cores
    )
) %>% mutate(Sample = args$sample_id)


log_info("Save as rds...")
saveRDS(res_df, file = outfile)

log_info("Finished!")
