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
        description = "Assessing relationship between proximity Malignant cells-Neurons and malignature signatures"
    )
    parser$add_argument("--input_dir", type = "character", help = "Directory with Seurat objects")
    parser$add_argument("--k_neighbors", type = "numeric", default = 15, help = "Number of neighbors (default=15)")
    parser$add_argument("--label", type = "character", help = "Variable in dataframe with labels")
    parser$add_argument("--lower_limit", type = "numeric", default = 0.10, help = "Quantile to use as lower limit (default=0.10)")
    parser$add_argument("--upper_limit", type = "numeric", default = 0.90, help = "Quantile to use as lower limit (default=0.90)")
    parser$add_argument("--n_iter", type = "numeric", default = 1000, help = "No. of itereations (default=1000)")
    parser$add_argument("--markers",
        type = "character", help = "Marker lists"
    )
    parser$add_argument("--n_cores",
        type = "numeric",
        default = 1, help = "Number of cores"
    )
    parser$add_argument("--nbin", type = "numeric", default = 24, help = "No. bins in AddModuleScore (default=24)")

    params <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    params <- list()
    params$log_level <- 5
    params$output_dir <- glue("{here::here()}/output/TESTING")
    params$k_neighbors <- 25
    params$label <- "BANKSY_snn_res.2_annot"
    params$lower_limit <- 0.10
    params$upper_limit <- 0.90
    params$n_iter <- 10
    params$n_cores <- 2
    params$markers <- "/cluster/projects/gaitigroup/Users/Joan/002_Project_GBM/000_misc/gene_lists/neftel_signatures.rds"
    params$nbin <- 24
    params$input_dir <- "data/VisiumHD/processed"
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
pacman::p_load(patchwork, Seurat)
coord_cols <- c("x_centroid", "y_centroid")

log_info("Load Neftel markers...")
markers <- readRDS(params$markers)
markers[["Invasive-high OPC/NPC1"]] <- GBMutils::load_invasive_signature("high")

markers <- list(
    Differentiated_like = union(markers$Neftel_MES, markers$Neftel_AC),
    Progenitor_like = union(markers$Neftel_OPC, markers$Neftel_NPC),
    `Invasive-high OPC/NPC1` = GBMutils::load_invasive_signature("high")
)

shuffling_ids <- paste0("cat_dist_sid_", seq_len(params$n_iter))

log_info("Load Seurat objects...")
seurat_obj_PT1 <- readRDS(file.path(params$input_dir, "Gaiti_Yiyan__6425_cortex_1_A1.rds"))
seurat_obj_PT3 <- readRDS(file.path(params$input_dir, "Gaiti_Yiyan__6425_cortex_3_D1.rds"))

log_info("Merge Seurat objects...")
seurat_obj <- merge(seurat_obj_PT1, seurat_obj_PT3, add.cell.ids = c("PT1", "PT3"))

log_info("Extract metadata...")
meta_df <- seurat_obj@meta.data

log_info("Extract Sample IDs...")
meta_df$sample <- str_split(rownames(meta_df), "_ID", simplify = TRUE)[, 1]

sample_IDs <- meta_df %>%
    pull(sample) %>%
    unique()


# For each sample
# 1. Compute pairwise distances between Malignant cells and Neurons
# 2. For each 'Malignant' cell, compute mean distance to k-nearest Neurons
mean_dist_knn_neurons_combi_df <- do.call(
    rbind,
    lapply(sample_IDs, function(current_sample, df, k_neighbors, label) {
        log_info(glue("Current Sample: {current_sample}..."))
        # Select malignant cells
        malign_cells <- df %>%
            filter(
                sample == current_sample,
                str_detect(!!sym(label), "Malignant")
            ) %>%
            select(x_centroid, y_centroid)
        # Select Neurons
        neurons <- df %>%
            filter(sample == current_sample, !!sym(label) == "Neuron") %>%
            select(x_centroid, y_centroid)

        log_info("Compute pair-wise distances between malignant cells and neurons...")
        dist_to_malign <- rdist::cdist(
            malign_cells, neurons,
            metric = "euclidean"
        ) %>%
            as.matrix()
        # Set cell IDs for malignant cells as rownames
        rownames(dist_to_malign) <- rownames(malign_cells)
        # Set cell IDs for Neurons as column names
        colnames(dist_to_malign) <- rownames(neurons)

        log_info("Compute mean distance of k shortest distances...")
        mean_dist_knn_neurons <- apply(dist_to_malign, 1,
            function(dist_vec, k_neighbors) {
                # Sorts distance vector from short - long
                # Use k shortest/smallest distances to compute mean
                return(mean(sort(dist_vec)[1:k_neighbors]))
            },
            k_neighbors = k_neighbors
        )
        log_info("Convert to dataframe...")
        mean_dist_knn_neurons_df <- data.frame(list(
            mean_dist = mean_dist_knn_neurons,
            sample = current_sample
        ))
        return(mean_dist_knn_neurons_df)
    },
    df = meta_df,
    k_neighbors = params$k_neighbors,
    label = params$label
    )
)

log_info("Define bins based on distance...")
# Determine quantiles for binnning
limits <- quantile(
    mean_dist_knn_neurons_combi_df$mean_dist,
    c(params$lower_limit, params$upper_limit)
)

# Define bin names
categories <- c(
    paste0("<", params$lower_limit * 100, "%"),
    paste0(
        params$lower_limit * 100, "% <= mean distance <= ",
        params$upper_limit * 100, "%"
    ),
    paste0(">", params$upper_limit * 100, "%")
)

log_info("Assign bin to malignant cells...")
mean_dist_knn_neurons_combi_df <- mean_dist_knn_neurons_combi_df %>%
    mutate(
        cat_dist = case_when(
            # d < lower limit
            mean_dist < limits[1] ~ categories[1],
            # d > upper limit
            mean_dist > limits[2] ~ categories[3],
            # lower limit <= d <= upper limit
            .default = categories[2],
        )
    )


log_info("Perform shuffling...")
# Shuffle labels N times to compute p-value
shuffled_labels <- do.call(
    cbind,
    lapply(seq_len(params$n_iter),
        function(ix, cat_vec) {
            sample(cat_vec)
        },
        cat_vec = mean_dist_knn_neurons_combi_df$cat_dist
    )
)
# Add labels for each columns
colnames(shuffled_labels) <- shuffling_ids

log_info("Combine observed + shuffled labels into single dataframe...")
cat_dist_labels <- cbind(
    mean_dist_knn_neurons_combi_df,
    data.frame(shuffled_labels)
)

log_info("Update metadata...")
meta_df <- (seurat_obj@meta.data %>% mutate(id = rownames(.))) %>%
    left_join(
        cat_dist_labels %>%
            mutate(id = rownames(.)),
        by = "id"
    )
rownames(meta_df) <- meta_df$id

log_info("Add labels of shuffles to a subset of seurat object (only keeping malignant cells)...")
seurat_obj@meta.data <- meta_df
seurat_obj <- subset(seurat_obj,
    cells = seurat_obj@meta.data %>%
        filter(str_detect(!!sym(params$label), "Malignant")) %>% dplyr::pull(id)
)

log_info("Save seurat object...")
saveRDS(seurat_obj,
    file = file.path(params$output_dir, glue("6425_combined_seurat_obj_meanshortest_k{params$k_neighbors}.rds"))
)

log_info("Save metadata...")
saveRDS(seurat_obj,
    file = file.path(params$output_dir, glue("6425_combined_seurat_obj_meta_meanshortest_k{params$k_neighbors}.rds"))
)
# Loop through observed labels ('cat_dist') and shuffled labels ('cat_dist_sid_')
vec_cols_oi <- c("cat_dist", shuffling_ids)

log_info("Pseudobulk and compute module scores...")
pseudobulk_res_df <- do.call(rbind, parallel::mclapply(vec_cols_oi,
    function(col_oi, obj, markers, nbin = 24) {
        pseudobulk_obj <- Seurat::AggregateExpression(obj,
            assays = "RNA",
            group.by = col_oi, slot = "counts", return.seurat = TRUE, verbose = FALSE
        )

        pseudobulk_obj <- NormalizeData(pseudobulk_obj, verbose = FALSE)
        pseudobulk_obj <- AddModuleScore(pseudobulk_obj, features = markers, nbin = nbin, assay = "RNA")

        df <- data.frame(pseudobulk_obj@meta.data)
        clusternames <- colnames(df)[str_starts(colnames(df), "Cluster")]
        colnames(df)[str_starts(colnames(df), "Cluster")] <- names(markers)

        df <- df %>%
            rename(cat_dist = !!sym(col_oi)) %>%
            select(all_of(names(markers)), cat_dist) %>%
            reshape2::melt(id.vars = "cat_dist", value.name = "AddModuleScore") %>%
            mutate(across(cat_dist, \(x) substr(x, 2, nchar(x))),
                ix = col_oi
            )
        return(df)
    },
    obj = seurat_obj, mc.cores = params$n_cores, markers = markers, mc.silent = TRUE, nbin = params$nbin
))


log_info("Save pseudobulk results...")
saveRDS(pseudobulk_res_df,
    file = file.path(params$output_dir, glue("6425_combined_pseudobulk_scores_k{params$k_neighbors}.rds"))
)

log_info("Finished!")
