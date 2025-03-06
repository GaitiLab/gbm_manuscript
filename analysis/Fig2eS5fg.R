# ---- Code to reproduce Figure 2e, S5f, g ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

pacman::p_load(
    GaitiLabUtils,
    GBMutils,
    glue,
    data.table,
    tidyverse,
    stringr,
    sf,
    rjson,
    Seurat,
    patchwork
)
logr <- init_logging()

# Required inputs
params <- list(
    output_dir = "output",
    plot_dir = "output/figures",
    k_neighbors = 30,
    n_iter = 10000,
    lower_limit = 0.10,
    upper_limit = 0.90,
    # Paths to Seurat objects generated using this manuscript's data, raw data and final metadata can be downloaded online, see publication
    seurat_obj_6425_A_path = "data/visiumhd/processed/6425_A/6425_A.rds",
    seurat_obj_6425_B_path = "data/visiumhd/processed/6425_B/6425_B.rds",
    degs_table_path = "misc/Table S2.xlsx", # can be downloaded online, see publication
    genelists_path = "misc/gene_signatures.xlsx",
    n_cores = 4
)

GaitiLabUtils::create_dir(params$plot_dir)
GaitiLabUtils::create_dir(params$output_dir)

# ---- Setup ---- #
malignant_signatures <- c(
    "Differentiated-like",
    "Progenitor-like",
    "Invasive-high OPC/NPC1"
)
cat_dist_simplified_labels <- setNames(
    c(
        "Nearest 10% malignant cells to neurons",
        "Furthest 10% malignant cells to neurons"
    ),
    c("closest", "furthest")
)

myfuns_nam <- list(close = min, far = max)
signif_level <- 0.05

# Create color palette for distance labels
distance_labels <- setNames(
    c(
        paste0("<", params$lower_limit * 100, "%"),
        paste0(
            params$lower_limit * 100,
            "% <= mean distance <= ",
            params$upper_limit * 100,
            "%"
        ),
        paste0(">", params$upper_limit * 100, "%")
    ),
    c("bottom", "between", "top")
)

# ---- Computing pseudobulk scores ---- #

# Load DEGs to extract invasive signature up
degs <- readxl::read_excel(
    params$params$degs_table_path,
    sheet = "DEGs",
    skip = 1
) %>%
    data.frame()

degs_up <- degs %>%
    filter(Direction == "Upregulated in PT OPC/NPC1-like cells") %>%
    pull(gene)

# Extract Neftel stignatures
gene_lists <- readxl::read_excel(params$genelists_path, skip = 1)
neftel_gene_list <- lapply(
    gene_lists %>% dplyr::select(starts_with("Neftel_")) %>% as.list(),
    function(gene_list) {
        gene_list[!is.na(gene_list)]
    }
)
# Collect malignant signatures
markers <- list(
    Differentiated_like = unique(c(
        neftel_gene_list$Neftel_MES2,
        neftel_gene_list$Neftel_MES1,
        neftel_gene_list$Neftel_AC
    )),
    Progenitor_like = unique(
        c(
            neftel_gene_list$Neftel_OPC,
            neftel_gene_list$Neftel_NP2,
            neftel_gene_list$Neftel_NPC1
        )
    ),
    `Invasive-high OPC/NPC1` = degs_up
)

shuffling_ids <- paste0("cat_dist_sid_", seq_len(params$n_iter))

log_info("Load Seurat objects...")
seurat_obj_6425_A <- readRDS(params$seurat_obj_6425_A_path)
seurat_obj_6425_B <- readRDS(params$seurat_obj_6425_B_path)

log_info("Merge Seurat objects...")
seurat_obj <- merge(
    seurat_obj_6425_A,
    seurat_obj_6425_B,
    add.cell.ids = c("6425_A", "6425_B")
)

log_info("Only subset 'malignant' cells...")
seurat_obj <- subset(seurat_obj, subset = Spatial_CellClass_L1 == "Malignant")

# Loop through observed labels ('cat_dist') and shuffled labels ('cat_dist_sid_')
vec_cols_oi <- c("cat_dist", shuffling_ids)

log_info("Pseudobulk and compute module scores...")
pseudobulk_scores_df_long <- do.call(
    rbind,
    parallel::mclapply(
        vec_cols_oi,
        function(col_oi, obj, markers) {
            pseudobulk_obj <- Seurat::AggregateExpression(
                obj,
                assays = "RNA",
                group.by = col_oi,
                slot = "counts",
                return.seurat = TRUE,
                verbose = FALSE
            )

            pseudobulk_obj <- NormalizeData(pseudobulk_obj, verbose = FALSE)
            pseudobulk_obj <- AddModuleScore(
                pseudobulk_obj,
                features = markers,
                assay = "RNA"
            )

            df <- data.frame(pseudobulk_obj@meta.data)
            clusternames <- colnames(df)[str_starts(colnames(df), "Cluster")]
            colnames(df)[str_starts(colnames(df), "Cluster")] <- names(markers)

            df <- df %>%
                rename(cat_dist = !!sym(col_oi)) %>%
                select(all_of(names(markers)), cat_dist) %>%
                reshape2::melt(
                    id.vars = "cat_dist",
                    value.name = "AddModuleScore"
                ) %>%
                mutate(
                    across(cat_dist, \(x) substr(x, 2, nchar(x))),
                    ix = col_oi
                )
            return(df)
        },
        obj = seurat_obj,
        mc.cores = params$n_cores,
        markers = markers,
        mc.silent = TRUE
    )
) %>%
    rename(signature = variable)

# Saving in case of reuse
saveRDS(
    pseudobulk_scores_df_long,
    file = file.path(
        params$output_dir,
        "6425_combined_pseudobulk_scores_k30.rds"
    )
)

# ---- Data wrangling ---- #
# Split into observed and shuffled
pseudobulk_scores_observed_df_long <- pseudobulk_scores_df_long %>%
    filter(ix == "cat_dist") %>%
    mutate(AddModuleScore = as.numeric(AddModuleScore)) %>%
    rename(AddModuleScore_obs = AddModuleScore) %>%
    select(-ix)

pseudobulk_scores_shuffled_df_long <- pseudobulk_scores_df_long %>%
    filter(ix != "cat_dist") %>%
    mutate(AddModuleScore = as.numeric(AddModuleScore)) %>%
    rename(AddModuleScore_shuffled = AddModuleScore)

# Combine observed and shuffled dataframes
stat_df_long <- pseudobulk_scores_shuffled_df_long %>%
    left_join(pseudobulk_scores_observed_df_long)

# Convert to wide format, compute basic statistics for permutations
stat_df_wide <- stat_df_long %>%
    pivot_wider(names_from = ix, values_from = AddModuleScore_shuffled) %>%
    rowwise() %>%
    mutate(
        perm_min = min(c_across(starts_with("cat_dist_sid"))),
        perm_max = max(c_across(starts_with("cat_dist_sid"))),
        perm_mean = mean(c_across(starts_with("cat_dist_sid"))),
        perm_median = median(c_across(starts_with("cat_dist_sid"))),
        perm_sd = sd(c_across(starts_with("cat_dist_sid"))),
        # Testing hypothesis H0 < H1 (alternative = 'greater')
        count_larger = sum(
            c_across(starts_with("cat_dist_sid_")) > AddModuleScore_obs
        ) +
            1,
        # Testing hypothesis H0 > H1 (alternative = 'less')
        count_smaller = sum(
            c_across(starts_with("cat_dist_sid_")) < AddModuleScore_obs
        ) +
            1,
        z_score = (AddModuleScore_obs - perm_mean) / perm_sd,
        interaction_type = ifelse(
            z_score > 0,
            "Co-localization",
            "Avoidance"
        ),
        p_value = ifelse(
            interaction_type == "Co-localization",
            # Co-localization
            count_larger / params$n_iter,
            count_smaller / params$n_iter
        ),
        p_value_coloc = count_larger / params$n_iter,
        p_value_avoid = count_smaller / params$n_iter,
        cat_dist_simplified = factor(
            ifelse(
                cat_dist == paste0("<", params$lower_limit * 100, "%"),
                cat_dist_simplified_labels[["closest"]],
                cat_dist_simplified_labels[["furthest"]]
            ),
            levels = cat_dist_simplified_labels
        )
    ) %>%
    # Only interested in bottom and top category
    filter(cat_dist != distance_labels[["between"]])

stat_df_wide <- stat_df_wide %>%
    mutate(
        is_significant = p_value < signif_level,
        signature = factor(
            str_replace_all(signature, "_", "-"),
            malignant_signatures
        )
    )

# ---- Create Figure 2e (barplot) ---- #
data_bar <- stat_df_wide %>%
    filter(cat_dist_simplified == cat_dist_simplified_labels[["closest"]]) %>%
    select(signature, z_score)

p_bar <- ggplot(
    data = data_bar,
    aes(
        fill = signature,
        x = signature,
        y = z_score
    )
) +
    geom_bar(
        stat = "identity",
        position = "dodge",
        show.legend = FALSE,
        fill = "grey",
        color = "black"
    ) +
    scale_fill_manual(values = palette) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    GBMutils::GBM_theme() +
    labs(
        title = cat_dist_simplified_labels[["closest"]],
        subtitle = glue(
            "k={params$k_neighbors} with no. iterations={prettyNum(params$n_iter, big.mark =',')}"
        ),
        y = "Deviation of signature score from the mean of 10,000 interactions (Z-score)",
        x = "Signature"
    ) +
    coord_flip()

log_info("Save figure as PDF...")
ggsave(
    p_bar,
    filename = glue(
        "Fig2e-malign_signature_vs_proximity_k{params$k_neighbors}_neurons_barplot_closest_to_neuron.pdf"
    ),
    height = 4,
    width = 10,
    path = params$plot_dir
)

# ---- Create Figure 2e (histogram) ---- #

# Split into shuffled and observed
stat_shuffled_df_subset_long <- stat_df_wide %>%
    filter(
        signature == "Invasive-high OPC/NPC1",
        cat_dist_simplified == cat_dist_simplified_labels[["closest"]]
    ) %>%
    select(cat_dist_simplified, signature, starts_with("cat_dist_sid")) %>%
    reshape2::melt(variable.name = "iteration_id", value.name = "score")

stat_obs_df_subset_wide <- stat_df_wide %>%
    select(
        cat_dist_simplified,
        signature,
        AddModuleScore_obs,
        perm_mean,
        z_score,
        p_value
    ) %>%
    distinct() %>%
    filter(
        signature == "Invasive-high OPC/NPC1",
        cat_dist_simplified == cat_dist_simplified_labels[["closest"]]
    )


p_histo <- ggplot(
    data = stat_shuffled_df_subset_long,
    aes(x = score)
) +
    geom_histogram(binwidth = 1e-3, show.legend = FALSE, fill = "grey") +
    geom_vline(
        data = stat_obs_df_subset_wide,
        aes(xintercept = AddModuleScore_obs),
        color = "blue"
    ) +
    geom_vline(
        data = stat_obs_df_subset_wide,
        aes(xintercept = perm_mean),
        color = "red",
        linetype = "dashed"
    ) +
    geom_text(
        data = stat_obs_df_subset_wide,
        aes(
            x = (AddModuleScore_obs - AddModuleScore_obs * 0.01),
            label = glue(
                "Observed Module score={round(AddModuleScore_obs, 5)}"
            ),
            y = 300
        ),
        colour = "blue",
        angle = 90,
        size = 2.5
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    GBMutils::GBM_theme() +
    labs(
        x = "Module Score",
        y = "No. of iterations",
        subtitle = glue(
            "k={params$k_neighbors} with no. iterations={prettyNum(params$n_iter, big.mark =',')}"
        ),
    ) +
    theme(aspect.ratio = 0.5)

log_info("Save figure as PDF...")
ggsave(
    p_histo,
    filename = glue(
        "Fig2e_histo.pdf"
    ),
    height = 4,
    path = params$plot_dir
)

# ---- Create Figure S5f ---- #

# Split into shuffled and observed
stat_shuffled_df_subset_long <- stat_df_wide %>%
    select(cat_dist_simplified, signature, starts_with("cat_dist_sid")) %>%
    reshape2::melt(variable.name = "iteration_id", value.name = "score") %>%
    filter(
        signature %in% c("Differentiated-like", "Progenitor-like"),
        cat_dist_simplified == "Nearest 10% malignant cells to neurons"
    ) %>%
    mutate(
        signature = factor(
            as.character(signature),
            levels = c("Differentiated-like", "Progenitor-like")
        )
    )

stat_obs_df_subset_wide <- stat_df_wide %>%
    select(
        cat_dist_simplified,
        signature,
        AddModuleScore_obs,
        perm_mean,
        z_score,
        p_value
    ) %>%
    distinct() %>%
    filter(
        signature %in% c("Differentiated-like", "Progenitor-like"),
        cat_dist_simplified == "Nearest 10% malignant cells to neurons"
    )


p_histo <- ggplot(
    data = stat_shuffled_df_subset_long,
    aes(x = score)
) +
    geom_histogram(binwidth = 1e-3, show.legend = FALSE, fill = "grey") +
    ggh4x::facet_grid2(
        signature ~ .,
        scales = "free_x",
        axes = "all",
        independent = "x"
    ) +
    geom_vline(
        data = stat_obs_df_subset_wide,
        aes(xintercept = AddModuleScore_obs),
        color = "blue"
    ) +
    geom_vline(
        data = stat_obs_df_subset_wide,
        aes(xintercept = perm_mean),
        color = "red",
        linetype = "dashed"
    ) +
    geom_text(
        data = stat_obs_df_subset_wide,
        aes(
            x = (AddModuleScore_obs - AddModuleScore_obs * 0.01),
            label = glue(
                "Observed Module score={round(AddModuleScore_obs, 5)}"
            ),
            y = 300
        ),
        colour = "blue",
        angle = 90,
        size = 1.5
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    GBMutils::GBM_theme() +
    labs(
        x = "Module Score",
        y = "No. of iterations",
        subtitle = glue(
            "k={params$k_neighbors} with no. iterations={prettyNum(params$n_iter, big.mark =',')}"
        ),
        title = "Nearest 10% malignant cells to neurons"
    ) +
    theme(
        aspect.ratio = 0.5,
        strip.text.x.top = element_text(size = rel(0.55)),
        strip.text.y.right = element_text(size = rel(0.35))
    )

log_info("Save figure as PDF...")
ggsave(
    p_histo,
    filename = glue(
        "FigS5f_malign_signature_vs_proximity_k{params$k_neighbors}_neurons_histo_DL_and_PL_closest.pdf"
    ),
    height = 8,
    path = params$plot_dir
)


# ---- Create Figure S5g ---- #

# Split into shuffled and observed
stat_shuffled_df_subset_long <- stat_df_wide %>%
    select(cat_dist_simplified, signature, starts_with("cat_dist_sid")) %>%
    reshape2::melt(variable.name = "iteration_id", value.name = "score") %>%
    filter(
        cat_dist_simplified == "Furthest 10% malignant cells to neurons"
    )

stat_obs_df_subset_wide <- stat_df_wide %>%
    select(
        cat_dist_simplified,
        signature,
        AddModuleScore_obs,
        perm_mean,
        z_score,
        p_value
    ) %>%
    distinct() %>%
    filter(
        cat_dist_simplified == "Furthest 10% malignant cells to neurons"
    )


p_histo <- ggplot(
    data = stat_shuffled_df_subset_long,
    aes(x = score)
) +
    geom_histogram(binwidth = 1e-3, show.legend = FALSE, fill = "grey") +
    ggh4x::facet_grid2(
        signature ~ .,
        scales = "free_x",
        axes = "all",
        independent = "x"
    ) +
    geom_vline(
        data = stat_obs_df_subset_wide,
        aes(xintercept = AddModuleScore_obs),
        color = "blue"
    ) +
    geom_vline(
        data = stat_obs_df_subset_wide,
        aes(xintercept = perm_mean),
        color = "red",
        linetype = "dashed"
    ) +
    geom_text(
        data = stat_obs_df_subset_wide,
        aes(
            x = (AddModuleScore_obs - AddModuleScore_obs * 0.01),
            label = glue(
                "Observed Module score={round(AddModuleScore_obs, 5)}"
            ),
            y = 300
        ),
        colour = "blue",
        angle = 90,
        size = 1.5
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    GBMutils::GBM_theme() +
    labs(
        x = "Module Score",
        y = "No. of iterations",
        subtitle = glue(
            "k={params$k_neighbors} with no. iterations={prettyNum(params$n_iter, big.mark =',')}"
        ),
        title = "Furthest 10% malignant cells to neurons"
    ) +
    theme(
        aspect.ratio = 0.5,
        strip.text.x.top = element_text(size = rel(0.55)),
        strip.text.y.right = element_text(size = rel(0.35))
    )

log_info("Save figure as PDF...")
ggsave(
    p_histo,
    filename = glue(
        "FigS5g_malign_signature_vs_proximity_k{params$k_neighbors}_neurons_histo_all_furthest.pdf"
    ),
    height = 8,
    path = params$plot_dir
)

p_bar <- ggplot(
    data = stat_df_wide %>%
        filter(cat_dist_simplified == cat_dist_simplified_labels[["furthest"]]),
    aes(
        fill = signature,
        x = signature,
        y = z_score
    )
) +
    geom_bar(
        stat = "identity",
        position = "dodge",
        show.legend = FALSE,
        fill = "grey",
        color = "black"
    ) +
    scale_fill_manual(values = palette) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    GBMutils::GBM_theme() +
    labs(
        title = cat_dist_simplified_labels[["furthest"]],
        subtitle = glue(
            "k={params$k_neighbors} with no. iterations={prettyNum(params$n_iter, big.mark =',')}"
        ),
        y = "Deviation of signature score from the mean of 10,000 interactions (Z-score)",
        x = "Signature"
    ) +
    coord_flip()

log_info("Save figure as PDF...")
ggsave(
    p_bar,
    filename = glue(
        "FigS5g_malign_signature_vs_proximity_k{params$k_neighbors}_neurons_barplot_furthest_from_neuron.pdf"
    ),
    height = 4,
    width = 10,
    path = params$plot_dir
)
