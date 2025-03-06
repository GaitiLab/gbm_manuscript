# ---- Code to reproduce Figure S7a-c ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# Load required packages
pacman::p_load(
    data.table,
    dplyr,
    amethyst,
    tibble,
    tidyr,
    future,
    furrr,
    purrr,
    ComplexHeatmap,
    colorRamp2,
    RColorBrewer,
    stringr
)

# Required inputs
params <- list(
    plot_dir = "output/submission/figures",
    path_to_amethys_obj = "path_to_amethyst_obj"
)

# Creating directories needed for outputs
GaitiLabUtils::create_dir(params$output_dir)
GaitiLabUtils::create_dir(params$plot_dir)

# --------------------------------- functions -------------------------------- #

# Function to cluster CNV profiles with k-means.
cluster_cnv_profiles <- function(cnv_matrix, k = 2) {
    # k-means cannot accept NA values.
    cnv_matrix[is.na(cnv_matrix)] <- 1

    # use two centroids to separate malignant and non-malignant
    clust <- kmeans(cnv_matrix, centers = k)

    clust_assignments <- clust$cluster
    clust_df <- data.frame(clust_assignments = clust_assignments)
    clust_df$clust_assignments <- factor(
        clust_df$clust_assignments,
        levels = c(1, 2)
    )

    return(clust_df)
}

# Function to normalize and log-transform each column
normalize_log_transform <- function(mat) {
    # scale counts
    scaled_matrix <- apply(mat, 2, function(x) {
        scaled_x <- (x / sum(x, na.rm = TRUE)) * 10000
        return(scaled_x)
    })

    log_matrix <- log2(1 + scaled_matrix)

    return(log_matrix)
}

#' CNV inference using amethyst object
cnvInference <- function(
    obj,
    reference_cells,
    query_cells,
    step_size = 25 * 1e6,
    num_threads = 4,
    plot_dir,
    nmin = 3,
    k = 2,
    cluster_cells = FALSE,
    save = FALSE,
    matrix_name = "cnv_matrix") {
    # Get CpG coverage matrix
    obj@genomeMatrices[["cov"]] <- makeWindows_cov(
        obj,
        type = "CG",
        genes = NULL,
        promoter = FALSE,
        stepsize = step_size,
        bed = NULL,
        index = "chr_cg",
        groupBy = NULL,
        threads = num_threads,
        futureType = "multicore",
        nmin = nmin,
        save = FALSE
    )

    sum_matrix <- obj@genomeMatrices[["cov"]]

    sum_matrix <- as.data.frame(normalize_log_transform(sum_matrix))

    # Select reference and query cells based on log-transformed total coverage
    ref_matrix <- sum_matrix %>%
        select(contains(paste0(reference_cells)))

    query_matrix <- sum_matrix %>%
        select(contains(paste0(query_cells)))

    # Compute log2 ratios between query and reference cells
    log2_ratios <- query_matrix

    for (cell in query_cells) {
        query_col <- query_matrix[[paste0(cell)]]

        # Compute log2 ratio between query cell and reference mean
        log2_ratio <- query_col - rowMeans(ref_matrix, na.rm = TRUE)

        log2_ratios[[paste0("ratio_", cell)]] <- log2_ratio
    }

    # Prepare matrix for heatmap (log2 ratios of query cells)
    ratio_matrix <- log2_ratios %>%
        select(contains("ratio"))

    # Reverse log transform
    ratio_matrix <- 2^ratio_matrix

    chromosomes <- gsub("_.*", "", rownames(ratio_matrix))
    ratio_matrix$chr <- chromosomes

    # Remove sex chromosomes
    ratio_matrix <- ratio_matrix %>%
        filter(!(chr %in% c("chrX", "chrY")))

    # Order chromosomes
    chromosome_numbers <- as.numeric(gsub("chr", "", ratio_matrix$chr))
    ratio_matrix$chr <- chromosome_numbers
    ratio_matrix <- ratio_matrix %>%
        arrange(chr) %>%
        select(-chr)

    windows <- rownames(ratio_matrix)

    # Transpose so rows are cells and columns are genomic windows
    cnv_matrix <- as.data.frame(t(ratio_matrix))
    colnames(cnv_matrix) <- windows

    if (cluster_cells == TRUE) {
        set.seed(111)
        # Cluster CNVs to classify tumor cells
        clust_df <- cluster_cnv_profiles(cnv_matrix, k = k)

        # Export cluster assignments
        write.csv(clust_df, file.path(plot_dir, "cnv_cluster_assignment.csv"))

        row_annot <- rowAnnotation(
            `Cluster` = clust_df$clust_assignments,
            show_annotation_name = FALSE,
            border = TRUE,
            col = list(`Cluster` = c("1" = "black", "2" = "grey"))
        )
    }

    if (save == TRUE) {
        # Export cnv matrix
        write.csv(cnv_matrix, file.path(plot_dir, paste0(matrix_name, ".csv")))
    }

    color.scheme <- rev(brewer.pal(9, "RdBu"))
    colors <- colorRamp2(
        c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2),
        color.scheme
    )

    chromosome_split <- factor(sort(chromosome_numbers), levels = c(seq(1, 22)))
    # Export chromosome split
    saveRDS(chromosome_split, file.path(plot_dir, "chromosome_split.rds"))

    chr_cols <- brewer.pal(22, "Set3")
    chr_cols <- c(chr_cols, chr_cols[1:10])
    names(chr_cols) <- unique(chromosome_split)
    col_annot <- HeatmapAnnotation(
        `Chromosome` = chromosome_split,
        show_annotation_name = FALSE,
        border = TRUE,
        col = list(`Chromosome` = chr_cols),
        height = unit(0.5, "cm")
    )

    if (cluster_cells == TRUE) {
        ht <- Heatmap(
            as.matrix(cnv_matrix),
            column_split = chromosome_split,
            left_annotation = row_annot,
            top_annotation = col_annot,
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            show_row_names = FALSE,
            clustering_method_rows = "ward.D2",
            show_column_names = FALSE,
            col = colors,
            name = "CNV",
            row_dend_width = unit(1, "cm"),
            border_gp = gpar(col = "black", lty = "solid"),
            column_gap = unit(0, "cm"),
            width = unit(30, "cm"),
            height = unit(15, "cm")
        )
    } else {
        ht <- Heatmap(
            as.matrix(cnv_matrix),
            column_split = chromosome_split,
            top_annotation = col_annot,
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            show_row_names = FALSE,
            clustering_method_rows = "ward.D2",
            show_column_names = FALSE,
            col = colors,
            name = "CNV",
            row_dend_width = unit(1, "cm"),
            border_gp = gpar(col = "black", lty = "solid"),
            column_gap = unit(0, "cm"),
            width = unit(30, "cm"),
            height = unit(15, "cm")
        )
    }

    return(ht)
}


# ---------------------------------------------------------------------------- #
#                                  Figure S7a                                  #
# ---------------------------------------------------------------------------- #

obj <- readRDS(params$path_to_amethys_obj)
obj <- subsetObject(obj, rownames(obj@metadata)[obj@metadata$Patient == "6425"])

ref_cells <- obj@metadata %>% filter(CellClass != "Malignant")
ref_cells <- rownames(ref_cells)

query_cells <- setdiff(rownames(obj@metadata), ref_cells)

cnv_results <- cnvInference(
    obj = obj,
    reference_cells = ref_cells,
    query_cells = query_cells,
    step_size = 25e6,
    num_threads = 15,
    plot_dir = params$plot_dir
)

pdf(
    file = file.path(params$plot_dir, "FigS7a_cnv_6425.pdf"),
    width = 15,
    height = 15
)
cnv_results
dev.off()

# ---------------------------------------------------------------------------- #
#                                  Figure S7b                                  #
# ---------------------------------------------------------------------------- #

obj <- readRDS(params$path_to_amethys_obj)
p <- dimFeature(obj, colorBy = CellClass, reduction = "umap") +
    scale_color_manual(
        values = c("#87B5B1", "#E7298A", "#8C6D31", "#CCB883", "#E6AB02")
    )
ggsave(p, filename = "FigS7b_cellclass_umap.pdf", path = params$plot_dir)

# ---------------------------------------------------------------------------- #
#                                  Figure S7c                                  #
# ---------------------------------------------------------------------------- #

cellclass500bwindows <- calcSmoothedWindows(
    obj,
    type = "CG",
    threads = 10,
    step = 500,
    smooth = 3,
    genome = "hg38",
    index = "chr_cg",
    groupBy = "CellClass",
    returnSumMatrix = FALSE,
    returnPctMatrix = TRUE
)

obj@genomeMatrices[["cg_cellclass_tracks"]] <- cellclass500bwindows

obj@genomeMatrices[["cg_cellclass_tracks"]] <- obj@genomeMatrices[[
    "cg_cellclass_tracks"
]][,
    !grepl("_t|_c", colnames(obj@genomeMatrices[["cg_cellclass_tracks"]])),
    with = FALSE
]
pdf(
    file = file.path(params$plot_dir, "FigS7c_markers_heatmap.pdf"),
    width = 15,
    height = 15
)
heatMap(
    obj,
    genes = c(
        "GFAP",
        "AQP4",
        "C1QA",
        "P2RY13",
        "NEUROD6",
        "LINGO1",
        "MOG",
        "PLP1"
    ),
    matrix = "cg_cellclass_tracks",
    nrow = 4,
    legend = FALSE,
    width = 500
)
dev.off()
