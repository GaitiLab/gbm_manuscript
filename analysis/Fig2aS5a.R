# ---- Code to reproduce Figure 2a & S5a ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(
    glue,
    data.table,
    tidyverse,
    stringr,
    GaitiLabUtils,
    GBMutils,
    ComplexHeatmap,
    readxl
)

# Required inputs
params <- list(
    plot_dir = "output/figures",
    interactions_path = "misc/Table S2.xlsx", # can be downloaded online, see publication
    condition_varname = "Region",
    group1 = "PT",
    group2 = "TC",
    remove_autocrine = TRUE,
    pval_type = "pval_adj",
    alpha = 0.05
)

# Create output directory
create_dir(params$plot_dir)

# Setup color palettes
color_group1 <- GBMutils::load_color_palette("Region")[[params$group1]]
color_group2 <- GBMutils::load_color_palette("Region")[[params$group2]]

# Setup celltype categories/labels
celltype_levels <- names(GBMutils::load_color_palette("CCI_CellClass_L2_2"))
malign <- celltype_levels[1:3]
tme_oi <- c(
    "Myeloid_Inflammatory",
    "Myeloid_Immunosuppressive",
    "GABAergic",
    "Glutamatergic",
    "Oligodendrocyte",
    "OPC"
)
celltype_levels <- c(malign, tme_oi)
label_levels <- c("Tumor", "TME")

# ---- Load data & Data wrangling ---- #
interactions_df <- readxl::read_excel(
    params$interactions_path,
    sheet = "All_predictions",
    skip = 1
) %>%
    # Only keep significant interactions based on (ajusted) Fisher' combined p-value
    filter(!!sym(params$pval_type) < params$alpha) %>%
    separate(
        source_target,
        c("source", "target"),
        sep = "__",
        remove = FALSE
    ) %>%
    # Remove direction by sorting source-target alphabetically
    rowwise() %>%
    mutate(
        source_target_undirected = paste0(
            sort(c(source, target)),
            collapse = "__"
        )
    ) %>%
    # Remove duplicate interactions, when they are found in both directions only keep 1
    distinct(
        !!sym(params$condition_varname),
        source_target_undirected,
        complex_interaction,
        .keep_all = FALSE
    ) %>%
    separate(source_target_undirected, c("source", "target"), sep = "__") %>%
    arrange(
        !!sym(params$condition_varname),
        source,
        target,
        .by_group = TRUE
    ) %>%
    # Count number of interactions per source-target per region
    group_by(!!sym(params$condition_varname), target, source) %>%
    summarise(n = n()) %>%
    ungroup()

# Convert to wide format based on levels in '{params$condition_varname}'
interactions_df_wide_condition <- interactions_df %>%
    # Replace missing values with 0
    pivot_wider(
        names_from = !!sym(params$condition_varname),
        values_from = n,
        values_fill = 0
    ) %>%
    dplyr::select(
        source,
        target,
        !!sym(params$group1),
        !!sym(params$group2)
    ) %>%
    mutate(diff_n = !!sym(params$group1) - !!sym(params$group2))

# Convert to wide format based on cell type labels with the difference in interactions as values
interactions_df_diff <- interactions_df_wide_condition %>%
    dplyr::select(source, target, diff_n) %>%
    arrange(target) %>%
    pivot_wider(names_from = target, values_from = diff_n, values_fill = 0) %>%
    column_to_rownames("source")

# Convert to matrix
mat <- data.matrix(interactions_df_diff)
colnames(mat) <- str_replace_all(colnames(mat), "_", "-")
rownames(mat) <- str_replace_all(rownames(mat), "_", "-")


# Make sure order of rows/columns are the same.
mat <- mat[
    intersect(str_replace_all(celltype_levels, "_", "-"), rownames(mat)),
    intersect(str_replace_all(celltype_levels, "_", "-"), colnames(mat))
]

# Fill whole matrix (mirror)
mat <- mat + t(mat)
mat[lower.tri(mat)] <- 0

if (params$remove_autocrine) {
    # Remove autocrine interactions (diagonal)
    diag(mat) <- 0
}

# Remove fully empty columns/rows
mat <- mat[, which(colSums(mat) != 0)]
mat <- mat[which(rowSums(mat) != 0), ]

# Setup heatmap annotation
label_helper <- function(name, malign, tme_oi) {
    if (name %in% stringr::str_replace_all(malign, "_", "-")) {
        return("Tumor")
    } else {
        return("TME")
    }
}
split_label_col <- factor(
    lapply(colnames(mat), label_helper, malign = malign, tme_oi = tme_oi) %>%
        unlist(),
    levels = label_levels
)
split_label_row <- factor(
    lapply(rownames(mat), label_helper, malign = malign, tme_oi = tme_oi) %>%
        unlist(),
    levels = label_levels
)

# Determine min/max values for heatmap
legend_max <- plyr::round_any(max(mat), 25, f = ceiling)
legend_min <- plyr::round_any(min(mat), 25, f = floor)

# Setup colors
color_fun <- circlize::colorRamp2(
    c(legend_min, 0, legend_max),
    c(color_group2, "white", color_group1)
)


# ---- Create Figure 2a ---- #
chosen_cell_func <- get_cell_function(
    matrix = mat,
    is_upper_tri = TRUE,
    add_annot = FALSE
)
hm <- create_hm(
    matrix = mat,
    col = color_fun,
    name = glue("{params$group1}-{params$group2}\ninteractions"),
    cell_fun = chosen_cell_func,
    column_title = "",
    row_title = "",
    column_title_rot = 0,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_height = 10,
    cell_width = 10,
    row_split = split_label_row,
    column_split = split_label_col,
    heatmap_legend_param = list(
        at = c(legend_min, legend_max),
        labels = c(
            glue("More in {params$group2}"),
            glue("More in {params$group1}")
        ),
        legend_height = unit(4, "lines")
    )
)

R.devices::suppressGraphics({
    save_hm(
        hm_obj = hm,
        output_file = file.path(params$plot_dir, "Fig2a.pdf")
    )
})

# ---- Create Figure S5a ---- #
chosen_cell_func <- get_cell_function(
    matrix = mat,
    is_upper_tri = TRUE,
    add_annot = TRUE
)

# Plot Heatmap and save.
hm <- create_hm(
    matrix = mat,
    col = color_fun,
    name = glue("{params$group1}-{params$group2}\ninteractions"),
    cell_fun = chosen_cell_func,
    column_title = "",
    row_title = "",
    column_title_rot = 0,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_height = 10,
    cell_width = 10,
    row_split = split_label_row,
    column_split = split_label_col,
    heatmap_legend_param = list(
        at = c(legend_min, legend_max),
        labels = c(
            glue("More in {params$group2}"),
            glue("More in {params$group1}")
        ),
        legend_height = unit(4, "lines")
    )
)

R.devices::suppressGraphics({
    save_hm(
        hm_obj = hm,
        output_file = file.path(params$plot_dir, "FigS5a.pdf")
    )
})
