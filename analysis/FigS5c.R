# ---- Code to reproduce Figure S5c ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(
    argparse,
    data.table,
    stringr,
    Seurat,
    Signac,
    log4r,
    infercnv,
    dplyr,
    GaitiLabUtils,
    GBMutils,
    ggplot2,
    ggrepel
)

logr <- GaitiLabUtils::init_logging()

# Required inputs
params <- list(
    input = "misc/InferCNV.csv",
    plot_dir = "output/figures"
)
GaitiLabUtils::create_dir(params$plot_dir)

# ---- Load data ---- #
scatter_df <- fread(params$input)

# Calculate axis limits
chr7_values <- scatter_df$Chr7
chr10_values <- scatter_df$Chr10

all_values <- c(chr7_values, chr10_values) - 1
axis_limit <- max(abs(all_values)) * 1.1
axis_max <- 1 + axis_limit
axis_min <- 1 - axis_limit

# Scatter plot
p <- ggplot(scatter_df, aes(x = Chr7, y = Chr10)) +
    # Shade the areas for Chr7 gain and Chr10 loss
    annotate(
        "rect",
        xmin = chr7gain_cutoff,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf,
        fill = "red",
        alpha = 0.1
    ) + # Chr7 gain
    annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = chr10loss_cutoff,
        fill = "blue",
        alpha = 0.1
    ) + # Chr10 loss
    # Solid lines at 1 for each axis
    geom_vline(
        xintercept = 1,
        linetype = "solid",
        color = "black",
        size = 0.5
    ) +
    geom_hline(
        yintercept = 1,
        linetype = "solid",
        color = "black",
        size = 0.5
    ) +
    # Scatter points
    geom_point(aes(color = Label, shape = Sample, alpha = Label), size = 3) +
    # Label malignant groups
    geom_text_repel(
        data = scatter_df[
            scatter_df$Chr7 > chr7gain_cutoff |
                scatter_df$Chr10 < chr10loss_cutoff,
        ],
        aes(label = Group),
        size = 3,
        color = "black",
        max.overlaps = Inf,
        box.padding = 0.3,
        segment.color = "grey50",
        segment.size = 0.5,
        segment.ncp = 3,
        segment.angle = 20,
        nudge_x = 0.02,
        nudge_y = 0.02
    ) +
    # Adjust transparency for normal cells
    scale_alpha_manual(values = c("Malignant" = 1, "Normal" = 0.25)) +
    scale_color_manual(values = c("Malignant" = "red", "Normal" = "black")) +
    scale_size_manual(values = c("Malignant" = 2, "Normal" = 1.5)) +
    # Set even axis limits
    scale_x_continuous(limits = c(axis_min, axis_max)) +
    scale_y_continuous(limits = c(axis_min, axis_max)) +
    labs(
        title = "CNV Signal for Chr7 and Chr10",
        x = "Chr7 Average CNV Signal",
        y = "Chr10 Average CNV Signal"
    ) +
    coord_fixed() + # Fix aspect ratio for even scaling
    GBM_theme() +
    theme(
        axis.text = element_text(size = 12)
    )
ggsave(
    plot = p,
    filename = "FigS5c_VisiumHC_CNV_signal.pdf",
    width = 8,
    height = 6,
    path = params$plot_dir
)
