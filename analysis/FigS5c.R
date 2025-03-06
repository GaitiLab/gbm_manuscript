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
    # path pointing to the output directory of inferCNV
    input_dir = "GBM_VisiumHD/results/InferCNV",
    plot_dir = "output/submission/figures"
)
GaitiLabUtils::create_dir(params$plot_dir)

# Find all CNV files
all_cnv_files <- list.files(
    params$input_dir,
    pattern = "CNV_signal.csv",
    full.names = TRUE,
    recursive = TRUE
)

# Read all files and combine into a single data frame
all_cnv_data <- lapply(all_cnv_files, fread) %>%
    bind_rows()

# Formatting data for plotting
all_cnv_data <- all_cnv_data %>%
    mutate(Cluster_name = paste0(Sample, "_", Group))

# Add a column to indicate whether to highlight
highlight_groups <- c(
    "6425_B_D1_g4",
    "6425_B_D1_g5",
    "6425_A_A1_g8",
    "6425_A_A1_g9",
    "6425_A_A1_g2",
    "6425_A_A1_g13",
    "6425_A_A1_g17",
    "6425_A_A1_g7",
    "6425_A_A1_g1",
    "6425_A_A1_g4",
    "6425_A_A1_g3",
    "6425_A_A1_g14",
    "6425_A_A1_g5",
    "6425_A_A1_g10",
    "6425_A_A1_g15"
)
all_cnv_data <- all_cnv_data %>%
    mutate(
        Label = case_when(
            Cluster_name %in% highlight_groups ~ "Malignant",
            str_detect(Cluster_name, "Malignant") ~ "Malignant",
            TRUE ~ "Normal"
        )
    ) %>%
    filter(!str_detect(Group, "Undetermined"))

# Order chromosomes for plotting
all_cnv_data$Chromosome <- factor(
    all_cnv_data$Chromosome,
    levels = c("Chr7", "Chr10")
)

# Set gain and loss cutoffs 2 MADs away from the median for each chromosome
chr7gain_cutoff <- median(all_cnv_data$Scaled_CNV[
    all_cnv_data$Chromosome == "Chr7"
]) +
    2 * mad(all_cnv_data$Scaled_CNV[all_cnv_data$Chromosome == "Chr7"])
chr10loss_cutoff <- median(all_cnv_data$Scaled_CNV[
    all_cnv_data$Chromosome == "Chr10"
]) -
    2 * mad(all_cnv_data$Scaled_CNV[all_cnv_data$Chromosome == "Chr10"])

# Plotting scatter plot of chr7 and chr10 with cutoffs
chr7_values <- all_cnv_data$Scaled_CNV[all_cnv_data$Chromosome == "Chr7"]
chr10_values <- all_cnv_data$Scaled_CNV[all_cnv_data$Chromosome == "Chr10"]
groups <- all_cnv_data$Cluster_name[all_cnv_data$Chromosome == "Chr7"]
Sample <- all_cnv_data$Sample[all_cnv_data$Chromosome == "Chr7"]

scatter_df <- data.frame(
    Group = groups,
    Sample = Sample,
    Chr7 = chr7_values,
    Chr10 = chr10_values,
    Label = all_cnv_data$Label[all_cnv_data$Chromosome == "Chr7"]
)

# Calculate axis limits
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
