# ---- Code to reproduce Figure 1i and S2o,p ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# Load required packages
pacman::p_load(
    data.table,
    stringr,
    dplyr,
    ggplot2,
    ggpubr,
    GBMutils,
    scales
)

# Set up logging
logr <- GaitiLabUtils::init_logging()

# Required inputs
params <- list(
    output_dir = "output/submission/",
    plot_dir = "output/figures",
    # TODO add where to download
    input_dir = "data/GBM_mIHC"
)

GaitiLabUtils::create_dir(params$output_dir)
GaitiLabUtils::create_dir(params$plot_dir)

# Get all the files in the input directory
input_files <- list.files(
    params$input_dir,
    full.names = TRUE,
    pattern = "\\.txt$",
    recursive = TRUE
)
sample_names <- unique(str_extract(basename(input_files), "^[^_]+"))

# working on a single sample
curr_sample_name <- "All"

### Parsing arguments
log_info("Input: ", params$input_dir)
log_info("Current sample: ", curr_sample_name)
log_info("Output directory: ", params$output_dir)

# ------------------------- Loading data ------------------------- #
log_info("Loading data...")
normal_region <- fread(file.path(
    params$input_dir,
    "GreenConsolidated_data.txt"
))
infiltrative_region <- fread(file.path(
    params$input_dir,
    "BlueConsolidated_data.txt"
))
tumor_region <- fread(file.path(params$input_dir, "RedConsolidated_data.txt"))

# Merge infiltrative_region and tumors together
normal_region$region <- "infiltrative"
infiltrative_region$region <- "infiltrative"
tumor_region$region <- "tumor bulk"
mIHC_df <- rbind(normal_region, infiltrative_region, tumor_region)

# color palette
color_palette <- c("infiltrative" = "#0173b2", "tumor bulk" = "#029e73")

# ------------------------- Labeling cells ------------------------- #
mIHC_df <- mIHC_df |>
    mutate(
        cell_type = case_when(
            str_detect(`Phenotype SOX2`, "\\+") &
                str_detect(`Phenotype ASCL1`, "\\+") ~
                "Invasive-high OPC/NPC1",
            str_detect(`Phenotype SOX2`, "\\+") &
                str_detect(`Phenotype OLIG2`, "\\+") ~
                "Progenitor_like",
            str_detect(`Phenotype SOX2`, "\\+") &
                str_detect(`Phenotype OLIG2`, "\\-") &
                str_detect(`Phenotype ASCL1`, "\\-") ~
                "Differentiated_like",
            str_detect(`Phenotype MAP2`, "\\+") &
                str_detect(`Phenotype SOX2`, "\\-") ~
                "Neuron",
            .default = "Other"
        )
    )
definition <- "SOX2+ & ASCL1+"

# count the number of cells in each category for each sample
summary <- mIHC_df |>
    group_by(cell_type, region) |>
    summarise(count = n())
write.csv(summary, file.path(params$output_dir, "cell_type_summary.csv"))

# ------------------------- Formatting data ------------------------- #
mIHC_df <- mIHC_df |>
    mutate(
        `Slide Name` = str_split(`Sample Name`, " ", simplify = TRUE)[, 1]
    ) |>
    # Exclude M4230251
    filter(`Slide Name` != "M4230251") |>
    filter(`Slide Name` != "M4210268")
mIHC_df$region <- factor(
    mIHC_df$region,
    levels = c("infiltrative", "tumor bulk")
)

# ------------------------- Overall analysis ------------------------- #
# function for df processing
df_aggregation <- function(df, column) {
    # mean and sem
    df <- df |>
        summarise(
            aggregation = mean(.data[[column]], na.rm = TRUE),
            error = sd(.data[[column]], na.rm = TRUE) / sqrt(n())
        )

    return(df)
}

df_filtering <- function(df, column) {
    df <- df |> filter(.data[[column]] > 0)
    return(df)
}

# ------------------------- SOX2+ percentage across regions ------------------------- #
log_info("Checking SOX2+ expression density across regions")
df <- mIHC_df |>
    filter(`Phenotype SOX2` != "") |>
    group_by(`Sample Name`, region, `Slide Name`) |>
    summarise(
        SOX2_pos = sum(`Phenotype SOX2` == "SOX2+") / n(),
        count = sum(`Phenotype SOX2` == "SOX2+"),
        total = n()
    ) |>
    ungroup() |>
    group_by(`Slide Name`, region)

df_FOV <- df_filtering(df, "SOX2_pos")
df_FOV <- df_FOV |>
    mutate(aggregation = SOX2_pos)

summary_data <- df_aggregation(df_FOV, "SOX2_pos")

ggplot(df_FOV, aes(x = region, y = aggregation, color = region)) +
    geom_jitter(width = 0.3, size = 0.8, alpha = 0.7) +
    geom_errorbar(
        data = summary_data,
        aes(x = region, ymin = aggregation - error, ymax = aggregation + error),
        width = 0.1,
        size = 0.6,
        color = "black"
    ) +
    geom_segment(
        data = summary_data,
        aes(
            x = as.numeric(region) - 0.2,
            xend = as.numeric(region) + 0.2,
            y = aggregation,
            yend = aggregation
        ),
        color = "black",
        linewidth = 1
    ) +
    labs(
        title = "Percentage of SOX2+ cells across regions",
        x = "Sample name",
        y = "Percentage of SOX2+ cells"
    ) +
    GBM_theme() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
        legend.position = "none"
    ) +
    scale_color_manual(values = color_palette) +
    scale_y_continuous(
        labels = percent,
        limits = c(
            min(0, min(summary_data$aggregation - summary_data$error)),
            max(df_FOV$aggregation) * 1.15
        )
    ) +
    facet_wrap(~`Slide Name`, nrow = 1) +
    geom_signif(
        comparisons = list(
            c("infiltrative", "tumor bulk")
        ),
        y_position = max(df_FOV$aggregation) * 1.05,
        tip_length = 0.03,
        size = 0.5,
        step_increase = 0.1,
        test = "wilcox.test",
        color = "black"
    )
ggsave(
    filename = "FigS2i.pdf",
    width = 6,
    height = 4,
    path = params$plot_dir
)

# ------------------------- Progenitor-like cells in SOX2+ cells ------------------------- #
log_info("Checking progenitor cell composition...")
df <- mIHC_df |>
    filter(`Phenotype SOX2` == "SOX2+") |>
    group_by(`Sample Name`, region, `Slide Name`) |>
    summarise(
        progenitor_pos = sum(cell_type == "Progenitor_like") / n(),
        count = sum(cell_type == "Progenitor_like"),
        total = n()
    ) |>
    ungroup() |>
    group_by(`Slide Name`, region)
df_FOV <- df_filtering(df, "progenitor_pos")
df_FOV <- df_FOV |>
    mutate(aggregation = progenitor_pos)
summary_data <- df_aggregation(df_FOV, "progenitor_pos")

p <- ggplot(df_FOV, aes(x = region, y = aggregation, color = region)) +
    geom_jitter(width = 0.3, size = 0.8, alpha = 0.7) +
    geom_errorbar(
        data = summary_data,
        aes(x = region, ymin = aggregation - error, ymax = aggregation + error),
        width = 0.1,
        size = 0.6,
        color = "black"
    ) +
    geom_segment(
        data = summary_data,
        aes(
            x = as.numeric(region) - 0.2,
            xend = as.numeric(region) + 0.2,
            y = aggregation,
            yend = aggregation
        ),
        color = "black",
        size = 1
    ) +
    labs(
        title = "Percentage of Progenitor-like cells (OLIG2+) in Malignant cells (SOX2+)",
        x = "Sample name",
        y = "Percentage of Progenitor-like cells"
    ) +
    GBM_theme() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
        legend.position = "none"
    ) +
    scale_color_manual(values = color_palette) +
    scale_y_continuous(
        labels = percent,
        limits = c(
            min(0, min(summary_data$aggregation - summary_data$error)),
            max(df_FOV$aggregation) * 1.15
        )
    ) +
    facet_wrap(~`Slide Name`, nrow = 1) +
    geom_signif(
        comparisons = list(
            c("infiltrative", "tumor bulk")
        ),
        y_position = max(df_FOV$aggregation) * 1.05,
        tip_length = 0.03,
        size = 0.5,
        step_increase = 0.1,
        test = "wilcox.test",
        color = "black"
    )
ggsave(
    plot = p,
    filename = "FigS2j.pdf",
    width = 6,
    height = 4,
    path = params$plot_dir
)

# ------------------------- Inv cell composition ------------------------- #
log_info("Checking for presence of inv high cells in SOX2+ cells...")
df <- mIHC_df |>
    filter(`Phenotype SOX2` == "SOX2+") |>
    group_by(`Sample Name`, region, `Slide Name`) |>
    summarise(
        inv_pos = sum(cell_type == "Invasive-high OPC/NPC1") / n(),
        count = sum(cell_type == "Invasive-high OPC/NPC1"),
        total = n()
    ) |>
    ungroup() |>
    group_by(`Slide Name`, region)
df_filtered <- df_filtering(df, "inv_pos")
df <- df_aggregation(df_filtered, "inv_pos")
df_tumor <- df |>
    filter(region == "tumor bulk") |>
    dplyr::select(`Slide Name`, aggregation)
write.csv(
    df,
    file.path(
        params$output_dir,
        "Fig1i_table.csv"
    )
)

# Prepare data for plotting
df_infiltrative <- df_filtered |>
    filter(region == "infiltrative") |>
    merge(df_tumor, by = "Slide Name") |>
    mutate(fold_change = inv_pos / aggregation) |>
    mutate(log2_fold_change = log2(fold_change))
df_infiltrative_summary <- df_infiltrative |>
    group_by(`Slide Name`, region, aggregation) |>
    summarise(
        mean_log2_fold_change = mean(log2_fold_change),
        sem = sd(log2_fold_change) / sqrt(n()),
        .groups = "drop"
    ) |>
    mutate(
        log2_fold_change_mad_min = mean_log2_fold_change - sem,
        log2_fold_change_mad_max = mean_log2_fold_change + sem,
        log2_fold_change = mean_log2_fold_change
    )

# Stats
results <-
    df_filtered |>
    merge(df_tumor, by = "Slide Name") |>
    mutate(fold_change = inv_pos / aggregation) |>
    group_by(`Slide Name`) |>
    summarize(
        p_value = tryCatch(
            {
                wilcox.test(
                    fold_change[region == "infiltrative"],
                    fold_change[region == "tumor bulk"],
                    exact = FALSE
                )$p.value
            },
            error = function(e) NA
        ),
        .groups = "drop"
    ) |>
    mutate(p_adj = p.adjust(p_value, method = "BH")) |>
    mutate(
        p_adj = ifelse(
            p_adj < 0.0001,
            formatC(p_adj, format = "e", digits = 2),
            round(p_adj, 4)
        ),
        p_value = ifelse(
            p_value < 0.0001,
            formatC(p_value, format = "e", digits = 2),
            round(p_value, 4)
        )
    )

# plot log2 fold change
p <- ggplot() +
    geom_bar(
        data = df_infiltrative_summary,
        aes(x = `Slide Name`, y = mean_log2_fold_change),
        stat = "identity",
        width = 0.6,
        fill = "#2B7095"
    ) +
    geom_errorbar(
        data = df_infiltrative_summary,
        aes(
            x = `Slide Name`,
            ymin = log2_fold_change_mad_min,
            ymax = log2_fold_change_mad_max
        ),
        width = 0.2,
        size = 0.6,
        color = "black"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_point(
        data = df_infiltrative,
        aes(x = `Slide Name`, y = log2_fold_change),
        position = position_jitter(width = 0.2),
        size = 1,
        color = "black"
    ) +
    geom_text(
        data = results,
        aes(
            x = `Slide Name`,
            y = max(df_infiltrative$log2_fold_change) * 1.05,
            label = p_value
        )
    ) +
    labs(
        title = "Log2 fold change of invasive-high OPC/NPC1 cells",
        x = "Sample name",
        y = "Log2 fold change"
    ) +
    GBM_theme() +
    theme(
        legend.position = "none" # Hide legend
    ) +
    coord_flip()
ggsave(
    plot = p,
    filename = "Fig1i.pdf",
    width = 6,
    height = 4,
    path = params$plot_dir
)
