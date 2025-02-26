# ---- Code to reproduce Figure 3a ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# Load needed packages
pacman::p_load(
    varhandle,
    ggplot2,
    ggrepel,
    mitch,
    fgsea,
    readr,
    msigdbr,
    escape,
    dittoSeq,
    tidyr,
    dplyr,
    stringr,
    data.table,
    ComplexHeatmap,
    colorRamp2,
    foreach,
    doParallel,
    tidyverse,
    viridis
)


params <- list(
    plot_dir = "output/submission/figures",
    # Gene level
    gene_auc_activator_mtx_path = "multiome_results/12_SCENIC_plus/outs/Plots/all_regions/Filtered_TFs/gene_auc_activator_mtx.csv",
    gene_TFs_to_label_path = "multiome_results/12_SCENIC_plus/outs/Plots/all_regions/Filtered_TFs/gene_TFs_to_label.csv",
    # Region level
    region_auc_activator_mtx_path = "multiome_results/12_SCENIC_plus/outs/Plots/all_regions/Filtered_TFs/region_auc_activator_mtx.csv",
    region_TFs_to_label_path = "multiome_results/12_SCENIC_plus/outs/Plots/all_regions/Filtered_TFs/region_TFs_to_label.csv",

    # Heatmap data rna and atac
    # TODO @Yiyan-YW these seem the only ones that are relevant??
    heatmap_data_rna_path = "multiome_results/12_SCENIC_plus/outs/RNA_heatmap.csv",
    heatmap_data_atac_path = "multiome_results/12_SCENIC_plus/outs/ATAC_heatmap.csv",
    TFs_to_plot_path = "multiome_results/12_SCENIC_plus/outs/TFs_to_plot.csv",
    genelists_csv_path = "misc/SuppTables/Table S2.xlsx"
)


heatmap_data_rna <- fread(
    "multiome_results/12_SCENIC_plus/outs/RNA_heatmap.csv"
) |>
    column_to_rownames("TF")
heatmap_data_atac <- fread(
    "multiome_results/12_SCENIC_plus/outs/ATAC_heatmap.csv"
) |>
    column_to_rownames("TF")
TFs_to_plot <- fread(
    "multiome_results/12_SCENIC_plus/outs/TFs_to_plot.csv",
    header = TRUE
)

GaitiLabUtils::create_dir(params$plot_dir)

OPC_NPC1_markers <- readxl::read_excel(params$genelists_csv_path, skip = 1) %>%
    filter(!is.na(Neftel_OPC), !is.na(Neftel_NPC1)) %>%
    dplyr::select(Neftel_OPC, Neftel_NPC1) %>%
    unlist() %>%
    unname() %>%
    unique()


generate_heatmap_with_labels <- function(
    auc_activator_mtx,
    TFs_to_label,
    output_file,
    row_order = NULL
) {
    # Load required libraries
    pacman::p_load(data.table, ComplexHeatmap, circlize, stringr)

    # Load matrix and labels data
    matrix_data <- fread(auc_activator_mtx)
    labels_data <- fread(TFs_to_label, header = TRUE)

    # Convert matrix data to data frame and set row names
    matrix_data <- as.data.frame(matrix_data)
    row.names(matrix_data) <- matrix_data$V1 # Assuming the first column contains the row names (e.g., gene names)
    matrix_data <- matrix_data[, -1] # Remove the first column after setting it as row names
    matrix_data <- matrix_data[grepl("\\+/\\+", rownames(matrix_data)), ] # Filter for activators only

    # Convert labels data to data frame
    labels_data <- as.data.frame(labels_data)
    colnames(labels_data) <- c("V1", "TF")
    labels_data <- labels_data[!duplicated(labels_data$TF), ]
    row.names(labels_data) <- labels_data$TF
    labels_data <- as.vector(labels_data[, -1])

    # Scale the matrix data (Z-score normalization across rows)
    scaled_matrix <- t(apply(matrix_data, 1, scale))
    colnames(scaled_matrix) <- colnames(matrix_data)

    # Set consistent row order if provided (using partial matching)
    if (!is.null(row_order)) {
        matched_rows <- row.names(scaled_matrix)[sapply(
            row_order,
            function(pattern) {
                match_idx <- grep(pattern, row.names(scaled_matrix))
                if (length(match_idx) > 0) {
                    return(match_idx[1]) # Return the first match if multiple found
                } else {
                    return(NA) # No match found
                }
            }
        )]
        matched_rows <- na.omit(matched_rows) # Remove unmatched rows
        scaled_matrix <- scaled_matrix[matched_rows, ] # Reorder rows
    }

    # Create a vector of labels and their positions
    manual_inv_makers <- c(
        "OLIG2",
        "MEOX2",
        "NKX6-2",
        "ASCL1",
        "SOX4",
        "TCF4",
        "ZEB1",
        "ZEB2",
        "E2F1",
        "ETV1",
        "HOXD3",
        "MEIS1",
        "SREBF2",
        "CPEB1",
        "E2F2",
        "HOXB3"
    )
    search_terms <- paste(
        c(OPC_NPC1_markers, manual_inv_makers),
        collapse = "|"
    )
    label_positions <- grep(search_terms, row.names(scaled_matrix))
    label_names <- str_split(
        row.names(scaled_matrix),
        "_direct_\\+/\\+_",
        simplify = TRUE
    )[, 1]
    label_names <- label_names[grepl(search_terms, label_names)]

    # Color palette for the heatmap
    if (grepl("region_heatmap", output_file)) {
        color_palette <- rocket(n = 256, direction = -1, begin = 0.2, end = 0.8)
    } else if (grepl("gene_heatmap", output_file)) {
        color_palette <- mako(n = 256, direction = -1, begin = 0.2, end = 0.8)
    } else {
        color_palette <- viridis(
            n = 256,
            direction = -1,
            begin = 0.2,
            end = 0.8
        )
    }

    # Define the heatmap annotation with lines pointing to the labels using anno_mark
    ha <- rowAnnotation(
        mark = anno_mark(
            at = label_positions,
            labels = label_names,
            side = "right",
            labels_gp = gpar(fontsize = 8),
            link_gp = gpar(col = "black")
        )
    )

    # Define the column order to preserve the current order
    column_order <- seq_len(ncol(scaled_matrix))

    # Save the heatmap to a PDF
    pdf(file = output_file, width = 4, height = 8)
    heatmap_width <- unit(1, "in") # Adjust the width of the heatmap

    # Create the heatmap
    if (!is.null(row_order)) {
        row_order <- seq_len(nrow(scaled_matrix))
        ht <- Heatmap(
            scaled_matrix,
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8),
            column_order = column_order,
            row_order = row_order,
            right_annotation = ha,
            show_row_names = FALSE,
            col = color_palette,
            width = heatmap_width,
            heatmap_legend_param = list(
                title = NULL,
                at = c(-1, 1),
                labels = c("Low TF activity", "High TF activity")
            )
        )
    } else {
        ht <- Heatmap(
            scaled_matrix,
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8),
            column_order = column_order,
            right_annotation = ha,
            show_row_names = FALSE,
            col = color_palette,
            width = heatmap_width,
            heatmap_legend_param = list(
                title = NULL,
                at = c(-1, 1),
                labels = c("Low TF activity", "High TF activity")
            )
        )
    }
    print(ht)

    # Close the PDF device
    dev.off()

    # Return the actual row order from the heatmap
    return(row.names(scaled_matrix)[row_order(ht)])
}

# Generate the gene heatmap using the same row order
row_order <- generate_heatmap_with_labels(
    auc_activator_mtx = params$gene_auc_activator_mtx_path,
    TFs_to_label = params$gene_TFs_to_label_path,
    output_file = file.path(params$plot_dir, "gene_heatmap_with_labels.pdf")
)
row_order <- str_split(row_order, "_direct_\\+/\\+_", simplify = TRUE)[, 1] # Extract gene names

# Generate the region heatmap and capture the row order
generate_heatmap_with_labels(
    auc_activator_mtx = params$region_auc_activator_mtx_path,
    TFs_to_label = params$TFs_to_label_path,
    output_file = file.path(
        params$plot_dir,
        "region_heatmap_with_labels.pdf"
    ),
    row_order = row_order # Use consistent row order
)
# TODO what is the above section for?? @Yiyan-YW (the generate_heatmap_with_labels); which figure does it match? seems like below is the actual Figure?

# --------------------------------------------- Plotting the heatmap --------------------------------------------- #
heatmap_data_rna <- fread(params$heatmap_data_rna_path) |>
    column_to_rownames("TF")
heatmap_data_atac <- fread(params$heatmap_data_atac_path) |>
    column_to_rownames("TF")
TFs_to_plot <- fread(params$TFs_to_plot, header = TRUE)

TF_names_to_plot <- TFs_to_plot %>%
    filter(Cell_type == "Invasive-high OPC/NPC1") %>%
    pull(TF)


pacman::p_load(ComplexHeatmap, circlize, RColorBrewer)

# Subset RNA and ATAC heatmap data
heatmap_data_rna_inv <- heatmap_data_rna[TF_names_to_plot, , drop = FALSE]
heatmap_data_atac_inv <- heatmap_data_atac[TF_names_to_plot, , drop = FALSE]

# Define discrete values and manual markers
discrete_values <- 0:5

# Extract the blue half and red half of the RdBu palette and create discrete color functions
full_rd_bu <- rev(brewer.pal(11, "RdBu"))
blue_half <- full_rd_bu[6:1] # "#F7F7F7" "#D1E5F0" "#92C5DE" "#4393C3" "#2166AC" "#053061"
red_half <- full_rd_bu[6:11] # "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B" "#67001F"
rna_discrete_col_fun <- structure(
    blue_half,
    names = as.character(discrete_values)
)
atac_discrete_col_fun <- structure(
    red_half,
    names = as.character(discrete_values)
)

# Label rows with manual markers
manual_inv_markers <- c(
    "OLIG2",
    "MEOX2",
    "NKX6-2",
    "ASCL1",
    "SOX4",
    "TCF4",
    "ZEB1",
    "ZEB2",
    "E2F1",
    "ETV1",
    "HOXD3",
    "MEIS1",
    "SREBF2",
    "CPEB1",
    "E2F2",
    "HOXB3"
)
search_terms <- paste(c(OPC_NPC1_markers, manual_inv_markers), collapse = "|")
label_positions <- grep(search_terms, rownames(heatmap_data_rna_inv))
label_names <- grep(search_terms, rownames(heatmap_data_rna_inv), value = TRUE)

# Row annotation with lines pointing to labels
ha <- rowAnnotation(
    mark = anno_mark(
        at = label_positions,
        labels = label_names,
        side = "right",
        labels_gp = gpar(fontsize = 8),
        link_gp = gpar(col = "black")
    )
)

# Fixed width for the heatmap
heatmap_width <- unit(1, "in")

# Save RNA Heatmap (Blue Half)
pdf(
    file.path(params$plot_dir, "Fig3a_RNA_heatmap_inv.pdf"),
    width = 4,
    height = 8
)
Heatmap(
    heatmap_data_rna_inv,
    name = "RNA",
    col = rna_discrete_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    right_annotation = ha,
    row_title = "TF",
    column_title = "CellType",
    width = heatmap_width,
    heatmap_legend_param = list(
        title = "Number of patients",
        at = discrete_values,
        labels = as.character(discrete_values)
    )
)
dev.off()

# Save ATAC Heatmap (Red Half)
pdf(
    file.path(params$plot_dir, "Fig3a_ATAC_heatmap_inv.pdf"),
    width = 4,
    height = 8
)
Heatmap(
    heatmap_data_atac_inv,
    name = "ATAC",
    col = atac_discrete_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    right_annotation = ha,
    row_title = "TF",
    column_title = "CellType",
    width = heatmap_width,
    heatmap_legend_param = list(
        title = "Number of patients",
        at = discrete_values,
        labels = as.character(discrete_values)
    )
)
dev.off()
