pacman::p_load(
  argparse,
  varhandle,
  log4r,
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

parser <- ArgumentParser(description = "Pipeline for SCENIC+ target gene enrichment analysis.")

parser$add_argument("-i", "--input",
  type = "character",
  default = "multiome_results/12_SCENIC_plus_snakemake/outs",
  help = "Input directory containing the eRegulon results."
)
parser$add_argument("-o", "--output_dir",
  type = "character",
  default = "multiome_results/12_SCENIC_plus_snakemake/outs",
  help = "Output directory to save results."
)
args <- parser$parse_args()

print("Parsing arguments...")
print("Input directory: ", args$input)
print("Output directory: ", args$output_dir)

generate_heatmap_with_labels <- function(auc_activator_mtx, TFs_to_label, genelists_csv, output_file, row_order = NULL) {
  # Load required libraries
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(stringr)
  
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
    matched_rows <- row.names(scaled_matrix)[sapply(row_order, function(pattern) {
      match_idx <- grep(pattern, row.names(scaled_matrix))
      if (length(match_idx) > 0) {
        return(match_idx[1]) # Return the first match if multiple found
      } else {
        return(NA) # No match found
      }
    })]
    matched_rows <- na.omit(matched_rows) # Remove unmatched rows
    scaled_matrix <- scaled_matrix[matched_rows, ] # Reorder rows
  }

  # Read Neftel OPC and NPC markers
  neftel_markers <- fread(genelists_csv)
  OPC_NPC1_markers <- unique(na.omit(unlist(neftel_markers[, .(Neftel_OPC, Neftel_NPC1)])))
  OPC_NPC1_markers <- OPC_NPC1_markers[OPC_NPC1_markers != ""] # Remove empty strings
  
  # Create a vector of labels and their positions
  manual_inv_makers <- c("OLIG2", "MEOX2", "NKX6-2", "ASCL1", "SOX4", "TCF4", "ZEB1", "ZEB2", "E2F1", "ETV1", "HOXD3", "MEIS1", "SREBF2", "CPEB1", "E2F2", "HOXB3")
  search_terms <- paste(c(OPC_NPC1_markers, manual_inv_makers), collapse = "|")
  label_positions <- grep(search_terms, row.names(scaled_matrix))
  label_names <- str_split(row.names(scaled_matrix), "_direct_\\+/\\+_", simplify = TRUE)[, 1]
  label_names <- label_names[grepl(search_terms, label_names)]
  
  # Color palette for the heatmap
  if (grepl("region_heatmap", output_file)) {
    color_palette <- rocket(n = 256, direction = -1, begin = 0.2, end = 0.8)
  } else if (grepl("gene_heatmap", output_file)) {
    color_palette <- mako(n = 256, direction = -1, begin = 0.2, end = 0.8)
  } else {
    color_palette <- viridis(n = 256, direction = -1, begin = 0.2, end = 0.8)
  }
  
  # Define the heatmap annotation with lines pointing to the labels using anno_mark
  ha <- rowAnnotation(mark = anno_mark(
    at = label_positions, # Positions of the labels
    labels = label_names, # Label names
    side = "right", # Place the labels on the right
    labels_gp = gpar(fontsize = 8), # Adjust font size of labels
    link_gp = gpar(col = "black") # Customize line appearance
  ))
  
  # Define the column order to preserve the current order
  column_order <- seq_len(ncol(scaled_matrix)) # Current column order
  
  # Save the heatmap to a PDF
  pdf(file = output_file, width = 4, height = 8)
  heatmap_width <- unit(1, "in") # Adjust the width of the heatmap

  # Create the heatmap
  if (!is.null(row_order)) {
    row_order <- seq_len(nrow(scaled_matrix)) # Default row order
    ht <- Heatmap(scaled_matrix,
            row_names_gp = gpar(fontsize = 8), # Adjust font size for row names
            column_names_gp = gpar(fontsize = 8), # Adjust font size for column names
            column_order = column_order, # Preserve the current column order
            row_order = row_order, # Apply row order if provided
            right_annotation = ha, # Add the row annotation with lines pointing to labels
            show_row_names = FALSE, # Hide original row names since we are adding custom labels with lines
            col = color_palette,
            width = heatmap_width,
            heatmap_legend_param = list(
              title = NULL, # No legend title
              at = c(-1, 1), # Positions for labels (arbitrary as numbers won't be shown)
              labels = c("Low TF activity", "High TF activity") # Custom labels
            )
    )
  } else {
    ht <- Heatmap(scaled_matrix,
            row_names_gp = gpar(fontsize = 8), # Adjust font size for row names
            column_names_gp = gpar(fontsize = 8), # Adjust font size for column names
            column_order = column_order, # Preserve the current column order
            right_annotation = ha, # Add the row annotation with lines pointing to labels
            show_row_names = FALSE, # Hide original row names since we are adding custom labels with lines
            col = color_palette,
            width = heatmap_width,
            heatmap_legend_param = list(
              title = NULL, # No legend title
              at = c(-1, 1), # Positions for labels (arbitrary as numbers won't be shown)
              labels = c("Low TF activity", "High TF activity") # Custom labels
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
  auc_activator_mtx = "multiome_results/12_SCENIC_plus_snakemake/outs/Plots/all_regions/Filtered_TFs/gene_auc_activator_mtx.csv",
  TFs_to_label = "multiome_results/12_SCENIC_plus_snakemake/outs/Plots/all_regions/Filtered_TFs/gene_TFs_to_label.csv",
  genelists_csv = "misc/data/genelists.csv",
  output_file = "multiome_results/12_SCENIC_plus_snakemake/outs/Plots/all_regions/Filtered_TFs/gene_heatmap_with_labels.pdf"
)
row_order <- str_split(row_order, "_direct_\\+/\\+_", simplify = TRUE)[, 1] # Extract gene names

# Generate the region heatmap and capture the row order
generate_heatmap_with_labels(
  auc_activator_mtx = "multiome_results/12_SCENIC_plus_snakemake/outs/Plots/all_regions/Filtered_TFs/region_auc_activator_mtx.csv",
  TFs_to_label = "multiome_results/12_SCENIC_plus_snakemake/outs/Plots/all_regions/Filtered_TFs/region_TFs_to_label.csv",
  genelists_csv = "misc/data/genelists.csv",
  output_file = "multiome_results/12_SCENIC_plus_snakemake/outs/Plots/all_regions/Filtered_TFs/region_heatmap_with_labels.pdf",
  row_order = row_order # Use consistent row order
)

# --------------------------------------------- Plotting the heatmap --------------------------------------------- #
heatmap_data_rna <- fread("multiome_results/12_SCENIC_plus_snakemake/outs/RNA_heatmap.csv") |>
  column_to_rownames("TF")
heatmap_data_atac <- fread("multiome_results/12_SCENIC_plus_snakemake/outs/ATAC_heatmap.csv") |>
  column_to_rownames("TF")
TFs_to_plot <- fread("multiome_results/12_SCENIC_plus_snakemake/outs/TFs_to_plot.csv", header = TRUE)

TF_names_to_plot <- TFs_to_plot %>%
  filter(Cell_type == "Invasive-high OPC/NPC1") %>%
  pull(TF)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Subset RNA and ATAC heatmap data
heatmap_data_rna_inv <- heatmap_data_rna[TF_names_to_plot, , drop = FALSE]
heatmap_data_atac_inv <- heatmap_data_atac[TF_names_to_plot, , drop = FALSE]

# Define discrete values and manual markers
discrete_values <- 0:5

# Extract the blue half and red half of the RdBu palette and create discrete color functions
full_rd_bu <- rev(brewer.pal(11, "RdBu"))
blue_half <- full_rd_bu[6:1] # "#F7F7F7" "#D1E5F0" "#92C5DE" "#4393C3" "#2166AC" "#053061"
red_half <- full_rd_bu[6:11] # "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B" "#67001F"
rna_discrete_col_fun <- structure(blue_half, names = as.character(discrete_values))
atac_discrete_col_fun <- structure(red_half, names = as.character(discrete_values))

# Label rows with manual markers
neftel_markers <- fread("misc/data/genelists.csv")
OPC_NPC1_markers <- unique(na.omit(unlist(neftel_markers[, .(Neftel_OPC, Neftel_NPC1)])))
OPC_NPC1_markers <- OPC_NPC1_markers[OPC_NPC1_markers != ""] # Remove empty strings
manual_inv_markers <- c("OLIG2", "MEOX2", "NKX6-2", "ASCL1", "SOX4", "TCF4", 
                        "ZEB1", "ZEB2", "E2F1", "ETV1", "HOXD3", "MEIS1", 
                        "SREBF2", "CPEB1", "E2F2", "HOXB3")
search_terms <- paste(c(OPC_NPC1_markers, manual_inv_markers), collapse = "|")
label_positions <- grep(search_terms, rownames(heatmap_data_rna_inv)) 
label_names <- grep(search_terms, rownames(heatmap_data_rna_inv), value = TRUE)

# Row annotation with lines pointing to labels
ha <- rowAnnotation(mark = anno_mark(
  at = label_positions, 
  labels = label_names,
  side = "right",
  labels_gp = gpar(fontsize = 8),
  link_gp = gpar(col = "black")
))

# Fixed width for the heatmap
heatmap_width <- unit(1, "in")

# Save RNA Heatmap (Blue Half)
pdf(file.path(args$output_dir, "Plots", "RNA_heatmap_inv.pdf"), width = 4, height = 8)
Heatmap(
  heatmap_data_rna_inv,
  name = "RNA",
  col = rna_discrete_col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
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
pdf(file.path(args$output_dir, "Plots", "ATAC_heatmap_inv.pdf"), width = 4, height = 8)
Heatmap(
  heatmap_data_atac_inv,
  name = "ATAC",
  col = atac_discrete_col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
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

