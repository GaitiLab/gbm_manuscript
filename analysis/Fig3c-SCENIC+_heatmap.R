library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(circlize)
library(viridis)

# Update paths to your files
region_auc_activator_mtx <- 'multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/Plots/all_regions/Filtered_TFs/region_auc_activator_mtx.csv'
TFs_to_label <- 'multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/Plots/all_regions/Filtered_TFs/TFs_to_label.csv'

# Load your data
matrix_data <- fread(region_auc_activator_mtx)
labels_data <- fread(TFs_to_label)

# Convert matrix data to data frame and set row names
matrix_data <- as.data.frame(matrix_data)
row.names(matrix_data) <- matrix_data$V1  # Assuming the first column contains the row names (e.g., gene names)
matrix_data <- matrix_data[, -1]  # Remove the first column after setting it as row names

# Convert labels data to data frame
labels_data <- as.data.frame(labels_data)
labels_data <- labels_data[-1, ]  # Remove the first row if it contains column names
row.names(labels_data) <- labels_data$V1
labels_data <- labels_data[, -1]

# Scale the matrix data (Z-score normalization across rows)
scaled_matrix <- t(apply(matrix_data, 1, scale))
colnames(scaled_matrix) <- colnames(matrix_data)

# Create a vector of labels and their positions
label_positions <- which(row.names(scaled_matrix) %in% labels_data$V2)
label_names <- labels_data$V2

# Define the heatmap annotation with lines pointing to the labels using anno_link
ha <- rowAnnotation(mark = anno_mark(
    at = label_positions,  # Positions of the labels
    labels = label_names,  # Label names
    side = "right",  # Place the labels on the right
    labels_gp = gpar(fontsize = 8),  # Adjust font size of labels
    link_gp = gpar(col = "black")  # Customize line appearance
))

# Save the heatmap to a PDF
output_file <- '/multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/Plots/all_regions/Filtered_TFs/heatmap_with_labels.pdf'
pdf(file = output_file, width = 4, height = 8)

# Create the heatmap
Heatmap(scaled_matrix, 
        name = "Scaled Data",
        row_names_gp = gpar(fontsize = 8),  # Adjust font size for row names
        column_names_gp = gpar(fontsize = 8),  # Adjust font size for column names
        right_annotation = ha,  # Add the row annotation with lines pointing to labels
        show_row_names = FALSE,  # Hide original row names since we are adding custom labels with lines
        col = rocket(n = 256, direction = -1, begin = 0.2, end = 0.8))  # Use rocket_r color scale

# Close the PDF device
dev.off()
