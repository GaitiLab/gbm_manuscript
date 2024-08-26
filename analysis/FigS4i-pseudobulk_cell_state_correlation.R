# Load seurat object
curr_seurat_data <- readRDS(args$input)
DefaultAssay(curr_seurat_data) <- "RNA"

curr_seurat_data <- subset(curr_seurat_data, subset = CellClass_L1 == "Malignant")
curr_seurat_data <- subset(curr_seurat_data, subset = Confident_Annotation == TRUE)

# For each patient, find the variable genes
Split_object <- SplitObject(curr_seurat_data, split.by = "Patient")
variable_genes_per_patient <- lapply(Split_object, function(patient_seurat) {
  patient_seurat <- FindVariableFeatures(patient_seurat, selection.method = "vst", nfeatures = 2000)
  VariableFeatures(patient_seurat)
})

# Find the consistently highly ranked genes
variable_genes_per_patient_ranked <- lapply(variable_genes_per_patient, function(genes) {
  variable_genes_per_patient_ranked <- data.frame(gene = genes, rank = seq_along(genes))
  variable_genes_per_patient_ranked
})
all_variable_genes_per_patient_ranked <- do.call(rbind, variable_genes_per_patient_ranked)

# Format the data
curr_seurat_data <- NormalizeData(curr_seurat_data) # Normalize the data
curr_seurat_data$cell_type_per_patient <- paste0(curr_seurat_data$Patient, "__", curr_seurat_data$CellClass_L5_2) # Get cell type per patient

number_of_genes = 2000
# Ranked by the number of times a gene is in the top 2000, if tied, ranked by the average rank, take the top 2000
all_variable_genes_per_patient_ranked_top_genes <- all_variable_genes_per_patient_ranked %>%
    group_by(gene) %>%
    summarise(n = n(), avg_rank = mean(rank)) %>%
    arrange(desc(n), avg_rank) %>%
    # Exclude MT genes, ribosomal genes, and non-coding genes
    filter(!grepl("^(MIR|LINC|ENSG|MT|RPL|RPS)|.*\\..*", gene)) %>%
    filter(!grepl("-AS", gene)) %>%
    head(number_of_genes)
genelist <- all_variable_genes_per_patient_ranked_top_genes$gene

# Get the average gene expression per cell type per patient
pseudobulked_cell_type_per_patient_RNA <- AggregateExpression(curr_seurat_data, 
                                                            group.by = "cell_type_per_patient", 
                                                            assay = "RNA",
                                                            return.seurat = TRUE
                                                            )
pseudobulked_cell_type_per_patient_RNA$CellClass_L5_2 <- sub(".*__", "", rownames(pseudobulked_cell_type_per_patient_RNA[[]]))

# Get the average gene expression per cell type
mean_cell_type_RNA <- as.data.frame(AverageExpression(pseudobulked_cell_type_per_patient_RNA,
                                                    group.by = "CellClass_L5_2",
                                                    assay = "RNA", 
                                                    slot = "data",
                                                    features = genelist
                                                    ))
colnames(mean_cell_type_RNA) <- sub("RNA.", "", colnames(mean_cell_type_RNA))

# Format column names
colnames(mean_cell_type_RNA) <- sub("Malignant_AC", "Malignant AC", colnames(mean_cell_type_RNA))
colnames(mean_cell_type_RNA) <- sub("Malignant_MES_AST", "Malignant MES-AST", colnames(mean_cell_type_RNA))
colnames(mean_cell_type_RNA) <- sub("Malignant_MES_INT", "Malignant MES-INT", colnames(mean_cell_type_RNA))
colnames(mean_cell_type_RNA) <- sub("Malignant_MES_HYP", "Malignant MES-HYP", colnames(mean_cell_type_RNA))
colnames(mean_cell_type_RNA) <- sub("Invasive.high.OPC.NPC1", "Invasive-high OPC/NPC1", colnames(mean_cell_type_RNA))
colnames(mean_cell_type_RNA) <- sub("Malignant_NPC1", "Malignant NPC1", colnames(mean_cell_type_RNA))
colnames(mean_cell_type_RNA) <- sub("Malignant_NPC2", "Malignant NPC2", colnames(mean_cell_type_RNA))
colnames(mean_cell_type_RNA) <- sub("Malignant_OPC", "Malignant OPC", colnames(mean_cell_type_RNA))

# Get the correlation matrix
correlations <- list()
# Loop through each column in mean_cell_type_RNA for cell type correlation
for (cell_type in colnames(mean_cell_type_RNA)) {
    # Calculate the correlation of each cell type with all others
    correlation_matrix <- cor(mean_cell_type_RNA[, cell_type], mean_cell_type_RNA, use = "complete.obs", method = "spearman")
    # Store the correlation matrix in the list
    correlations[[cell_type]] <- correlation_matrix
}
correlation_df <- do.call(rbind, correlations)
rownames(correlation_df) <- colnames(mean_cell_type_RNA)
# Save the correlation matrix
write.csv(correlation_df, paste0(args$output_dir, "/correlation_matrix_", curr_param, "_genes.csv"))

# Heatmap of correlation matrix
correlation_df <- correlation_df[grepl("Malignant|Invasive", rownames(correlation_df)), 
                                grepl("Malignant|Invasive", colnames(correlation_df))]
scaled_correlation_df <- (correlation_df - min(correlation_df)) / (max(correlation_df) - min(correlation_df))

# Set the color palette
col_fun <- colorRamp2(c(0, 0.5, 1), c("#2367AC", "white", "#B21F2C")) # blue-white-red

# Plot the heatmap
genes <- paste0("variable_genes_", length(genelist))
log_info("Analysis for ", genes)

pdf_output_path <- paste0(args$output_dir, "/correlation_heatmap_", genes, "_psuedobulked_patient_mean.pdf")
pdf(pdf_output_path, width = 8.5, height = 8)
ht <- Heatmap(scaled_correlation_df, 
            show_row_names = TRUE, 
            show_column_names = TRUE,
            col = col_fun,
            row_dend_reorder = TRUE,
            column_dend_reorder = TRUE,
            column_title = "Cell types", 
            row_title = "Cell types")

draw(ht)
dev.off()
