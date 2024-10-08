# ---------------------------------------- Functions ---------------------------------------- #
score_signature <- function(seurat_obj, signature_genes, name) {
  assay <- "RNA"

  # Original
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(signature_genes),
    name = name,
    assay = assay
  )

  # Perturbed
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(signature_genes),
    name = paste0(name, "_perturbed"),
    assay = paste0(assay, "_perturbed")
  )
  return(seurat_obj)
}

validate_and_get_columns <- function(curr_seurat_obj, signature, annotation_column) {
  # Validate that a signature name is provided and is not empty
  if (is.null(signature) || signature == "") {
    stop("A non-empty 'signature' name must be provided.")
  }
  
  # REVIEW: This is a temporary solution to handle the different signatures
  # Find the correct columns based on the signature
  if (signature == "invasive_signature") {
    original_col <- sym("original_inv")
    perturbed_col <- sym("perturbed_inv")
  } else if (signature == "opc_signature") {
    original_col <- sym("original_opc")
    perturbed_col <- sym("perturbed_opc")
  } else if (signature == "synaptic_signaling") {
    original_col <- sym("synaptic_signaling1")
    perturbed_col <- sym("synaptic_signaling_perturbed1")
  } else if (signature == "synapse") {
    original_col <- sym("synapse1")
    perturbed_col <- sym("synapse_perturbed1")
  } else if (signature == "postsynapse") {
    original_col <- sym("postsynapse1")
    perturbed_col <- sym("postsynapse_perturbed1")
  } else {
    stop("Invalid signature name provided.")
  }

  # Validate that the annotation column is provided and is not empty
  if (is.null(annotation_column) || annotation_column == "") {
    warning("No 'annotation_column' provided.")
    annotation_col <- NULL
  } else if (!annotation_column %in% colnames(curr_seurat_obj[[]])) {
    stop("The provided 'annotation_column' does not exist in the Seurat object.")
  } else {
    annotation_col <- sym(annotation_column)
  }

  # y-axis label based on the signature
  y_label_dict <- c("invasive_signature" = "Invasive signature score",
                    "opc_signature" = "(OPC - COP)",
                    "synaptic_signaling" = "Synaptic signaling score",
                    "synapse" = "Synapse score",
                    "postsynapse" = "Postsynapse score")
  y_label <- y_label_dict[signature]

  list(original_col = original_col, perturbed_col = perturbed_col, annotation_col = annotation_col, y_label = y_label)
}

plot_signature_changes_individually <- function(curr_seurat_obj, cell_type, TF_to_perturb, perturbation_dir, color_palette, plot_title = NULL, signature = NULL, annotation_column = NULL, same_y_axis = FALSE) {
  # Validate and get columns
  cols <- validate_and_get_columns(curr_seurat_obj, signature, annotation_column)
  original_col <- cols$original_col
  perturbed_col <- cols$perturbed_col
  annotation_col <- cols$annotation_col
  y_label <- cols$y_label

  if (same_y_axis) {
    # Calculate the range for the y-axis
    signature_changes <- curr_seurat_obj[[]] %>% 
      mutate(Cell_id = rownames(.)) %>%
      mutate(Sample = sub(".*___", "", Cell_id)) %>%
      group_by(Sample, !!annotation_col) %>%
      summarize(
        original = median(!!original_col),
        perturbed = median(!!perturbed_col),
        .groups = 'drop'  # Ensure the result is ungrouped
      )
    y_range <- range(c(signature_changes$original, signature_changes$perturbed))
  } else {
    # Calculate the range for the y-axis
    signature_changes <- curr_seurat_obj[[]] %>% 
      filter(!!annotation_col == cell_type) %>%
      mutate(Cell_id = rownames(.)) %>%
      mutate(Sample = sub(".*___", "", Cell_id)) %>%
      group_by(Sample) %>%
      summarize(
        original = median(!!original_col),
        perturbed = median(!!perturbed_col),
        .groups = 'drop'  # Ensure the result is ungrouped
      )
    y_range <- range(c(signature_changes$original, signature_changes$perturbed))
  }

  ### Mean per sample ###
  # Filter for the specified cell type
  signature_changes <- curr_seurat_obj[[]] %>% 
    filter(!!annotation_col == cell_type) %>%
    mutate(Cell_id = rownames(.)) %>%
    mutate(Sample = sub(".*___", "", Cell_id)) %>%
    group_by(Sample) %>%
    summarize(
      original = median(!!original_col),
      perturbed = median(!!perturbed_col),
      .groups = 'drop'  # Ensure the result is ungrouped
    )
  write.csv(signature_changes, paste0(perturbation_dir, "/", signature, "_genes_change_sample_", gsub("[^[:alnum:]]", "_", cell_type), "_", gsub("[^[:alnum:]]", "_", TF_to_perturb), "_", gsub("[^[:alnum:]]", "_", annotation_column), ".csv"))
  
  # Pivot data for plotting
  signature_changes_long <- signature_changes %>%
    pivot_longer(cols = c(original, perturbed), names_to = "Type", values_to = "Value")

  # Calculate mean changes per sample
  mean_changes <- signature_changes %>%
    summarize(mean_changes = mean(perturbed - original)) %>%
    pull(mean_changes)
  
  # Generate the plot
  p <- ggplot(signature_changes_long, aes(x = Type, y = Value)) +
    geom_boxplot(aes(fill = Type), width = 0.3) +
    geom_jitter(width = 0) +
    geom_line(aes(group = Sample), linetype = "dotted") +
    labs(title = ifelse(is.null(plot_title), paste("Changes in mean invasive signature expression for ", cell_type, " cells"), plot_title),
         x = "",
         y = y_label) +
    ylim(floor(y_range[1] * 100)/100-0.01, ceiling(y_range[2] * 100)/100) +
    GBM_theme() +
    theme(axis.title = element_text(face = "plain")) +
    stat_compare_means(paired = TRUE) +
    scale_fill_manual(values = color_palette) +
    geom_text(aes(x = 1.5, y = floor(y_range[1] * 100)/100-0.01, label = paste("Mean sample median change:", round(mean_changes, 3))), vjust = 1)
  
  # Save the plot with a unique name based on the signature
  file_name <- paste0(perturbation_dir, "/", signature, "_genes_change_sample_", gsub("[^[:alnum:]]", "_", cell_type), "_", gsub("[^[:alnum:]]", "_", TF_to_perturb), "_", gsub("[^[:alnum:]]", "_", annotation_column), ".pdf")
  ggsave(file_name, plot = p, width = 4, height = 5)
}


plot_transition_probabilities_per_sample <- function(curr_seurat_obj, TF_to_perturb, perturbation_dir, annotation_column) {
  # Ensure the annotation column is a valid column in the Seurat object
  if (!annotation_column %in% colnames(curr_seurat_obj[[]])) {
    stop("The provided 'annotation_column' does not exist in the Seurat object.")
  } else {
    annotation_col <- sym(annotation_column)
  }
  
  # Prepare the data with transition information
  all_cell <- curr_seurat_obj[[]] %>% 
    select(original_opc, perturbed_opc, !!annotation_col, Sample) %>%
    mutate(transition = case_when(
      original_opc > perturbed_opc ~ "Differentiation",
      original_opc < perturbed_opc ~ "De-differentiation"
    ))

  # Calculate the transition probability for each cell type
  transition_prob <- all_cell %>%
    group_by(!!annotation_col, transition, Sample) %>%
    summarize(
      count = n(),
      .groups = 'drop'
    )

  # Pivot the data to have one row per cell type with columns for Differentiation and De-differentiation counts
  transition_prob_per_cell_type <- transition_prob %>%
    pivot_wider(names_from = transition, values_from = count, values_fill = list(count = 0)) %>%
    mutate(
      diff_prob = Differentiation / (Differentiation + `De-differentiation`),
      dediff_prob = `De-differentiation` / (Differentiation + `De-differentiation`)
    )

  # Transform the data to long format for plotting
  long_data <- transition_prob_per_cell_type %>%
    select(!!annotation_col, Sample, diff_prob, dediff_prob) %>%
    pivot_longer(cols = c(diff_prob, dediff_prob), names_to = "Transition_Type", values_to = "Probability") %>%
    mutate(Transition_Type = case_when(
      Transition_Type == "diff_prob" ~ "Differentiation",
      Transition_Type == "dediff_prob" ~ "De-differentiation"
    )) %>%
    mutate(!!annotation_col := case_when(
      !!annotation_col == "Invasive-high OPC/NPC1" ~ "Invasive-high OPC/NPC1",
      !!annotation_col == "Malignant_OPC" ~ "Malignant OPC",
      !!annotation_col == "Malignant_NPC1" ~ "Malignant NPC1",
      !!annotation_col == "Malignant_NPC2" ~ "Malignant NPC2",
      !!annotation_col == "Malignant_AC" ~ "Malignant AC",
      !!annotation_col == "Malignant_MES_AST" ~ "Malignant MES-AST",
      !!annotation_col == "Malignant_MES_HYP" ~ "Malignant MES-HYP",
      !!annotation_col == "Malignant_MES_INT" ~ "Malignant MES-INT",
      !!annotation_col == "Progenitor_like" ~ "Progenitor-like",
      !!annotation_col == "Differentiated_like" ~ "Differentiated-like",
      TRUE ~ as.character(!!annotation_col) # Retain original value if no match
    ))

  # Calculate mean and standard error of the mean (SEM) for each cell type and transition type
  summary_data <- long_data %>%
    group_by(!!annotation_col, Transition_Type) %>%
    summarize(
      mean_prob = mean(Probability),
      sem_prob = sd(Probability) / sqrt(n()),
      .groups = 'drop'
    )
  write.csv(summary_data, paste0(perturbation_dir, "/transition_probability_per_sample_mean_", TF_to_perturb, "_", annotation_column, ".csv"))
  
  # Rank the cell types based on the differentiation probability
  rank_by_diff_prob <- summary_data %>%
    filter(Transition_Type == "Differentiation") %>%
    arrange(desc(mean_prob)) %>%
    pull(!!annotation_col)
  summary_data[[annotation_column]] <- factor(summary_data[[annotation_column]], levels = rank_by_diff_prob)

  # OPC_dev color palette
  color_palette <- c("Differentiation" = "#F0E9BA", "De-differentiation" = "#ED8B22")

  # Plotting
  p <- ggplot(summary_data, aes(x = !!annotation_col, y = mean_prob, fill = Transition_Type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      geom_errorbar(aes(ymin = mean_prob - sem_prob, ymax = mean_prob + sem_prob), 
                position = position_dodge(width = 0.7), width = 0.1, 
                color = "black", size = 0.5) +
    scale_fill_manual(values = color_palette) +
    labs(title = paste0("Differentiation Probability Across Cell Types for ", TF_to_perturb, " KO"),
         x = "Cell Type",
         y = "Transition Probability") +
    GBM_theme()
   
  # Save the plot
  ggsave(paste0(perturbation_dir, "/transition_probability_per_sample_", TF_to_perturb, "_", annotation_column, ".pdf"), plot = p, 
         width = 1.5 * length(unique(summary_data[[annotation_column]])) + 3, height = 4)
 
}

# ---------------------------------------- Data processing ---------------------------------------- #
# Load the Seurat object
merged_obj <- readRDS(args$merged_obj)
metadata <- merged_obj[[]] %>% 
  as.data.frame() %>%
  filter(CellClass_L1 == "Malignant" & Confident_Annotation == TRUE & Platform == "Multiome") %>%
  mutate(
    barcode = sub(".*__", "", CellID),
    SCENICPLUS_label = paste0(barcode, "___", Sample)
  )
cells_to_keep <- metadata$SCENICPLUS_label

# Load the perturbation data
original_matrix <- fread(paste0(perturbation_dir, "/original_matrix.csv")) %>%
  as.data.frame() %>% 
  column_to_rownames("V1") %>%
  filter(rownames(.) %in% cells_to_keep)

perturbation_matrix <- fread(paste0(perturbation_dir, "/perturbed_matrix_", args$TF_to_perturb, ".csv")) %>%
  as.data.frame() %>%
  column_to_rownames("V1") %>%
  filter(rownames(.) %in% cells_to_keep)

# Load the DE genes
DE_genes <- fread(args$DE_genes, header = TRUE)
Inv_high <- DE_genes$inv_up
Inv_low <- DE_genes %>% 
  filter(inv_down != "") %>%
  pull(inv_down)

# Create a Seurat object for the original and perturbed data
# rotate the matrix to have cells as rows and genes as columns
original_matrix <- t(original_matrix)
perturbation_matrix <- t(perturbation_matrix)

# Create the Seurat object
curr_seurat_obj <- CreateSeuratObject(
  counts = original_matrix,
  assay = "RNA"
)
perturbed_assay <- CreateAssayObject(
  counts = perturbation_matrix
)
curr_seurat_obj[["RNA_perturbed"]] <- perturbed_assay

# reorder the metadata to match the order of the cells in the Seurat object
metadata_update <- metadata[match(colnames(curr_seurat_obj), metadata$SCENICPLUS_label), ]
rownames(metadata_update) <- metadata_update$SCENICPLUS_label
curr_seurat_obj <- AddMetaData(curr_seurat_obj, metadata_update)

# Remove cells with extreme expression values
curr_seurat_obj <- subset(curr_seurat_obj, subset = nCount_RNA_perturbed > 200 & nCount_RNA_perturbed < median(nCount_RNA_perturbed) + 3 * mad(nCount_RNA_perturbed))

# Score the signature genes
curr_seurat_obj <- score_signature(curr_seurat_obj, Inv_high, "invasive_up")
curr_seurat_obj <- score_signature(curr_seurat_obj, Inv_low, "invasive_down")

# Plot the changes in invasive signature genes
curr_seurat_obj@meta.data <- curr_seurat_obj[[]] %>% 
  mutate(original_inv = invasive_up1 - invasive_down1,
         perturbed_inv = invasive_up_perturbed1 - invasive_down_perturbed1)

plot_signature_changes_individually(curr_seurat_obj = curr_seurat_obj, 
                                    cell_type = "Invasive-high OPC/NPC1", 
                                    TF_to_perturb = args$TF_to_perturb, 
                                    perturbation_dir = perturbation_dir, 
                                    color_palette = color_palette, 
                                    plot_title = "Invasive-high OPC/NPC1",
                                    signature = "invasive_signature",
                                    annotation_column = "CCI_CellClass_L2_2")

plot_signature_changes_individually(curr_seurat_obj = curr_seurat_obj, 
                                    cell_type = "Differentiated_like", 
                                    TF_to_perturb = args$TF_to_perturb, 
                                    perturbation_dir = perturbation_dir, 
                                    color_palette = color_palette, 
                                    plot_title = "Differentiated_like",
                                    signature = "invasive_signature",
                                    annotation_column = "CCI_CellClass_L2_2")

# Plot the changes in OPC signature genes
# Add module score for OPC and COP cells
dev_gene_list <- fread("opc_vs_cop.csv", header = TRUE)

# genes with positive log2FC are OPC, negative is COP
opc_genes <- dev_gene_list %>% 
  filter(log2fc < 0) %>%
  pull(gene)

cop_genes <- dev_gene_list %>%
  filter(log2fc > 0) %>%
  pull(gene)

# Add module score for OPC and COP cells
curr_seurat_obj <- score_signature(curr_seurat_obj, opc_genes, "opc_dev")
curr_seurat_obj <- score_signature(curr_seurat_obj, cop_genes, "cop_dev")

curr_seurat_obj@meta.data <- curr_seurat_obj[[]] %>% 
  mutate(original_opc = opc_dev1 - cop_dev1,
         perturbed_opc = opc_dev_perturbed1 - cop_dev_perturbed1)

# Plot the changes in OPC signature genes
plot_signature_changes_individually(curr_seurat_obj = curr_seurat_obj, 
                                    cell_type = "Invasive-high OPC/NPC1", 
                                    TF_to_perturb = args$TF_to_perturb, 
                                    perturbation_dir = perturbation_dir, 
                                    color_palette = color_palette,
                                    same_y_axis = TRUE, 
                                    plot_title = "Invasive-high OPC/NPC1", 
                                    signature = "opc_signature",
                                    annotation_column = "CCI_CellClass_L2_2")

plot_signature_changes_individually(curr_seurat_obj = curr_seurat_obj, 
                                    cell_type = "Differentiated_like", 
                                    TF_to_perturb = args$TF_to_perturb, 
                                    perturbation_dir = perturbation_dir, 
                                    color_palette = color_palette, 
                                    same_y_axis = TRUE, 
                                    plot_title = "Differentiated_like",
                                    signature = "opc_signature",
                                    annotation_column = "CCI_CellClass_L2_2")

# ---------------------------------------- Transition probability ---------------------------------------- #
# Sample level transition probability
plot_transition_probabilities_per_sample(curr_seurat_obj = curr_seurat_obj, 
                                        TF_to_perturb = args$TF_to_perturb, 
                                        perturbation_dir = perturbation_dir, 
                                        annotation_column = "CellClass_L5_2")

plot_transition_probabilities_per_sample(curr_seurat_obj = curr_seurat_obj,
                                        TF_to_perturb = args$TF_to_perturb,
                                        perturbation_dir = perturbation_dir,
                                        annotation_column = "CCI_CellClass_L2_2")