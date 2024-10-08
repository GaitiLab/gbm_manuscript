# Load the Seurat object
curr_seurat_data <- readRDS(paste0(args$output_dir, "/GBM3D.rds"))

plot_cor_TF_and_inv_per_cell <- function(curr_seurat_data, tf_name, chromVar_names_path, plot_dir) {
  # Get invasive signature score at single cell level
  invasive_score <- curr_seurat_data[[]] %>%
    select(invasive_score) %>%
    rownames_to_column(var = "Cell_id")

  # ChromVar score for the given transcription factor
  chromVar_matrix <- GetAssayData(object = curr_seurat_data, assay = "ChromVar")
  chromVar_matrix <- chromVar_matrix %>% as.matrix() %>% t()

  # Split the matrix into two parts
  chromVar_matrix_1 <- chromVar_matrix[,1:(ncol(chromVar_matrix)/2)]
  chromVar_matrix_2 <- chromVar_matrix[,(ncol(chromVar_matrix)/2 + 1):ncol(chromVar_matrix)]

  # Check if the ChromVar object has the TF names of interest
  if(length(grep("TCF12", colnames(chromVar_matrix_2))) == 0){
    log_info("Need to correct the ChromVar TF names...")
    # Function to find the best match where row name is a substring of a longer name in chromVar_names
    match_names <- function(target_names, full_names) {
      sapply(target_names, function(name) {
        matches <- grep(paste0("^", name, "(_|$)"), full_names, value = TRUE)  # Find matches using regex
        if(length(matches) > 0) matches[1] else NA  # Return the first match or NA if none
      })
    }

    # ChromVar TF names
    chromVar_names <- fread(chromVar_names_path)$chromVar_names
    
    # Apply the matching function to get matched names from chromVar_names for each rowname
    matched_names <- match_names(colnames(chromVar_matrix_2), chromVar_names)
    valid_matches <- matched_names[!is.na(matched_names)]
    chromVar_matrix_2 <- chromVar_matrix_2[,colnames(chromVar_matrix_2) %in% names(valid_matches)]
    colnames(chromVar_matrix_2) <- valid_matches[match(colnames(chromVar_matrix_2), names(valid_matches))]
    names(colnames(chromVar_matrix_2)) <- NULL
  }

  # Match the order of the ChromVar matrix
  colnames(chromVar_matrix_1) <- gsub("-", "_", colnames(chromVar_matrix_1)) # replace - with underscore
  common_names <- intersect(colnames(chromVar_matrix_1), colnames(chromVar_matrix_2))
  chromVar_matrix_1 <- chromVar_matrix_1[,common_names]
  chromVar_matrix_2 <- chromVar_matrix_2[,common_names]

  # Add the two matrices
  chromVar_matrix <- chromVar_matrix_1 + chromVar_matrix_2

  chromVar_matrix <- chromVar_matrix %>% 
    as.data.frame() %>%
    rownames_to_column(var = "Cell_id") %>%
    left_join(invasive_score, by = c("Cell_id" = "Cell_id"))

  # Check if the TF is present in the ChromVar matrix
  tf_name_matrix <- grep(tf_name, colnames(chromVar_matrix), value = TRUE)
  if(length(tf_name_matrix) == 0){
    print(paste0("Transcription factor ", tf_name, " not found in ChromVar matrix"))
  } else {
    # In case there is more than one TF
    for (i in 1:length(tf_name_matrix)){
      df <- chromVar_matrix[,c("Cell_id", tf_name_matrix[i], "invasive_score")]
      colnames(df) <- c("Cell_id", "ChromVar", "invasive_score")

      # Plot
      ggplot(df, aes(x = invasive_score, y = ChromVar)) +
        geom_hex(bins = 100) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +
        stat_cor(method = "spearman") +
        labs(x = "Invasive signature score", 
              y = paste0(tf_name, " motif binding site"),
              subtitle = tf_name_matrix[i]) +
        GBM_theme() +
        theme(legend.key.height = unit(0.25, 'cm'),  # Adjust height of legend key
              legend.key.width = unit(1, 'cm')) +   # Adjust width of legend key
        scale_fill_gradient(limits = c(0, 100))

      ggsave(paste0(plot_dir, "/Correlation_invasive_signature_ChromVar_per_cell_", tf_name_matrix[i], ".pdf"), width = 4, height = 5)
    }
  }
}

# Plot for each TF
for (tf_name in TFs_to_plot){
  plot_cor_TF_and_inv_per_cell(curr_seurat_data, tf_name, "misc/data/chromVar_names.csv", plot_dir)
}
