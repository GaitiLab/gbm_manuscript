# TF with high deviation from ChromVAR in inv high and high correlation to gene expression
# Identifying deviant TF motifs
seGroupMotif <- getGroupSE(ArchRProj = archr_proj, useMatrix = "MotifMatrix", groupBy = "Seurat_CCI_CellClass_L2_2")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

# Identify the maximum delta in z-score between Invasive-high OPC/NPC1 and Progenitor_like
rowData(seZ)$Inv_delta <- assay(seZ)[, "Invasive-high OPC/NPC1"] - assay(seZ)[, "Progenitor_like"]

# Identifying correlated TF motifs and TF expression
corGEM_MM <- correlateMatrices(
    ArchRProj = archr_proj,
    useMatrix1 = "GeneExpressionMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "LSI_Combined"
)

# Add maximum delta deviation to the correlation data frame
corGEM_MM$Inv_delta <- as.numeric(rowData(seZ)[match(corGEM_MM$MotifMatrix_name, rowData(seZ)$name), "Inv_delta"])

# Filter out duplicated TFs
corGEM_MM <- corGEM_MM[order(abs(corGEM_MM$cor), decreasing = TRUE), ]
corGEM_MM <- corGEM_MM[which(!duplicated(gsub("\\-.*","",corGEM_MM$MotifMatrix_name))), ]

# TF to label in the plot
high_SCENIC_TF <- motifs_df$TF[motifs_df$Subgroup == "Invasive-high OPC/NPC1"]
high_SCENIC_TF <- paste0(high_SCENIC_TF, "_")

low_SCENIC_TF <- motifs_df$TF[motifs_df$Subgroup == "Differentiated_like"]
low_SCENIC_TF <- paste0(low_SCENIC_TF, "_")
low_SCENIC_TF <- sort(low_SCENIC_TF)

# NPC/OPC markers
NPC_OPC_markers <- c("SOX2", "SOX4", "SOX11", "OLIG1", "OLIG2")
NPC_OPC_markers <- paste0(NPC_OPC_markers, "_")

# AP-1 family TFs
AP1_TF <- c("FOS", "FOSB", "FOSL1", "FOSL2", "FOSL2", "JUN", "JUNB", "JUND")

# Motifs to highlight in ChromVAR
Positive_TF <- corGEM_MM$MotifMatrix_name[corGEM_MM$cor > 0 & corGEM_MM$Inv_delta > quantile(corGEM_MM$Inv_delta, 0.9)]
Positive_TF <- Positive_TF[!is.na(Positive_TF)] # remove NA
Positive_TF_name <- substr(Positive_TF, 1, regexpr("_", Positive_TF) - 1) # remove number after _
Positive_TF_name <- paste0(Positive_TF_name, "_")

# Find high_SCENIC_TF intersect with Positive_TF
motifs_to_plot <- high_SCENIC_TF[high_SCENIC_TF %in% Positive_TF_name]

# Find TFs in each category
high_SCENIC_TF_to_plot <- high_SCENIC_TF[!sapply(high_SCENIC_TF, function(x) any(sapply(motifs_to_plot, function(y) grepl(y, x))))] # Identified in SCENIC but not in ChromVAR
high_SCENIC_TF_to_plot <- corGEM_MM$MotifMatrix_name[sapply(corGEM_MM$MotifMatrix_name, function(x) any(sapply(high_SCENIC_TF_to_plot, function(y) grepl(y, x))))]
motifs_to_plot <- Positive_TF[sapply(Positive_TF, function(x) any(sapply(motifs_to_plot, function(y) grepl(y, x))))] # Identified in both SCENIC and ChromVAR

# Label the TFs
corGEM_MM$TFRegulator <- "Not significant"
corGEM_MM$TFRegulator[corGEM_MM$MotifMatrix_name %in% high_SCENIC_TF_to_plot] <- "Candidate identified in SCENIC+"
corGEM_MM$TFRegulator[corGEM_MM$MotifMatrix_name %in% motifs_to_plot] <- "Putative regulator"
corGEM_MM$TFRegulator[sapply(corGEM_MM$MotifMatrix_name, function(x) any(sapply(NPC_OPC_markers, function(y) grepl(y, x))))] <- "NPC/OPC markers"
corGEM_MM$TFRegulator[sapply(corGEM_MM$MotifMatrix_name, function(x) any(sapply(AP1_TF, function(y) grepl(y, x))))] <- "AP-1 family TFs"

# Fix MotifMatrix_name - delete everything after MotifMatrix_name
corGEM_MM$MotifMatrix_name <- substr(corGEM_MM$MotifMatrix_name, 1, regexpr("_", corGEM_MM$MotifMatrix_name) - 1)

# Visualize the correlation between TF motif and gene expression
ggplot(data.frame(corGEM_MM), aes(cor, Inv_delta, color = TFRegulator)) +
  geom_point(data = data.frame(corGEM_MM[corGEM_MM$TFRegulator == "Not significant",]), aes(cor, Inv_delta, color = TFRegulator), size = 1, alpha = 0.75) +
  geom_point(data = data.frame(corGEM_MM[corGEM_MM$TFRegulator != "Not significant",]), aes(cor, Inv_delta, color = TFRegulator), size = 2) +
  GBM_theme() +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = quantile(corGEM_MM$Inv_delta, 0.95), lty = "dashed", color = "darkgrey") +
    scale_color_manual(values = c("Not significant"="darkgrey",
                                "Putative regulator"="#2B7095",
                                "Candidate identified in SCENIC+" = "#EDAE49FF",
                                "NPC/OPC markers"="#7C9EB5",
                                "AP-1 family TFs"="#C05E00")) +
  geom_label_repel(data = data.frame(corGEM_MM[corGEM_MM$TFRegulator != "Not significant",]), aes(label = MotifMatrix_name), size = 3, nudge_x = 0.1, nudge_y = 0.1, max.overlaps=10) +
  labs(
    y = "TF motif accessibility difference between \ninvasive-high OPC/NPC1 and progenitor-like (Δz-score)",
    x = "Correlation of TF motif accessibility and TF expression"
  ) +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(min(corGEM_MM$Inv_delta)*1.05, max(corGEM_MM$Inv_delta)*1.05)
  ) +
  scale_x_continuous(
    expand = c(0,0), 
    limits = c(-1, 1)
  ) +
  theme(legend.position = "bottom",
        legend.title = element_blank())
ggsave(paste0(plot_dir, "/corGEM_MM_TF_Regulator_95.pdf"), width = 7, height = 7.5)

# ---------------------------- Scatter plot for each TF ---------------------------- #
p <- readRDS(paste0(plot_dir, "/TF_deviations.rds"))

motifs_to_plot <- c("TCF12", "TCF4", "TCF3", "ASCL1", "ZEB1", "KLF12", "OLIG1", "OLIG2", "SOX2", "SOX4", "SOX10", "SOX11", "FOSL1", "FOSL2")

df_neuronal <- NULL
for (i in 1:length(p)) {
  p1 <- p[[i]]

  # get the name of the TF
  curr_TF <- names(p)[i]
  curr_TF <- strsplit(curr_TF, "_")[[1]][1]
  curr_TF <- strsplit(curr_TF, "z:")[[1]][2]
  log_info("Currently analyzing: ", curr_TF)

  if (!(curr_TF %in% motifs_to_plot)) next

  # add the data to the dataframe
  df_neuronal <- rbind(df_neuronal, p1$data %>% mutate(TF = curr_TF))
}

# Load the expression data for the TFs
exp_mat <- readRDS(paste0(plot_dir, "/TF_exp_mtx.csv"))
exp_mat$archr_cellnames <- sub("-1.*", "", exp_mat$archr_cellnames)

for (TF_name in unique(df_neuronal$TF)) {
  log_info("Currently analyzing: ", TF_name)
  curr_TF_chromVAR <- df_neuronal[df_neuronal$TF == TF_name,] %>%
    rownames_to_column(var = "archr_cellnames") 
  curr_TF_chromVAR$archr_cellnames <- sub("-1.*", "", curr_TF_chromVAR$archr_cellnames)
  curr_TF_exp <- exp_mat %>% 
    filter(archr_cellnames %in% curr_TF_chromVAR$archr_cellnames) %>%
    select(archr_cellnames, TF_name)
  curr_df <- merge(curr_TF_chromVAR, curr_TF_exp, by = "archr_cellnames") %>%
    mutate(Sample = sub("#.*", "", archr_cellnames),
          Barcode = sub(".*#", "", archr_cellnames))

  # Calculate the average by sample and cell type
  curr_df <- curr_df %>%
    group_by(Sample, x) %>%
    summarise(mean_gex = mean(get(TF_name)),
              mean_acc = mean(y),
              n_cell = n()) %>%
    ungroup() %>%
    filter(n_cell > 10) # filter out low sample

  # Aggregate at cell type level
  curr_df <- curr_df %>%
    group_by(x) %>%
    summarise(mean_cell_type_gex = mean(mean_gex),
              mean_cell_type_acc = mean(mean_acc),
              total_n = sum(n_cell))

  # Rename cell type
  curr_df <- curr_df %>%
    mutate(x = case_when(
      x == "Malignant_NPC1" ~ "NPC1-like",
      x == "Malignant_NPC2" ~ "NPC2-like",
      x == "Malignant_OPC" ~ "OPC-like",
      x == "Invasive-high_OPC_NPC1" ~ "Invasive-high OPC/NPC1",
      x == "Malignant_AC" ~ "AC-like",
      x == "Malignant_MES_AST" ~ "MES-AC-like",
      x == "Malignant_MES_INT" ~ "MES-INT-like",
      x == "Malignant_MES_HYP" ~ "MES-HYP-like")
    )

  color_palette <- c("NPC1-like" = "#86CCE8",
                    "NPC2-like" = "#ADD7E5",
                    "OPC-like" = "#AECFA5",
                    "AC-like" = "#DFB63B",
                    "MES-AC-like" = "#F37F72",
                    "MES-INT-like" = "#B02425",
                    "MES-HYP-like" = "#F1634B",
                    "Invasive-high OPC/NPC1" = "#2B7095")

  # Plot with legend
  ggplot(curr_df, aes(x = mean_cell_type_gex, y = mean_cell_type_acc)) +
    geom_point(aes(size = total_n, color = "#000000", fill = x), shape = 21, stroke = 1) + # Set stroke and shape here
    geom_text_repel(aes(label = x), size = 5) +
    scale_color_manual(values = color_palette) + # Use custom color palette
    scale_fill_manual(values = color_palette) + # Use custom color palette for fill
    labs(x = "Scaled mean gene expression", y = "Motif accessibility (mean z-score)") +
    GBM_theme() +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    ggtitle(TF_name)
  ggsave(paste0(plot_dir, "/Neuronal_motif_deviation_score_scatterplot_legend_", TF_name, ".pdf"), width = 5, height = 5.5)
}

