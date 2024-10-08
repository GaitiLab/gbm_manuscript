motif_results_paths <- list.files(paste0(args$output_dir, "/homer_motif/Top1000"), pattern = "knownResults.txt", full.names = TRUE, recursive = TRUE)

# Function to perform BH correction
bh_correction <- function(p_values) {
  p.adjust(p_values, method = "BH")
}

motif_results <- setNames(lapply(motif_results_paths, function(x) {
  # Read the data
  data <- fread(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data$cell_type <- basename(dirname(x))
  data$FDR <- bh_correction(data$`P-value`)
  return(data)
}), sapply(motif_results_paths, function(x) {
  # Extract cell type from the path
  basename(dirname(x))
}))
all_motif_results <- do.call(rbindlist, list(motif_results, fill = TRUE))

# Format motif name
all_motif_results$motif <- gsub("\\(.*", "", all_motif_results$`Motif Name`)
all_motif_results <- all_motif_results %>%
  mutate(motif = toupper(motif),
         target_sequence_with_motif = sub("%", "", `% of Target Sequences with Motif`),
         background_sequence_with_motif = sub("%", "", `% of Background Sequences with Motif`)) %>%
  mutate(target_sequence_with_motif = as.numeric(target_sequence_with_motif),
         background_sequence_with_motif = as.numeric(background_sequence_with_motif)) %>%
  mutate(pct_diff = target_sequence_with_motif - background_sequence_with_motif) %>%
  mutate(`P-value` = as.numeric(`P-value`),
         FDR = as.numeric(FDR))

# Specify motifs of interest
Motif_of_interest <- c("ASCL1", "TCF3", "TCF4", "TCF12", "ZEB1", "KLF12", "NFIA", "SOX10")
Differentation_motif <- c("SOX10")

# OPC vs COP rank plot
OPC_motif_results <- all_motif_results_diff %>%
  filter(cell_type == "OPC") %>%
  mutate(p_value_OPC = FDR,
         pct_target_OPC = target_sequence_with_motif,
         motif = `Motif Name`) %>%
  select(motif, p_value_OPC, pct_target_OPC)

COP_motif_results <- all_motif_results_diff %>%
  filter(cell_type == "COP") %>%
  mutate(p_value_COP = FDR,
          pct_target_COP = target_sequence_with_motif, 
          motif = `Motif Name`) %>%
  select(motif, p_value_COP, pct_target_COP)

# Creating diff_motif_results with necessary columns and calculations
diff_motif_results <- merge(OPC_motif_results, COP_motif_results, by = "motif") %>%
  mutate(cell_type_pct_diff = pct_target_OPC - pct_target_COP,
         motif = gsub("\\(.*", "", motif)) %>%
  mutate(motif = toupper(motif)) %>%
  arrange(desc(cell_type_pct_diff)) %>%
  mutate(rank = row_number()) %>%
  group_by(motif) %>%
  filter(rank == min(rank)) %>%
  ungroup() %>%
  mutate(rank = row_number()) %>%
  mutate(significant_status = case_when(
    p_value_OPC < 0.05 & p_value_COP > 0.05 & cell_type_pct_diff > 0 ~ "Significant in OPC",
    p_value_OPC > 0.05 & p_value_COP < 0.05 & cell_type_pct_diff < 0 ~ "Significant in COP",
    p_value_OPC < 0.05 & p_value_COP < 0.05 ~ "Significant in both",
    TRUE ~ "Not significant"
  )) %>%
  mutate(fill_color = ifelse(motif %in% Differentation_motif, "Differentiation", 
                             ifelse(motif %in% Motif_of_interest, "Stemness", "Other"))) %>%
  mutate(p_value = ifelse(cell_type_pct_diff > 0, p_value_OPC, p_value_COP))

# Plotting the data
ggplot(diff_motif_results, aes(x = rank, y = cell_type_pct_diff)) + 
  geom_point(aes(fill = significant_status, size = -log10(p_value), alpha = significant_status), shape = 21, stroke = 0.5) +
  geom_text_repel(
    data = subset(diff_motif_results, fill_color == "Stemness"),
    aes(label = motif),
    size = 3,
    nudge_x = 20,
    nudge_y = 4, 
    segment.size = 0.2,
    force = 4,
    segment.linetype = "dashed",
    color = "#CD2027"
  ) +
  geom_text_repel(
    data = subset(diff_motif_results, fill_color == "Differentiation"),
    aes(label = motif),
    size = 3,
    nudge_x = -40,
    nudge_y = 4, 
    segment.size = 0.2,
    force = 4,
    segment.linetype = "dashed",
    color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  scale_size_continuous(range = c(0, 5)) +
  scale_fill_manual(values = c("Significant in OPC" = "#ED8B22", "Significant in COP" = "#F0E9BA", "Significant in both" = "black", "Not significant" = "grey")) +
  scale_alpha_manual(values = c("Significant in OPC" = 1, "Significant in COP" = 1, "Significant in both" = 1, "Not significant" = 0.05)) +
  labs(
    x = "Rank",
    y = "Motif enrichment differences \nin OPC vs. COP",
    fill = "Fill Color",
    color = "Significant in OPC",
    alpha = "Significance"
  ) +
  # Make sure the y-axis symmetry
  ylim(c(-max(abs(diff_motif_results$cell_type_pct_diff), na.rm = TRUE), max(abs(diff_motif_results$cell_type_pct_diff), na.rm = TRUE))) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_text(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )
ggsave(paste0(plot_dir, "/Differential_motif_enrichment_OPC_vs_COP.pdf"), width = 7, height = 4, dpi = 300)