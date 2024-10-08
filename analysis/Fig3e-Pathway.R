Inv_high_TF <- c("TCF12", "TCF4", "ASCL1", "ZEB1", "TCF3", "KLF12")
Inv_high_TF_gene_list <- list()
for (TF_name in Inv_high_TF) {
  Inv_high_TF_gene_list <- append(Inv_high_TF_gene_list, target_gene_list[[TF_name]])
}
Inv_high_TF_gene_list <- unique(unlist(Inv_high_TF_gene_list))

run_ora(
  msigdb_cat_list = msigdb_cat_list,
  all_gene_sets = all_gene_sets,
  fg_gene_list = Inv_high_TF_gene_list,
  bg_gene_list = all_target_genes,
  output_name = "Inv_high_TF",
  out_dir = ORA_dir
)

# Color
low_color = "#FDBB84FF"
high_color = "#B30000FF"

# log_info("Plotting ORA results...")
all_sig_results <- all_sig_results %>% 
  mutate(gene_fraction = overlap / size)

# Only keep top 20 for each msigdb_cat for visualization
all_sig_results <- all_sig_results %>%
  group_by(msigdb_cat) %>%
  top_n(20, -padj) %>%
  filter(msigdb_cat == "C2_CP:REACTOME") %>%
  filter(!pathway %in% terms_to_exclude)
all_sig_results <- all_sig_results %>% 
  mutate(pathway = str_replace_all(pathway, "_", " ")) %>%
  mutate(pathway = str_replace_all(pathway, "REACTOME ", "")) %>%
  mutate(pathway = tolower(pathway)) %>%
  mutate(pathway = str_to_title(pathway))

# Plot
sig_result_plots <- ggplot(all_sig_results, aes(y = reorder(pathway, gene_fraction), x = gene_fraction)) +
  geom_segment(aes(yend = reorder(pathway, gene_fraction), xend = 0), color = "grey") +  # Add lines for the lollipop stick
  geom_point(aes(color = -log10(padj), size = 2)) +  # Points as the lollipop heads
  facet_wrap(. ~ msigdb_cat, scales = "free", ncol = 1) +
  scale_color_gradient(low = low_color, high = high_color) +
  theme_minimal() +
  labs(title = "Pathway Enrichment (hypergeometric test) p-adj<0.05", 
       x = "Overlapping gene fraction in the geneset",
       y = "Pathway") +
  default_theme()
ggsave(paste0(ORA_dir, "/ORA_Inv_high_TF_all_formatted.pdf"),
    width = 8, height = length(unique(all_sig_results$msigdb_cat))*5, units = "in", bg = "white")

# ----------------------------- NOTCH network ----------------------------- #
# Load the data
ORA_dir <- paste0(args$output_dir, "/ORA")
ORA_results <- fread(paste0(ORA_dir, "/ORA_Inv_high_TF_KLF12_included.csv"))

NOTCH_target <- ORA_results %>%
  filter(pathway == "REACTOME_NOTCH_HLH_TRANSCRIPTION_PATHWAY") %>%
  pull(overlapGenes)
NOTCH_target <- strsplit(NOTCH_target, split = "\\|")[[1]]

# Inv_high TFs
Inv_high_TF <- c("TCF12", "TCF4", "ASCL1", "ZEB1", "TCF3", "KLF12")
Inv_high_TF_gene_list <- list()
for (TF in Inv_high_TF) {
  Inv_high_TF_gene_list[[TF]] <- target_gene_list[[TF]]
}

# Create an edge list for the network
edges <- data.frame()
for (TF in Inv_high_TF) {
  targets <- Inv_high_TF_gene_list[[TF]]
  # Only include targets that are in the NOTCH target genes
  common_targets <- intersect(targets, NOTCH_target)
  if (length(common_targets) > 0) {
    temp <- data.frame(from = TF, to = common_targets)
    edges <- rbind(edges, temp)
  }
}

# Create the graph object
g <- graph_from_data_frame(edges, directed = TRUE)

# Set vertex colors: TFs colored as #2B7095 and others as default
V(g)$color <- ifelse(V(g)$name %in% Inv_high_TF, "#2B7095", "lightgray")

# Create a custom layout
layout <- layout_with_fr(g)  # Initial layout with Fruchterman-Reingold

# Set y-coordinates to separate TFs and target genes
y_positions <- ifelse(V(g)$name %in% Inv_high_TF, 1, -1)

# Spread TFs and target genes apart horizontally
x_positions_TF <- seq(-3, 3, length.out = sum(V(g)$name %in% Inv_high_TF))
x_positions_genes <- seq(-6, 6, length.out = sum(V(g)$name %in% NOTCH_target))

# Assign x and y positions
layout[V(g)$name %in% Inv_high_TF, ] <- cbind(x_positions_TF, 2)
layout[V(g)$name %in% NOTCH_target, ] <- cbind(x_positions_genes, -2)

# Plot the network with increased space
plot(g, 
     layout = layout, 
     vertex.label.color = "black", 
     vertex.size = 100, 
     edge.arrow.size = 0.5, 
     main = "TF to NOTCH Target Gene Network",
     rescale = FALSE,
     ylim = c(-2, 2),
     xlim = c(-6, 6))

# save the plot
plot_path <- "TF_to_NOTCH_target_gene_network.pdf"
pdf(plot_path, width = 12, height = 6)
plot(g, 
     layout = layout, 
     vertex.label.color = "black", 
     vertex.size = 100, 
     edge.arrow.size = 0.5, 
     main = "TF to NOTCH Target Gene Network",
     rescale = FALSE,
     ylim = c(-2, 2),
     xlim = c(-6.2, 6.2))
dev.off()