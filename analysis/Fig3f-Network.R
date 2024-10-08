library(igraph)
library(readr)

# Loading eRegulon metadata 
eRegulons <- fread(paste0(args$input, "/eRegulons_simplified.csv"))

# Identify background genes = all genes predicted by SCENIC+ as targets of any TF
all_target_genes <- unique(eRegulons$Gene)

# Get the list of target genes for each TF
target_gene_per_TF <- eRegulons %>%
  select(TF, Gene) %>%
  distinct() %>%
  pivot_wider(names_from = TF, values_from = Gene) %>%
  t() %>%
  as.data.frame()
colnames(target_gene_per_TF) <- "Genes"

# Iterate through each row of the data frame to create a list of target genes for each TF
target_gene_list <- list()
for (i in 1:nrow(target_gene_per_TF)) {
  TF_name <- rownames(target_gene_per_TF)[i]
  target_gene_list[[TF_name]] <- unlist(target_gene_per_TF$Genes[i])
}

# Load the data
gene_list <- read_csv("multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/df/SCENICplus_activator_gene_list.csv")

# Inv_high TFs
TF_priority_list <- read_csv("multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/df/all_regions/TFs_to_plot.csv")
TF_priority <- TF_priority_list$TF[TF_priority_list$Cell_type == "Invasive-high OPC/NPC1"]

# Filter the gene list to keep columns that are in the TF priority list
TF_piority_existed <- intersect(TF_priority, colnames(gene_list))
gene_list <- gene_list[, c(TF_piority_existed)]

# Create an empty list to store edges
edges <- vector("list", length = 0)

# Loop over each column and row to create edges and keep only TFs that are in the priority list
for (gene in names(gene_list)) {
    for (target in gene_list[[gene]]) {
        if (!is.na(target) && (target %in% TF_priority) && (gene != target)) {  # Check if the target is in the priority TF list
            edges <- append(edges, list(c(gene, target)))
        }
    }
}

# Create a dataframe from the list of edges
edge_df <- do.call(rbind, edges)
gene_network <- graph_from_edgelist(as.matrix(edge_df), directed = TRUE)

# Load the list of differentially expressed TFs
opc_diff_tfs_path <- "misc/data/OPC_diff_TFs.csv"
opc_diff_tfs <- read_csv(opc_diff_tfs_path)
opc_diff_tfs <- opc_diff_tfs$x  # Assuming 'x' is the column with the TF names

# Calculate the degree of each node
degree_values <- degree(gene_network, mode = "out")

# Normalize degrees for size scaling
max_degree <- max(degree_values)
normalized_degrees <- degree_values / max_degree

# Set base size and scale factor for the nodes
base_size <- 8
scale_factor <- 15

# Set node and edge attributes
V(gene_network)$size <- base_size + (normalized_degrees * scale_factor)  # Bigger dots for more connected nodes
V(gene_network)$color <- ifelse(V(gene_network)$name %in% opc_diff_tfs, 
                                "#9BB98CCC",  # OPC_diff_TFs in green
                                "#EFD496CC")   # Other nodes in yellow
E(gene_network)$color <- "gray"

# Set node labels
V(gene_network)$label <- V(gene_network)$name
V(gene_network)$label.cex <- 0.5  # Control label size

# Choose a layout; here we use layout_with_fr for the Fruchterman-Reingold layout
layout <- layout_with_fr(gene_network)

# Save the plot
plot_path <- "multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/Plots/PT,TE/gene_network.pdf"
pdf(plot_path, width = 10, height = 10)
# Plot the network with custom attributes
network <- plot(gene_network, layout=layout.fruchterman.reingold, main="fruchterman.reingold",
     vertex.label.color = "black", 
     vertex.label.family = "sans",
     vertex.shape = "circle",
     edge.arrow.size = 0.5,
     edge.curved = 0.1)  # Slight curve for edges

# Add legend
legend("topright",  # Location of the legend in the plot
       legend=c("TFs reported in OPC differentiation", "Candidate identified in SCENIC+"),  # Labels
       col=c("#9BB98CCC", "#EFD496CC"),  # Colors corresponding to the labels
       pch=21,  # Type of point, 21 is a filled circle which matches nodes
       pt.bg=c("#9BB98CCC", "#EFD496CC"),  # Background color for points
       cex=0.8,  # Text size
       bty="n")  # No box around the legend

# Add legend for node sizes
size_legend_labels <- c(min(degree_values), median(degree_values), max(degree_values))
size_legend_values <- c(base_size, base_size + (0.5 * scale_factor), base_size + scale_factor)

legend("bottomright",  # Location of the legend in the plot
       legend=size_legend_labels,  # Labels
       pt.cex=size_legend_values / 5,  # Scale point sizes to match the sizes in the plot
       pch=21,  # Type of point, 21 is a filled circle which matches nodes
       cex=0.8,  # Text size
       bty="n",  # No box around the legend
       title="Node Degree (Number of TF targets)",
       y.intersp=2)  # Increase vertical spacing between legend items

network
dev.off()
   
# ----------------------------------------- Premutation test ----------------------------------------- #
# plot distribution of node degrees + visualize OPC TFs
network_df <- data.frame(
  node = V(gene_network)$name,
  degree = degree_values,
  OPC_TF = V(gene_network)$name %in% opc_diff_tfs
)

# rank the TFs by degree
network_df <- network_df %>%
  arrange(desc(degree))

# Calculate the node degree of top 10% of TFs
ranked_degrees <- sort(degree_values, decreasing = TRUE)
top_10_percent <- ranked_degrees[1:round(length(ranked_degrees) * 0.1)]
master_regulator_number <- length(top_10_percent)
mean_opc_degree <- mean(network_df$degree[network_df$OPC_TF])

ggplot(network_df, aes(x = degree, fill = OPC_TF)) +
  geom_histogram(binwidth = 1, alpha = 0.8) +
  geom_vline(xintercept = mean_opc_degree, linetype = "dashed", color = "red") +
  labs(title = "Distribution of node degrees",
       x = "Node degree",
       y = "Number of TFs",
       fill = "OPC TF") +
  default_theme() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values = c(`TRUE` = "#9BB98C", `FALSE` = "#EFD496"))
ggsave("multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/Plots/PT,TE/node_degree_distribution.pdf", width = 6, height = 4)

# Define master regulators as top 80% of TFs by degree
degree_cutoff <- quantile(network_df$degree, 0.8)

# Identify actual master regulators
network_df <- network_df %>%
  mutate(Master_Regulator = degree >= degree_cutoff)

# Observed count of actual OPC TFs among master regulators
observed_count <- sum(network_df$OPC_TF & network_df$Master_Regulator)

# Permutation test
set.seed(123)  # For reproducibility
n_permutations <- 10000
perm_counts <- replicate(n_permutations, {
  # Shuffle the OPC_TF labels
  shuffled_OPC_TF <- sample(network_df$OPC_TF)
  sum(shuffled_OPC_TF & network_df$Master_Regulator)
})

# Calculate p-value
p_value <- mean(perm_counts >= observed_count)

# Visualize the permutation test results
ggplot(data.frame(perm_counts), aes(x = perm_counts)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black") +
  geom_vline(xintercept = observed_count, linetype = "dashed", color = "red") +
  labs(title = "Permutation Test for Master Regulators",
       x = "Number of OPC TFs among Master Regulators",
       y = "Frequency") +
  default_theme() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
ggsave("/multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/Plots/PT,TE/permutation_test.pdf", width = 4, height = 3)
