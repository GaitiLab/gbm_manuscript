library(igraph)
library(readr)

# ---------------------------------------- Venn diagram for overlaping TFs for lineage-specific genes ---------------------------------------- #
# Load the list of lineage-specific genelist
lineage_genes <- fread("misc/data/Marker_genes.csv", header = TRUE)
lineage_gene_list <- split(lineage_genes$Gene, lineage_genes$class)

invasive_high_TFs <- read_csv(
    "multiome_results/12_SCENIC_plus/outs/TFs_to_plot.csv"
)
invasive_high_TFs <- invasive_high_TFs$TF[
    invasive_high_TFs$Cell_type == "Invasive-high OPC/NPC1"
]
lineage_gene_list[["Invasive-high TFs"]] <- invasive_high_TFs

lineage_gene_list <- lapply(
    lineage_gene_list,
    function(x) x[x %in% invasive_high_TFs]
)

output_file <- paste0(args$output_dir, "/Plots/venn_lineage.pdf")
pdf(output_file, width = 5, height = 5)
plot(euler(lineage_gene_list, shape = "ellipse"), quantities = TRUE)
dev.off()

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
gene_list <- read_csv(
    "multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/df/SCENICplus_activator_gene_list.csv"
)

# Inv_high TFs
TF_priority_list <- read_csv(
    "multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/df/all_regions/TFs_to_plot.csv"
)
TF_priority <- TF_priority_list$TF[
    TF_priority_list$Cell_type == "Invasive-high OPC/NPC1"
]

# Filter the gene list to keep columns that are in the TF priority list
TF_piority_existed <- intersect(TF_priority, colnames(gene_list))
gene_list <- gene_list[, c(TF_piority_existed)]

# Create an empty list to store edges
edges <- vector("list", length = 0)

# Loop over each column and row to create edges and keep only TFs that are in the priority list
for (gene in names(gene_list)) {
    for (target in gene_list[[gene]]) {
        if (!is.na(target) && (target %in% TF_priority) && (gene != target)) {
            # Check if the target is in the priority TF list
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
opc_diff_tfs <- opc_diff_tfs$x # Assuming 'x' is the column with the TF names

# Calculate the degree of each node
degree_values <- degree(gene_network, mode = "out")

# Normalize degrees for size scaling
max_degree <- max(degree_values)
normalized_degrees <- degree_values / max_degree

# Set base size and scale factor for the nodes
base_size <- 8
scale_factor <- 15

# Set node and edge attributes
V(gene_network)$size <- base_size + (normalized_degrees * scale_factor) # Bigger dots for more connected nodes
V(gene_network)$color <- ifelse(
    V(gene_network)$name %in% opc_diff_tfs,
    "#9BB98CCC", # OPC_diff_TFs in green
    "#EFD496CC"
) # Other nodes in yellow
E(gene_network)$color <- "gray"

# Set node labels
V(gene_network)$label <- V(gene_network)$name
V(gene_network)$label.cex <- 0.5 # Control label size

# Choose a layout; here we use layout_with_fr for the Fruchterman-Reingold layout
layout <- layout_with_fr(gene_network)

# Save the plot
plot_path <- "multiome_results/12_SCENIC_plus/ATAC/Integrated_confident/Plots/PT,TE/gene_network.pdf"
pdf(plot_path, width = 10, height = 10)
# Plot the network with custom attributes
network <- plot(
    gene_network,
    layout = layout.fruchterman.reingold,
    main = "fruchterman.reingold",
    vertex.label.color = "black",
    vertex.label.family = "sans",
    vertex.shape = "circle",
    edge.arrow.size = 0.5,
    edge.curved = 0.1
) # Slight curve for edges

# Add legend
legend(
    "topright", # Location of the legend in the plot
    legend = c(
        "TFs reported in OPC differentiation",
        "Candidate identified in SCENIC+"
    ), # Labels
    col = c("#9BB98CCC", "#EFD496CC"), # Colors corresponding to the labels
    pch = 21, # Type of point, 21 is a filled circle which matches nodes
    pt.bg = c("#9BB98CCC", "#EFD496CC"), # Background color for points
    cex = 0.8, # Text size
    bty = "n"
) # No box around the legend

# Add legend for node sizes
size_legend_labels <- c(
    min(degree_values),
    median(degree_values),
    max(degree_values)
)
size_legend_values <- c(
    base_size,
    base_size + (0.5 * scale_factor),
    base_size + scale_factor
)

legend(
    "bottomright", # Location of the legend in the plot
    legend = size_legend_labels, # Labels
    pt.cex = size_legend_values / 5, # Scale point sizes to match the sizes in the plot
    pch = 21, # Type of point, 21 is a filled circle which matches nodes
    cex = 0.8, # Text size
    bty = "n", # No box around the legend
    title = "Node Degree (Number of TF targets)",
    y.intersp = 2
) # Increase vertical spacing between legend items

network
dev.off()
