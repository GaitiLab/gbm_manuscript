# ---- Code to reproduce Figure 4c ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #
# Load required packages
pacman::p_load(igraph, readr, data.table, tidyverse, eulerr)

# Required inputs
params <- list(
    gene_lists_path = "misc/venn_uncommitted_opc_TFs.csv",
    plot_dir = "output/figures",
)

GaitiLabUtils::create_dir(params$plot_dir)

# Load gene lists
gene_lists_df <- read.csv(gene_lists_path, stringsAsFactors = FALSE)
gene_lists <- split(gene_lists_df$Gene, gene_lists_df$Category)

# Plot the Venn diagram
output_file <- file.path(params$plot_dir, "Fig4c_venn_lineage.pdf")
pdf(output_file, width = 5, height = 5)
plot(euler(lineage_gene_list, shape = "ellipse"), quantities = TRUE)
dev.off()
