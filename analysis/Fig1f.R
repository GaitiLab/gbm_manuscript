# ---- Code to reproduce Figure 1f ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(
    data.table,
    patchwork,
    ggplot2,
    dplyr,
    tidyr,
    tibble,
    Seurat,
    DESeq2,
    SingleCellExperiment,
    apeglm,
    Matrix.utils
)

# Required inputs
params <- list(
    # Path to Seurat object generated using this manuscript's data, raw data and final metadata can be downloaded online, see publication
    seurat_obj_path = "",
    output_dir = "output",
    plot_dir = "output/figures",
    protein_coding_genes_list_path = "misc/ensembl_protein_coding_genes.csv",
    genelists_path = "misc/gene_signatures.xlsx"
)

# Creating directories needed for outputs
GaitiLabUtils::create_dir(params$output_dir)
GaitiLabUtils::create_dir(params$plot_dir)

# ---- Load data & Data wrangling ---- #
seurat_object <- readRDS(params$seurat_obj_path)

# Prepare groups for comparison
# Sample level region comparisons
seurat_object$Region <- ifelse(seurat_object$Region == "PT", "PT", "Rest")

seurat_object <- subset(
    seurat_object,
    # Keep Stringent
    subset = is_malignant_confident == TRUE
)

seurat_object <- subset(
    seurat_object,
    # Only select matched patients
    subset = Patient %in% c("6237", "6245", "6467", "6419", "6514", "6425")
)

seurat_object <- subset(
    seurat_object,
    # Subset for cell types of interest
    subset = CellClass_L3 %in% c("Malignant_OPC", "Malignant_NPC1")
)

# Convert to SingleCellExperiment
counts <- seurat_object@assays$RNA@counts

# Protein coding genes only
protein_coding_genes <- read.csv(params$protein_coding_genes_list_path)
counts <- counts[rownames(counts) %in% protein_coding_genes$hgnc_symbol, ]

metadata <- seurat_object[[]]
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

# Aggregate counts
groups <- colData(sce)[, c("Sample")] # Group single cell data by sample to perform pseudobulking
aggr_counts <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum")
aggr_counts <- t(aggr_counts)

# Prepare metadata
metadata <- colData(sce) %>%
    as.data.frame() %>%
    dplyr::select(c("Sample", "Patient", "Region"))

# Remove duplicates
metadata <- metadata[!duplicated(metadata), ]
# Rename rows
rownames(metadata) <- metadata$Sample
# Cell counts
t <- table(colData(sce)$Sample)
df <- as.data.frame(t)
colnames(df) <- c("Sample", "cell_count")

## Join data frames
df <- plyr::join(df, metadata, by = "Sample")

# Update rownames of metadata to match colnames of count matrix, as needed later for DE
index <- match(colnames(aggr_counts), df$Sample)
df <- df[index, , drop = FALSE]

print(df)

# Converting col metadata to factors
df$Patient <- factor(df$Patient, levels = unique(df$Patient))
df$Region <- factor(df$Region, levels = unique(df$Region))

# Create DESeq2 object
deseq <- DESeqDataSetFromMatrix(
    aggr_counts,
    colData = df,
    design = ~ Patient + Region
)

# Set to min(#sample in group 1, #sample in group 2).
keep <- rowSums(counts(deseq) >= 5) >= 12
deseq <- deseq[keep, ]

genes_used <- rownames(deseq)
genes_used <- data.frame(gene = genes_used)
write.csv(genes_used, file.path(params$output_dir, "Fig1f_genes_used.csv"))

# Run DE analysis
deseq <- DESeq(deseq)

pdf(
    file = file.path(params$plot_dir, "Fig1f_deseq_disp_ests.pdf"),
    width = 10,
    height = 10
)
plotDispEsts(deseq)
dev.off()

contrast <- c("Region", "PT", "Rest")
res <- results(deseq, contrast = contrast, alpha = 0.05)

res <- lfcShrink(deseq, contrast = contrast, type = "ashr", res = res)

res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()

# ---------------------------------- volcano --------------------------------- #
de_results <- as.data.frame(res_tbl)

# export results for Fig1g script
write.csv(de_results, file.path(params$output_dir, "Fig1f_de_results.csv"))

# Thresholds
fc <- log2(1.5)
p_valadj <- 0.05

de_results$diffexpressed <- "NO"
de_results$diffexpressed[
    de_results$log2FoldChange > fc & de_results$padj < p_valadj
] <- "UP"
de_results$diffexpressed[
    de_results$log2FoldChange < -fc & de_results$padj < p_valadj
] <- "DOWN"

rownames(de_results) <- de_results$gene

# Label genes
gene_lists <- readxl::read_excel(params$genelists_path, skip = 1)

# Venkataramani 2022
up <- gene_lists %>%
    filter(!is.na(Venkataramani_invasivity_up)) %>%
    pull(Venkataramani_invasivity_up)
custom_labs_up <- up[
    up %in% (de_results$gene[de_results$diffexpressed == "UP"])
]

# Hai 2024
down <- gene_lists %>%
    filter(!is.na(Hai_connectivity_up)) %>%
    pull(Hai_connectivity_up)

custom_labs_down <- down[
    down %in% (de_results$gene[de_results$diffexpressed == "DOWN"])
]

df <- de_results
df$diffexpressed <- case_when(
    df$diffexpressed == "NO" ~ "NO",
    rownames(df) %in% custom_labs_up ~ "Invasivity (Venkataramani et al.)",
    rownames(df) %in% custom_labs_down ~ "Connectivity (Hai et al.)",
    .default = "YES"
)

df$diffexpressed <- factor(
    df$diffexpressed,
    levels = c(
        "NO",
        "YES",
        "Invasivity (Venkataramani et al.)",
        "Connectivity (Hai et al.)"
    )
)
df$sort <- ifelse(
    df$diffexpressed %in%
        c("Invasivity (Venkataramani et al.)", "Connectivity (Hai et al.)"),
    1,
    0
)
df$alpha <- ifelse(
    df$diffexpressed %in%
        c("Invasivity (Venkataramani et al.)", "Connectivity (Hai et al.)"),
    1,
    0.05
)
df$label_1 <- case_when(
    df$diffexpressed == "Invasivity (Venkataramani et al.)" ~ rownames(df),
    TRUE ~ NA
)
df$label_2 <- case_when(
    df$gene %in% c("CHI3L1", "GAP43", "ANXA2") ~ rownames(df),
    TRUE ~ NA
)

df <- df %>% arrange(sort)
p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(
        aes(colour = diffexpressed, alpha = alpha),
        size = 3.5,
        stroke = NA
    ) +
    scale_colour_manual(values = c("black", "black", "darkorange", "red4")) +
    geom_label_repel(
        aes(label = label_1),
        size = 3,
        na.rm = TRUE,
        colour = "darkorange",
        nudge_x = 1,
        nudge_y = 1
    ) +
    geom_label_repel(
        aes(label = label_2),
        size = 3,
        na.rm = TRUE,
        colour = "red4",
        nudge_x = 0.1,
        nudge_y = 0.5
    ) +
    geom_vline(xintercept = fc, linetype = "dashed") +
    geom_vline(xintercept = -fc, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    (theme_foundation(base_size = 14, base_family = "Helvetica") +
        theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2),
                hjust = 0.5
            ),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1.1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = rel(1)),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size = unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face = "plain"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(
                colour = "#ffffff",
                fill = "#ffffff"
            ),
            strip.text = element_text(face = "bold", size = rel(1.1))
        )) +
    ylim(c(0, 17)) +
    xlim(c(-3, 3)) +
    xlab(bquote(~ Log[2] ~ FoldChange)) +
    ylab(bquote(~ -Log[10] ~ italic(FDR)))
ggsave(
    plot = p,
    filename = "Fig1f_volcano.pdf",
    path = params$plot_dir,
    height = 10,
    width = 7
)
