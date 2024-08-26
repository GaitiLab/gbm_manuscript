# Code to reproduce Fig 1i

if (!("pacman" %in% rownames(installed.packages()))) {
  install.packages("pacman")
}

pacman::p_load(data.table,
               patchwork,
               ggplot2,
               dplyr,
               tidyr,
               tibble,
               Seurat,
               log4r,
               DESeq2,
               SingleCellExperiment,
               apeglm,
               Matrix.utils)

# Enter paths here
path_to_seurat_object <- ""
output_dir <- ""

# Creating directories needed for outputs
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_dir <- paste0(output_dir, "/path_to_output_directory")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ----------------------------- Setting up logger ---------------------------- #
log.file <- paste0(args$output, "/log.txt")
console.appender <- console_appender(layout = default_log_layout())
file.appender <- file_appender(log.file, append = TRUE, 
                               layout = default_log_layout())
logr <- log4r::logger(threshold = 1, 
                      appenders = list(console.appender, file.appender))

object.appender <- file_appender(log.file, append = TRUE, 
                                 layout = bare_log_layout())
obj.logr <- log4r::logger(threshold = 1,
                          appenders = list(object.appender))

log_info <- function(...) {
  log4r::info(logr, paste0(...))
}

# Load object
log_info("Loading Seurat object")
seurat_object <- readRDS(path_to_seurat_object)

# Prepare groups for comparison
seurat_object$Region <- ifelse(seurat_object$Region == "PT", "PT", "Rest") # Sample level region comparisons
seurat_object <- subset(seurat_object, subset = is_malignant_confident == TRUE) # high confidence malignant cells
seurat_object <- subset(seurat_object, subset = Patient %in% c("6237", "6245", "6467", "6419")) # Patients with matched PT and Tumor samples that have malignant cells
seurat_object <- subset(seurat_object, subset = CellClass_L3 %in% c("Malignant_OPC", "Malignant_NPC1")) # Subset for cell types of interest

# Convert to SingleCellExperiment
log_info("Converting to SingleCellExperiment")
counts <- seurat_object@assays$RNA@counts

metadata <- seurat_object[[]] 
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Aggregate counts 
log_info("Aggregating counts by group of interest")
groups <- colData(sce)[, c("Sample")] # Group single cell data by sample to perform pseudobulking
aggr_counts <- aggregate.Matrix(t(counts(sce)),
                                groupings = groups, fun = "sum")
aggr_counts <- t(aggr_counts)

# Prepare metadata
log_info("Preparing metadata")
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  select(c("Sample", "Patient", "Region")) 

# Remove duplicates
metadata <- metadata[!duplicated(metadata), ]
# Rename rows
rownames(metadata) <- metadata$Sample
# Cell counts
t <- table(colData(sce)$Sample)
df <- as.data.frame(t)
colnames(df) <- c("Sample", "cell_count")

## Join data frames
df <- plyr::join(df, metadata, 
                 by = "Sample")

## Update rownames of metadata to match colnames of count matrix, as needed later for DE
index <- match(colnames(aggr_counts), df$Sample)
df <- df[index, , drop = FALSE]

write.csv(df, file.path(plot_dir, "samples_tested.csv"))

# Converting col metadata to factors
df$Patient <- factor(df$Patient, levels = unique(df$Patient))
df$Region <- factor(df$Region, levels = unique(df$Region))

# Create DESeq2 object
log_info("Creating DESeq2 object")
deseq <- DESeqDataSetFromMatrix(aggr_counts,
                                colData = df,
                                design = ~ Patient + Region)

log_info("Filtering low count genes")
keep <- rowSums(counts(deseq) >= 5) >= 6 # Set to min(#sample in group 1, #sample in group 2). 
deseq <- deseq[keep, ]

genes_used <- rownames(deseq)
genes_used <- data.frame(gene = genes_used)
write.csv(genes_used, file.path(plot_dir, "genes_used.csv"))

log_info("Running DESeq2")
# Run DE analysis
deseq <- DESeq(deseq)

pdf(file = paste0(plot_dir, "/deseq_disp_ests.pdf"), width = 10, height = 10)
plotDispEsts(deseq)
dev.off()

contrast <- c("Region", "PT", "Rest")
res <- results(deseq, 
  contrast = contrast,
  alpha = 0.05)

res <- lfcShrink(deseq,
  contrast = contrast,
  type = "ashr",
  res = res)

res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()

# ---------------------------------- volcano --------------------------------- #
de_results <- as.data.frame(res_tbl)

# Thresholds
fc <- log2(2)
p_valadj <- 0.05

de_results$diffexpressed <- "NO"
de_results$diffexpressed[de_results$avg_logFC > fc & de_results$FDR < p_valadj] <- "UP"
de_results$diffexpressed[de_results$avg_logFC < -fc & de_results$FDR < p_valadj] <- "DOWN"
de_results$de_label <- NA
de_results$de_label[de_results$diffexpressed != "NO"] <- de_results$Gene[de_results$diffexpressed != "NO"]

rownames(de_results) <- de_results$Gene

# Label genes
up <- read.csv("gene_lists/invasivity.csv") # Venkataramani 2022
up <- up %>% filter(direction == "Anticorrelated")
up <- up$Gene
custom_labs_up <- up[up %in% (de_results$Gene[de_results$diffexpressed == "UP"])]

down <- read.csv("gene_lists/hai_connectivity.csv") # Hai 2024
down <- down %>% filter(direction == "Up")
down <- down$Gene
custom_labs_down <- down[down %in% (de_results$Gene[de_results$diffexpressed == "DOWN"])]

custom_labs <- c(custom_labs_up)

df <- de_results
df$diffexpressed <- case_when(
  df$diffexpressed == "NO" ~ "NO",
  rownames(df) %in% custom_labs_up ~ "Invasivity (Venkataramani et al.)",
  rownames(df) %in% custom_labs_down ~ "Connectivity (Hai et al.)",
  .default = "YES"
)

df$diffexpressed <- factor(df$diffexpressed, levels = c("NO", "YES", "Invasivity (Venkataramani et al.)", "Connectivity (Hai et al.)"))
df$sort <- ifelse(df$diffexpressed %in% c("Invasivity (Venkataramani et al.)", "Connectivity (Hai et al.)"), 1, 0)
df$alpha <- ifelse(df$diffexpressed %in% c("Invasivity (Venkataramani et al.)", "Connectivity (Hai et al.)"), 1, 0.05)
df$label <- ifelse(df$diffexpressed == "Invasivity (Venkataramani et al.)", rownames(df), NA)

df <- df %>% arrange(sort)
p <- ggplot(df, aes(x = avg_logFC, y = -log10(FDR))) +
  geom_point(aes(colour = diffexpressed, alpha = alpha), size = 3.5, stroke = NA) +
  scale_colour_manual(values = c("black", "black", "darkorange", "red4")) +
  geom_label_repel(aes(label = label), size = 3, colour = "darkorange", na.rm = TRUE, nudge_x = 2, nudge_y = 1) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  (theme_foundation(base_size = 14, base_family = "Helvetica") 
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
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
            strip.background = element_rect(colour = "#ffffff", fill = "#ffffff"),
            strip.text = element_text(face = "bold", size = rel(1.1))
        )) +
  ylim(c(0,25)) +
  xlim(c(-4,4)) +
  xlab(bquote(~Log[2] ~ FoldChange)) +
  ylab(bquote(~-Log[10] ~ italic(FDR)))
ggsave(filename = "volcano.pdf", path = plot_dir, height = 10, width = 7)