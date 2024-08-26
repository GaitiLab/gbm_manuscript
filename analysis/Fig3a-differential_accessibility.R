params <- list(
    input = "/multiome_results/10_ArchR_new",
    archr_threads = 8,
    genome_version = "hg38",
    celltype_column = "CellClass_L5_2",
    patient_column = "Patient",
    region_column = "Region",
    confidence_column = "Confident_Annotation",
    p_cutoff = 0.1,
    include_types = "Malignant_NPC1,Malignant_NPC2,Malignant_OPC,Invasive-high_OPC_NPC1,Malignant_AC,Malignant_MES_AST,Malignant_MES_INT,Malignant_MES_HYP",
    exclude_types = NULL,
    include_regions = "PT,TE",
    topVarGenes = 10000,
    folder_name = "CellClass_L5_2_archr_top10000"
)

# Load the ArchR library
addArchRThreads(threads = as.numeric(params$archr_threads))
addArchRGenome(params$genome_version)

# Load the ArchR project
archr_proj <- loadArchRProject(params$input)
cell_type_column <- paste0("Seurat_", params$celltype_column)
patient_column <- paste0("Seurat_", params$patient_column)
region_column <- paste0("Seurat_", params$region_column)
confidence_column <- paste0("Seurat_", params$confidence_column)

# Filter cells based on cell types, regions, and confidence
if (!is.null(params$include_types)) {
    include_types <- str_split(params$include_types, ",")[[1]]
    archr_proj <- archr_proj[which(getCellColData(archr_proj, select = cell_type_column, drop = TRUE) %in% include_types)]
} else if (!is.null(params$exclude_types)) {
    exclude_types <- str_split(params$exclude_types, ",")[[1]]
    archr_proj <- archr_proj[which(!(getCellColData(archr_proj, select = cell_type_column, drop = TRUE) %in% exclude_types))]
}

archr_proj <- archr_proj[which(getCellColData(archr_proj, select = confidence_column, drop = TRUE) == "TRUE")]

if (!is.null(params$include_regions)) {
  include_regions <- str_split(params$include_regions, ",")[[1]]
  archr_proj <- archr_proj[which(getCellColData(archr_proj, select = region_column, drop = TRUE) %in% include_regions)]
}

# Define groups
df <- getCellColData(archr_proj, select = c(cell_type_column, patient_column, region_column), drop = TRUE)
low <- c("Malignant_NPC1", "Malignant_NPC2", "Malignant_OPC")
high <- "Invasive-high_OPC_NPC1"

# Differential accessibility analysis
archr_proj_low_and_high <- archr_proj[which(getCellColData(archr_proj, select = cell_type_column, drop = TRUE) %in% c(low, high))]
archr_proj_low <- archr_proj_low_and_high[which(getCellColData(archr_proj_low_and_high, select = cell_type_column, drop = TRUE) %in% low)]
archr_proj_high <- archr_proj_low_and_high[which(getCellColData(archr_proj_low_and_high, select = cell_type_column, drop = TRUE) %in% high)]

curr_results <- getMarkerFeatures(ArchRProj = archr_proj_low_and_high,
                                  useMatrix = "GeneScoreMatrix",
                                  groupBy = cell_type_column,
                                  useGroups = high,
                                  bgdGroups = low,
                                  bias = c("TSSEnrichment", "log10(nFrags)"),
                                  testMethod = "wilcoxon",
                                  maxCells = 100000,
                                  logFile = createLogFile(name = "getMarkerFeatures", logDir = log_dir))

marker_list <- data.table(name = curr_results@elementMetadata@listData$name,
                          chr = as.character(curr_results@elementMetadata@listData$seqnames),
                          start = curr_results@elementMetadata@listData$start,
                          end = curr_results@elementMetadata@listData$end,
                          log2FC = curr_results@assays@data$Log2FC[[1]],
                          Mean = curr_results@assays@data$Mean[[1]],
                          FDR = curr_results@assays@data$FDR[[1]],
                          Pval = curr_results@assays@data$Pval[[1]],
                          MeanDiff = curr_results@assays@data$MeanDiff[[1]],
                          AUC = curr_results@assays@data$AUC[[1]],
                          MeanBGD = curr_results@assays@data$MeanBGD[[1]])

# Get highly variable genes
curr_cells <- archr_proj_low_and_high$cellNames
gene_variances <- rowIQRs(gene_score_mat[,curr_cells], useNames = TRUE)
gene_variances <- gene_variances[order(gene_variances, decreasing = TRUE)]

top_genes <- names(gene_variances)[1:params$topVarGenes]
marker_list <- marker_list[!is.na(FDR),]
marker_list <- marker_list[name %in% top_genes,]
marker_list$FDR <- p.adjust(marker_list$Pval, method = "fdr")

marker_list <- marker_list %>% mutate(
    comparison = paste0(high, "_vs_All"),
    cell_type = paste0(high, "_vs_All"),
)

# Volcano plot for differentially accessible genes and NOTCH signaling pathway genes
marker_list$label <- ifelse(marker_list$FDR < 0.05 & abs(as.numeric(marker_list$log2FC)) > 0.38, "adj_p<0.05 & log2FC>0.38", "adj_p>=0.05 or log2FC<0.38")
marker_list <- marker_list[!is.na(FDR),]

# Specify the log2FC and FDR cutoffs
curr.log2FC <- 0.38
curr.fdr <-0.05

# Align significant and insignificant genes
marker_list <- marker_list %>%
    mutate(diffaccessible = case_when(
    log2FC > curr.log2FC & FDR <= curr.fdr ~ "Accessible in Invasive-high OPC/NPC1",
    log2FC < -curr.log2FC & FDR <= curr.fdr ~ "Accessible in Progenitor-like (NPC + OPC)",
    TRUE ~ "Not significant"
    ))

# change the order of the factor
marker_list$diffaccessible <- factor(marker_list$diffaccessible, levels = c("Accessible in Invasive-high OPC/NPC1", "Not significant", "Accessible in Progenitor-like (NPC + OPC)"))

# Generate label for plot annotation
marker_list$dalabel <- NA
marker_list$dalabel[marker_list$diffaccessible != "Not significant"] <- marker_list$name[marker_list$diffaccessible != "Not significant"]

# Overlay the volcano plot with the NOTCH signaling pathway genes
specific_gene <- all_gene_sets %>%
                    filter(gs_name == "WP_NOTCH_SIGNALING") %>%
                    pull(gene_symbol)
                    marker_list <- marker_list %>%
                    mutate(highlight = case_when(
                        name %in% specific_gene & diffaccessible == "Accessible in Invasive-high OPC/NPC1" ~ "yes",
                        TRUE ~ "no"
                    ))

# Create the volcano plot with highlighted genes
marker_list %>% filter(highlight == "yes")
volc <- ggplot(data = marker_list, aes(x = log2FC, y = -log10(FDR), label = name)) +
    geom_point(aes(color = highlight), size = 1, alpha = 0.25, stroke = NA) +
    geom_point(data = marker_list[marker_list$highlight != "no", ], aes(color = highlight), size = 1.5) +
    geom_vline(xintercept = curr.log2FC, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = -curr.log2FC, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(curr.fdr), linetype = "dashed", alpha = 0.5) +
    expand_limits(x = c(-ceiling(max(abs(marker_list$log2FC))), ceiling(max(abs(marker_list$log2FC))))) +
    geom_label_repel(data = marker_list[marker_list$highlight != "no", ],
                    fill = "white",
                    label.padding = 0.1,
                    box.padding = 0.5,
                    size = 3,
                    min.segment.length = 0,
                    aes(color = highlight),
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    force = 2, 
                    nudge_x = 0.3,
                    nudge_y = 0.3) +
    scale_color_manual(values = color_palette, guide = guide_legend(title = "")) +
    guides(color = guide_legend(title = "")) +
    labs(
        x = expression(Log[2] * " Fold Change"),
        y = expression(-Log[10] * " FDR")) +
    theme_minimal() +
    theme(legend.key = element_blank(),
        legend.position = "bottom",
        plot.title.position = "plot",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = paste0("DiffAccess_", "all", "_low_vs_high_highlight_specific_gene.pdf"),
        path = out_dir,
        units = "in", width = 5, height = 6)