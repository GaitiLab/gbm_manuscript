# Plot genome browser tracks with Seurat
log_info("Creating genome browser tracks with Seurat...")
integrated_seurat <- readRDS(args$seurat_object)
metadata <- integrated_seurat[[]] %>%
  filter(Platform == "Multiome") %>%
  filter(CellClass_L1 == "Malignant") %>%
  filter(Confident_Annotation == "TRUE") %>%
  select(CellID, CellClass_L5_2)

# Load seurat object for visualization
curr_seurat_data <- readRDS("merged.rds")
curr_seurat_data <- subset(curr_seurat_data, subset = CellClass_L1 == "Malignant")

# Add metadata
curr_seurat_data <- subset(curr_seurat_data, subset = CellID %in% metadata$CellID) # Filter cells
metadata <- metadata[match(curr_seurat_data$CellID, metadata$CellID), ] # Match metadata rows
metadata <- metadata %>%
  mutate(CellClass_L5_2 = case_when(
    CellClass_L5_2 == "Malignant_NPC1" ~ "Malignant NPC1",
    CellClass_L5_2 == "Malignant_NPC2" ~ "Malignant NPC2",
    CellClass_L5_2 == "Malignant_OPC" ~ "Malignant OPC",
    CellClass_L5_2 == "Malignant_AC" ~ "Malignant AC",
    CellClass_L5_2 == "Malignant_MES_AST" ~ "Malignant MES-AST",
    CellClass_L5_2 == "Malignant_MES_INT" ~ "Malignant MES-INT",
    CellClass_L5_2 == "Malignant_MES_HYP" ~ "Malignant MES-HYP",
    CellClass_L5_2 == "Invasive-high OPC/NPC1" ~ "Invasive-high OPC/NPC1"
  ))
curr_seurat_data <- AddMetaData(curr_seurat_data, metadata$CellClass_L5_2, col.name = "CellClass_L5_2")

# Define custom color palette
color_palette <- c("Malignant NPC1" = "#86CCE8",
                  "Malignant NPC2" = "#ADD7E5",
                  "Malignant OPC" = "#AECFA5",
                  "Malignant AC" = "#DFB63B",
                  "Malignant MES-AST" = "#F37F72",
                  "Malignant MES-INT" = "#B02425",
                  "Malignant MES-HYP" = "#F1634B",
                  "Invasive-high OPC/NPC1" = "#2B7095")

# Set ident
Idents(curr_seurat_data) <- "CellClass_L5_2"
DefaultAssay(curr_seurat_data) <- "ATAC"

curr_seurat_data$CellClass_L5_2 <- factor(curr_seurat_data$CellClass_L5_2, levels = c("Invasive-high OPC/NPC1", "Malignant NPC1", "Malignant NPC2", "Malignant OPC", "Malignant AC", "Malignant MES-AST", "Malignant MES-INT", "Malignant MES-HYP"))

# CoveragePlot
cov_plot <- CoveragePlot(
  object = curr_seurat_data,
  region = "chr19-39497000-39509000",
  group.by = "CellClass_L5_2",
  annotation = FALSE,
  peaks = FALSE) +
  scale_fill_manual(values = color_palette)
ggsave(
  filename = paste0(args$output_dir, "/CoveragePlot_DLL3.pdf"),
  plot = cov_plot,
  width = 10,
  height = 10
)

# GenePlot
gene_plot <- AnnotationPlot(
  object = curr_seurat_data,
  region = c("chr19-39497000-39509000")
)
ggsave(
  filename = paste0(args$output_dir, "/AnnotationPlot_DLL3.pdf"),
  plot = gene_plot,
  width = 10,
  height = 10
)

# expression
expr_plot <- ExpressionPlot(
  object = curr_seurat_data,
  features = "DLL3",
  assay = "RNA",
  slot = "data") +
  scale_fill_manual(values = color_palette)
ggsave(
  filename = paste0(args$output_dir, "/ExpressionPlot_DLL3.pdf"),
  plot = expr_plot,
  width = 10,
  height = 10
)

combined <- CombineTracks(
  plotlist = list(cov_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(8, 1),
  widths = c(5, 1)
)
ggsave(
  filename = paste0(args$output_dir, "/Combined_DLL3.pdf"),
  plot = combined,
  width = 7,
  height = 5
)