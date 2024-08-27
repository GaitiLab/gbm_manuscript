# Code to reproduce Fig S4g

if(!("pacman" %in% rownames(installed.packages()))){
  install.packages("pacman")
}

pacman::p_load(data.table,
               ggplot2,
               dplyr,
               ggpubr,
               Seurat,
               reticulate,
               leiden,
               GBMutils,
               glmGamPoi)

region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

# Enter paths here
path_to_seurat_object <- ""
output_dir <- ""

# Creating directories needed for outputs
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_dir <- paste0(output_dir, "/path_to_output_directory")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---------------------------------------------------------------------------- #
#                                 Preprocessing                                #
# ---------------------------------------------------------------------------- #

path <- "/Krishna_2023" # path to Krishna 2023 data. Data can be downloaded on GEO

dirs <- list.dirs(path, full.names = TRUE, recursive = FALSE)

objs <- lapply(dirs, function(dir) {
    matrix <- Read10X(dir)
    so <- CreateSeuratObject(counts = matrix)

    so$orig.ident <- basename(dir)
    print(basename(dir))
    so$Patient <- paste0("Patient", substr(basename(dir), 4, 4))

    # Following preprocessing steps reported in Krishna 2023
    so[["mito_pct"]] <- PercentageFeatureSet(so, pattern = "^MT-")
    so <- subset(so, subset = (nFeature_RNA > 500 & nFeature_RNA < 10000))
    so <- subset(so, subset = mito_pct < 20)

    so <- SCTransform(so, vars.to.regress = c("nCount_RNA", "mito_pct"), vst.flavor = "v2", method = "glmGamPoi")

    print(so)
    return (so)
})

combined_obj <- merge(objs[[1]], objs[-1])
print(combined_obj)

saveRDS(combined_obj, "/krishna2023.rds")

# Take this object, run infercnv using oligo from our data set to annotate malignant cells. Then, subtype and perform DE, gene set scoring

so <- readRDS("/krishna2023.rds")

so$CNV <- case_when(
  so$Patient == "Patient3" & so$has_dupli_chr5 == TRUE ~ "Malignant",
  so$Patient == "Patient3" & so$has_dupli_chr5 == FALSE ~ "Non-malignant",
  so$has_dupli_chr7 == TRUE & so$has_loss_chr10 ~ "Malignant",
  .default = "Non-malignant"
)

DefaultAssay(so) <- "RNA"
so <- NormalizeData(so)
so <- FindVariableFeatures(so, nfeatures = 3000, assay = "RNA")
so <- ScaleData(so)
so <- RunPCA(so, assay = "RNA")
so <- FindNeighbors(so, reduction = "pca", dims = 1:30)
so <- RunUMAP(so, dims = 1:30)
p <- DimPlot(so, reduction = "umap", group.by = "CNV")
ggsave(filename = "cnv_umap.pdf", path = plot_dir)

# Neftel subtyping and score invasive sig

neftel_markers <- read.csv("/Neftel_gene_list.csv")

gene_list <- list()
for(i in 1:nrow(neftel_markers)) {
  gene_list[[toString(neftel_markers[i,"Cell"])]] <- append(gene_list[[toString(neftel_markers[i,"Cell"])]],toString(neftel_markers[i, "Gene"]))
}

cell_type <- c()
for(i in 1:length(gene_list)){
  cell_type <- c(cell_type, paste0("Neftel_", names(gene_list)[i]))
}
names(gene_list) <- cell_type

# Invasive signature
degs <- read.csv("/de_results.csv")
degs_up <- degs %>% 
  filter(log2FoldChange > 1 & padj < 0.05)
degs_dn <- degs %>% 
  filter(log2FoldChange < -1 & padj < 0.05)

gene_sets <- c(gene_list, list(degs_up = degs_up$gene, degs_dn = degs_dn$gene))

so$cell_id <- rownames(so[[]])

dfs <- lapply(unique(so$orig.ident), function(sample) {
  curr_so <- subset(so, subset = orig.ident == sample)

  curr_so <- AddModuleScore(curr_so, features = gene_sets, name = "Geneset", assay = "SCT")
  curr_df <- curr_so[[]] %>% 
    select(all_of(c(paste0("Geneset", seq(1,8)), "cell_id")))
  colnames(curr_df) <- c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2", "Inv_up", "Inv_dn", "cell_id")

  neftel_df <- curr_df[,seq(1,6)]
  cell_subtypes <- ifelse(apply(X = neftel_df, MARGIN = 1, FUN = max, na.rm = T) > 0, colnames(neftel_df)[max.col(neftel_df, ties.method = "first")], "Undetermined")
  curr_df$Neftel_subtype <- cell_subtypes
  return(curr_df)
})

df <- rbindlist(dfs)
df <- as.data.frame(df)
rownames(df) <- df$cell_id
df <- df %>% select(-cell_id)

print(head(df))

so <- AddMetaData(so, metadata = df)

print(table(so$orig.ident, so$Neftel_subtype))

so <- subset(so, subset = Neftel_subtype %in% c("OPC", "NPC1"))
so$Invasive_sig <- so$Inv_up - so$Inv_dn
so$Type <- substr(so$orig.ident, 1, 3)

print(table(so$orig.ident, so$Neftel_subtype))

df <- so[[]] %>% filter(Patient != "Patient3") # Removed due to low number of OPC/NPC1-like cells (only 2 in HFC sample)

df$orig.ident <- factor(df$orig.ident, levels = c("HFC1", "LFC1", "HFC2", "LFC2"))
p <- ggplot(df, aes(x = orig.ident, y = Invasive_sig)) +
  geom_boxplot(aes(fill = Type)) +
  scale_fill_manual(values = c("#336699", "#86BBD8")) +
  GBM_theme() +
  stat_compare_means(comparisons = list(c("HFC1", "LFC1"), c("HFC2", "LFC2"))) +
  xlab("Sample") +
  ylab("Invasive Signature Score") +
  ggtitle("Invasive Signature in HFC vs LFC OPC/NPC1-like cells (Krishna2023)")
ggsave(filename = "invasive_sig_scored_opc_npc1.pdf", path = plot_dir, width = 12, height = 12)

