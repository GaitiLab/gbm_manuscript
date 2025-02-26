# ---- Code to reproduce Figure S4g ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #
pacman::p_load(
    data.table,
    ggplot2,
    dplyr,
    ggpubr,
    Seurat,
    reticulate,
    leiden,
    GBMutils,
    glmGamPoi
)

region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

# Enter paths here
params <- list(
    seurat_obj_path = "",
    output_dir = "output/submission",
    plot_dir = "output/submission/figures",
    # path to Krishna 2023 data. Data can be downloaded on GEO
    krishna2023_dir_path = "/Krishna_2023",
    path_to_genelists = "misc/SuppTables/Table S2.xlsx",
    path_to_de_genes_table = "misc/SuppTables/Table S3.xlsx"
)

# Creating directories needed for outputs
GaitiLabUtils::create_dir(params$output_dir)
GaitiLabUtils::create_dir(params$plot_dir)


# ---------------------------------------------------------------------------- #
#                                 Preprocessing                                #
# ---------------------------------------------------------------------------- #

dirs <- list.dirs(
    params$krishna2023_dir_path,
    full.names = TRUE,
    recursive = FALSE
)

objs <- lapply(dirs, function(dir) {
    matrix <- Read10X(dir)
    seurat_obj <- CreateSeuratObject(counts = matrix)

    seurat_obj$orig.ident <- basename(dir)
    print(basename(dir))
    seurat_obj$Patient <- paste0("Patient", substr(basename(dir), 4, 4))

    # Following preprocessing steps reported in Krishna 2023
    seurat_obj[["mito_pct"]] <- PercentageFeatureSet(
        seurat_obj,
        pattern = "^MT-"
    )
    seurat_obj <- subset(
        seurat_obj,
        subset = (nFeature_RNA > 500 & nFeature_RNA < 10000)
    )
    seurat_obj <- subset(seurat_obj, subset = mito_pct < 20)

    seurat_obj <- SCTransform(
        seurat_obj,
        vars.to.regress = c("nCount_RNA", "mito_pct"),
        vst.flavor = "v2",
        method = "glmGamPoi"
    )

    print(seurat_obj)
    return(seurat_obj)
})

combined_obj <- merge(objs[[1]], objs[-1])
print(combined_obj)

saveRDS(combined_obj, file.path(params$output_dir, "FigS4g_krishna2023.rds"))

# Take this object, run infercnv using oligo from our data set to annotate malignant cells. Then, subtype and perform DE, gene set scoring

seurat_obj <- readRDS(file.path(params$output_dir, "FigS4g_krishna2023.rds"))

seurat_obj$CNV <- case_when(
    seurat_obj$Patient == "Patient3" & seurat_obj$has_dupli_chr5 == TRUE ~
        "Malignant",
    seurat_obj$Patient == "Patient3" & seurat_obj$has_dupli_chr5 == FALSE ~
        "Non-malignant",
    seurat_obj$has_dupli_chr7 == TRUE & seurat_obj$has_loss_chr10 ~ "Malignant",
    .default = "Non-malignant"
)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000, assay = "RNA")
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, assay = "RNA")
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
p <- DimPlot(seurat_obj, reduction = "umap", group.by = "CNV")
ggsave(filename = "cnv_umap.pdf", path = plot_dir)

# Neftel subtyping and score invasive sig

gene_lists <- readxl::read_excel(params$path_to_genelists, skip = 1)
gene_list <- lapply(
    gene_lists %>% dplyr::select(starts_with("Neftel_")) %>% as.list(),
    function(gene_list) {
        gene_list[!is.na(gene_list)]
    }
)

# Invasive signature
# TODO @Bensonwu02 is this correct
# degs <- read.csv("/de_results.csv")
# degs_up <- degs %>%
#     filter(log2FoldChange > 1 & padj < 0.05)
# degs_dn <- degs %>%
#     filter(log2FoldChange < -1 & padj < 0.05)

degs <- readxl::read_excel(
    params$path_to_de_genes_table,
    sheet = "DEGs",
    skip = 1
) %>%
    data.frame()
degs_up <- degs %>%
    filter(Direction == "Upregulated in PT OPC/NPC1-like cells")
degs_dn <- degs %>%
    filter(
        Direction == "Upregulated in tumor bulk (TE+TC) OPC/NPC1-like cells"
    )

gene_sets <- c(gene_list, list(degs_up = degs_up$gene, degs_dn = degs_dn$gene))

seurat_obj$cell_id <- rownames(seurat_obj[[]])

dfs <- lapply(unique(seurat_obj$orig.ident), function(sample) {
    curr_so <- subset(seurat_obj, subset = orig.ident == sample)

    curr_so <- AddModuleScore(
        curr_so,
        features = gene_sets,
        name = "Geneset",
        assay = "SCT"
    )
    curr_df <- curr_so[[]] %>%
        select(all_of(c(paste0("Geneset", seq(1, 8)), "cell_id")))
    colnames(curr_df) <- c(
        "MES2",
        "MES1",
        "AC",
        "OPC",
        "NPC1",
        "NPC2",
        "Inv_up",
        "Inv_dn",
        "cell_id"
    )

    neftel_df <- curr_df[, seq(1, 6)]
    cell_subtypes <- ifelse(
        apply(X = neftel_df, MARGIN = 1, FUN = max, na.rm = T) > 0,
        colnames(neftel_df)[max.col(neftel_df, ties.method = "first")],
        "Undetermined"
    )
    curr_df$Neftel_subtype <- cell_subtypes
    return(curr_df)
})

df <- rbindlist(dfs)
df <- as.data.frame(df)
rownames(df) <- df$cell_id
df <- df %>% select(-cell_id)

print(head(df))

seurat_obj <- AddMetaData(seurat_obj, metadata = df)

print(table(seurat_obj$orig.ident, seurat_obj$Neftel_subtype))

seurat_obj <- subset(seurat_obj, subset = Neftel_subtype %in% c("OPC", "NPC1"))
seurat_obj$Invasive_sig <- seurat_obj$Inv_up - seurat_obj$Inv_dn
seurat_obj$Type <- substr(seurat_obj$orig.ident, 1, 3)

print(table(seurat_obj$orig.ident, seurat_obj$Neftel_subtype))

df <- seurat_obj[[]] %>% filter(Patient != "Patient3") # Removed due to low number of OPC/NPC1-like cells (only 2 in HFC sample)

df$orig.ident <- factor(
    df$orig.ident,
    levels = c("HFC1", "LFC1", "HFC2", "LFC2")
)
p <- ggplot(df, aes(x = orig.ident, y = Invasive_sig)) +
    geom_boxplot(aes(fill = Type)) +
    scale_fill_manual(values = c("#336699", "#86BBD8")) +
    GBM_theme() +
    stat_compare_means(
        comparisons = list(c("HFC1", "LFC1"), c("HFC2", "LFC2"))
    ) +
    xlab("Sample") +
    ylab("Invasive Signature Score") +
    ggtitle(
        "Invasive Signature in HFC vs LFC OPC/NPC1-like cells (Krishna2023)"
    )
ggsave(
    filename = "FigS4g_invasive_sig_scored_opc_npc1.pdf",
    path = params$plot_dir,
    width = 12,
    height = 12
)
