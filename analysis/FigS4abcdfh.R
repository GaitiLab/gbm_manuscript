# ---- Code to reproduce Figure S4 a,b,c,d,f,h ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# Load needed packages
pacman::p_load(
    data.table,
    ggplot2,
    Seurat,
    tidyr,
    dplyr,
    GBMutils,
    annotables,
    VennDiagram,
    GaitiLabUtils
)

region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

# Enter paths here
params <- list(
    seurat_obj_path = "",
    output_dir = "output/submission",
    plot_dir = "output/submission/figures",
    fgseares_c2_cp_path = "misc/Table S2.xlsx",
    degs_table_path = "misc/Table S2.xlsx",
    genelists_path = "misc/gene_signatures.xlsx",
    gbmap_core_obj_path = "gbmap_core.rds"
)

GaitiLabUtils::create_dir(params$output_dir)
GaitiLabUtils::create_dir(params$plot_dir)

# ---------------------------------------------------------------------------- #
#                                  Figure S4a                                  #
# ---------------------------------------------------------------------------- #

res <- readxl::read_excel(
    params$fgseares_c2_cp_path,
    sheet = "C2_CP pathways in inv high",
    skip = 1
    )
    
pathways <- c(
    "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
    "REACTOME_NEURONAL_SYSTEM",
    "KEGG_NOTCH_SIGNALING_PATHWAY",
    "WP_NOTCH_SIGNALING_PATHWAY",
    "REACTOME_SODIUM_CALCIUM_EXCHANGERS",
    "REACTOME_ION_CHANNEL_TRANSPORT",
    "WP_FATTY_ACID_BIOSYNTHESIS",
    "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES",
    "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
    "REACTOME_SIGNALING_BY_INTERLEUKINS",
    "WP_VEGFA_VEGFR2_SIGNALING",
    "WP_COMPLEMENT_SYSTEM",
    "PID_HIF1_TFPATHWAY",
    "REACTOME_INNATE_IMMUNE_SYSTEM",
    "REACTOME_GLYCOLYSIS"
)

plot_df <- res %>% filter(pathway %in% pathways)

plot_df$pathway <- case_when(
    plot_df$pathway == "REACTOME_CHOLESTEROL_BIOSYNTHESIS" ~
        "Reactome \n Cholesterol Biosynthesis",
    plot_df$pathway == "REACTOME_NEURONAL_SYSTEM" ~
        "Reactome \n Neuronal System",
    plot_df$pathway == "REACTOME_TRANSMISSION_ACROSS_CHEMICAL_SYNAPSES" ~
        "Reactome \n Transmission Across Chemical Synapses",
    plot_df$pathway == "KEGG_NOTCH_SIGNALING_PATHWAY" ~
        "KEGG \n Notch Signaling Pathway",
    plot_df$pathway ==
        "KEGG_MEDICUS_REFERENCE_CA2_ENTRY_VOLTAGE_GATED_CA2_CHANNEL" ~
        "KEGG \n Ca2+ Entry Voltage Gated Ca2+ Channel",
    plot_df$pathway == "REACTOME_ION_CHANNEL_TRANSPORT" ~
        "Reactome \n Ion Channel Transport",
    plot_df$pathway == "REACTOME_SODIUM_CALCIUM_EXCHANGERS" ~
        "Reactome \n Sodium Calcium Exchangers",
    plot_df$pathway == "WP_NOTCH_SIGNALING_PATHWAY" ~
        "WP \n Notch Signaling Pathway",
    plot_df$pathway == "WP_FATTY_ACID_BIOSYNTHESIS" ~
        "WP \n Fatty Acid Biosynthesis",
    plot_df$pathway == "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM" ~
        "Reactome \n Cytokine Signaling in Immune System",
    plot_df$pathway == "REACTOME_SIGNALING_BY_INTERLEUKINS" ~
        "Reactome \n Signaling by Interleukins",
    plot_df$pathway == "WP_VEGFA_VEGFR2_SIGNALING" ~
        "WP \n VEGFA VEGFR2 Signaling",
    plot_df$pathway == "WP_COMPLEMENT_SYSTEM" ~ "WP \n Complement System",
    plot_df$pathway == "PID_HIF1_TFPATHWAY" ~ "PID \n HIF1 TF Pathway",
    plot_df$pathway == "REACTOME_INNATE_IMMUNE_SYSTEM" ~
        "Reactome \n Innate Immune System",
    plot_df$pathway == "REACTOME_GLYCOLYSIS" ~ "Reactome \n Glycolysis"
)

p <- ggplot(plot_df, aes(x = NES, y = reorder(pathway, NES))) +
    geom_bar(aes(fill = NES), stat = "identity", alpha = 0.8) +
    scale_fill_gradient2(low = "#0571B0", high = "#CA0020", mid = "#F7F7F7") +
    GBM_theme() +
    labs(
        title = "C2:CP Pathways Enriched in PT OPC/NPC1-like",
        x = "NES",
        y = "",
        colour = "-log10(FDR)"
    ) +
    theme(
        legend.position = "right",
        legend.direction = "vertical",
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.key.size = unit(1, "cm")
    )

ggsave(
    filename = "FigS4a_c2_cp_select_pathways.pdf",
    path = params$plot_dir,
    height = 12,
    width = 12
)

# ---------------------------------------------------------------------------- #
#                                  Figure S4b                                  #
# ---------------------------------------------------------------------------- #
# Load object
seurat_obj <- readRDS(params$seurat_obj_path)

seurat_obj <- subset(
    seurat_obj,
    subset = Patient %in% c("6237", "6245", "6419", "6467", "6514", "6425")
)
seurat_obj <- subset(
    seurat_obj,
    subset = CellClass_L5_2 %in%
        c(
            "Invasive-high OPC/NPC1",
            "Invasive-low OPC/NPC1",
            "Malignant_OPC",
            "Malignant_NPC1"
        )
)
seurat_obj$Region <- ifelse(seurat_obj$Region == "PT", "PT", "Tumour")
seurat_obj$Bulk_group <- paste0(seurat_obj$Patient, "_", seurat_obj$Region)
genes <- c(degs_up$gene, degs_dn$gene)

degs <- as.data.frame(AverageExpression(
    seurat_obj,
    features = genes,
    assays = "RNA",
    group.by = "Bulk_group"
))

degs <- t(degs)
degs <- scale(degs)
degs <- t(degs)

degs <- as.matrix(degs)

sample_names <- colnames(degs)
sample_cat <- ifelse(grepl("PT", sample_names), "PT", "Tumour")
sample_cat_factor <- factor(sample_cat, levels = c("Tumour", "PT"))

patient_names <- sapply(colnames(degs), function(x) {
    return(substr(x, 5, 8))
})
patient_names_factor <- factor(patient_names)
patient_cols <- rainbowJam(n = length(unique(patient_names)))
names(patient_cols) <- levels(patient_names_factor)

gene_cat <- ifelse(rownames(degs) %in% degs_up$gene, "Up in PT", "Up in Tumour")
gene_cat_factor <- factor(gene_cat)

col_annotation <- HeatmapAnnotation(
    df = data.frame(
        Region = sample_cat[order(sample_cat_factor)],
        Patient = patient_names[order(sample_cat_factor)]
    ),
    col = list(
        Region = c("PT" = "#0173b2", "Tumour" = "#029e73"),
        Patient = patient_cols
    ),
    which = "col",
    annotation_legend_param = list(
        title_gp = gpar(fontsize = 14),
        labels_gp = gpar(fontsize = 12)
    )
)

color.scheme <- rev(brewer.pal(9, "RdBu"))
similarity_colors <- colorRamp2(
    c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
    color.scheme
)
sorted_indices <- order(sample_cat_factor)
degs_sorted <- degs[, sorted_indices]

genes_to_label <- c(
    "NEU4",
    "DLL3",
    "XYLT1",
    "LRRN1",
    "BCAN",
    "HES5",
    "HES6",
    "EGFR",
    "TNK2",
    "CSPG5",
    "SEZ6L",
    "MFNG",
    "HES5",
    "GRIK1",
    "VEGFA",
    "CD44",
    "FOSL2",
    "C3",
    "GLIS3",
    "FN1",
    "CHI3L1",
    "OSMR",
    "GAP43",
    "ADM"
)
index <- which(rownames(degs_sorted) %in% genes_to_label)
gene_labs <- rowAnnotation(
    foo = anno_mark(at = index, labels = genes_to_label, lines_gp = gpar())
)

ht <- Heatmap(
    degs_sorted,
    row_split = gene_cat,
    column_split = factor(
        sample_cat[sorted_indices],
        levels = c("Tumour", "PT")
    ),
    cluster_column_slices = FALSE,
    name = "Scaled Expression",
    width = ncol(degs_sorted) * unit(1, "cm"),
    height = nrow(degs_sorted) * unit(0.05, "cm"),
    top_annotation = col_annotation,
    right_annotation = gene_labs,
    column_title = "Invasive Signature",
    row_title = "Genes",
    column_names_gp = gpar(fontsize = 14),
    row_title_gp = gpar(fontsize = 14),
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = similarity_colors,
    heatmap_legend_param = list(
        title_gp = gpar(fontsize = 14),
        labels_gp = gpar(fontsize = 12)
    )
)

pdf(
    file = file.path(params$plot_dir, "FigS4b_invasive_sig_heatmap.pdf"),
    width = 15,
    height = 15
)
pushViewport(viewport(gp = gpar(fontfamily = "Helvetica")))
draw(ht, newpage = FALSE)
popViewport()
dev.off()

# ---------------------------------------------------------------------------- #
#                                  Figure S4c                                  #
# ---------------------------------------------------------------------------- #
degs <- readxl::read_excel(
    params$degs_table_path,
    sheet = "DEGs",
    skip = 1
) %>%
    data.frame()
degs_up <- degs %>%
    filter(Direction == "Upregulated in PT OPC/NPC1-like cells")

gene_lists <- readxl::read_excel(params$genelists_path, skip = 1)

neftel_gene_list <- lapply(
    gene_lists %>% dplyr::select(starts_with("Neftel_")) %>% as.list(),
    function(gene_list) {
        gene_list[!is.na(gene_list)]
    }
)

opc <- neftel_gene_list$Neftel_OPC
npc1 <- neftel_gene_list$Neftel_NPC1

cols <- c(
    "Neftel - NPC1" = "skyblue",
    "Neftel - OPC" = "#afd0a6",
    "Invasive signature - Up" = "#2B7095"
)

pdf(
    file.path(params$plot_dir, "FigS4c_inv_up_vs_neftel_venn.pdf"),
    width = 10,
    height = 10
)
venn <- venn.diagram(
    x = list(
        "Invasive signature - Up" = up_sig$gene,
        "Neftel - OPC" = opc,
        "Neftel - NPC1" = npc1
    ),
    category.names = c(
        "Invasive signature - Up",
        "Neftel - OPC",
        "Neftel - NPC1"
    ),
    fill = cols,
    output = TRUE,
    filename = NULL
)
grid.draw(venn)
dev.off()

# ---------------------------------------------------------------------------- #
#                                  Figure S4f                                  #
# ---------------------------------------------------------------------------- #
degs <- readxl::read_excel(
    params$degs_table_path,
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

degs_list <- list(degs_up$gene, degs_dn$gene)

gbmap <- readRDS(params$gbmap_core_obj_path) # path to gbmap obj
gbmap_genes <- rownames(gbmap@assays$RNA@counts)

# Convert gene symbols
genes <- grch38[match(gbmap_genes, grch38$ensgene), "symbol"]
df <- data.frame(gbmap_gene = gbmap_genes, symbol = genes)
df <- df %>%
    mutate(gene = ifelse(symbol == "", gbmap_gene, symbol))
rownames(gbmap@assays$RNA@counts) <- df$gene
rownames(gbmap@assays$RNA@data) <- df$gene

gbmap <- subset(
    gbmap,
    subset = annotation_level_1 == "Neoplastic" &
        sector %in% c("Core", "Periphery")
)
print(colnames(gbmap[[]]))
print(table(gbmap$author, gbmap$sector))
print(gbmap)

plots <- list()
dfs <- list()
for (dataset in unique(gbmap$author)) {
    print(dataset)
    seurat_obj <- subset(gbmap, subset = author == dataset)
    print(seurat_obj)
    seurat_obj <- AddModuleScore(
        seurat_obj,
        features = list(degs_up$gene, degs_dn$gene),
        name = "DEG"
    )
    seurat_obj$inv_sig <- seurat_obj$DEG1 - seurat_obj$DEG2

    seurat_obj$sector <- factor(
        seurat_obj$sector,
        levels = c("Periphery", "Core")
    )

    dfs[[dataset]] <- seurat_obj[[]]
}
gbmap <- as.data.frame(rbindlist(dfs))

df <- gbmap %>% filter(annotation_level_3 %in% c("NPC-like", "OPC-like"))

df <- df %>%
    group_by_at(c("donor_id", "sector", "author")) %>%
    summarise(mean_inv_sig = mean(inv_sig))

patients <- as.data.frame(table(df$donor_id)) %>% filter(Freq == 2)
patients <- patients$Var1

df <- df %>% filter(donor_id %in% patients)

color_blind_friendly_colors <- c(
    "Wu2020" = "#e69f00",
    "Yu2020" = "#d55e00",
    "Darmanis2017" = "#cc79a7",
    "Core" = '#029e73',
    "Periphery" = "#0173b2"
)

p <- ggplot(df, aes(x = sector, y = mean_inv_sig, fill = sector)) +
    geom_boxplot(aes(fill = sector), alpha = 0.8, width = 0.3, size = 1) +
    geom_line(
        aes(group = donor_id),
        colour = "black",
        alpha = 0.3,
        linewidth = 0.8
    ) +
    geom_point(
        aes(group = donor_id, fill = author, shape = 21),
        stroke = 0.8,
        size = 3,
        colour = "black"
    ) +
    scale_shape_identity() +
    stat_compare_means(
        method = "wilcox.test",
        paired = T,
        comparisons = list(c("Core", "Periphery")),
        position = position_dodge(width = 0.8)
    ) +
    scale_fill_manual(values = color_blind_friendly_colors) +
    scale_colour_manual(
        values = c("Core" = '#029e73', "Periphery" = "#0173b2")
    ) +
    ylab("Invasive Signature Score") +
    xlab("Region") +
    GBM_theme() +
    ggtitle(
        "Mean Invasive Signature in OPC/NPC-like cells in Yu2020, Wu2020, Darmanis2017"
    ) +
    guides(
        fill = guide_legend(
            override.aes = list(
                shape = 21,
                size = 3,
                stroke = 0.8,
                colour = "black"
            )
        )
    )
ggsave(
    p,
    filename = "FigS4f_inv_sig_gbmap_paired.pdf",
    path = params$plot_dir,
    height = 10,
    width = 10
)

# ---------------------------------------------------------------------------- #
#                                  Figure S4h                                  #
# ---------------------------------------------------------------------------- #
# Load object
seurat_obj <- readRDS(params$seurat_obj_path)

df <- seurat_obj[[]]
df$Region <- ifelse(df$Region == "PT", "PT", "Tumour")
df$Region <- factor(df$Region, levels = c("PT", "Tumour"))

df <- df %>%
    filter(is_malignant_confident == TRUE) %>%
    filter(CellClass_L2 %in% c("Malignant_OPC", "Malignant_NPC")) %>%
    select(all_of(c("Sample", "Region", "CellClass_L5_2", "Patient"))) %>%
    group_by_at(c("Sample", "Region", "CellClass_L5_2", "Patient")) %>%
    summarise(count = n(), .groups = "drop") %>%
    tidyr::complete(Sample, CellClass_L5_2, fill = list(count = 0)) %>%
    group_by(Sample) %>%
    tidyr::fill(Region, Patient, .direction = "downup") %>%
    mutate(total = sum(count)) %>%
    mutate(Fraction = count / total) %>%
    select(all_of(c(
        "Sample",
        "Region",
        "CellClass_L5_2",
        "Fraction",
        "Patient"
    ))) %>%
    filter(CellClass_L5_2 == "Invasive-high OPC/NPC1")
df$Region <- factor(df$Region, levels = c("PT", "Tumour"))

df2 <- seurat_obj[[]] %>%
    filter(
        is_malignant_confident == TRUE &
            Patient %in% c("6245", "6237", "6419", "6467", "6425", "6514")
    ) %>%
    filter(CellClass_L2 %in% c("Malignant_OPC", "Malignant_NPC")) %>%
    select(all_of(c("Region", "CellClass_L5_2", "Patient")))
df2$Region <- ifelse(df2$Region == "PT", "PT", "Tumour")
df2$Region <- factor(df2$Region, levels = c("PT", "Tumour"))

df2 <- df2 %>%
    group_by_at(c("Region", "CellClass_L5_2", "Patient")) %>%
    summarise(count = n(), .groups = "drop") %>%
    tidyr::complete(Patient, CellClass_L5_2, fill = list(count = 0)) %>%
    group_by(Patient) %>%
    tidyr::fill(Region, .direction = "downup") %>%
    group_by_at(c("Patient", "Region")) %>%
    mutate(total = sum(count)) %>%
    mutate(Fraction = count / total) %>%
    select(all_of(c("Region", "CellClass_L5_2", "Fraction", "Patient"))) %>%
    filter(CellClass_L5_2 == "Invasive-high OPC/NPC1")

p <- ggplot(df, aes(x = Region, y = Fraction, fill = Region)) +
    geom_boxplot(alpha = 0.9, outlier.shape = NA) +
    geom_point(
        aes(group = Patient),
        size = 2.5,
        position = position_dodge(width = 0.75)
    ) +
    geom_line(
        data = df2,
        aes(x = Region, y = Fraction, group = Patient),
        colour = "black",
        alpha = 0.3,
        linewidth = 0.8
    ) +
    xlab("Region") +
    ylab("Proportion") +
    ggtitle("Proportion of invasive-high OPC/NPC1-like cells across regions") +
    scale_fill_manual(values = c(region_cols)) +
    GBM_theme() +
    theme(
        legend.text = element_text(size = 18),
        legend.key.size = unit(1, "cm")
    ) +
    geom_signif_lmm(
        data_df = df,
        response = "Fraction",
        condition = "Region",
        latent_vars = c("Patient"),
        comparisons = list(c("PT", "Tumour")),
    )
ggsave(
    filename = "FigS4h_invasive_high_props.pdf",
    path = params$plot_dir,
    width = 5,
    height = 7
)
