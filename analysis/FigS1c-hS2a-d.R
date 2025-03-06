# ---- Code to reproduce Figure S1c-h S2a-d ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

pacman::p_load(
    data.table,
    ggplot2,
    Seurat,
    tidyr,
    dplyr,
    patchwork,
    scales,
    here,
    GBMutils,
    GaitiLabUtils
)

# Required inputs
params <- list.files(
    seurat_obj_path = "",
    plot_dir = "output/submission/figures",
    scvi_latent_rep_path = "snrnaseq_umap_coords.csv"
)


# ------------------ Creating directories needed for outputs ----------------- #
GaitiLabUtils::create_dir(params$plot_dir)


# Loading Seurat object
seurat_obj <- readRDS(params$path_to_seurat_object)

print(seurat_obj)

DefaultAssay(seurat_obj) <- "RNA"

# load color palette
region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

color_palette <- load_color_palette(name = "CellClass_L2")
to_remove <- c(
    "Malignant_MES",
    "Malignant_NPC",
    "Malignant_OPC",
    "Malignant_AC"
)
nm_cols <- c(
    color_palette[setdiff(names(color_palette), malignant_l2)],
    "Myeloid" = "#787335"
)
malignant_cols <- c(
    "Malignant_AC" = "#e0b73d",
    "Malignant_MES_INT" = "firebrick",
    "Malignant_MES_HYP" = "tomato",
    "Malignant_MES_AST" = "salmon",
    "Malignant_NPC1" = "skyblue",
    "Malignant_NPC2" = "lightblue",
    "Malignant_OPC" = "#afd0a6"
)
prog_like_cols <- c(
    "Malignant_NPC1" = "skyblue",
    "Malignant_NPC2" = "lightblue",
    "Malignant_OPC" = "#afd0a6",
    "Invasive-high OPC/NPC1" = "#2B7095"
)
color_palette <- c(nm_cols, malignant_cols)

# ---------------------------------------------------------------------------- #
#                                  Figure S1c                                  #
# ---------------------------------------------------------------------------- #
# Import scVI embeddings
dat <- read.csv(params$scvi_latent_rep_path, row.names = 1) # path to scvi latent rep
mat <- as.matrix(dat)
seurat_obj[["scVIumap"]] <- CreateDimReducObject(
    embeddings = mat,
    key = "XscVIumap_",
    assay = "RNA"
)

# malignant PT highlighted umap
so_umap$umap <- ifelse(
    so_umap$is_malignant_confident == TRUE & so_umap$Region == "PT",
    "PT-Malignant",
    "Other"
)
p <- DimPlot(
    so_umap,
    reduction = "scVIumap",
    group.by = "umap",
    cols = c("PT-Malignant" = "black", "Other" = "grey")
)
ggsave(
    p,
    filename = "S1c_scvi_pt_malignant_umap.pdf",
    path = plot_dir,
    height = 7,
    width = 9
)

# EGFR expression of malignant PT vs non-malignant
metadata <- seurat_obj@meta.data %>%
    mutate(
        vio = case_when(
            seurat_obj$is_malignant_confident == TRUE &
                seurat_obj$Region == "PT" ~
                "PT-Malignant",
            seurat_obj$is_malignant_confident == TRUE &
                seurat_obj$Region == "TE" ~
                "TE-Malignant",
            seurat_obj$is_malignant_confident == TRUE &
                seurat_obj$Region == "TC" ~
                "TC-Malignant",
            seurat_obj$CellClass_L1 != "Malignant" ~ "Non-malignant"
        )
    )
seurat_obj@meta.data <- metadata
seurat_obj$vio <- factor(
    seurat_obj$vio,
    levels = c("Non-malignant", "PT-Malignant", "TE-Malignant", "TC-Malignant")
)
p <- VlnPlot(
    seurat_obj,
    features = c("EGFR", "PTPRZ1"),
    pt.size = 0,
    group.by = "vio",
    cols = c(
        "PT-Malignant" = "#7BB7AD",
        "TE-Malignant" = "#7BB7AD",
        "TC-Malignant" = "#7BB7AD",
        "Non-malignant" = "grey"
    )
)
ggsave(
    plot = p,
    filename = "S1c_PT_EGFR_PTPRZ1_expression.pdf",
    path = params$plot_dir,
    height = 5,
    width = 7
)


# ---------------------------------------------------------------------------- #
#                                  Figure S1d                                  #
# ---------------------------------------------------------------------------- #
df <- seurat_obj[[]]
df <- df %>%
    filter(CellClass_L1 != "Malignant") %>%
    select(all_of(c("Region", "Sample", "CellClass_L1"))) %>%
    group_by_at(c("Region", "Sample", "CellClass_L1")) %>%
    summarise(n = n())
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))

nm_prop <- seurat_obj[[]] %>%
    filter(CellClass_L1 != "Malignant") %>%
    select(all_of(c("Region", "Sample", "CellClass_L1"))) %>%
    group_by_at(c("Region", "Sample", "CellClass_L1")) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by_at(c("Region", "Sample")) %>%
    mutate(total = sum(n), prop = n / total) %>%
    filter(CellClass_L1 == "Oligodendrocyte") %>%
    arrange(prop)
df$Sample <- factor(df$Sample, levels = nm_prop$Sample)

p <- ggplot(df, aes(x = Sample, y = n, fill = factor(CellClass_L1))) +
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    scale_fill_manual(values = nm_cols) +
    labs(y = "Percent") +
    scale_y_continuous(labels = scales::percent) +
    facet_grid(. ~ Region, scales = "free_x", space = "free_x") +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 50),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 50),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 50),
        strip.background = element_blank(),
        panel.spacing = unit(1, "cm"),
        legend.text = element_text(size = 35),
        legend.key.size = unit(3, "cm"),
        legend.position = "bottom",
        legend.title = element_blank()
    )
ggsave(
    plot = p,
    filename = "S1d_non_malig_cell_type_prop_by_sample.pdf",
    path = params$plot_dir,
    height = 25,
    width = 25
)

# ---------------------------------------------------------------------------- #
#                                  Figure S1e                                  #
# ---------------------------------------------------------------------------- #
Idents(seurat_obj) <- "CellClass_L1"
Idents(seurat_obj) <- factor(
    Idents(seurat_obj),
    levels = c(
        "Astrocyte",
        "Endothelial",
        "Malignant",
        "Myeloid",
        "Neuron",
        "Oligodendrocyte",
        "OPC",
        "Pericyte",
        "T_cell"
    )
)

markers <- c(
    "SLC1A2",
    "ADGRV1",
    "FLT1",
    "ABCB1",
    "EGFR",
    "PTPRZ1",
    "DOCK8",
    "APBB1IP",
    "CNTNAP2",
    "SYT1",
    "PLP1",
    "MBP",
    "PCDH15",
    "CA10",
    "IL7R",
    "SKAP1"
)

p <- DotPlot(
    seurat_obj,
    features = markers,
    dot.scale = 8,
    assay = "RNA",
    scale = TRUE
) +
    RotatedAxis() +
    coord_flip() +
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal"
    )
ggsave(
    plot = p,
    filename = "S1e_markers_dotplot.pdf",
    path = params$plot_dir,
    height = 7,
    width = 5
)

# ---------------------------------------------------------------------------- #
#                                  Figure S1f                                  #
# ---------------------------------------------------------------------------- #

df <- seurat_obj[[]] %>%
    filter(CellClass_L1 == "Neuron") %>%
    filter(CellClass_L2 != "Neuron") %>%
    group_by_at("CellClass_L2") %>%
    summarise(n = n())

p <- ggplot(df, aes(x = "", y = n, fill = factor(CellClass_L2))) +
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    scale_fill_manual(values = c("#D9C28F", "#B8A56C")) +
    labs(y = "Proportion") +
    scale_y_continuous(labels = scales::percent) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        panel.spacing = unit(1, "cm"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom",
        legend.title = element_blank()
    )
ggsave(
    plot = p,
    filename = "S1f_neuron_subtype_distribution.pdf",
    path = plot_dir,
    width = 3
)

# Neuron subtype markers
neurons <- subset(seurat_obj, subset = CellClass_L1 == "Neuron")
neurons <- subset(neurons, subset = CellClass_L2 != "Neuron")

Idents(neurons) <- "CellClass_L2"
Idents(neurons) <- factor(
    Idents(neurons),
    levels = c("Glutamatergic", "GABAergic")
)

markers <- c("SLC17A7", "SLC17A6", "GRIN1", "GAD2", "GAD1", "SLC32A1")

p <- DotPlot(
    neurons,
    features = markers,
    dot.scale = 10,
    assay = "RNA",
    scale = TRUE
) +
    coord_flip() +
    RotatedAxis()
ggsave(
    plot = p,
    filename = "S1f_neuron_markers.pdf",
    path = params$plot_dir,
    height = 6,
    width = 4
)

rm(neurons)

# ---------------------------------------------------------------------------- #
#                                  Figure S1g                                  #
# ---------------------------------------------------------------------------- #

df <- seurat_obj[[]] %>%
    filter(CellClass_L1 == "Myeloid") %>%
    select(all_of(c("Region", "CCI_CellClass_L2"))) %>%
    mutate(Region = ifelse(Region == "PT", "PT", "Tumor")) %>%
    group_by_at(c("Region", "CCI_CellClass_L2")) %>%
    summarise(n = n())
df$Region <- factor(df$Region, levels = c("PT", "Tumor"))

p <- ggplot(df, aes(x = Region, y = n, fill = factor(CCI_CellClass_L2))) +
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    scale_fill_manual(values = c("#000000", "#D3D3D3")) +
    labs(y = "Proportion") +
    scale_y_continuous(labels = scales::percent) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        panel.spacing = unit(1, "cm"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom",
        legend.title = element_blank()
    )
ggsave(
    plot = p,
    filename = "S1g_myeloid_subtype_distribution.pdf",
    path = params$plot_dir,
    width = 5
)

# Myeloid subtype markers
myeloid <- subset(seurat_obj, subset = CellClass_L1 == "Myeloid")

Idents(myeloid) <- "CCI_CellClass_L2"
Idents(myeloid) <- factor(
    Idents(myeloid),
    levels = c("Myeloid_Immunosuppressive", "Myeloid_Inflammatory")
)

markers <- c("C1QA", "TGFBI", "CD163", "CD83", "IL1B", "CCL3")

p <- DotPlot(
    myeloid,
    features = markers,
    dot.scale = 10,
    assay = "RNA",
    scale = TRUE
) +
    coord_flip() +
    RotatedAxis()
ggsave(
    plot = p,
    filename = "S1g_myeloid_markers.pdf",
    path = params$plot_dir,
    height = 6,
    width = 4
)

rm(myeloid)
# ---------------------------------------------------------------------------- #
#                                  Figure S1h                                  #
# ---------------------------------------------------------------------------- #

df <- seurat_obj[[]] %>%
    filter(is_malignant_confident == TRUE | CellClass_L1 != "Malignant") %>%
    select(all_of(c("Sample", "Region", "CellClass_L1", "Patient")))
df$CellClass_L1 <- ifelse(
    df$CellClass_L1 == "Malignant",
    "Malignant",
    "Non-Malignant"
)
df <- df %>%
    group_by_at(c("Sample", "Region", "CellClass_L1", "Patient")) %>%
    summarise(count = n()) %>%
    group_by(Sample) %>%
    mutate(total = sum(count)) %>%
    mutate(Fraction = count / total) %>%
    filter(CellClass_L1 == "Malignant")
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))

p <- ggplot(df, aes(x = Region, y = Fraction, fill = Region)) +
    geom_boxplot(alpha = 0.9) +
    geom_point(size = 2.5, position = position_dodge(width = 0.75)) +
    xlab("Region") +
    ylab("Relative Abundance of Malignant Cells") +
    ggtitle("Tumour Cell Purity") +
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
        comparisons = comps,
        step_increase = c(0, 0.1, 0.1)
    )
ggsave(
    plot = p,
    filename = "S1h_tumour_cell_purity.pdf",
    path = params$plot_dir,
    width = 10,
    height = 10
)

# ---------------------------------------------------------------------------- #
#                                  Figure S2a                                  #
# ---------------------------------------------------------------------------- #
conf_cells <- seurat_obj[[]] %>% filter(is_malignant_confident == TRUE)
conf_cell_id <- rownames(conf_cells)
df <- read.csv("neftel_metamodule_scores.csv", row.names = 1)
df <- df %>% filter(Cell_id %in% conf_cell_id)

df$Type <- NA
df$Type[grep("NPC1", df$Neftel.modules)] <- "NPC"
df$Type[grep("NPC2", df$Neftel.modules)] <- "NPC"
df$Type[grep("MES1", df$Neftel.modules)] <- "MES"
df$Type[grep("MES2", df$Neftel.modules)] <- "MES"
df$Type[grep("OPC", df$Neftel.modules)] <- "OPC"
df$Type[grep("AC", df$Neftel.modules)] <- "AC"

NPC_num <- paste0("NPC-like (n=", length(df$Type[df$Type == "NPC"]), ")")
MES_num <- paste0("MES-like (n=", length(df$Type[df$Type == "MES"]), ")")
OPC_num <- paste0("OPC-like (n=", length(df$Type[df$Type == "OPC"]), ")")
AC_num <- paste0("AC-like (n=", length(df$Type[df$Type == "AC"]), ")")

annotations <- data.frame(
    xpos = c(-Inf, -Inf, Inf, Inf),
    ypos = c(-Inf, Inf, -Inf, Inf),
    annotateText = c(AC_num, OPC_num, MES_num, NPC_num),
    hjustvar = c(-0.3, -0.3, 1.2, 1.2),
    vjustvar = c(-1, 2, -1, 2)
)

a <- ggplot(df, aes(x = logneg, y = D)) +
    scale_x_continuous(limits = c(-3, 3)) +
    scale_y_continuous(limits = c(-3, 3)) +
    geom_point(
        aes(col = NPC),
        data = subset(df, Type == "NPC"),
        size = 0.5,
        alpha = 1.0,
        stroke = NA
    ) +
    scale_color_gradient2(
        "NPC",
        low = "cadetblue1",
        mid = "cadetblue3",
        high = "navy"
    ) +
    new_scale("color") +
    geom_point(
        aes(col = OPC),
        data = subset(df, Type == "OPC"),
        size = 0.5,
        alpha = 1.0,
        stroke = NA
    ) +
    scale_color_gradient2(
        "OPC",
        low = "darkseagreen1",
        mid = "darkseagreen3",
        high = "darkolivegreen4"
    ) +
    new_scale("color") +
    geom_point(
        aes(col = AC),
        data = subset(df, Type == "AC"),
        size = 0.5,
        alpha = 1.0,
        stroke = NA
    ) +
    scale_color_gradient2(
        "AC",
        low = "darkgoldenrod1",
        mid = "darkgoldenrod2",
        high = "darkgoldenrod4"
    ) +
    new_scale("color") +
    geom_point(
        aes(col = MES),
        data = subset(df, Type == "MES"),
        size = 0.5,
        alpha = 1.0,
        stroke = NA
    ) +
    scale_color_gradient2(
        "MES",
        low = "lavenderblush1",
        mid = "coral2",
        high = "darkred"
    ) +
    new_scale("color") +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    coord_cartesian(xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0)) +
    ylab("Relative meta-module score\n[log2(|(NPC+OPC) - (AC+MES)|+1)]") +
    xlab("Relative meta-module score\n[log2(|(NPC+MES) - (OPC+AC)|+1)]") +
    theme(
        plot.title = element_text(
            hjust = 0.5,
            color = "black",
            family = "Helvetica"
        ),
        axis.title.y = element_text(
            color = "black",
            size = 18,
            face = "bold",
            margin = margin(t = 0, r = 10, b = 0, l = 0),
            family = "Helvetica"
        ),
        axis.title.x = element_text(
            color = "black",
            size = 18,
            face = "bold",
            margin = margin(t = 0, r = 10, b = 0, l = 0),
            family = "Helvetica"
        ),
        axis.text.x = element_text(size = 13, family = "Helvetica"),
        axis.text.y = element_text(size = 13, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(5, 5, 5, 5),
        legend.text = element_text(face = "bold", family = "Helvetica"),
        legend.spacing.x = unit(0.2, "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    geom_text(
        data = annotations,
        aes(
            x = xpos,
            y = ypos,
            hjust = hjustvar,
            vjust = vjustvar,
            label = annotateText
        ),
        size = 5,
        family = "Helvetica"
    )

ggsave(
    plot = a,
    filename = "S2a_Neftel_quadrant_plot.pdf",
    path = params$plot_dir,
    width = 15,
    height = 15
)
# ---------------------------------------------------------------------------- #
#                                  Figure S2b                                  #
# ---------------------------------------------------------------------------- #

df <- seurat_obj[[]] %>%
    filter(is_malignant_confident == TRUE) %>%
    select(all_of(c("Region", "Sample", "CellClass_L3", "Patient"))) %>%
    group_by_at(c("Region", "Sample", "CellClass_L3", "Patient")) %>%
    summarise(n = n())
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))
write.csv(df, file.path(plot_dir, "malignant_cell_num_by_sample.csv"))

m_prop <- seurat_obj[[]] %>%
    filter(is_malignant_confident == TRUE) %>%
    select(all_of(c("Region", "Sample", "CellClass_L3"))) %>%
    group_by_at(c("Region", "Sample", "CellClass_L3")) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by_at(c("Region", "Sample")) %>%
    mutate(total = sum(n), prop = n / total) %>%
    filter(CellClass_L3 == "Malignant_OPC") %>%
    arrange(prop)
df$Sample <- factor(df$Sample, levels = m_prop$Sample)

p <- ggplot(df, aes(x = Sample, y = n, fill = factor(CellClass_L3))) +
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    scale_fill_manual(values = malignant_cols) +
    labs(y = "Percent") +
    scale_y_continuous(labels = scales::percent) +
    facet_grid(. ~ Region, scales = "free_x", space = "free_x") +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 50),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 50),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 50),
        strip.background = element_blank(),
        panel.spacing = unit(1, "cm"),
        legend.text = element_text(size = 35),
        legend.key.size = unit(3, "cm"),
        legend.position = "bottom",
        legend.title = element_blank()
    )
ggsave(
    plot = p,
    filename = "S2b_cell_state_prop_by_sample.pdf",
    path = params$plot_dir,
    height = 25,
    width = 25
)

# ---------------------------------------------------------------------------- #
#                                  Figure S2c                                  #
# ---------------------------------------------------------------------------- #
# Neftel malignant state markers

malig_dotplot <- subset(seurat_obj, subset = is_malignant_confident == TRUE)

Idents(malig_dotplot) <- "CellClass_L3"
Idents(malig_dotplot) <- factor(
    Idents(malig_dotplot),
    levels = c(
        "Malignant_OPC",
        "Malignant_NPC1",
        "Malignant_NPC2",
        "Malignant_MES_INT",
        "Malignant_MES_AST",
        "Malignant_AC",
        "Malignant_MES_HYP"
    )
)

markers <- c(
    "MKI67",
    "TNR",
    "BCAN",
    "OLIG2",
    "OLIG1",
    "DLL1",
    "DLL3",
    "SOX4",
    "SOX11",
    "TUBB3",
    "CD24",
    "VIM",
    "CD44",
    "CLU",
    "APOE",
    "AQP4",
    "GFAP",
    "NDRG1",
    "PKM",
    "ADM",
    "ATF3"
)

p <- DotPlot(
    malig_dotplot,
    features = markers,
    dot.scale = 8,
    assay = "RNA",
    scale = TRUE
) +
    RotatedAxis()
ggsave(
    plot = p,
    filename = "S2c_malignant_markers_dotplot.pdf",
    path = plot_dir,
    height = 5,
    width = 8
)

rm(malig_dotplot)

# ---------------------------------------------------------------------------- #
#                                  Figure S2d                                  #
# ---------------------------------------------------------------------------- #
df <- seurat_obj[[]] %>%
    filter(is_malignant_confident == TRUE) %>%
    group_by_at(c("Sample", "Region", "CellClass_L3", "Patient")) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(total = sum(count)) %>%
    mutate(Fraction = count / total) %>%
    select(all_of(c("Sample", "Region", "CellClass_L3", "Fraction", "Patient")))

df$Region <- ifelse(df$Region == "PT", "PT", "Tumour")
df$Region <- factor(df$Region, levels = c("PT", "Tumour"))

plots <- list()
for (subtype in unique(df$CellClass_L3)) {
    curr_data <- df %>% filter(CellClass_L3 == subtype)

    plots[[subtype]] <- ggplot(
        curr_data,
        aes(x = Region, y = Fraction, fill = Region)
    ) +
        geom_boxplot(alpha = 0.9, outlier.shape = NA) +
        geom_point(
            aes(group = Patient),
            size = 2.5,
            position = position_dodge(width = 0.75)
        ) +
        xlab("Region") +
        ylab("Proportion of Malignant Cells") +
        scale_fill_manual(values = c(region_cols[1], region_cols[3])) +
        GBM_theme() +
        theme(
            legend.text = element_text(size = 18),
            legend.key.size = unit(1, "cm")
        ) +
        ggtitle(subtype) +
        geom_signif_lmm(
            data_df = curr_data,
            response = "Fraction",
            condition = "Region",
            latent_vars = c("Patient"),
            comparisons = list(c("PT", "Tumour"))
        )
}
plots <- wrap_plots(plots, ncol = 3, nrow = 3, guides = "collect")
ggsave(
    plots,
    filename = "S2d_malignant_sample_cell_type_props.pdf",
    path = params$plot_dir,
    width = 10,
    height = 14
)

# ---------------------------------------------------------------------------- #
#                                  Figure S2e                                  #
# ---------------------------------------------------------------------------- #
gene_lists <- readxl::read_excel(params$genelists_path, skip = 1)

# Developmental, injury response (Richards 2021)
richards_dev <- gene_lists %>%
    filter(!is.na(Richards_Developmental)) %>%
    pull(Richards_Developmental)

richards_injury <- gene_lists %>%
    filter(!is.na(Richards_Injury_Response)) %>%
    pull(Richards_Injury_Response)

richards <- list(richards_dev, richards_injury)
names(richards) <- c("Richards_Developmental", "Richards_Injury_Response")

# Hallmark hypoxia
hypoxia_markers <- gene_lists %>%
    filter(!is.na(Hallmark_Hypoxia)) %>%
    pull(Hallmark_Hypoxia)

hypoxia_markers <- list(hypoxia_markers)
names(hypoxia_markers) <- "Hypoxia"

# Hallmark TNF alpha
tnf_markers <- gene_lists %>%
    filter(!is.na(Hallmark_TNF_alpha)) %>%
    pull(Hallmark_TNF_alpha)

tnf_markers <- list(tnf_markers)
names(tnf_markers) <- "TNF_alpha"

# Combine gene signatures

genesets <- c(
    richards_gene_list,
    hypoxia_markers,
    tnf_markers,
)

seurat_obj <- AddModuleScore(seurat_obj, features = genesets, name = "geneset")

score_df <- seurat_obj[[]][, grepl("geneset", colnames(seurat_obj[[]]))]
colnames(score_df) <- auc_colnames
score_df$Richards <- score_df$Richards_Developmental -
    score_df$Richards_Injury_Response
score_df <- score_df[, -which(colnames(score_df) == "Richards_Developmental")]
score_df <- score_df[, -which(colnames(score_df) == "Richards_Injury_Response")]

df <- as.data.table(seurat_obj@meta.data)[, .(
    Region,
    Sample,
    Patient,
    CellClass_L3,
    CellClass_L2
)]
df <- cbind(df, score_df)

boxplots <- list()
dotplot_sample_level <- list()
progenitor_subset <- list()
diff_subset <- list()
for (gene_sig in colnames(score_df)) {
    # Sample level
    curr_data <- df[, .(Region, Sample, Patient, Program = get(gene_sig))]
    curr_data$Region <- factor(curr_data$Region, levels = c("PT", "TE", "TC"))
    df_list <- lapply(unique(curr_data$Patient), function(x) {
        df <- curr_data %>%
            filter(Patient == x)
        df$Program <- scale(df$Program)
        return(df)
    })
    curr_data <- rbindlist(df_list)
    curr_data <- curr_data %>%
        group_by_at(c("Sample", "Region")) %>%
        summarise(Program = mean(Program))
    comps <- list(c("PT", "TE"), c("TE", "TC"), c("PT", "TC"))

    dotplot_sample_level[[gene_sig]] <- ggplot(
        curr_data,
        aes(x = Region, y = Program)
    ) +
        geom_boxplot(aes(fill = Region), outlier.shape = NA) +
        scale_fill_manual(values = c(region_cols)) +
        geom_point(
            aes(group = Sample),
            size = 4,
            position = position_dodge(width = 0.75)
        ) +
        GBM_theme() +
        xlab("") +
        ylab("Mean Normalized Module Score") +
        ggtitle(title) +
        stat_compare_means(
            label = "p.format",
            method = "wilcox.test",
            comparisons = comps
        )
}

plots <- wrap_plots(
    dotplot_sample_level,
    ncol = 4,
    nrow = 5,
    guides = "collect"
) &
    theme(legend.position = "bottom")
ggsave(
    plots,
    filename = "S2e_mod_score_gene_sig_boxplots_sample_level.pdf",
    path = params$plot_dir,
    height = 25,
    width = 18
)

# ---------------------------------------------------------------------------- #
#                                  Figure S2f                                  #
# ---------------------------------------------------------------------------- #

df <- seurat_obj[[]] %>%
    filter(CellClass_L1 != "Malignant") %>%
    select(all_of(c("Sample", "Region", "CellClass_L1", "Patient"))) %>%
    filter(
        CellClass_L1 %in% c("Neuron", "Oligodendrocyte", "Myeloid", "T_cell")
    ) %>%
    group_by_at(c("Sample", "Region", "CellClass_L1", "Patient")) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(total = sum(count)) %>%
    mutate(Fraction = count / total) %>%
    select(all_of(c("Sample", "Region", "CellClass_L1", "Fraction", "Patient")))
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))

plots <- list()
for (celltype in unique(df$CellClass_L1)) {
    curr_data <- df %>% filter(CellClass_L1 == celltype)
    plots[[celltype]] <- ggplot(
        curr_data,
        aes(x = Region, y = Fraction, fill = Region)
    ) +
        geom_boxplot(alpha = 0.9, outlier.shape = NA) +
        geom_point(
            aes(group = Sample),
            size = 2.5,
            position = position_dodge(width = 0.75)
        ) +
        ggtitle(celltype) +
        xlab("Region") +
        ylab("Proportion of Non-Malignant Cells") +
        scale_fill_manual(values = c(region_cols)) +
        GBM_theme() +
        theme(
            legend.text = element_text(size = 18),
            legend.key.size = unit(1, "cm")
        ) +
        geom_signif_lmm(
            data_df = curr_data,
            response = "Fraction",
            condition = "Region",
            latent_vars = c("Patient"),
            comparisons = comps,
            step_increase = c(0, 0.1, 0.1)
        )
}
plots <- wrap_plots(plots, ncol = 2, nrow = 2, guides = "collect")
ggsave(
    plots,
    filename = "S2f_nonmalignant_sample_cell_type_props.pdf",
    path = params$plot_dir,
    height = 9,
    width = 9
)
