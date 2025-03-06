# ---- Code to reproduce Figure 4a,b and 6a,b,c,g ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #
pacman::p_load(
    data.table,
    patchwork,
    ggplot2,
    dplyr,
    tidyr,
    ggpubr,
    Seurat,
    harmony,
    destiny,
    scales,
    sceasy
)

region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

# Required inputs
params <- list(
    # Path to Seurat object generated using this manuscript's data, raw data and final metadata can be downloaded online, see publication
    seurat_obj_path = "",
    plot_dir = "output/figures",
    # Download data from https://github.com/linnarsson-lab/developing-human-brain
    braun_h5ad_path = "developing_opc.h5ad",
    braun_subsampled_h5ad_path = "developing_brain_subsampled.h5ad",
    genelists_path = "misc/gene_signatures.xlsx",
    degs_table_path = "misc/Table S2.xlsx" # can be downloaded online, see publication
)

# Creating directories needed for outputs
GaitiLabUtils::create_dir(params$output_dir)
GaitiLabUtils::create_dir(params$plot_dir)

# ---------------------------------------------------------------------------- #
#                                 Convert h5ad                                 #
# ---------------------------------------------------------------------------- #
# Extract Oligo and Pre-OPC from Braun et al. dataset
seurat_obj <- sceasy::convertFormat(
    params$braun_h5ad_path,
    from = "anndata",
    to = "seurat",
    outFile = "filename.rds"
)
print(seurat_obj)
print(colnames(seurat_obj[[]]))

seurat_obj$CellClass_L3 <- case_when(
    seurat_obj$Clusters %in% c("609", "610") ~ "COP",
    seurat_obj$Clusters %in% paste0(seq(611, 616)) ~ "OPC",
    seurat_obj$Clusters == "19" ~ "Pre-OPC"
)
seurat_obj$CellClass_L3 <- factor(
    seurat_obj$CellClass_L3,
    levels = c("Pre-OPC", "OPC", "COP")
)

# Randomly subsampled 5k cells for each cellclass in Braun et al. dataset
dev_brain <- sceasy::convertFormat(
    params$braun_subsampled_h5ad_path,
    from = "anndata",
    to = "seurat",
    outFile = "filename.rds"
)

# ---------------------------------------------------------------------------- #
#                                   Figure 4a                                  #
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

dev_brain <- AddModuleScore(
    dev_brain,
    features = list(degs_up$gene, degs_dn$gene),
    name = "INV"
)
dev_brain$invasive_signature <- dev_brain$INV1 - dev_brain$INV2

df <- dev_brain[[]] %>%
    group_by(CellClass) %>%
    summarise(median_inv_sig = median(invasive_signature)) %>%
    arrange(-median_inv_sig)

plot_df <- dev_brain[[]]

# Glial lineage
plot_df_glial <- plot_df %>%
    filter(CellClass %in% c("Radial glia", "Glioblast", "Oligo"))
plot_df_glial$CellClass <- factor(
    plot_df_glial$CellClass,
    levels = c("Radial glia", "Glioblast", "Oligo")
)
plot_df_glial <- plot_df_glial %>%
    group_by(CellClass) %>%
    summarise(
        mean_score = mean(invasive_signature),
        stdev_score = sd(invasive_signature)
    ) %>%
    mutate(ymin = mean_score - stdev_score, ymax = mean_score + stdev_score)
plot_df_glial$CellClass <- factor(
    plot_df_glial$CellClass,
    levels = c("Radial glia", "Glioblast", "Oligo")
)
p <- ggplot(plot_df_glial, aes(x = CellClass, y = mean_score, group = 1)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "black") +
    geom_point(color = "black") +
    geom_line(color = "black") +
    GBM_theme() +
    ylim(c(-0.2, 0.35))
ggsave(
    filename = "Fig4a_inv_sig_across_glial_lineage_cont.pdf",
    path = params$plot_dir,
    width = 7,
    height = 7
)

# Neuronal lineage
plot_df_neuronal <- plot_df %>%
    filter(
        CellClass %in% c("Radial glia", "Neuronal IPC", "Neuroblast", "Neuron")
    )
plot_df_neuronal$CellClass <- factor(
    plot_df_neuronal$CellClass,
    levels = c("Radial glia", "Neuronal IPC", "Neuroblast", "Neuron")
)
plot_df_neuronal <- plot_df_neuronal %>%
    group_by(CellClass) %>%
    summarise(
        mean_score = mean(invasive_signature),
        stdev_score = sd(invasive_signature)
    ) %>%
    mutate(ymin = mean_score - stdev_score, ymax = mean_score + stdev_score)
plot_df_neuronal$CellClass <- factor(
    plot_df_neuronal$CellClass,
    levels = c("Radial glia", "Neuronal IPC", "Neuroblast", "Neuron")
)
p <- ggplot(plot_df_neuronal, aes(x = CellClass, y = mean_score, group = 1)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "black") +
    geom_point(color = "black") +
    geom_line(color = "black") +
    GBM_theme() +
    ylim(c(-0.2, 0.35))
ggsave(
    plot = p,
    filename = "Fig4a_inv_sig_across_neuronal_lineage_cont.pdf",
    path = params$plot_dir,
    width = 7,
    height = 7
)

# ---------------------------------------------------------------------------- #
#                            Integrate with harmony                            #
# ---------------------------------------------------------------------------- #
print(head(seurat_obj[[]]))

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("CellCycle"))
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

seurat_obj <- RunHarmony(
    seurat_obj,
    group.by.vars = c("Platform", "Sample"),
    reduction.use = "pca"
)

harmony_embeddings <- Embeddings(seurat_obj, reduction = "harmony")

# ---------------------------------------------------------------------------- #
#                               DPT with harmony                               #
# ---------------------------------------------------------------------------- #

dm <- DiffusionMap(data = harmony_embeddings)
dm_coords <- eigenvectors(dm)
dm_coords <- dm_coords[, c(1, 2)]
colnames(dm_coords) <- c("DC1", "DC2")
df <- cbind(dm_coords, seurat_obj[[]] %>% dplyr::select(CellClass_L3))

df$rank <- rank(-df$DC1)
index <- 1:length(df$DC1)
dpt <- DPT(dm, tips = index[df$rank == 1])
df$dpt <- dpt$dpt

write.csv(df, file.path(params$output_dir, "diffusion_comps.csv"))

# ---------------------------------------------------------------------------- #
#                                 Figure S6e,f                                 #
# ---------------------------------------------------------------------------- #

# Load gene sets
gene_lists <- readxl::read_excel(params$genelists_path, skip = 1)

neftel_gene_list <- lapply(
    gene_lists %>% dplyr::select(starts_with("Neftel_")) %>% as.list(),
    function(gene_list) {
        gene_list[!is.na(gene_list)]
    }
)

opc <- neftel_gene_list$Neftel_OPC

synapse <- gene_lists %>%
    filter(!is.na(GOCC_SYNAPTIC_MEMBRANE)) %>%
    pull(GOCC_SYNAPTIC_MEMBRANE)

synaptic_signaling <- gene_lists %>%
    filter(!is.na(GOBP_SYNAPTIC_SIGNALING)) %>%
    pull(GOBP_SYNAPTIC_SIGNALING)

gene_list <- list(degs_up$gene, degs_dn$gene, opc, synapse, synaptic_signaling)

# Score cells in normal opc lineage
seurat_obj <- AddModuleScore(seurat_obj, features = gene_list, name = "DEG")
seurat_obj$invasive_signature <- seurat_obj$DEG1 - seurat_obj$DEG2

cols <- c("#09283CFF", "#F18B00FF", "#F2EBBBFF")

p <- ggplot(seurat_obj[[]], aes(x = CellClass_L3, y = invasive_signature)) +
    geom_boxplot(aes(fill = CellClass_L3)) +
    GBM_theme() +
    stat_compare_means(
        comparisons = list(
            c("Pre-OPC", "OPC"),
            c("OPC", "COP"),
            c("Pre-OPC", "COP")
        )
    ) +
    scale_fill_manual(values = cols)
ggsave(
    plot = p,
    filename = "FigS6e_inv_sig_in_opc_lineage.pdf",
    path = params$plot_dir,
    width = 4
)

p <- ggplot(seurat_obj[[]], aes(x = CellClass_L3, y = DEG3)) +
    geom_boxplot(aes(fill = CellClass_L3)) +
    GBM_theme() +
    stat_compare_means(
        comparisons = list(
            c("Pre-OPC", "OPC"),
            c("OPC", "COP"),
            c("Pre-OPC", "COP")
        )
    ) +
    scale_fill_manual(values = cols)
ggsave(
    plot = p,
    filename = "FigS6e_neftel_opc_in_opc_lineage.pdf",
    path = params$plot_dir,
    width = 4
)

p <- ggplot(seurat_obj[[]], aes(x = CellClass_L3, y = DEG4)) +
    geom_boxplot(aes(fill = CellClass_L3)) +
    GBM_theme() +
    stat_compare_means(
        comparisons = list(
            c("Pre-OPC", "OPC"),
            c("OPC", "COP"),
            c("Pre-OPC", "COP")
        )
    ) +
    scale_fill_manual(values = cols)
ggsave(
    plot = p,
    filename = "FigS6f_synapse_in_opc_lineage.pdf",
    path = params$plot_dir,
    width = 4
)

p <- ggplot(seurat_obj[[]], aes(x = CellClass_L3, y = DEG5)) +
    geom_boxplot(aes(fill = CellClass_L3)) +
    GBM_theme() +
    stat_compare_means(
        comparisons = list(
            c("Pre-OPC", "OPC"),
            c("OPC", "COP"),
            c("Pre-OPC", "COP")
        )
    ) +
    scale_fill_manual(values = cols)
ggsave(
    plot = p,
    filename = "FigS6f_synaptic_signaling_in_opc_lineage.pdf",
    path = params$plot_dir,
    width = 4
)

# ---------------------------------------------------------------------------- #
#                                  Figure S6i                                  #
# ---------------------------------------------------------------------------- #

p <- DotPlot(
    seurat_obj,
    features = c("ZEB1", "TCF4", "SOX10", "PLP1"),
    group.by = "CellClass_L3",
    cols = c("lightgrey", "red3")
) +
    theme(axis.text.x = element_text(angle = 90))
ggsave(
    plot = p,
    filename = "FigS6i.pdf",
    path = params$plot_dir,
    height = 6,
    width = 4
)

# ---------------------------------------------------------------------------- #
#                                   Figure 4b                                  #
# ---------------------------------------------------------------------------- #

diff_comp <- read.csv(file.path(params$output_dir, "diffusion_comps.csv"))
pseudotime <- diff_comp %>% dplyr::select(dpt)
colnames(pseudotime) <- "pseudotime"
df <- cbind(seurat_obj[[]], pseudotime)

min_max_scale <- function(v) {
    return((v - min(v)) / (max(v) - min(v)))
}

df$invasive_signature <- min_max_scale(df$invasive_signature)
df$DEG3 <- min_max_scale(df$DEG3)


p <- ggplot(df, aes(x = pseudotime)) +
    geom_smooth(
        aes(y = invasive_signature),
        method = "gam",
        formula = y ~ s(x, k = 5, bs = "cs"),
        colour = "darkred",
        se = FALSE
    ) +
    geom_smooth(
        aes(y = DEG3),
        method = "gam",
        formula = y ~ s(x, k = 5, bs = "cs"),
        colour = "darkseagreen3",
        se = FALSE
    ) +
    xlab("Pseudotime") +
    ylab("Module Score") +
    theme_classic() +
    GBM_theme() +
    theme(legend.position = "right", legend.direction = "vertical")

# Density
densities <- do.call(
    rbind,
    lapply(split(df, df$CellClass_L3), function(group_data) {
        dens <- density(
            group_data$pseudotime,
            bw = 2,
            from = min(df$pseudotime),
            to = max(df$pseudotime)
        )
        data.frame(
            x = dens$x,
            density = dens$y,
            group = group_data$CellClass_L3[1]
        )
    })
)

densities <- densities %>%
    group_by(group) %>%
    mutate(density = density / max(density)) %>%
    ungroup()

densities$group <- factor(densities$group, levels = c("Pre-OPC", "OPC", "COP"))
p2 <- ggplot(densities, aes(x = x, y = group)) +
    geom_tile(aes(fill = group, alpha = density)) +
    scale_fill_manual(values = cols) +
    scale_alpha(range = c(0, 1)) +
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
    ) +
    labs(
        x = "Pseudotime",
        fill = "Group",
        alpha = "Density"
    )

p <- p / p2 + plot_layout(heights = c(4, 1))
ggsave(
    plot = p,
    filename = "Fig4b.pdf",
    path = params$plot_dir,
    height = 8,
    width = 9
)

# ---------------------------------------------------------------------------- #
#                                  Figure S6d                                  #
# ---------------------------------------------------------------------------- #

dcs <- diff_comp %>% dplyr::select(all_of(c("DC1", "DC2")))
df <- cbind(df, dcs)
df <- df %>% arrange(DC1)
p1 <- ggplot(df, aes(x = DC1, y = DC2)) +
    geom_point(
        aes(colour = CellClass_L3),
        stroke = NA,
        size = 4,
        alpha = 0.75
    ) +
    scale_colour_manual(values = cols) +
    GBM_theme() +
    guides(
        colour = guide_legend(
            title = "Cell Class",
            title.theme = element_text(size = 15),
            label.theme = element_text(size = 10),
            override.aes = list(size = 10)
        )
    )

viridis_cols <- viridis_pal()(3)
p2 <- ggplot(df, aes(x = DC1, y = DC2)) +
    geom_point(aes(colour = pseudotime), stroke = NA, size = 4, alpha = 0.75) +
    scale_colour_gradient2(
        low = viridis_cols[1],
        mid = viridis_cols[2],
        high = viridis_cols[3],
        midpoint = 2.5
    ) +
    GBM_theme() +
    guides(
        colour = guide_colorbar(
            title = "Pseudotime",
            title.theme = element_text(size = 15),
            label.theme = element_text(size = 10),
            override.aes = list(size = 10),
            barwidth = 5,
            barheight = 2,
            reverse = FALSE
        )
    )

df <- df %>% arrange(invasive_signature)
magma_cols <- viridis_pal(option = "A")(3)
p3 <- ggplot(df, aes(x = DC1, y = DC2)) +
    geom_point(aes(colour = invasive_signature), stroke = NA, size = 4) +
    scale_colour_gradient2(
        low = magma_cols[1],
        mid = magma_cols[2],
        high = magma_cols[3],
        midpoint = 0.5
    ) +
    GBM_theme() +
    guides(
        colour = guide_colorbar(
            title = "Scaled Invasive Signature",
            title.theme = element_text(size = 15),
            label.theme = element_text(size = 10),
            override.aes = list(size = 10),
            barwidth = 8,
            barheight = 2,
            reverse = FALSE
        )
    )
p <- p1 /
    p2 /
    p3
ggsave(
    plot = p,
    filename = "FigS6d.pdf",
    path = params$plot_dir,
    height = 15,
    width = 7
)
