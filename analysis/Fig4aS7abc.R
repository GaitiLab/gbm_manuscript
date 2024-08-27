# Code to reproduce Fig 4a and S7a,b,c

if(!("pacman" %in% rownames(installed.packages()))){
  install.packages("pacman")
}

pacman::p_load(data.table,
               patchwork,
               ggplot2,
               dplyr,
               tidyr,
               ggpubr,
               Seurat,
               harmony,
               destiny,
               scales,
               sceasy)

region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

# Enter paths here
path_to_seurat_object <- ""
output_dir <- ""

# Creating directories needed for outputs
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_dir <- paste0(output_dir, "/path_to_output_directory")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---------------------------------------------------------------------------- #
#                                 Convert h5ad                                 #
# ---------------------------------------------------------------------------- #
h5ad_file <- "/cluster/projects/gaitigroup/Users/Benson/Parsebio/scanpyobjects/developing_opc.h5ad"
so <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
                       outFile='filename.rds')
print(so)
print(colnames(so[[]]))

so$CellClass_L3 <- case_when(
    so$Clusters %in% c("609", "610") ~ "COP",
    so$Clusters %in% paste0(seq(611, 616)) ~ "OPC",
    so$Clusters == "19" ~ "Pre-OPC"
)
so$CellClass_L3 <- factor(so$CellClass_L3, levels = c("Pre-OPC", "OPC", "COP"))

# ---------------------------------------------------------------------------- #
#                            Integrate with harmony                            #
# ---------------------------------------------------------------------------- #
print(head(so[[]]))

so <- NormalizeData(so)
so <- FindVariableFeatures(so, nfeatures = 2000)
so <- ScaleData(so, vars.to.regress = c("CellCycle"))
so <- RunPCA(so, npcs = 30)

so <- RunHarmony(so, group.by.vars = c("Platform", "Sample"), reduction.use = "pca")

harmony_embeddings <- Embeddings(so, reduction = "harmony")

# ---------------------------------------------------------------------------- #
#                               DPT with harmony                               #
# ---------------------------------------------------------------------------- #

dm <- DiffusionMap(data = harmony_embeddings)
dm_coords <- eigenvectors(dm)
dm_coords <- dm_coords[, c(1,2)]
colnames(dm_coords) <- c("DC1", "DC2")
df <- cbind(dm_coords, so[[]] %>% select(CellClass_L3))

df$rank <- rank(-df$DC1)
index <- 1:length(df$DC1)
dpt <- DPT(dm, tips = index[df$rank == 1])
df$dpt <- dpt$dpt

write.csv(df, file.path(plot_dir, "diffusion_comps.csv"))

# ---------------------------------------------------------------------------- #
#                                 Figure S7b,c                                 #
# ---------------------------------------------------------------------------- #

# Load gene sets
degs <- read.csv("/de_results.csv")
degs_up <- degs %>% 
  filter(log2FoldChange > 1 & padj < 0.05)
degs_dn <- degs %>% 
  filter(log2FoldChange < -1 & padj < 0.05)

neftel_markers <- read.csv(file.path("/Neftel_gene_list.csv")) # Neftel 2019

neftel_gene_list <- list()
for(i in 1:nrow(neftel_markers)) {
neftel_gene_list[[toString(neftel_markers[i,"Cell"])]] <- append(neftel_gene_list[[toString(neftel_markers[i,"Cell"])]],toString(neftel_markers[i, "Gene"]))
}
cell_type <- c()
for(i in 1:length(neftel_gene_list)){
cell_type <- c(cell_type, paste0("Neftel_", names(neftel_gene_list)[i]))
}
names(neftel_gene_list) <- cell_type

opc <- neftel_gene_list$Neftel_OPC

postsynapse <- gmtPathways("/GOCC_POSTSYNAPTIC_MEMBRANE.v2023.2.Hs.gmt")
postsynapse <- postsynapse$GOCC_POSTSYNAPTIC_MEMBRANE
synapse <- gmtPathways("/GOCC_SYNAPTIC_MEMBRANE.v2023.2.Hs.gmt")
synapse <- synapse$GOCC_SYNAPTIC_MEMBRANE
synaptic_signaling <- gmtPathways("/GOBP_SYNAPTIC_SIGNALING.v2023.2.Hs.gmt")
synaptic_signaling <- synaptic_signaling$GOBP_SYNAPTIC_SIGNALING

gene_list <- list(degs_up$gene, degs_dn$gene, opc, postsynapse, synapse, synaptic_signaling)

# Score cells in normal opc lineage
so <- AddModuleScore(so, features = gene_list, name = "DEG")
so$invasive_signature <- so$DEG1 - so$DEG2

cols <- c("#09283CFF", "#F18B00FF", "#F2EBBBFF")

p <- ggplot(so[[]], aes(x = CellClass_L3, y = invasive_signature)) +
    geom_boxplot(aes(fill = CellClass_L3)) +
    GBM_theme() + 
    stat_compare_means(comparisons = list(c("Pre-OPC", "OPC"), c("OPC", "COP"), c("Pre-OPC", "COP"))) +
    scale_fill_manual(values = cols)
ggsave(filename = "inv_sig_in_opc_lineage.pdf", path = plot_dir, width = 4)

p <- ggplot(so[[]], aes(x = CellClass_L3, y = DEG3)) +
    geom_boxplot(aes(fill = CellClass_L3)) +
    GBM_theme() + 
    stat_compare_means(comparisons = list(c("Pre-OPC", "OPC"), c("OPC", "COP"), c("Pre-OPC", "COP"))) +
    scale_fill_manual(values = cols)
ggsave(filename = "neftel_opc_in_opc_lineage.pdf", path = plot_dir, width = 4)

p <- ggplot(so[[]], aes(x = CellClass_L3, y = DEG4)) +
    geom_boxplot(aes(fill = CellClass_L3)) +
    GBM_theme() + 
    stat_compare_means(comparisons = list(c("Pre-OPC", "OPC"), c("OPC", "COP"), c("Pre-OPC", "COP"))) +
    scale_fill_manual(values = cols)
ggsave(filename = "postsynapse_in_opc_lineage.pdf", path = plot_dir, width = 4)

p <- ggplot(so[[]], aes(x = CellClass_L3, y = DEG5)) +
    geom_boxplot(aes(fill = CellClass_L3)) +
    GBM_theme() + 
    stat_compare_means(comparisons = list(c("Pre-OPC", "OPC"), c("OPC", "COP"), c("Pre-OPC", "COP"))) +
    scale_fill_manual(values = cols)
ggsave(filename = "synapse_in_opc_lineage.pdf", path = plot_dir, width = 4)

p <- ggplot(so[[]], aes(x = CellClass_L3, y = DEG6)) +
    geom_boxplot(aes(fill = CellClass_L3)) +
    GBM_theme() + 
    stat_compare_means(comparisons = list(c("Pre-OPC", "OPC"), c("OPC", "COP"), c("Pre-OPC", "COP"))) +
    scale_fill_manual(values = cols)
ggsave(filename = "synaptic_signaling_in_opc_lineage.pdf", path = plot_dir, width = 4)

# ---------------------------------------------------------------------------- #
#                                   Figure 4a                                  #
# ---------------------------------------------------------------------------- #

diff_comp <- read.csv("/diffusion_comps.csv")
pseudotime <- diff_comp %>% select(dpt)
colnames(pseudotime) <- "pseudotime"
df <- cbind(so[[]], pseudotime)

min_max_scale <- function(v) {
  return((v - min(v)) / (max(v) - min(v)))
}

df$invasive_signature <- min_max_scale(df$invasive_signature)
df$DEG6 <- min_max_scale(df$DEG6)
df$DEG7 <- min_max_scale(df$DEG7)
df$DEG3 <- min_max_scale(df$DEG3)
df$DEG4 <- min_max_scale(df$DEG4)

p <- ggplot(df, aes(x = pseudotime)) +
  geom_smooth(aes(y = invasive_signature), method = "loess", se = TRUE, colour = "darkred") +
  geom_smooth(aes(y = DEG3), method = "loess", se = TRUE, colour = "darkseagreen3") + 
  xlab("Pseudotime") +
  ylab("Module Score") +
  GBM_theme() +
  theme(legend.position = "right", 
    legend.direction = "vertical")

p2 <- ggplot(df, aes(x = CellClass_L3, y = pseudotime)) +
  geom_boxplot(aes(fill = CellClass_L3)) +
  scale_fill_manual(values = cols) +
  ylab("Pseudotime") +
  coord_flip() +
  GBM_theme() +
  theme(legend.position = "right", 
    legend.direction = "vertical")

p / p2 + plot_layout(heights = c(5, 2))
ggsave(filename = "pseudotime.pdf", path = plot_dir, height = 8, width = 9)

# ---------------------------------------------------------------------------- #
#                                  Figure S7a                                  #
# ---------------------------------------------------------------------------- #

dcs <- diff_comp %>% select(all_of(c("DC1", "DC2")))
df <- cbind(df, dcs)
df <- df %>% arrange(DC1)
p1 <- ggplot(df, aes(x = DC1, y = DC2)) +
  geom_point(aes(colour = CellClass_L3), stroke = NA, size = 4, alpha = 0.75) +
  scale_colour_manual(values = cols) +
  GBM_theme() +
  guides(colour = guide_legend(title = "Cell Class",
                          title.theme = element_text(size = 15),
                          label.theme = element_text(size = 10),
                          override.aes = list(size = 10)))

viridis_cols <- viridis_pal()(3)
p2 <- ggplot(df, aes(x = DC1, y = DC2)) +
  geom_point(aes(colour = pseudotime), stroke = NA, size = 4, alpha = 0.75) +
  scale_colour_gradient2(low = viridis_cols[1], mid = viridis_cols[2], high = viridis_cols[3], midpoint = 2.5) +
  GBM_theme() +
  guides(colour = guide_colorbar(title = "Pseudotime",
                            title.theme = element_text(size = 15),
                            label.theme = element_text(size = 10),
                            override.aes = list(size = 10),
                            barwidth = 5,
                            barheight = 2,
                            reverse = FALSE))

df <- df %>% arrange(invasive_signature)
magma_cols <- viridis_pal(option = "A")(3)
p3 <- ggplot(df, aes(x = DC1, y = DC2)) +
  geom_point(aes(colour = invasive_signature), stroke = NA, size = 4) +
  scale_colour_gradient2(low = magma_cols[1], mid = magma_cols[2], high = magma_cols[3], midpoint = 0.5) +
  GBM_theme() +
  guides(colour = guide_colorbar(title = "Scaled Invasive Signature",
                            title.theme = element_text(size = 15),
                            label.theme = element_text(size = 10),
                            override.aes = list(size = 10),
                            barwidth = 8,
                            barheight = 2,
                            reverse = FALSE))
p1 /
p2 /
p3
ggsave(filename = "dpt.pdf", path = plot_dir, height = 15, width = 7)