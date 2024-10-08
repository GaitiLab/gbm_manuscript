# Code to reproduce Figure 1 b,c,d,e

if(!("pacman" %in% rownames(installed.packages()))){
  install.packages("pacman")
}

pacman::p_load(data.table,
               ggplot2,
               Seurat,
               tidyr,
               dplyr,
               scales,
               here,
               leiden,
               reticulate,
               igraph,
               tibble,
               GBMutils)

# Enter paths here
path_to_seurat_object <- ""
output_dir <- ""

# Creating directories needed for outputs 
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_dir <- paste0(output_dir, "/path_to_output_directory")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Loading Seurat object

so <- readRDS(path_to_seurat_object)

print(so)

DefaultAssay(so) <- "RNA"


# load color palette
region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

color_palette <- load_color_palette(name = "CellClass_L2")
to_remove <- c("Malignant_MES", "Malignant_NPC", "Malignant_OPC", "Malignant_AC")
nm_cols <- color_palette[setdiff(names(color_palette), to_remove)]
malignant_cols <- c("Malignant_AC" = "#e0b73d", "Malignant_MES_INT" = "firebrick", "Malignant_MES_HYP" = "tomato", "Malignant_MES_AST" = "salmon",  "Malignant_NPC1" = "skyblue", "Malignant_NPC2" = "lightblue", "Malignant_OPC" = "#afd0a6")
prog_like_cols <- c("Malignant_NPC1" = "skyblue", "Malignant_NPC2" = "lightblue", "Malignant_OPC" = "#afd0a6", "Invasive-high OPC/NPC1" = "#2B7095")
color_palette <- c(nm_cols, malignant_cols)

# ---------------------------------------------------------------------------- #
#                                   Figure 1b                                  #
# ---------------------------------------------------------------------------- #
# Import scVI embeddings 
dat <- read.csv("/scvi_latent_rep.csv", row.names = 1) # path to scvi latent rep
mat <- as.matrix(dat)
so[["scVI"]] <- CreateDimReducObject(embeddings = mat, key = "XscVI_", assay = "RNA")

so <- FindNeighbors(so, k.param = 20, reduction = "scVI", dims = 1:30, graph.name = c("scVI_nn", "scVI_snn"))
so <- FindClusters(so, 
                   verbose = FALSE,
                   resolution = 0.3, 
                   algorithm = 4,
                   method = "igraph", 
                   graph.name = "scVI_snn")

so <- RunUMAP(so,
              dims = 1:30,
              reduction = "scVI",
              reduction.name = "scVIumap",
              reduction.key = "scVIUMAP_",
              n.neighbors = 20)

so[["scVI_clusters0.3"]] <- so[["seurat_clusters"]]

# Annotation umap 
so$umap <- ifelse(so$CellClass_L1 == "Malignant", "Malignant", so$CellClass_L2)
p <- DimPlot(so, reduction = "scVIumap", group.by = "umap", cols = c("Malignant" = "#7BB7AD", nm_cols))
ggsave(p, filename = "scvi_l2_umap.pdf", path = plot_dir, height = 7, width = 9)

# ---------------------------------------------------------------------------- #
#                                   Figure 1c                                  #
# ---------------------------------------------------------------------------- #

df <- so[[]] %>% 
  mutate(CellClass_L0 = case_when(
    CellClass_L1 != "Malignant" ~ "Non-Malignant",
    .default = "Malignant"
  )) %>% 
  filter(CellClass_L0 == "Non-Malignant" | is_malignant_confident == TRUE) %>% 
  select(all_of(c("Region", "Sample", "CellClass_L0"))) %>% 
  group_by_at(c("Region", "Sample", "CellClass_L0")) %>% 
  summarise(n = n())
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))
df$CellClass_L0 <- factor(df$CellClass_L0, levels = c("Malignant", "Non-Malignant"))

nm_prop <- so[[]] %>% 
  mutate(CellClass_L0 = case_when(
    CellClass_L1 != "Malignant" ~ "Non-Malignant",
    .default = "Malignant"
  )) %>% 
  filter(CellClass_L0 == "Non-Malignant" | is_malignant_confident == TRUE) %>% 
  select(all_of(c("Region", "Sample", "CellClass_L0"))) %>% 
  group_by_at(c("Region", "Sample", "CellClass_L0")) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by_at(c("Region", "Sample")) %>% 
  mutate(total = sum(n), prop = n / total) %>% 
  filter(CellClass_L0 == "Non-Malignant") %>% 
  arrange(prop)
df$Sample <- factor(df$Sample, levels = nm_prop$Sample)
print(nm_prop)
print(df)

write.csv(df, file.path(plot_dir, "malig_non_malig_abundance.csv"))

p <- ggplot(df, aes(x = Sample, y = n, fill = factor(CellClass_L0))) +
      geom_bar(position = "fill", stat = "identity") + 
      theme_classic() +
      scale_fill_manual(values = c("#87B5B1", "#AA9F8B")) +      
      labs(y = "Relative proportion") +
      scale_y_continuous(labels=scales::percent) +
      facet_grid(.~Region, scales = "free_x", space = "free_x") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(size = 50),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 50),
            axis.ticks.x = element_blank(),
            axis.line = element_line(size = 2),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            strip.text = element_text(size = 50),
            strip.background = element_blank(),
            panel.spacing = unit(1, "cm"),
            legend.text = element_text(size = 35),
            legend.key.size = unit(3, "cm"),
            legend.position = "bottom",
            legend.title = element_blank()) 
ggsave(filename = "malig_vs_non_malig_prop_by_sample.pdf", path = plot_dir, height = 25, width = 25)

# ---------------------------------------------------------------------------- #
#                                  Figure 1d-e                                 #
# ---------------------------------------------------------------------------- #
all_metadata <- so[[]]

all_metadata$Region <- factor(all_metadata$Region, levels = c("PT", "TE", "TC"))

df <- all_metadata %>%
  filter(CellClass_L2 != "NoSignals" & (is_malignant_confident == TRUE | CellClass_L1 != "Malignant")) %>%
  mutate(CellClass_L0 = case_when(
    CellClass_L1 == "Malignant" ~ "Malignant",
    TRUE ~ "Non-malignant"
  )) %>% 
  group_by(Sample, CellClass_L3, CellClass_L0, Region) %>%
  summarise(n = n()) %>%
  group_by(Sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  mutate(proportion = ifelse(CellClass_L0 == "Malignant", abs(proportion), -abs(proportion)))

df$proportion <- abs(df$proportion)
df$Region <- ifelse(df$Region == "PT", "PT", "Tumour")

df <- df %>% 
  group_by_at(c("Sample", "Region", "CellClass_L0")) %>% 
  mutate(L0_total = sum(n), L0_prop = n / L0_total) %>% 
  group_by_at(c("Region", "CellClass_L3")) %>% 
  mutate(l3_n = n())
df <- df %>% 
  group_by_at(c("CellClass_L3", "Region", "CellClass_L0", "l3_n")) %>% 
  summarise(mean_prop = mean(L0_prop), sd_prop = sd(L0_prop)) %>% 
  mutate(sem_prop = sd_prop / l3_n)

m <- df %>% 
  filter(CellClass_L0 == "Malignant") %>% 
  select(all_of(c("CellClass_L3", "mean_prop", "Region", "sem_prop", "sd_prop", "l3_n"))) %>% 
  pivot_wider(names_from = "Region", values_from = c("mean_prop", "sem_prop", "sd_prop", "l3_n")) %>% 
  mutate(diff = mean_prop_PT - mean_prop_Tumour, se_diff = sqrt((sem_prop_PT^2/l3_n_PT) + (sem_prop_Tumour^2/l3_n_Tumour))) %>% 
  arrange(diff)
m$CellClass_L3 <- factor(m$CellClass_L3, levels = m$CellClass_L3)

ggplot(m, aes(x = mean_prop_PT, y = mean_prop_Tumour)) +
  geom_point(aes(fill = CellClass_L3), stroke = 1, colour = "black", shape = 21, size = 5) +
  scale_fill_manual(values = malignant_cols) +
  ylim(c(0, 0.5)) +
  xlim(c(0, 0.5)) +
  ylab("Relative proportion of malignant cells in Tumor (TE+TC)") +
  xlab("Relative proportion of malignant cells in PT") +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  theme(
    legend.direction = "vertical", 
    legend.position = "right"
  )
ggsave(filename = "m_abundance_scatter.pdf", path = plot_dir, height = 6, width = 8)

nm <- df %>% 
  filter(CellClass_L0 == "Non-malignant") %>% 
  select(all_of(c("CellClass_L3", "mean_prop", "Region", "sem_prop", "sd_prop", "l3_n"))) %>% 
  pivot_wider(names_from = "Region", values_from = c("mean_prop", "sem_prop", "sd_prop", "l3_n")) %>% 
  mutate(diff = mean_prop_PT - mean_prop_Tumour, se_diff = sqrt((sem_prop_PT^2/l3_n_PT) + (sem_prop_Tumour^2/l3_n_Tumour))) %>% 
  arrange(diff)
nm$CellClass_L3 <- factor(nm$CellClass_L3, levels = nm$CellClass_L3)

ggplot(nm, aes(x = mean_prop_PT, y = mean_prop_Tumour)) +
  geom_point(aes(fill = CellClass_L3), stroke = 1, colour = "black", shape = 21, size = 5) +
  scale_fill_manual(values = nm_cols) +
  ylim(c(0, 0.5)) +
  xlim(c(0, 0.5)) +
  ylab("Relative proportion of malignant cells in Tumor (TE+TC)") +
  xlab("Relative proportion of malignant cells in PT") +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  theme(
    legend.direction = "vertical", 
    legend.position = "right"
  )
ggsave(filename = "nm_abundance_scatter.pdf", path = plot_dir, height = 6, width = 8)
