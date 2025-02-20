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
malignant_l2 <- c("Malignant_MES", "Malignant_NPC", "Malignant_OPC", "Malignant_AC")
nm_cols <- c(color_palette[setdiff(names(color_palette), malignant_l2)], "Myeloid" = "#787335")
malignant_cols <- c("Malignant_AC" = "#e0b73d", "Malignant_MES_INT" = "firebrick", "Malignant_MES_HYP" = "tomato", "Malignant_MES_AST" = "salmon",  "Malignant_NPC1" = "skyblue", "Malignant_NPC2" = "lightblue", "Malignant_OPC" = "#afd0a6")
prog_like_cols <- c("Malignant_NPC1" = "skyblue", "Malignant_NPC2" = "lightblue", "Malignant_OPC" = "#afd0a6", "Invasive-high OPC/NPC1" = "#2B7095")
color_palette <- c(nm_cols, malignant_cols)

# ---------------------------------------------------------------------------- #
#                                   Figure 1b                                  #
# ---------------------------------------------------------------------------- #
# Import umap embeddings 
dat <- read.csv("snrnaseq_umap_coords.csv", row.names = 1) # path to scvi latent rep
mat <- as.matrix(dat)
so[["scVIumap"]] <- CreateDimReducObject(embeddings = mat, key = "XscVIumap_", assay = "RNA")

# Annotation umap 
p <- DimPlot(so, reduction = "scVIumap", group.by = "CellClass_L1", cols = c("Malignant" = "#7BB7AD", nm_cols))
ggsave(p, filename = "scvi_umap.pdf", path = plot_dir, height = 7, width = 9)

# ---------------------------------------------------------------------------- #
#                                   Figure 1c                                  #
# ---------------------------------------------------------------------------- #

df <- so[[]] %>% 
  mutate(CellClass_L0 = case_when(
    CellClass_L1 != "Malignant" ~ "Non-Malignant",
    .default = "Malignant"
  )) %>% 
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
  select(all_of(c("Region", "Sample", "CellClass_L0"))) %>% 
  group_by_at(c("Region", "Sample", "CellClass_L0")) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by_at(c("Region", "Sample")) %>% 
  mutate(total = sum(n), prop = n / total) %>% 
  filter(CellClass_L0 == "Non-Malignant") %>% 
  arrange(prop)
df$Sample <- factor(df$Sample, levels = nm_prop$Sample)

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
            axis.line = element_line(size = 1),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            strip.text = element_text(size = 50),
            strip.background = element_blank(),
            panel.spacing = unit(1, "cm"),
            legend.text = element_text(size = 35),
            legend.key.size = unit(1, "cm"),
            legend.position = "bottom",
            legend.title = element_blank()) 
ggsave(filename = "malig_vs_non_malig_prop_by_sample.pdf", path = plot_dir, height = 25, width = 25)

# ---------------------------------------------------------------------------- #
#                                  Figure 1d-e                                 #
# ---------------------------------------------------------------------------- #
all_metadata <- so[[]]

m <- all_metadata %>% 
  filter(is_malignant_confident == TRUE) %>% 
  mutate(Region = ifelse(Region == "PT", "PT", "Tumor")) %>% 
  group_by_at(c("Sample", "CellClass_L3", "Region")) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by(Sample) %>% 
  mutate(sample_total = sum(n), percent = (n / sample_total) * 100) %>% 
  group_by_at(c("CellClass_L3", "Region")) %>% 
  summarise(mean_percent = mean(percent)) %>% 
  pivot_wider(names_from = "Region", values_from = "mean_percent") 
m$CellClass_L3 <- factor(m$CellClass_L3, levels = m$CellClass_L3)

ggplot(m, aes(x = PT, y = Tumor)) +
  geom_point(aes(fill = CellClass_L3), stroke = 1, colour = "black", shape = 21, size = 5) +
  scale_fill_manual(values = malignant_cols) +
  ylim(c(0, 50)) +
  xlim(c(0, 50)) +
  ylab("Relative percent of malignant cells in Tumor (TE+TC)") +
  xlab("Relative percent of malignant cells in PT") +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  theme(
    legend.direction = "vertical", 
    legend.position = "right"
  )
ggsave(filename = "m_abundance_scatter.pdf", path = plot_dir, height = 6, width = 8)

nm <- all_metadata %>% 
  filter(CellClass_L1 != "Malignant") %>% 
  mutate(Region = ifelse(Region == "PT", "PT", "Tumor")) %>% 
  group_by_at(c("Sample", "CellClass_L1", "Region")) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by(Sample) %>% 
  mutate(sample_total = sum(n), percent = (n / sample_total) * 100) %>% 
  group_by_at(c("CellClass_L1", "Region")) %>% 
  summarise(mean_percent = mean(percent)) %>% 
  pivot_wider(names_from = "Region", values_from = "mean_percent") 

nm$CellClass_L1 <- factor(nm$CellClass_L1, levels = nm$CellClass_L1)

ggplot(nm, aes(x = PT, y = Tumor)) +
  geom_point(aes(fill = CellClass_L1), stroke = 1, colour = "black", shape = 21, size = 5) +
  scale_fill_manual(values = nm_cols) +
  ylim(c(0, 65)) +
  xlim(c(0, 65)) +
  ylab("Relative percent of malignant cells in Tumor (TE+TC)") +
  xlab("Relative percent of malignant cells in PT") +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  theme(
    legend.direction = "vertical", 
    legend.position = "right"
  )
ggsave(filename = "nm_abundance_scatter.pdf", path = plot_dir, height = 6, width = 8)
