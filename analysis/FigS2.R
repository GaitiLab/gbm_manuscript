# Code to reproduce Figure S2

if(!("pacman" %in% rownames(installed.packages()))){
  install.packages("pacman")
}

pacman::p_load(data.table,
               ggplot2,
               Seurat,
               tidyr,
               dplyr,
               patchwork,
               scales,
               here,
               GBMutils)

source(here::here("R/utils/gbm_project.R"))               

# Enter paths here
path_to_seurat_object <- ""
output_dir <- ""

# ------------------ Creating directories needed for outputs ----------------- #
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
#                                  Figure S2a                                  #
# ---------------------------------------------------------------------------- #

so$umap <- ifelse(so$is_malignant_confident == TRUE & so$Region == "PT", "PT-Malignant", "Other")
p <- DimPlot(so, reduction = "scVIumap", group.by = "umap", cols = c("PT-Malignant" = "black", "Other" = "grey"))
ggsave(p, filename = "scvi_pt_malignant_umap.pdf", path = plot_dir, height = 7, width = 9)

# EGFR expression of malignant PT vs non-malignant
so_subset <- subset(so, subset = is_malignant_confident == TRUE | CellClass_L1 != "Malignant")
metadata <- so_subset@meta.data %>% 
  mutate(vio = case_when(
  so_subset$is_malignant_confident == TRUE & so_subset$Region == "PT" ~ "PT-Malignant",
  so_subset$is_malignant_confident == TRUE & so_subset$Region == "TE" ~ "TE-Malignant",
  so_subset$is_malignant_confident == TRUE & so_subset$Region == "TC" ~ "TC-Malignant",
  so_subset$CellClass_L1 != "Malignant" ~ "Non-malignant"
))
so_subset@meta.data <- metadata
so_subset$vio <- factor(so_subset$vio, levels = c("Non-malignant", "PT-Malignant", "TE-Malignant", "TC-Malignant"))
p <- VlnPlot(so_subset, features = c("EGFR", "PTPRZ1"), pt.size = 0, group.by = "vio", cols = c("PT-Malignant" = "#7BB7AD", "TE-Malignant" = "#7BB7AD","TC-Malignant" = "#7BB7AD", "Non-malignant" = "grey"))
ggsave(filename = "PT_EGFR_PTPRZ1_expression.pdf", path = plot_dir, height = 5, width = 7)

# ---------------------------------------------------------------------------- #
#                                  Figure S2b                                  #
# ---------------------------------------------------------------------------- #
Idents(so) <- "CellClass_L1"
Idents(so) <- factor(Idents(so), levels = c("Astrocyte", "Endothelial", "Malignant", "Myeloid", "Neuron", "Oligodendrocyte", "OPC", "Pericyte", "T_cell"))

markers <- c("SLC1A2", "ADGRV1", "FLT1", "ABCB1", "DOCK8", "APBB1IP", "CNTNAP2", "SYT1", "PLP1",
    "MBP", "PCDH15", "CA10", "DCN", "EBF1", "IL7R", "SKAP1")

p <- DotPlot(so, features = markers, dot.scale = 8, assay = "RNA", scale = TRUE) +
    RotatedAxis()
ggsave(filename = "markers_dotplot.pdf", path = plot_dir, height = 5, width = 8)

# ---------------------------------------------------------------------------- #
#                                  Figure S2c                                  #
# ---------------------------------------------------------------------------- #
df <- so[[]] %>% 
  filter(CellClass_L1 != "Malignant") %>% 
  select(all_of(c("Region", "Sample", "CellClass_L2"))) %>% 
  group_by_at(c("Region", "Sample", "CellClass_L2")) %>% 
  summarise(n = n())
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))
write.csv(df, file.path(plot_dir, "non_malig_abundance_by_sample.csv"))

nm_prop <- so[[]] %>% 
  filter(CellClass_L1 != "Malignant") %>% 
  select(all_of(c("Region", "Sample", "CellClass_L2"))) %>% 
  group_by_at(c("Region", "Sample", "CellClass_L2")) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by_at(c("Region", "Sample")) %>% 
  mutate(total = sum(n), prop = n / total) %>% 
  filter(CellClass_L2 == "Oligodendrocyte") %>% 
  arrange(prop)
df$Sample <- factor(df$Sample, levels = nm_prop$Sample)

p <- ggplot(df, aes(x = Sample, y = n, fill = factor(CellClass_L2))) +
      geom_bar(position = "fill", stat = "identity") + 
      theme_classic() +
      scale_fill_manual(values = nm_cols) +      
      labs(y = "Proportion (%)") +
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
ggsave(filename = "non_malig_cell_type_prop_by_sample.pdf", path = plot_dir, height = 25, width = 25)


# ---------------------------------------------------------------------------- #
#                                  Figure S2d                                  #
# ---------------------------------------------------------------------------- #

df <- so[[]] %>% 
  filter(is_malignant_confident == TRUE | CellClass_L1 != "Malignant") %>% 
  select(all_of(c("Sample", "Region", "CellClass_L1", "Patient")))
df$CellClass_L1 <- ifelse(df$CellClass_L1 == "Malignant", "Malignant", "Non-Malignant")
df <- df %>% 
  group_by_at(c("Sample", "Region", "CellClass_L1", "Patient")) %>% 
  summarise(count = n()) %>% 
  group_by(Sample) %>% 
  mutate(total = sum(count)) %>% 
  mutate(Fraction = count/total) %>% 
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
  theme(legend.text = element_text(size = 18),
        legend.key.size = unit(1, "cm")) +
  geom_signif_lmm(
    data_df = df,
    response = "Fraction",
    condition = "Region",
    latent_vars = c("Patient"),
    comparisons = comps,
    step_increase = c(0, 0.1, 0.1)
  )
ggsave(filename = "tumour_cell_purity.pdf", path = plot_dir, width = 10, height = 10)

# ---------------------------------------------------------------------------- #
#                                  Figure S2e                                  #
# ---------------------------------------------------------------------------- #
df <- so[[]] %>% 
  filter(is_malignant_confident == TRUE) %>% 
  select(all_of(c("Region", "Sample", "CellClass_L3", "Patient"))) %>% 
  group_by_at(c("Region", "Sample", "CellClass_L3", "Patient")) %>% 
  summarise(n = n())
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))
write.csv(df, file.path(plot_dir, "malignant_cell_num_by_sample.csv"))

m_prop <- so[[]] %>% 
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
      labs(y = "Proportion (%)") +
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
ggsave(filename = "cell_type_prop_by_sample.pdf", path = plot_dir, height = 25, width = 25)

# ---------------------------------------------------------------------------- #
#                                  Figure S2f                                  #
# ---------------------------------------------------------------------------- #

so[["Progenitor_Differentiated"]] <- ifelse(so$is_malignant_confident == TRUE & so$CellClass_L2 %in% c("Malignant_MES", "Malignant_AC"), "Differentiated_like",
                                            ifelse(so$is_malignant_confident == TRUE & so$CellClass_L3 %in% c("Malignant_NPC1", "Malignant_OPC"), "Progenitor_like", 
                                            so$CellClass_L3))
comps <- list(c("PT", "TE"), c("TE", "TC"), c("PT", "TC"))

df <- so[[]] %>% 
  filter(is_malignant_confident == TRUE) %>% 
  group_by_at(c("Sample", "Region", "CellClass_L3", "Patient")) %>% 
  summarise(count = n()) %>% 
  group_by(Sample) %>% 
  mutate(total = sum(count)) %>% 
  mutate(Fraction = count/total) %>% 
  select(all_of(c("Sample", "Region", "CellClass_L3", "Fraction", "Patient"))) 

df$Region <- ifelse(df$Region == "PT", "PT", "Tumour")
df$Region <- factor(df$Region, levels = c("PT", "Tumour"))

plots <- list()
for (subtype in unique(df$CellClass_L3)) {
  curr_data <- df %>% filter(CellClass_L3 == subtype)

  plots[[subtype]] <- ggplot(curr_data, aes(x = Region, y = Fraction, fill = Region)) +
    geom_boxplot(alpha = 0.9, outlier.shape = NA) + 
    geom_point(aes(group = Patient), size = 2.5, position = position_dodge(width = 0.75)) +
    xlab("Region") +
    ylab("Proportion of Malignant Cells") +
    scale_fill_manual(values = c(region_cols[1], region_cols[3])) +
    GBM_theme() +
    theme(legend.text = element_text(size = 18),
        legend.key.size = unit(1, "cm")) +
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
ggsave(plots, filename = "malignant_sample_cell_type_props.pdf", path = plot_dir, width = 10, height = 14)

# ---------------------------------------------------------------------------- #
#                                  Figure S2h                                  #
# ---------------------------------------------------------------------------- #

df <- so[[]] %>% 
  filter(CellClass_L1 != "Malignant") %>% 
  select(all_of(c("Sample", "Region", "CellClass_L1", "Patient"))) %>% 
  filter(CellClass_L1 %in% c("Neuron", "Oligodendrocyte", "Myeloid", "T_cell")) %>% 
  group_by_at(c("Sample", "Region", "CellClass_L1", "Patient")) %>% 
  summarise(count = n()) %>% 
  group_by(Sample) %>% 
  mutate(total = sum(count)) %>% 
  mutate(Fraction = count/total) %>% 
  select(all_of(c("Sample", "Region", "CellClass_L1", "Fraction", "Patient"))) 
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))

plots <- list()
for (celltype in unique(df$CellClass_L1)) {
  curr_data <- df %>% filter(CellClass_L1 == celltype)
  plots[[celltype]] <- ggplot(curr_data, aes(x = Region, y = Fraction, fill = Region)) +
    geom_boxplot(alpha = 0.9, outlier.shape = NA) + 
    geom_point(aes(group = Sample), size = 2.5, position = position_dodge(width = 0.75)) +
    ggtitle(celltype) +
    xlab("Region") +
    ylab("Proportion of Non-Malignant Cells") +
    scale_fill_manual(values = c(region_cols)) +
    GBM_theme() +
    theme(legend.text = element_text(size = 18),
        legend.key.size = unit(1, "cm")) +
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
ggsave(plots, filename = "nonmalignant_sample_cell_type_props.pdf", path = plot_dir, height = 9, width = 9)



