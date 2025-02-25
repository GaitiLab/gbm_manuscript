# Load seurat object
curr_seurat_data <- readRDS(args$input)
DefaultAssay(curr_seurat_data) <- "RNA"

# Invasive signature
inv_sig <- read.csv("misc/data/inv_sig.csv") %>% column_to_rownames("X")
inv_up_list <- list(inv_sig$inv_up)
inv_down_list <- list(inv_sig$inv_down)

# Add invasive signature score to the Seurat object
curr_seurat_data <- AddModuleScore(curr_seurat_data, features = inv_up_list, name = "inv_up", nbin = 15)
curr_seurat_data <- AddModuleScore(curr_seurat_data, features = inv_down_list, name = "inv_down", nbin = 15)
curr_seurat_data$invasive_score <- curr_seurat_data$inv_up1 - curr_seurat_data$inv_down1

# Get metadata and calculate mean invasive score
metadata <- curr_seurat_data[[]]
metadata_malignant <- metadata %>% 
  select(Sample, invasivity, invasive_score, CellClass_L1, Region, Patient) %>%
  filter(CellClass_L1 == "Malignant") %>%
  mutate(Region = case_when(
    Region == "PT" ~ "Peri-tumoral",
    Region == "TE" ~ "Tumor edge",
    Region == "TC" ~ "Tumor core",
  ))

color_palette <- c("Peri-tumoral" = "#0173b2", "Tumor edge" = "#de8f05", "Tumor core" = "#029e73")

metadata_malignant <- metadata_malignant %>%
  group_by(Sample, Region) %>%
  summarise(across(c(invasivity, invasive_score), mean, na.rm = TRUE))

p <- ggscatter(metadata_malignant, x = "invasivity", y = "invasive_score", 
               xlab = "Mean invasivity signature module score", ylab = "Mean Neuronal signature module score",
                add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman") +
     geom_point(aes(color = Region)) +
     scale_color_manual(values = color_palette)
ggsave(paste0(args$output_dir, "/invasivity_invasive_score_correlation_malignant.pdf"), width=12, height=12, units = "cm")
