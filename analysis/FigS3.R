# Code to reproduce Figure S3

if(!("pacman" %in% rownames(installed.packages()))){
  install.packages("pacman")
}

pacman::p_load(data.table,
               patchwork,
               ggplot2,
               dplyr,
               tidyr,
               ComplexHeatmap,
               grid,
               colorRamp2,
               fgsea,
               GBMutils,
               scales)

region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

# Enter paths here
output_dir <- ""

# Creating directories needed for outputs
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_dir <- paste0(output_dir, "/path_to_output_directory")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---------------------------------------------------------------------------- #
#                                  Figure S3a                                  #
# ---------------------------------------------------------------------------- #

# Load gene sets

invasivity <- read.csv("/invasivity.csv") # Venkataramani 2022 
garofano <- read.csv("/garofano.csv") # Garofano 2021
neftel_markers <- read.csv(file.path("/Neftel_gene_list.csv")) # Neftel 2019

invasivity_up <- invasivity %>% 
  filter(direction == "Anticorrelated")
invasivity_up <- list(invasivity_up$Gene)
names(invasivity_up) <- "Invasivity_up"

garofano <- as.vector(garofano)

neftel_gene_list <- list()
for(i in 1:nrow(neftel_markers)) {
neftel_gene_list[[toString(neftel_markers[i,"Cell"])]] <- append(neftel_gene_list[[toString(neftel_markers[i,"Cell"])]],toString(neftel_markers[i, "Gene"]))
}
cell_type <- c()
for(i in 1:length(neftel_gene_list)){
cell_type <- c(cell_type, paste0("Neftel_", names(neftel_gene_list)[i]))
}
names(neftel_gene_list) <- cell_type

markers <- read.csv("consensus_mps_markers.csv", row.names = 1) # load consensus factor markers
markers <- as.vector(markers)
names(markers) <- paste0("Factor", seq(5,1))

gene_lists <- c(neftel_gene_list, garofano, invasivity_up, markers)

jaccard_similarity <- function(vec1, vec2) {
  intersection <- length(intersect(vec1, vec2))
  union <- length(unique(c(vec1, vec2)))
  return(intersection / union)
}

n <- length(gene_lists)

similarity_matrix <- matrix(NA, n, n)
rownames(similarity_matrix) <- colnames(similarity_matrix) <- names(gene_lists)

for (i in 1:n) {
  for (j in 1:n) {
    similarity_matrix[i, j] <- jaccard_similarity(gene_lists[[i]], gene_lists[[j]])
  }
}

pdf(file = paste0(plot_dir, "/jaccard_matrix.pdf"), width = 10, height = 10)
pushViewport(viewport(gp = gpar(fontfamily = "Helvetica")))
heatmap <- Heatmap(similarity_matrix, 
        name = "Jaccard Similarity", 
        na_col = "white", 
        height = nrow(similarity_matrix)*unit(1, "cm"),
        width = ncol(similarity_matrix)*unit(1, "cm"),
        col = colorRamp2(c(0, 0.375, 0.75), c("white", "orange", "red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE, 
        show_heatmap_legend = TRUE,
        show_row_names = TRUE, 
        show_column_names = TRUE)
draw(heatmap, newpage = FALSE)
popViewport()
dev.off()

# ---------------------------------------------------------------------------- #
#                                  Figure S3b                                  #
# ---------------------------------------------------------------------------- #

metadata <- read.csv("metacell_metadata.csv", row.names = 1)
metadata$Region <- factor(metadata$Region, levels = c("PT", "TE", "TC"))

df <- metadata
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))
df <- df %>% 
  select(all_of(c("Region", "Patient", paste0("Factor", seq(5,1))))) %>% 
  pivot_longer(-c("Region", "Patient"), names_to = "Factor", values_to = "Activation")
plots <- list()
for (factor in unique(df$Factor)) {
  curr_df <- df %>% filter(Factor == factor)
  plots[[factor]] <- ggplot(curr_df, aes(x = Region, y = Activation, fill = Region)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA) +
  scale_fill_manual(values = c(region_cols)) +
  GBM_theme() +
  theme(legend.position = "none") + 
  ylab("Factor Activation") +
  ggtitle(factor) + 
  geom_signif_lmm(
    data_df = curr_df,
    response = "Activation",
    condition = "Region",
    latent_vars = c("Patient"),
    comparisons = list(c("PT", "TE"), c("TE", "TC"), c("PT", "TC")),
    step_increase = c(0, 0.1, 0.1)
  )
}
plots <- wrap_plots(plots, ncol = 5, nrow = 1, guides = "collect")
ggsave(plots, filename = "Factors_boxplot.pdf", path = plot_dir, width = 12, height = 5)

# ---------------------------------------------------------------------------- #
#                                  Figure S3c                                  #
# ---------------------------------------------------------------------------- #
markers <- as.vector(as.data.frame(markers))

factors <- read.csv("/consensus_mps.csv") # load conensus factors 

ora_dir <- file.path(plot_dir, "/ora/")
if(!dir.exists(ora_dir)) dir.create(ora_dir, recursive = TRUE)

msigdb_cat_list <- list(c("H"),
                        c("C2", "CGP"),
                        c("C2", "CP"),
                        c("C5", "GO:BP"),
                        c("C5", "GO:CC"),
                        c("C5", "GO:MF"))

msigdb_data <- fread("msigdb.csv") # path to msigdb gene sets

universe_genes <- rownames(factors)

for(curr_factor in names(markers)){
  
  all_sig_results <- c()
  for(j in seq(1, length(msigdb_cat_list))){
    curr_cat <- msigdb_cat_list[[j]]
    
    if(length(curr_cat) == 2){
      curr_gene_sets <- msigdb_data[gs_cat == curr_cat[1] & gs_subcat == curr_cat[2]]
    }else{
      curr_gene_sets <- msigdb_data[gs_cat == curr_cat[1]]
    }
    
    curr_msigdbr_list <- split(x = curr_gene_sets$gene_symbol, f = curr_gene_sets$gs_name)
    
    curr_ora_results <- fora(curr_msigdbr_list, markers[[curr_factor]], universe_genes)
    sig_ora_results <- curr_ora_results[padj < 0.1]
    
    if(nrow(sig_ora_results) > 10){
      setorder(sig_ora_results, padj)
      sig_ora_results <- sig_ora_results[1:10]
    }
    
    sig_ora_results$msigdb_cat <- paste0(curr_cat, collapse = "_")
    all_sig_results <- rbindlist(list(all_sig_results, sig_ora_results))
  }
  
    all_sig_results <- all_sig_results %>% arrange(padj)
    fwrite(all_sig_results, file = paste0(ora_dir, "/", curr_factor, ".csv"))
}

ora_res <- list.files("/ora", full.names = TRUE)
ora_res <- ora_res[grepl(".csv", ora_res)]

plots <- list()
for (path in ora_res) {
  res <- read.csv(path)
  res <- res[1:5,]

  res$pct_overlap <- res$overlap / res$size

  factor <- gsub(".csv", "", basename(path))
  title <- case_when(
    factor == "Factor5" ~ "Factor 5 - Hypoxia", # ordering is modified to align with factor-factor similarity heatmap.
    factor == "Factor4" ~ "Factor 4 - Neuronal",
    factor == "Factor3" ~ "Factor 3 - Injury",
    factor == "Factor2" ~ "Factor 2 - Cilia",
    factor == "Factor1" ~ "Factor 1 - Cell Cycle",
  )

  if (factor == "Factor5") {
    res$pathway <- case_when(
      res$pathway == "MENSE_HYPOXIA_UP" ~ "Mense et al. \n Hypoxia",
      res$pathway == "HALLMARK_MTORC1_SIGNALING" ~ "Hallmark \n MTORC1 Signaling",
      res$pathway == "HALLMARK_HYPOXIA" ~ "Hallmark \n Hypoxia",
      res$pathway == "ELVIDGE_HYPOXIA_BY_DMOG_UP" ~ "Elvidge et al. \n Hypoxia by DMOG",
      res$pathway == "ELVIDGE_HIF1A_AND_HIF2A_TARGETS_DN" ~ "Elvidge et al. \n HIF1A and HIF2A Targets",
    )
  } else if (factor == "Factor4") {
    res$pathway <- case_when(
      res$pathway == "GOCC_SYNAPTIC_MEMBRANE" ~ "GO:CC \n Synaptic Membrane",
      res$pathway == "GOCC_POSTSYNAPTIC_MEMBRANE" ~ "GO:CC \n Postsynaptic Membrane",
      res$pathway == "BENPORATH_ES_WITH_H3K27ME3" ~ "Benporath et al. \n ES with H3K27Me3",
      res$pathway == "GOCC_INTRINSIC_COMPONENT_OF_POSTSYNAPTIC_MEMBRANE" ~ "GO:CC \n Component of Postsynaptic Membrane",
      res$pathway == "GOCC_CATION_CHANNEL_COMPLEX" ~ "GO:CC \n Cation Channel Complex",
    )
  } else if (factor == "Factor3") {
    res$pathway <- case_when(
      res$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" ~ "Hallmark \n TNFa signaling via NFkB",
      res$pathway == "BROCKE_APOPTOSIS_REVERSED_BY_IL6" ~ "Brocke et al. \n Apoptosis reversed by IL6",
      res$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" ~ "Hallmark \n EMT",
      res$pathway == "GOCC_ACTOMYOSIN" ~ "GO:CC \n Actomyosin",
      res$pathway == "MCDOWELL_ACUTE_LUNG_INJURY_UP" ~ "McDowell et al. \n Acute Lung Injury",
    )
  } else if (factor == "Factor2") {
    res$pathway <- case_when(
      res$pathway == "GOBP_CILIUM_MOVEMENT" ~ "GO:BP \n Cilium Movement",
      res$pathway == "GOCC_CILIUM" ~ "GO:CC \n Cilium",
      res$pathway == "GOCC_9PLUS2_MOTILE_CILIUM" ~ "GO:CC \n 9+2 Motile Cilium",
      res$pathway == "GOCC_MOTILE_CILIUM" ~ "GO:CC \n Motile Cilium",
      res$pathway == "GOBP_MICROTUBULE_BASED_MOVEMENT" ~ "GO:BP \n Microtubule-based Movement",
    )
  } else {
    res$pathway <- case_when(
      res$pathway == "FISCHER_DREAM_TARGETS" ~ "Fischer et al. \n Dream targets",
      res$pathway == "MARSON_BOUND_BY_E2F4_UNSTIMULATED" ~ "Marson et al. \n Bound by E2F4 Unstimulated",
      res$pathway == "GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP" ~ "Gobert et al. \n Oligodendrocyte Differentiation Up",
      res$pathway == "FLORIO_NEOCORTEX_BASAL_RADIAL_GLIA_DN" ~ "Elvidge et al. \n Neocortex Apical Basal Radial Glia",
      res$pathway == "FISCHER_G2_M_CELL_CYCLE" ~ "Fischer et al. \n G2M Cell Cycle",
    )
  }
  p <- ggplot(res, aes(x = -log10(padj), y = reorder(pathway, -log10(padj)))) +
    geom_segment(aes(xend = 0, yend = pathway), color = "royalblue") +
    geom_point(aes(size = pct_overlap), color = "royalblue") +
    scale_size_continuous(range = c(4, 10)) +
    GBM_theme() +
    labs(title = title,
        x = "-log10(FDR)",
        y = "",
        size = "Gene set overlap") +
    theme(legend.position = "left",
          legend.direction = "vertical")
  p <- grid.arrange(p, widths = unit(9, "in"))

  plots[[factor]] <- p

  ggsave(plot = p, filename = paste0(factor, "_ora_plot.pdf"), path = ora_dir, width = 10, height = 5)
}

plots <- wrap_plots(plots, ncol = 3, nrow = 2, guides = "collect")
ggsave(plot = plots, filename = "factors_ora.pdf", path = plot_dir, width = 28, height = 10)

# ---------------------------------------------------------------------------- #
#                                  Figure S4d                                  #
# ---------------------------------------------------------------------------- #

colnames(metadata) <- gsub("Neftel1", "MES2", colnames(metadata))
colnames(metadata) <- gsub("Neftel2", "MES1", colnames(metadata))
colnames(metadata) <- gsub("Neftel3", "AC", colnames(metadata))
colnames(metadata) <- gsub("Neftel4", "OPC", colnames(metadata))
colnames(metadata) <- gsub("Neftel5", "NPC1", colnames(metadata))
colnames(metadata) <- gsub("Neftel6", "NPC2", colnames(metadata))

dfs <- lapply(unique(metadata$Patient), function(patient) {
  curr_df <- metadata %>% filter(Patient == patient)

  curr_df$MES1 <- scale(curr_df$MES1)
  curr_df$MES2 <- scale(curr_df$MES2)
  curr_df$AC <- scale(curr_df$AC)
  curr_df$OPC <- scale(curr_df$OPC)
  curr_df$NPC1 <- scale(curr_df$NPC1)
  curr_df$NPC2 <- scale(curr_df$NPC2)
  curr_df$Factor1 <- scale(curr_df$Factor1)
  curr_df$Factor2 <- scale(curr_df$Factor2)
  curr_df$Factor3 <- scale(curr_df$Factor3)
  curr_df$Factor4 <- scale(curr_df$Factor4)
  curr_df$Factor5 <- scale(curr_df$Factor5)

  return(curr_df)
})
plot_df <- rbindlist(dfs)

for (state in c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2")) {

  cor_test <- cor.test(plot_df[[state]], plot_df[["Factor4"]], method = "pearson")
  print(paste0(state, " corr:", cor_test$estimate, ", pval:", cor_test$p.value))

  p <- ggplot(plot_df, aes(x = Factor4, y = get(state))) +
    geom_point(aes(colour = Region), size = 3, alpha = 0.2, stroke = NA) +
    stat_density_2d(aes(fill = Region, alpha = ..level..), geom = "polygon", contour_var = "ndensity", bins = 5) +
    scale_alpha_continuous(range = c(0.1, 0.5)) +
    scale_fill_manual(values = c(region_cols)) +
    scale_colour_manual(values = c(region_cols)) +
    ylab(paste0("Scaled ", state, " score")) +
    xlab("Scaled Factor 4 (Neuronal) activation") +
    GBM_theme()

  ggsave(filename = paste0(state, "_score_vs_Factor4_scaled.pdf"), path = plot_dir)
}
