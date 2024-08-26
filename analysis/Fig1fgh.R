# Code to reproduce Fig 1 f,g,h

if(!("pacman" %in% rownames(installed.packages()))){
  install.packages("pacman")
}

pacman::p_load(data.table,
               patchwork,
               ggplot2,
               dplyr,
               ComplexHeatmap,
               grid,
               colorRamp2,
               dendsort,
               lsa, 
               tibble,
               purrr,
               GBMutils
               )

region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

# Enter paths here
path_to_seurat_object <- ""
output_dir <- ""

# Creating directories needed for outputs
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_dir <- paste0(output_dir, "/path_to_output_directory")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---------------------------------------------------------------------------- #
#                                  Figure f                                    #
# ---------------------------------------------------------------------------- #

# consensus by platform

multiome_factors <- read.table("/multiome_metacells.spectra.k_11.dt_0_5.consensus.txt") # path to cNMF output for multiome

multiome_geps <- lapply(seq(1,length(rownames(multiome_factors))), function(gep_num) {
  gep <- as.vector(multiome_factors[gep_num,])
})

names(multiome_geps) <- paste0("MultiomeFactor", seq(1, length(multiome_geps)))

parse_factors <- read.table("/parse_metacells.spectra.k_10.dt_0_5.consensus.txt") # path to cNMF output for evercode

parse_geps <- lapply(seq(1,length(rownames(parse_factors))), function(gep_num) {
  gep <- as.vector(parse_factors[gep_num,])
})

names(parse_geps) <- paste0("ParseFactor", seq(1, length(parse_geps)))

filtered_geps <- list()  
for (i in seq_along(parse_geps)) {
  for (j in seq_along(multiome_geps)) {
    parse_gep <- parse_geps[[i]]
    multiome_gep <- multiome_geps[[j]]
    if (cosine(as.numeric(parse_gep), as.numeric(multiome_gep)) >= 0.6) { # Filter by greater than or equal to cosine similarity of 0.6
      filtered_geps[[names(parse_geps)[i]]] <- parse_geps[[i]]
      filtered_geps[[names(multiome_geps)[j]]] <- multiome_geps[[j]]
    }
  }
}
geps <- filtered_geps


sim_mtx <- matrix(nrow = length(geps), ncol = length(geps),
                         dimnames = list(names(geps), names(geps)))

# Compute the cosine similarity index for each pair of GEPs
for (i in seq_along(geps)) {
  for (j in seq_along(geps)) {
    v_i <- geps[[i]]
    v_j <- geps[[j]]
    sim_mtx[i, j] <- cosine(as.numeric(v_i), as.numeric(v_j))
  }
}

clustering <- dendsort(hclust(dist(sim_mtx)), type = "average")
clusters <- cutree(clustering, k = 5)

similarity_colors <- colorRamp2(c(0.5, 0.6, 0.7), c("white", "#C0DAD7", "#87B5B1"))

# Platform annotation
stack <- stack(clusters)
stack$Platform <- case_when(
  grepl("Parse", stack$ind) ~ "ParseBio",
  grepl("Multiome", stack$ind) ~ "Multiome"
)
platform_colours <- c("#5773CCFF", "#FFB900FF")
platform_colours <- c("#002A32", "#C4A29E")
names(platform_colours) <- unique(stack$Platform)

annotation_top <- HeatmapAnnotation(`Platform` = stack$Platform[match(colnames(sim_mtx), stack$ind)],
                                             col = list(`Platform` = platform_colours),
                                             show_annotation_name = c(TRUE),
                                             annotation_height = unit(0.5, "cm"),
                                             show_legend = c(TRUE),
                                             gp = gpar(fontsize = 1)
)

pdf(file = paste0(plot_dir, "/geps.pdf"), width = 10, height = 10)
Heatmap(sim_mtx, show_row_names = FALSE, show_column_names = FALSE,
                       height = nrow(sim_mtx)*unit(2, "cm"),
                       width = ncol(sim_mtx)*unit(2, "cm"),
                       col = similarity_colors,
                       name = "Factor-Factor Similarity",
                       column_title = "", row_title = "",
                       cluster_rows = clustering,
                       cluster_columns = clustering,
                       top_annotation = annotation_top,
                       show_row_dend = TRUE,
                       show_column_dend = TRUE)
decorate_heatmap_body("Factor-Factor Similarity", {
  grid.rect(gp = gpar(col = "black", lwd = 2))
})
dev.off()

gep_clusters <- split(names(geps), clusters)
consensus_mps <- lapply(seq(1, length(gep_clusters)), function(i) {
  cluster <- gep_clusters[[i]]
  factors <- geps[cluster]
  genes <- names(factors[[1]])
  consensus <- map(factors, ~as.numeric(as.character(.))) %>% 
    reduce(`+`) %>% 
    `/`(length(factors))
  names(consensus) <- genes
  return (consensus)
})

consensus_mps <- as.data.frame(consensus_mps)
colnames(consensus_mps) <- paste0("Factor", seq(length(colnames(consensus_mps)), 1)) # To match order that consensus factors appear in factor-factor heatmap, naming order is reversed

# write
write.csv(consensus_mps, "/consensus_mps.csv")

# top 100 markers

markers <- sapply(colnames(consensus_mps), function(factor) {
  print(factor)
  curr_mp <- consensus_mps %>% 
    select(factor) %>% 
    arrange(-!!sym(factor))
  return (rownames(curr_mp)[1:100])
})

write.csv(markers, "/consensus_mps_markers.csv")

# ---------------------------------------------------------------------------- #
#                                  Figure g,h                                  #
# ---------------------------------------------------------------------------- #

# function for lmm test

#' perform a nested t-test using LMMs.
#' 
#' @param data_df data frame containing data
#' @param response name of column in data_df that corresponds to the response variable
#' @param condition name of column in data_df that corresponds to condition of interest
#' @param latent_vars vector containing columns in data_df that correspond to variables that describe the nested structure of the data e.g. Patient
#' 
#' @return p-value derived from likelihood ratio test of model incorporating condition vs null model
#' @export
LMM_test <- function(data_df, response, condition, latent_vars) {
    
    require(lme4)
    
    full_formula <- as.formula(
        paste(response, "~", condition, "+", paste("(1 |", latent_vars, ")", collapse = "+"))
    )
    
    null_formula <- as.formula(
        paste(response, "~ 1 +", paste("(1 |", latent_vars, ")", collapse = "+"))
    )
    
    # Fit the models
    model <- lmer(full_formula, data = data_df)
    model_null <- lmer(null_formula, data = data_df)
    
    # Perform the likelihood ratio test
    res <- anova(model_null, model)
    
    # Return the p-value from the Chi-square test
    return(res$`Pr(>Chisq)`[2])
}


#' geom signif with lmm test. See geom_signif from ggsignif for more information
#' 
#' @param data_df data frame containing data
#' @param response name of column in data_df that corresponds to the response variable
#' @param condition name of column in data_df that corresponds to condition of interest
#' @param latent_vars vector containing columns in data_df that correspond to variables that describe the nested structure of the data e.g. Patient
#' @param comparisons list of vectors of length two that define the comparisons to make
#' @param y_position position of brackets 
#' @param tip_length length of bracket tips
#' @param size size of brackets
#' @param step_increase distance between brackets
#' 
#' @return geom signif
#' @export
geom_signif_lmm <- function(data_df, response, condition, latent_vars,
    comparisons, y_position = NULL, tip_length = 0.03, size = 0.5, step_increase = 0) {
    
    require(lme4)
    require(ggsignif)
    require(dplyr)

    p_vals <- sapply(seq_along(comparisons), function(comp) {
        
        curr_comparison <- comparisons[[comp]]
        curr_df <- data_df %>% 
            filter(!!sym(condition) %in% curr_comparison)
        
        return(round(LMM_test(curr_df, response, condition, latent_vars), digits = 5))
    })

    return(
        geom_signif(
            comparisons = comparisons,
            map_signif_level = FALSE,
            annotations = p_vals,
            y_position = y_position,
            tip_length = tip_length,
            size = size,
            step_increase = step_increase
        )
    )
}

metadata <- read.csv("/metacell_metadata.csv", row.names = 1) # path to metacell metadata
metadata$Region <- factor(metadata$Region, levels = c("PT", "TE", "TC"))

df <- metadata
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))
df <- df %>% filter(CellClass_L3 %in% c("Malignant_OPC", "Malignant_NPC1"))
p <- ggplot(df, aes(x = Region ,y = Factor4)) +
  geom_boxplot(aes(fill = Region), alpha = 0.9) +
  scale_fill_manual(values = c(region_cols)) +
  GBM_theme() +
  ylab("Factor activation") +
  ggtitle("Activation of Factor 4 (Neuronal) in OPC/NPC1-like metacells") +
  geom_signif_lmm(
    data_df = df,
    response = "Factor4",
    condition = "Region",
    latent_vars = c("Patient"),
    comparisons = list(c("PT", "TE"), c("TE", "TC"), c("PT", "TC")),
    step_increase = c(0, 0.1, 0.1)
  )
ggsave(filename = "OPC_NPC1_ONLY_Factor_4.pdf", path = plot_dir, height = 10, width = 10)

df <- metadata
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))
df <- df %>% filter(CellClass_L3 %in% c("Malignant_OPC", "Malignant_NPC1"))
p <- ggplot(df, aes(x = Region ,y = Factor5)) +
  geom_boxplot(aes(fill = Region), alpha = 0.9) +
  scale_fill_manual(values = c(region_cols)) +
  GBM_theme() +
  ylab("Factor activation") +
  ggtitle("Activation of Factor 5 (Hypoxia) in OPC/NPC1-like metacells") +
  geom_signif_lmm(
    data_df = df,
    response = "Factor5",
    condition = "Region",
    latent_vars = c("Patient"),
    comparisons = list(c("PT", "TE"), c("TE", "TC"), c("PT", "TC")),
    step_increase = c(0, 0.1, 0.1)
  )
ggsave(filename = "OPC_NPC1_ONLY_Factor_5.pdf", path = plot_dir, height = 10, width = 10)
