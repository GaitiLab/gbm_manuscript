# Code to reproduce Figure 1a

library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyr)
library(colorRamp2)

cell_border_fun <- function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "black", lwd = 1, fill = NA))
}

df <- read.csv("Table S1.csv")
df <- df %>% 
    group_by_at(c("Patient.ID", "Sex","IDH", "CDKN2A", "TP53", "ATRX", "BRAF","Region.annotation")) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = "Region.annotation", values_from = "n")

col_annot_df <- df %>% ungroup() %>% select(all_of(c("Sex")))

col_annot <- HeatmapAnnotation(
    df = col_annot_df,
    which = "col",
    col = list(
        Sex = c("Male" = "white", "Female" = "black")
    ),
    gp = gpar(col = "black")
)

heatmap_df <- df %>% ungroup() %>%  select(all_of(c("IDH", "CDKN2A", "TP53", "ATRX", "BRAF")))
heatmap_df <- as.matrix(heatmap_df)
rownames(heatmap_df) <- df$Patient.ID
ht1 <- Heatmap(t(heatmap_df),
    name = "Molecular features",
    col = c("WT" = "white", "MUT" = "grey", "DEL" = "black"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    border = TRUE,
    cell_fun = cell_border_fun,
    show_column_names = FALSE,
    top_annotation = col_annot)

heatmap_df2 <- df %>% ungroup() %>% select(all_of(c("PT", "TE", "TC")))
heatmap_df2 <- as.matrix(heatmap_df2)
rownames(heatmap_df2) <- df$Patient.ID
heatmap_df2 <- ifelse(is.na(heatmap_df2), 0, heatmap_df2)
col_scheme <- brewer.pal(5,"Blues")
col_scheme[1] <- "white"

ht2 <- Heatmap(t(heatmap_df2),
    name = "Number of samples from region",
    col = colorRamp2(c(0, 1, 2, 3, 4), col_scheme),
    cluster_rows = FALSE, cluster_columns = FALSE,
    border = TRUE,
    cell_fun = cell_border_fun)

heatmap <- ht1 %v% ht2 

pdf(file = "/Users/bensonwu/Downloads/cohort_metadata_heatmap_top_up.pdf", width = 10, height = 10)
draw(heatmap, ht_gap = unit(0.5, "cm"))
dev.off()
