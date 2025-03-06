# ---- Code to reproduce Figure 1a ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(dplyr, ComplexHeatmap, RColorBrewer, tidyr, circlize)

# Helper function for drawing Heatmaps
cell_border_fun <- function(j, i, x, y, width, height, fill) {
    grid.rect(
        x = x,
        y = y,
        width = width,
        height = height,
        gp = gpar(col = "black", lwd = 1, fill = NA)
    )
}

# Required inputs
params <- list(
    plot_dir = "output/figures",
    input_table_path = "misc/Table S1.xlsx" # can be downloaded online, see publication
)

GaitiLabUtils::create_dir(params$plot_dir)

# ---- Load data & Data wrangling ---- #
mutations_oi <- c("IDH", "CDKN2A", "TP53", "ATRX", "BRAF")
regions_oi <- c("PT", "TE", "TC")

df <- readxl::read_excel(
    params$input_table_path,
    sheet = "Patient and sample info",
    skip = 2
) %>%
    tidyr::fill(everything()) %>%
    mutate(
        IDH = ifelse(
            stringr::str_detect(`Molecular alterations`, "IDH wildtype"),
            "WT",
            "MUT"
        ),
        CDKN2A = ifelse(
            stringr::str_detect(
                `Molecular alterations`,
                "CDKN2A homozygous deletion"
            ),
            "DEL",
            "WT"
        ),
        TP53 = ifelse(
            stringr::str_detect(
                `Molecular alterations`,
                "TP53 p\\.C277Wfs\\*27  mutation"
            ) |
                stringr::str_detect(`Molecular alterations`, "P53 Mutated"),
            "MUT",
            "WT"
        ),
        ATRX = ifelse(
            stringr::str_detect(`Molecular alterations`, "ATRX retained"),
            "WT",
            "MUT"
        ),
        BRAF = ifelse(
            stringr::str_detect(`Molecular alterations`, "BRAF mutation"),
            "MUT",
            "WT"
        ),
    ) %>%
    group_by_at(c(
        "Patient ID",
        "Sex",
        "IDH",
        "CDKN2A",
        "TP53",
        "ATRX",
        "BRAF",
        "Region annotation"
    )) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = "Region annotation", values_from = "n")

col_annot_df <- df %>%
    ungroup() %>%
    dplyr::select(all_of(c("Sex"))) %>%
    pull()


# ---- Create Figure 1a ---- #

col_annot <- HeatmapAnnotation(
    Sex = col_annot_df,
    which = "col",
    col = list(
        Sex = c("Male" = "white", "Female" = "black")
    ),
    gp = gpar(col = "black")
)

heatmap_df <- df %>%
    ungroup() %>%
    dplyr::select(all_of(mutations_oi))
heatmap_df <- as.matrix(heatmap_df)
rownames(heatmap_df) <- df %>% pull(`Patient ID`)
ht1 <- Heatmap(
    t(heatmap_df),
    name = "Molecular features",
    col = c("WT" = "white", "MUT" = "grey", "DEL" = "black"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    border = TRUE,
    cell_fun = cell_border_fun,
    show_column_names = FALSE,
    top_annotation = col_annot,
)

heatmap_df2 <- df %>%
    ungroup() %>%
    dplyr::select(all_of(regions_oi))
heatmap_df2 <- as.matrix(heatmap_df2)
rownames(heatmap_df2) <- df %>% pull(`Patient ID`)
heatmap_df2 <- ifelse(is.na(heatmap_df2), 0, heatmap_df2)
col_scheme <- brewer.pal(5, "Blues")
col_scheme[1] <- "white"

ht2 <- Heatmap(
    t(heatmap_df2),
    name = "Number of samples from region",
    col = circlize::colorRamp2(c(0, 1, 2, 3, 4), col_scheme),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    border = TRUE,
    cell_fun = cell_border_fun
)

heatmap <- ht1 %v% ht2

pdf(
    file = file.path(params$plot_dir, "Fig1a.pdf"),
    width = 10,
    height = 10
)
draw(heatmap, ht_gap = unit(0.5, "cm"))
dev.off()
