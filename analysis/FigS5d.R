# ---- Code to reproduce Figure S5d ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #
pacman::p_load(
    GaitiLabUtils,
    GBMutils,
    glue,
    data.table,
    tidyverse,
    stringr,
    sf,
    rjson,
    Seurat
)
logr <- GaitiLabUtils::init_logging()

# Required inputs
params <- list(
    # path pointing to the parent directory with sample directories (masks + seurat obj)
    input_dir = "data/visiumhd/processed",
    plot_dir = "output/submission/figures"
)

sample_ids <- c(
    "6425_A",
    "6425_B"
)
palette_L1 <- GBMutils::load_color_palette("Spatial_CellClass_L1")

sample_6425b_genes <- list(
    OPC = c(
        "XYLT1",
        "LUZP2",
        "MEGF11",
        "PCDH15",
        "DSCAM",
        "SMOC1",
        "TNR",
        "LHFPL3"
    ),
    Oligodendrocyte = c(
        "PLP1",
        "TF",
        "TMEM144",
        "MOBP",
        "MBP",
        "DOCK5",
        "EDIL3",
        "RNF220",
        "ST18",
        "CTNNA3"
    ),
    Astrocyte = c(
        "SPARCL1",
        "NEBL",
        "TPD52L1",
        "GPC5",
        "TRPM3",
        "HPSE2",
        "SLC1A2",
        "COL5A3",
        "ADGRV1",
        "RYR3",
        "ATP1A2",
        "RANBP3L"
    ),
    Myeloid = c(
        "ADAM28",
        "SLC11A1",
        "APBB1IP",
        "ATP8B4",
        "MSR1",
        "DOCK8",
        "PLXDC2",
        "ST6GAL1",
        "PALD1"
    ),
    BV = c(
        "IGFBP7",
        "VWF",
        "ARHGAP29",
        "ITGA1",
        "EBF1",
        "MECOM",
        "COL4A1",
        "ADGRL4",
        "RBMS3",
        "COBLL1",
        "ABCB1",
        "FLT1",
        "ADGRF5",
        "PTPRB"
    ),
    Neuron = c(
        "SYT1",
        "GALNTL6",
        "TENM2",
        "CNTN5",
        "KCNC2",
        "RYR2",
        "CDH18",
        "CNTNAP2",
        "NRG3",
        "RALYL",
        "RBFOX1",
        "DLGAP2",
        "CCSER1",
        "RIMS2"
    )
)

sample_6425a_genes <- list(
    BV = c(
        "VWF",
        "ABCB1",
        "ADGRL4",
        "FLT1",
        "ARHGAP29",
        "EBF1",
        "RBMS3",
        "MECOM",
        "ITGA1",
        "ADGRF5",
        "PTPRB",
        "COBLL1",
        "IGFBP7",
        "COL4A1"
    ),
    Neuron = c(
        "SYT1",
        "RYR2",
        "KCNC2",
        "RALYL",
        "DLGAP2",
        "RIMS2",
        "CDH18",
        "RBFOX1"
    )
)

gene_lists <- list(`6425_B` = sample_6425b_genes, `6425_A` = sample_6425a_genes)


current_sample_params <- list()
for (current_sample_id in sample_ids) {
    log_info(glue("Current sample: {current_sample_id}..."))
    current_sample_params$seurat_obj_path <- file.path(
        params$input_dir,
        current_sample_id,
        paste0(current_sample_id, ".rds")
    )

    current_gene_list <- gene_lists[[current_sample_id]]

    log_info("Load Seurat object...")
    seurat_obj <- readRDS(current_sample_params$seurat_obj_path)
    DefaultAssay(seurat_obj) <- "RNA"
    Idents(seurat_obj) <- "Spatial_CellClass_L1"

    # Use standardized data
    seurat_obj <- ScaleData(
        seurat_obj,
        features = unname(unlist(current_gene_list)),
        assay = "RNA"
    )
    avg_exp <- AverageExpression(
        seurat_obj,
        features = unname(unlist(current_gene_list)),
        layer = "scale.data",
        assay = "RNA"
    )[["RNA"]] %>%
        data.frame() %>%
        rownames_to_column("gene")

    log_info(
        "Annotate genes in VisiumHD using the cell type markers from scRNAseq..."
    )

    gene_list_as_df <- do.call(
        rbind,
        lapply(seq_len(length(names(current_gene_list))), function(ix) {
            data.frame(
                label = names(current_gene_list)[ix],
                gene = current_gene_list[[names(current_gene_list)[ix]]]
            )
        })
    )

    df <- avg_exp %>%
        left_join(
            gene_list_as_df
        ) %>%
        filter(!is.na(label)) %>%
        arrange(label)

    # Order to match publication
    rownames(df) <- df$gene
    df <- df[gene_list_as_df$gene, ]
    df <- df %>% remove_rownames()

    top_annot_df <- df %>% select(gene, label)

    mat <- df %>%
        select(
            -c(
                label
            )
        ) %>%
        t()

    colnames(mat) <- mat[1, ]
    mat <- mat[2:nrow(mat), ]
    mat <- mat %>% data.matrix()
    old_rownames <- rownames(mat)

    mat <- apply(mat, 2, as.numeric, simplify = TRUE)
    rownames(mat) <- old_rownames

    # Setup legend for ComplexHeatmap
    legend_min <- plyr::round_any(min(mat), 0.5, round)
    legend_mid <- 0
    legend_max <- plyr::round_any(max(mat), 0.5, round)
    limit <- max(abs(c(legend_min, legend_max)))

    # Setup colors
    col_fun <- circlize::colorRamp2(
        # Values
        c(-limit, legend_mid, limit),
        # Colors
        GBMutils::load_color_palette("Divergent_1")
    )

    lgd <- ComplexHeatmap::Legend(
        col_fun = col_fun,
        title = "Mean scaled\nexpression",
        direction = "horizontal",
    )
    lgd_list <- list(lgd)

    # Setup gene annotation categories (labels from annotation file)
    top_annot <- ComplexHeatmap::HeatmapAnnotation(
        Label = top_annot_df %>% pull(label),
        col = list(
            Label = palette_L1[!is.na(names(palette_L1))]
        )
    )

    # Determine order of columns (cluster genes within each cell type)
    dend <- ComplexHeatmap::cluster_within_group(
        mat,
        as.character(top_annot_df$label)
    )
    # Use the same order to order the rows
    mat <- mat[
        c(
            unique(as.character(top_annot_df$label)[order.dendrogram(dend)]),
            "Malignant",
            "Undetermined"
        ),
    ]

    log_info("Create heatmap...")
    hm <- GaitiLabUtils::create_hm(
        mat,
        cell_width = 5,
        cell_height = 5,
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        column_split = top_annot_df %>% pull(label) %>% unique() %>% length(),
        # Font sizes
        row_names_gp = grid::gpar(fontsize = 12),
        column_names_gp = grid::gpar(fontsize = 8),
        column_title_gp = grid::gpar(fontsize = 12),
        # Dendrograms
        show_row_dend = FALSE,
        show_column_dend = FALSE,

        # Clustering
        cluster_rows = FALSE,
        cluster_columns = dend,
        show_heatmap_legend = FALSE,
        row_title = "Cluster",
        name = "Mean scaled expression",

        # Annotations
        top_annotation = top_annot
    )
    # Save heatmap
    GaitiLabUtils::save_hm(
        hm_obj = hm,
        output_file = file.path(
            params$plot_dir,
            glue(
                "FigS5d-{current_sample_id}_expression_heatmaps.pdf"
            )
        ),
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom",
        merge_legend = TRUE,
        heatmap_legend_list = lgd_list
    )
}
