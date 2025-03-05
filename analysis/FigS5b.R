# ---- Code to reproduce Figure S5b ---- #

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
    readxl,
    ggplot2,
    ggrastr,
    ggtext
)

# Required inputs
params <- list(
    plot_dir = "output/submission/figures",
    # scRNAseq Seurat object
    seurat_obj_path = "",
    interactions_path = "misc/SuppTables/Table S4.xlsx",
    condition_varname = "Region",
    pval_type = "pval_adj",
    condition_oi = "PT",
    pair_y = "Glutamatergic__Invasive-high OPC/NPC1",
    pair_x = "Glutamatergic__Progenitor_like",
    neuron_subtype = "Glutamatergic",
    alpha = 0.05,
    # Ref db can be downloaded from https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline/blob/main/assets/interactions_db/ref_db.rds
    ref_db = "misc/internal/ref_db.rds"
)

GaitiLabUtils::create_dir(params$plot_dir)

label_color <- setNames(
    c(
        "background",
        c("#2B7095", "#C05E00")
    ),
    c("black,", "NRXN-NLGN interactions", "Glutamate production/uptake")
)

# ---- Load data & data wrangling ---- #

# Setup for visualization
ref_db <- readRDS(params$ref_db)
genes_in_db <- ref_db %>%
    dplyr::select(ligand_complex, receptor_complex) %>%
    unlist() %>%
    str_split(., ":") %>%
    unname(.) %>%
    unlist() %>%
    unique()

# Load metadata
seurat_obj <- readRDS(params$seurat_obj_path)

meta <- seurat_obj@meta.data %>%
    dplyr::select(Sample, !!sym(params$condition_varname)) %>%
    distinct() %>%
    remove_rownames()

# Format expression
groupbyvar <- c("Sample", "CCI_CellClass_L2_2")

# If Seurat v4
# avg_expr <- AverageExpression(
#     seurat_obj,
#     assays = "RNA",
#     group.by = groupbyvar,
#     slot = "counts"
# )[["RNA"]]

# If Seurat v5
avg_exp <- AverageExpression(
    seurat_obj,
    assays = "RNA",
    layer = "counts",
    group.by = groupbyvar
)[["RNA"]]


column_names <- colnames(avg_exp)

# Remove prefix 'g-' that is added by Seurat
new_colnames <- str_sub_all(column_names, 2) %>%
    unlist() %>%
    # col1 = cell type, col2 = Sample
    str_split(., "_", simplify = TRUE) %>%
    apply(., 1, rev) %>%
    t() %>%
    data.frame() %>%
    unite(colnames, X1, X2, sep = ":") %>%
    pull()
colnames(avg_exp) <- new_colnames

# Format gene expression matrix
avg_exp_long <- data.frame(data.matrix(avg_exp), check.names = FALSE) %>%
    rownames_to_column("gene") %>%
    # Convert to long-format
    pivot_longer(
        names_sep = ":",
        names_to = c("label", "Sample"),
        cols = all_of(colnames(avg_exp)),
        values_to = "avg_expression"
    ) %>%
    # Ensure that formatting of `Sample` matches with `Sample` in metadata for merging
    mutate(Sample = str_replace_all(Sample, "-", "_")) %>%
    # Adding regional info
    left_join(meta) %>%
    group_by(!!sym(params$condition_varname), label, gene) %>%
    summarise(mean = mean(avg_expression)) %>%
    mutate(
        Region = factor(
            !!sym(params$condition_varname),
            levels = names(GBMutils::load_color_palette("Region"))
        )
    ) %>%
    ungroup()

exp_df <- avg_exp_long %>%
    # Taking log-mean
    mutate(mean = log2(mean + 1)) %>%
    filter(
        label %in%
            c(
                params$neuron_subtype,
                "Progenitor-like",
                "Invasive-high OPC/NPC1"
            ),
        # Only keep genes that are in database
        gene %in% genes_in_db
    ) %>%
    pivot_wider(names_from = label, values_from = mean, values_fill = NA) %>%
    mutate(diff = !!sym("Invasive-high OPC/NPC1") - !!sym("Progenitor-like"))


interactions_of_interest_df <- readxl::read_excel(
    params$interactions_path,
    sheet = "All_predictions",
    skip = 1
) %>%
    separate(
        complex_interaction,
        into = c("ligand_complex", "receptor_complex"),
        sep = "__",
        remove = FALSE
    ) %>%
    filter(
        # Only keep interactions detected for Neuron - Progenitor-like OR Neuron - Invasive-high OPC/NPC1
        source_target %in% c(params$pair_x, params$pair_y),
        # Only keep interactions with a (adjusted) Fisher combined p-value < params$alpha
        !!sym(params$pval_type) < params$alpha,
        # Keep interactions that are detected in PT region
        !!sym(params$condition_varname) == params$condition_oi,
        # Only keep certain interactions
        # NRXN-NLGN / NRXN-NLGN interactions
        ((str_detect(ligand_complex, "NRXN") &
            str_detect(receptor_complex, "NLGN")) |
            (str_detect(ligand_complex, "NLGN") &
                str_detect(receptor_complex, "NRXN")) |
            # Glutamate
            (str_detect(complex_interaction, "GLS|GRIA|GRIK")))
    ) %>%
    # Remove duplicate interactions
    distinct(complex_interaction, .keep_all = TRUE) %>%
    dplyr::select(ligand_complex, receptor_complex)

# Extract genes from interactions of interest
genes_oi <- interactions_of_interest_df %>%
    dplyr::select(ligand_complex, receptor_complex) %>%
    unlist() %>%
    unname() %>%
    str_split(., "\\:") %>%
    unlist() %>%
    unique()

df_region <- exp_df %>%
    # Only keep the interactions from region PT
    filter(
        !!sym(params$condition_varname) == params$condition_oi,
    ) %>%
    mutate(
        # Label genes that are of interest NRXN, NLGN and glumate-related genes
        is_gene_oi = case_when(
            gene %in% genes_oi & str_detect(gene, "NRXN|NLGN") ~
                "NRXN-NLGN interactions",
            gene %in% genes_oi & (str_detect(gene, "GLS|GRIA|GRIK")) ~
                "Glutamate production/uptake",
            .default = "background"
        )
    )


# ---- Create Figure S5b ---- #
min_limit <- plyr::round_any(min(df_region$diff), 0.5, f = floor)
max_limit <- plyr::round_any(max(df_region$diff), 0.5, f = ceiling)
global_max_lim <- max(abs(c(min_limit, max_limit)))

p <- ggplot(data = df_region, aes(x = diff, y = !!sym(params$neuron_subtype))) +
    geom_point_rast(size = 2, alpha = .1, color = "black", stroke = NA) +
    labs(
        y = paste0("log<sub>2</sub>(mean exp + 1) in ", params$neuron_subtype),
        x = paste0(
            paste0(
                "log<sub>2</sub>(mean exp + 1) in ",
                "Invasive-high OPC/NPC1"
            ),
            " - ",
            paste0("log<sub>2</sub>(mean exp + 1) in ", "Progenitor-like")
        ),
        subtitle = glue(
            "No. of genes {prettyNum(nrow(df_region), big.mark=',')}"
        )
    ) +
    GBMutils::GBM_theme() +
    theme(
        aspect.ratio = 1,
        plot.subtitle = element_text(size = ggplot2::rel(1.3)),
        axis.title.x = ggtext::element_markdown(
            size = ggplot2::rel(1.2),
            face = "plain"
        ),
        axis.title.y = ggtext::element_markdown(
            size = ggplot2::rel(1.2),
            face = "plain"
        ),
        axis.text = ggplot2::element_text(size = ggplot2::rel(1.2)),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))
    ) +
    # Labeling
    ggnewscale::new_scale("colour") +
    scale_color_manual(
        values = label_color,
        guide = guide_legend(
            override.aes = list(size = 7, alpha = 1),
            title = "Involved in",
            nrow = 2
        )
    ) +
    geom_point_rast(
        data = df_region %>% filter(is_gene_oi != "background"),
        aes(color = is_gene_oi),
        alpha = 1,
        show.legend = TRUE,
        size = 2,
        stroke = NA
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggrepel::geom_label_repel(
        data = df_region %>% filter(is_gene_oi != "background"),
        aes(label = gene, color = is_gene_oi),
        alpha = 1,
        show.legend = FALSE,
        min.segment.length = 0.1,
        size = 6
    ) +
    coord_fixed(xlim = c(-global_max_lim, global_max_lim))

ggsave(
    plot = p,
    filename = "FigS5b.pdf",
    width = 12,
    height = 12,
    path = params$plot_dir
)
