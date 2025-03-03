# ---- Code to reproduce Figure S3 ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

pacman::p_load(
    data.table,
    patchwork,
    ggplot2,
    dplyr,
    tidyr,
    ComplexHeatmap,
    grid,
    colorRamp2,
    fgsea,
    GBMutils,
    scales,
    GaitiLabUtils,
    Seurat
)

region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

# Required inputs
params <- list(
    output_dir = "output/submission",
    plot_dir = "output/submission/figures",
    seurat_obj_path = "",
    protein_coding_genes_list_path = "misc/ensembl_protein_coding_genes.csv",
    multiome_factors_spectra_path = "misc/multiome_top_up_rna.spectra.k_7.dt_0_5.consensus.txt",
    parsebio_factors_spectra_path = "misc/parsebio_top_up_rna.spectra.k_7.dt_0_5.consensus.txt",
    usage_mtx_path = "misc/usage_mtx.csv",
    genelists_path = "misc/SuppTables/Table S2.xlsx",
)

# Creating directories needed for outputs
GaitiLabUtils::create_dir(params$output_dir)
GaitiLabUtils::create_dir(params$plot_dir)

# export data for cNMF
seurat_obj <- readRDS(params$seurat_obj_path)
seurat_obj <- subset(seurat_obj, subset = is_malignant_confident == TRUE)

split_list <- SplitObject(seurat_obj, split.by = "Patient")
var_features <- c()
for (x in split_list) {
    DefaultAssay(x) <- "RNA"
    features <- FindVariableFeatures(x, verbose = FALSE, nfeatures = 2000)
    features <- features@assays$RNA@var.features
    var_features <- c(var_features, features)
}
var_features <- data.frame(features = var_features) %>%
    group_by(features) %>%
    summarise(count = n()) %>%
    arrange(-count)
protein_coding_genes <- read.csv(params$protein_coding_genes_list_path)
var_features <- var_features[
    var_features$features %in% protein_coding_genes$hgnc_symbol,
]

var_features <- var_features$features[1:2000]

parse <- subset(seurat_obj, subset = Platform == "ParseBio")
rna_parse_counts <- parse@assays$RNA@counts[var_features, ]
rna_parse_counts <- as.data.table(
    t(as.matrix(rna_parse_counts)),
    keep.rownames = "ID"
)
fwrite(
    rna_parse_counts,
    file.path(params$output_dir, "malignant_RNA_counts_parse.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

rna_parse_data <- parse@assays$RNA@data[var_features, ]
rna_parse_data <- as.data.table(
    t(as.matrix(rna_parse_data)),
    keep.rownames = "ID"
)
fwrite(
    rna_parse_data,
    file.path(params$output_dir, "malignant_RNA_data_parse.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

multiome <- subset(seurat_obj, subset = Platform == "Multiome")
rna_multiome_counts <- multiome@assays$RNA@counts[var_features, ]
rna_multiome_counts <- as.data.table(
    t(as.matrix(rna_multiome_counts)),
    keep.rownames = "ID"
)
fwrite(
    rna_multiome_counts,
    file.path(params$output_dir, "malignant_RNA_counts_multiome.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

rna_multiome_data <- multiome@assays$RNA@data[var_features, ]
rna_multiome_data <- as.data.table(
    t(as.matrix(rna_multiome_data)),
    keep.rownames = "ID"
)
fwrite(
    rna_multiome_data,
    file.path(params$output_dir, "malignant_RNA_data_multiome.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ---------------------------------------------------------------------------- #
#                                    Fig S3a                                   #
# ---------------------------------------------------------------------------- #
# consensus by platform

# Multiome
multiome_factors <- read.table(params$multiome_factors_spectra_path)

multiome_geps <- lapply(
    seq(1, length(rownames(multiome_factors))),
    function(gep_num) {
        gep <- as.vector(multiome_factors[gep_num, ])
    }
)

names(multiome_geps) <- paste0("MultiomeFactor", seq(1, length(multiome_geps)))

# Parsebio
parse_factors <- read.table(params$parsebio_factors_spectra_path)

parse_geps <- lapply(
    seq(1, length(rownames(parse_factors))),
    function(gep_num) {
        gep <- as.vector(parse_factors[gep_num, ])
    }
)

names(parse_geps) <- paste0("ParseFactor", seq(1, length(parse_geps)))

filtered_geps <- list()
for (i in seq_along(parse_geps)) {
    for (j in seq_along(multiome_geps)) {
        parse_gep <- parse_geps[[i]]
        multiome_gep <- multiome_geps[[j]]
        print(cosine(as.numeric(parse_gep), as.numeric(multiome_gep)))
        if (cosine(as.numeric(parse_gep), as.numeric(multiome_gep)) >= 0.5) {
            filtered_geps[[names(parse_geps)[i]]] <- parse_geps[[i]]
            filtered_geps[[names(multiome_geps)[j]]] <- multiome_geps[[j]]
        }
    }
}
geps <- filtered_geps

sim_mtx <- matrix(
    nrow = length(geps),
    ncol = length(geps),
    dimnames = list(names(geps), names(geps))
)

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

similarity_colors <- colorRamp2(
    c(0.5, 0.65, 0.8),
    c("white", "#C0DAD7", "#87B5B1")
)

# Platform annotation
stack <- stack(clusters)
stack$Platform <- case_when(
    grepl("Parse", stack$ind) ~ "ParseBio",
    grepl("Multiome", stack$ind) ~ "Multiome"
)
platform_colours <- c("#5773CCFF", "#FFB900FF")
platform_colours <- c("#002A32", "#C4A29E")
names(platform_colours) <- unique(stack$Platform)

consensus_factor_cols <- brewer.pal(5, "Set3")
names(consensus_factor_cols) <- unique(stack$values)

annotation_top <- HeatmapAnnotation(
    `Platform` = stack$Platform[match(colnames(sim_mtx), stack$ind)],
    `Consensus Factor` = stack$values[match(colnames(sim_mtx), stack$ind)],
    col = list(
        `Platform` = platform_colours,
        `Consensus Factor` = consensus_factor_cols
    ),
    show_annotation_name = c(TRUE, FALSE),
    annotation_height = unit(0.5, "cm"),
    show_legend = c(TRUE, FALSE),
    gp = gpar(fontsize = 1)
)

annotation_left <- rowAnnotation(
    `Consensus Factor` = stack$values[match(colnames(sim_mtx), stack$ind)],
    col = list(`Consensus Factor` = consensus_factor_cols),
    show_annotation_name = c(FALSE),
    annotation_height = unit(0.5, "cm"),
    show_legend = c(FALSE),
    gp = gpar(fontsize = 1)
)


pdf(
    file = file.path(params$plot_dir, "FigS3a_geps.pdf"),
    width = 10,
    height = 10
)
Heatmap(
    sim_mtx,
    show_row_names = FALSE,
    show_column_names = FALSE,
    height = nrow(sim_mtx) * unit(1, "cm"),
    width = ncol(sim_mtx) * unit(1, "cm"),
    col = similarity_colors,
    name = "Factor-Factor Similarity",
    column_title = "",
    row_title = "",
    cluster_rows = clustering,
    cluster_columns = clustering,
    top_annotation = annotation_top,
    left_annotation = annotation_left,
    show_row_dend = TRUE,
    show_column_dend = TRUE
)
decorate_heatmap_body("Factor-Factor Similarity", {
    grid.rect(gp = gpar(col = "black", lwd = 2))
})
dev.off()

gep_clusters <- split(names(geps), clusters)
consensus_mps <- lapply(seq(1, length(gep_clusters)), function(i) {
    cluster <- gep_clusters[[i]]
    factors <- geps[cluster]
    genes <- names(factors[[1]])
    consensus <- map(factors, ~ as.numeric(as.character(.))) %>%
        reduce(`+`) %>%
        `/`(length(factors))
    names(consensus) <- genes
    return(consensus)
})

consensus_mps <- as.data.frame(consensus_mps)
colnames(consensus_mps) <- paste0(
    "Factor",
    seq(length(colnames(consensus_mps)), 1)
)
write.csv(
    consensus_mps,
    file.path(params$output_dir, "FigS3a_consensus_mps.csv")
)

# top 100 markers

markers <- sapply(colnames(consensus_mps), function(factor) {
    print(factor)
    curr_mp <- consensus_mps %>%
        select(factor) %>%
        arrange(-!!sym(factor))
    return(rownames(curr_mp)[1:100])
})

write.csv(
    markers,
    file.path(params$output_dir, "FigS3a_consensus_mps_markers.csv")
)

# ---------------------------------------------------------------------------- #
#                                    Fig S3b                                   #
# ---------------------------------------------------------------------------- #
markers <- as.vector(as.data.frame(markers))

ora_dir <- file.path(params$plot_dir, "ora")
GaitiLabUtils::create_dir(ora_dir)
msigdb_cat_list <- list(
    c("H"),
    c("C2", "CGP"),
    c("C2", "CP"),
    c("C5", "GO:BP"),
    c("C5", "GO:CC"),
    c("C5", "GO:MF")
)

# TODO @Bensonwu02 'msigdb' hasn't been specified before.
msigdb_data <- fread(msigdb)

universe_genes <- colnames(multiome_factors)

for (curr_factor in names(markers)) {
    log_info(paste0("ORA for ", curr_factor))

    all_sig_results <- c()
    for (j in seq(1, length(msigdb_cat_list))) {
        curr_cat <- msigdb_cat_list[[j]]

        if (length(curr_cat) == 2) {
            curr_gene_sets <- msigdb_data[
                gs_cat == curr_cat[1] & gs_subcat == curr_cat[2]
            ]
        } else {
            curr_gene_sets <- msigdb_data[gs_cat == curr_cat[1]]
        }

        curr_msigdbr_list <- split(
            x = curr_gene_sets$gene_symbol,
            f = curr_gene_sets$gs_name
        )

        curr_ora_results <- fora(
            curr_msigdbr_list,
            markers[[curr_factor]],
            universe_genes
        )
        sig_ora_results <- curr_ora_results[padj < 0.5]

        if (nrow(sig_ora_results) > 10) {
            setorder(sig_ora_results, padj)
            sig_ora_results <- sig_ora_results[1:10]
        }

        sig_ora_results$msigdb_cat <- paste0(curr_cat, collapse = "_")
        all_sig_results <- rbindlist(list(all_sig_results, sig_ora_results))
    }

    if (nrow(all_sig_results) == 0) {
        sig_result_plots <- NULL
    } else {
        if (length(unique(all_sig_results$msigdb_cat)) == 1) {
            sig_result_plots <- ggplot(
                all_sig_results,
                aes(
                    y = reorder(pathway, -padj),
                    x = -log10(padj),
                    color = -log10(padj)
                )
            ) +
                geom_point(aes(size = size)) +
                scale_color_gradient(low = "orange", high = "red") +
                labs(title = "Pathway Enrichment", y = "Pathway") +
                theme_classic()

            ggsave(
                filename = paste0(curr_factor, ".pdf"),
                plot = sig_result_plots,
                path = ora_dir,
                width = 8,
                height = 15
            )
        } else {
            plots <- list()
            for (cat in unique(all_sig_results$msigdb_cat)) {
                df <- all_sig_results %>%
                    filter(msigdb_cat == cat)

                plots[[cat]] <- ggplot(
                    df,
                    aes(
                        y = reorder(pathway, -padj),
                        x = -log10(padj),
                        color = -log10(padj)
                    )
                ) +
                    geom_point(aes(size = size)) +
                    scale_color_gradient(low = "orange", high = "red") +
                    labs(title = "Pathway Enrichment", y = "Pathway") +
                    ggtitle(cat) +
                    theme_classic()
            }
            plots <- wrap_plots(plots, ncol = 3, nrow = 3)
            ggsave(
                plots,
                filename = paste0("fora_indiv_", curr_factor, ".pdf"),
                path = ora_dir,
                height = 15,
                width = 25
            )

            all_sig_results <- all_sig_results %>% arrange(padj)
        }

        fwrite(
            all_sig_results,
            file = paste0(ora_dir, "/", curr_factor, ".csv")
        )
    }
}

# TODO @Bensonwu02 changed this to explicit path as you specified before
ora_res <- list.files(ora_dir, full.names = TRUE)
# ora_res <- list.files("/ora", full.names = TRUE)

ora_res <- ora_res[grepl("^Factor.*\\.csv$", basename(ora_res))]

plots <- list()
for (path in ora_res) {
    res <- read.csv(path)
    res$pct_overlap <- (res$overlap / res$size) * 100

    factor <- gsub(".csv", "", basename(path))
    title <- case_when(
        factor == "Factor5" ~ "Factor 5 (Hypoxia)",
        factor == "Factor4" ~ "Factor 4 (Neural Crest)",
        factor == "Factor3" ~ "Factor 3 (Neuronal)",
        factor == "Factor2" ~ "Factor 2 (Cilia)",
        factor == "Factor1" ~ "Factor 1 (Cell Cycle)",
    )

    if (factor == "Factor5") {
        res$pathway <- case_when(
            res$pathway == "HALLMARK_HYPOXIA" ~ "Hallmark \n Hypoxia",
            res$pathway == "MENSE_HYPOXIA_UP" ~ "Mense et al. \n Hypoxia Up",
            res$pathway == "ELVIDGE_HYPOXIA_UP" ~
                "Elvidge et al. \n Hypoxia Up",
            res$pathway == "HALLMARK_MTORC1_SIGNALING" ~
                "Hallmark \n MTORC1 Signaling"
        )
        res <- res[
            res$pathway %in%
                c(
                    "Hallmark \n Hypoxia",
                    "Mense et al. \n Hypoxia Up",
                    "Elvidge et al. \n Hypoxia Up",
                    "Hallmark \n MTORC1 Signaling"
                ),
        ]
    } else if (factor == "Factor4") {
        res$pathway <- case_when(
            res$pathway == "LEE_NEURAL_CREST_STEM_CELL_UP" ~
                "Lee et al. \n Neural Crest Stem Cell Up",
            res$pathway == "VERHAAK_GLIOBLASTOMA_CLASSICAL" ~
                "Verhaak et al. \n Glioblastoma Classical"
        )
        res <- res[
            res$pathway %in%
                c(
                    "Lee et al. \n Neural Crest Stem Cell Up",
                    "Verhaak et al. \n Glioblastoma Classical"
                ),
        ]
    } else if (factor == "Factor3") {
        res$pathway <- case_when(
            res$pathway == "GOCC_SYNAPSE" ~ "GO:CC \n Synapse",
            res$pathway == "GOCC_GLUTAMATERGIC_SYNAPSE" ~
                "GO:CC \n Glutamatergic Synapse",
            res$pathway == "GOCC_SYNAPTIC_MEMBRANE" ~
                "GO:CC \n Synaptic Membrane",
            res$pathway == "GOBP_SYNAPSE_ORGANIZATION" ~
                "GO:BP \n Synapse Organization"
        )
        res <- res[
            res$pathway %in%
                c(
                    "GO:CC \n Synapse",
                    "GO:CC \n Glutamatergic Synapse",
                    "GO:CC \n Synaptic Membrane",
                    "GO:BP \n Synapse Organization"
                ),
        ]
    } else if (factor == "Factor2") {
        res$pathway <- case_when(
            res$pathway == "LIM_MAMMARY_LUMINAL_MATURE_DN" ~
                "Lim et al. \n Mammary Luminal Mature Down",
            res$pathway == "GOBP_MICROTUBULE_BASED_MOVEMENT" ~
                "GO:BP \n Microtubule Based Movement",
            res$pathway == "GOBP_CILIUM_ORGANIZATION" ~
                "GO:BP \n Cilium Organization",
            res$pathway == "GOBP_CILIUM_MOVEMENT" ~ "GO:BP \n Cilium Movement"
        )
        res <- res[
            res$pathway %in%
                c(
                    "Lim et al. \n Mammary Luminal Mature Down",
                    "GO:BP \n Microtubule Based Movement",
                    "GO:BP \n Cilium Organization",
                    "GO:BP \n Cilium Movement"
                ),
        ]
    } else {
        res$pathway <- case_when(
            res$pathway == "FISCHER_DREAM_TARGETS" ~
                "Fischer et al. \n Dream targets",
            res$pathway == "MARSON_BOUND_BY_E2F4_UNSTIMULATED" ~
                "Marson et al. \n Bound by E2F4 Unstimulated",
            res$pathway == "GOBP_CELL_CYCLE" ~ "GO:BP \n Cell Cycle",
            res$pathway == "HALLMARK_G2M_CHECKPOINT" ~
                "Hallmark \n G2M Checkpoint",
        )
        res <- res[
            res$pathway %in%
                c(
                    "Fischer et al. \n Dream targets",
                    "Marson et al. \n Bound by E2F4 Unstimulated",
                    "GO:BP \n Cell Cycle",
                    "Hallmark \n G2M Checkpoint"
                ),
        ]
    }
    p <- ggplot(
        res,
        aes(x = -log10(padj), y = reorder(pathway, -log10(padj)))
    ) +
        geom_bar(stat = "identity", color = "grey") +
        GBM_theme() +
        labs(title = title, x = "-log10(FDR)", y = "")

    plots[[factor]] <- p
}

plots <- wrap_plots(plots, ncol = 1, nrow = 5, guides = "collect")
ggsave(
    plot = plots,
    filename = "FigS3b_factors_ora.pdf",
    path = params$plot_dir,
    width = 10,
    height = 20
)

# ---------------------------------------------------------------------------- #
#                                  Figure S3c                                  #
# ---------------------------------------------------------------------------- #
usage <- read.csv(params$usage_mtx_path, row.names = 1)

rna <- rbind(rna_multiome_data, rna_parse_data)

colnames(usage) <- paste0("Factor", seq(length(colnames(usage)), 1)) # Ordering is reversed to match order in which factors appear in factor-factor heatmap
rownames(usage) <- rownames(rna)

seurat_obj[["nmf"]] <- CreateAssayObject(t(usage))
seurat_obj <- AddMetaData(seurat_obj, metadata = as.data.frame(usage))

# Extract Neftel stignatures
gene_lists <- readxl::read_excel(params$genelists_path, skip = 1)
gene_list <- lapply(
    gene_lists %>% dplyr::select(starts_with("Neftel_")) %>% as.list(),
    function(gene_list) {
        gene_list[!is.na(gene_list)]
    }
)

seurat_obj <- AddModuleScore(seurat_obj, features = gene_list, name = "Neftel")

metadata <- seurat_obj[[]]
metadata$Region <- factor(metadata$Region, levels = c("PT", "TE", "TC"))

df <- metadata
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))
df <- df %>%
    select(c(
        "Region",
        "Patient",
        paste0("Factor", seq(1, length(colnames(usage))))
    )) %>%
    pivot_longer(
        -c("Region", "Patient"),
        names_to = "Factor",
        values_to = "Activation"
    )

plots <- list()
for (factor in unique(df$Factor)) {
    curr_df <- df %>% filter(Factor == factor)
    plots[[factor]] <- ggplot(
        curr_df,
        aes(x = Region, y = Activation, fill = Region)
    ) +
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
ggsave(
    plots,
    filename = "FigS3c_Factors_boxplot_lmm.pdf",
    path = params$plot_dir,
    width = 12,
    height = 5
)

df <- metadata %>%
    mutate(proliferating = ifelse(Factor1 > 0, "Cycling", "NonCycling")) %>%
    group_by_at(c("Region", "proliferating")) %>%
    summarise(num_cycling = n())
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))
p <- ggplot(df, aes(x = Region, y = num_cycling, fill = proliferating)) +
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    scale_fill_manual(values = c("#6E7E85", "#1C0F13"))
ggsave(filename = "FigS3c_factor1_cycling_cells.pdf", path = params$plot_dir)

# ---------------------------------------------------------------------------- #
#                                 Figure S3d-e                                 #
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
corr_df <- list()
states <- c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2")

for (factor in paste0("Factor", seq(1, length(colnames(usage))))) {
    cors <- c()

    for (state in states) {
        cor_test <- cor.test(
            plot_df[[state]],
            plot_df[[factor]],
            method = "pearson"
        )
        print(paste0(
            state,
            " corr:",
            cor_test$estimate,
            ", pval:",
            cor_test$p.value
        ))
        cors <- c(cors, cor_test$estimate)

        p <- ggplot(plot_df, aes(x = get(factor), y = get(state))) +
            stat_density_2d(
                aes(fill = Region, alpha = ..level..),
                geom = "polygon",
                contour_var = "ndensity",
                bins = 5
            ) +
            annotate(
                "text",
                label = paste0(
                    state,
                    " corr:",
                    round(cor_test$estimate, 2),
                    ", pval:",
                    format(cor_test$p.value, scientific = TRUE, digits = 3)
                ),
                x = 0,
                y = 3
            ) +
            scale_alpha_continuous(range = c(0.1, 0.5)) +
            scale_fill_manual(values = c(region_cols)) +
            scale_colour_manual(values = c(region_cols)) +
            ylab(paste0("Scaled ", state, " score")) +
            xlab(paste0("Scaled ", factor, " activation")) +
            GBM_theme()
        ggsave(
            filename = paste0(
                "FigS3e_",
                state,
                "_score_vs_",
                factor,
                "_scaled.pdf"
            ),
            path = params$plot_dir
        )
    }

    names(cors) <- states
    corr_df[[factor]] <- cors
}

corr_df <- as.data.frame(corr_df)

pdf(
    file = file.path(params$plot_dir, "FigS3d_factor_state_correlation.pdf"),
    width = 12,
    height = 10
)
pushViewport(viewport(gp = gpar(fontfamily = "Helvetica")))
heatmap <- Heatmap(
    as.matrix(corr_df),
    name = "Factor-State correlation",
    height = nrow(corr_df) * unit(2, "cm"),
    width = ncol(corr_df) * unit(5, "cm"),
    col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), rev(brewer.pal(5, "RdBu"))),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_heatmap_legend = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE
)
draw(heatmap, newpage = FALSE)
popViewport()
dev.off()

# ---------------------------------------------------------------------------- #
#                                 Figure S3f-j                                 #
# ---------------------------------------------------------------------------- #

df <- metadata
df$Region <- factor(df$Region, levels = c("PT", "TE", "TC"))
df <- df %>%
    filter(CellClass_L3 %in% c("Malignant_OPC", "Malignant_NPC1")) %>%
    select(c(
        "Region",
        "Patient",
        paste0("Factor", seq(1, length(colnames(usage))))
    )) %>%
    pivot_longer(
        -c("Region", "Patient"),
        values_to = "Activation",
        names_to = "Factor"
    )

df$Factor <- case_when(
    df$Factor == "Factor5" ~ "Factor 5 (Hypoxia)",
    df$Factor == "Factor4" ~ "Factor 4 (Neural Crest)",
    df$Factor == "Factor3" ~ "Factor 3 (Neuronal)",
    df$Factor == "Factor2" ~ "Factor 2 (Cilia)",
    df$Factor == "Factor1" ~ "Factor 1 (Cell Cycle)",
)

for (factor in unique(df$Factor)) {
    curr_factor <- df %>%
        filter(Factor == factor)

    p <- ggplot(curr_factor, aes(x = Region, y = Activation)) +
        geom_boxplot(aes(fill = Region), alpha = 0.9, outlier.shape = NA) +
        scale_fill_manual(values = region_cols) +
        GBM_theme() +
        ylab("Factor activation") +
        ggtitle(paste0("Activation of ", factor, " in OPC/NPC1-like cells"))

    p <- p +
        geom_signif_lmm(
            data_df = curr_factor,
            response = "Activation",
            condition = "Region",
            latent_vars = c("Patient"),
            comparisons = list(c("PT", "TE"), c("TE", "TC"), c("PT", "TC")),
            step_increase = c(0, 0.1, 0.1)
        )
    ggsave(
        filename = paste0("FigS3f_OPC_NPC1_", factor, ".pdf"),
        path = params$plot_dir,
        height = 8,
        width = 8
    )
}


df <- metadata %>%
    mutate(proliferating = ifelse(Factor1 > 0, "Cycling", "NonCycling")) %>%
    mutate(Region = ifelse(Region == "PT", "PT", "Tumor")) %>%
    filter(CellClass_L3 %in% c("Malignant_OPC", "Malignant_NPC1")) %>%
    group_by_at(c("Region", "proliferating", "Patient")) %>%
    summarise(num_cycling = n()) %>%
    group_by_at(c("Region", "Patient")) %>%
    mutate(total = sum(num_cycling), percent = (num_cycling / total) * 100) %>%
    filter(proliferating == "Cycling")
df$Region <- factor(df$Region, levels = c("PT", "Tumor"))
p <- ggplot(df, aes(x = Region, y = percent)) +
    geom_boxplot(aes(fill = Region)) +
    geom_point(size = 2) +
    theme_classic() +
    scale_fill_manual(values = c(region_cols[1], "#7F7F7F")) +
    stat_compare_means(comparisons = list(c("PT", "Tumor")))
ggsave(
    filename = "FigS3f_OPC_NPC1_ONLY_factor1_cycling_cells.pdf",
    path = params$plot_dir,
    height = 7,
    width = 5
)
