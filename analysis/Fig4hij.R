# ---- Code to reproduce Figure 4h,i,j ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(
    data.table,
    ggplot2,
    dplyr,
    amethyst,
    umap,
    tibble,
    tidyr,
    future,
    furrr,
    purrr,
    stringr,
    RColorBrewer,
    Matrix,
    patchwork,
    irlba,
    leiden,
    fgsea,
    lme4,
    mitch,
    ggrepel
)

# Required inputs
params <- list(
    output_dir = "output/submission",
    plot_dir = "output/submission/figures",
    amethyst_obj_path = "path_to_amethyst_obj",
    c2_gmt_path = "c2.all.v2023.1.Hs.symbols.gmt", # download from msigdb
    genelists_path = "misc/gene_signatures.xlsx"
)

GaitiLabUtils::create_dir(params$output_dir)
GaitiLabUtils::create_dir(params$plot_dir)

set.seed(123)

# --------------------------------- functions -------------------------------- #

#' Function to perform DMR at genes. Modified from amethyst R package
#'
#' @param obj Amethyst object
#' @param gene_matrix gene x cell methylation matrix in genomeMatrices in obj
#' @param genes genes to analyze
#' @param group_by grouping variable, by default cluster_id
#' @param method wilcox or lmm
#' @param filter_cov whether to filter cells based on coverage at gene promoter
#' @param min_cov min coverage to keep cell if using filter_cov
#' @param cov_matrix gene x cell coverage matrix
#' @param min_cells minimum number of cells with coverage in group for testing a gene
#' @param threads number of threads for multithreading
#'
#' @return dataframe with results
findGeneMarkers <- function(
    obj,
    matrix = "gene_promoter_cg",
    genes = rownames(obj@genomeMatrices[["gene_cg"]]),
    group_by = "cluster_id",
    min_cells = 5,
    threads = 1,
    verbose = FALSE
) {
    options(scipen = 3)

    print(cov_matrix)

    genematrix <- as.matrix(obj@genomeMatrices[[matrix]])
    print(rownames(genematrix)[1:10])
    print(colnames(genematrix)[1:10])

    obj@metadata[["grouping_var"]] <- obj@metadata[[group_by]]
    grouping_var <- group_by
    print(grouping_var)

    membership <- obj@metadata |> dplyr::select("grouping_var")

    # need to ensure that only cells in the data matrix are being used
    membership$cell_id <- rownames(membership)
    membership <- membership |> dplyr::filter(cell_id %in% colnames(genematrix))
    membership$cell_id <- NULL
    print(head(membership))

    # Set up multithreading
    if (threads > 1) {
        future::plan(future::multicore, workers = threads)
    }

    results <- furrr::future_map(
        .x = genes,
        .f = function(gene) {
            gene_results <- list() # Initialize outside the loop
            for (id in unique(membership$grouping_var)) {
                members <- rownames(
                    membership |> dplyr::filter(grouping_var == id)
                )
                nonmembers <- rownames(
                    membership |> dplyr::filter(grouping_var != id)
                )
                tryCatch(
                    {
                        members_data <- genematrix[gene, members]
                        nonmembers_data <- genematrix[gene, nonmembers]

                        if (
                            sum(!is.na(members_data)) < min_cells ||
                                sum(!is.na(nonmembers_data)) < min_cells
                        ) {
                            stop("Not enough cells covered for testing.")
                        }

                        # test differential methylation
                        gene_results[[id]] <- data.frame(
                            "p.val" = stats::wilcox.test(
                                x = members_data,
                                y = nonmembers_data
                            )$p.value,
                            "gene" = gene,
                            "grouping_var" = id,
                            mean_1 = mean(members_data, na.rm = TRUE),
                            mean_2 = mean(nonmembers_data, na.rm = TRUE)
                        ) |>
                            dplyr::mutate(
                                logFC = log2(mean_2 / mean_1),
                                Delta = mean_2 - mean_1,
                                direction = ifelse(
                                    mean_1 > mean_2,
                                    "hypermethylated",
                                    "hypomethylated"
                                )
                            )
                        print(gene_results)
                    },
                    error = function(e) {
                        cat(
                            "Error processing gene:",
                            gene,
                            "and cluster:",
                            id,
                            ". Error: ",
                            e$message,
                            "\n"
                        )
                        gene_results[[id]] <- NA
                    }
                )
            }
            gene_results
        },
        .progress = verbose
    )

    if (threads > 1) {
        future::plan(future::sequential)
        gc()
    }

    results <- do.call(rbind, lapply(results, function(x) do.call(rbind, x)))

    results <- results |>
        dplyr::group_by(grouping_var) |>
        dplyr::mutate(p.adj = stats::p.adjust(p.val, method = "bonferroni")) |>
        dplyr::select(p.val, p.adj, everything())
    return(results)
}

# ---------------------------------------------------------------------------- #
#                                   Figure 4h                                  #
# ---------------------------------------------------------------------------- #
amethyst_obj <- readRDS(params$amethyst_obj_path)
amethyst_obj <- subsetObject(
    amethyst_obj,
    rownames(amethyst_obj@metadata)[
        amethyst_obj@metadata$CellClass == "Malignant"
    ]
)

amethyst_obj@genomeMatrices[["gene_promoter_cg"]] <- makeWindows(
    amethyst_obj,
    genes = unique(amethyst_obj@ref$gene_name[
        amethyst_obj@ref$gene_type == "protein_coding"
    ]),
    promoter = TRUE,
    type = "CG",
    metric = "percent",
    threads = 12,
    index = "chr_cg"
)

amethyst_obj@metadata$cluster_id <- amethyst_obj@metadata$Region
print(table(amethyst_obj@metadata$Region))
print(table(amethyst_obj@metadata$CNV))

saveRDS(amethyst_obj, file.path(params$output_dir, "combined_obj_new.rds"))

amethyst_obj <- readRDS(file.path(
    params$output_dir,
    "combined_obj_current.rds"
))


amethyst_obj@metadata$cluster_id <- amethyst_obj@metadata$Region

markers <- findGeneMarkers(
    amethyst_obj,
    matrix = "gene_promoter_cg",
    genes = rownames(amethyst_obj@genomeMatrices[["gene_promoter_cg"]]),
    threads = 10,
    method = "wilcox",
    filter_cov = FALSE
)
write.csv(
    markers,
    file.path(params$output_dir, "pt_te_dmr_gene_promoter_cg.csv")
)

## GSEA
curr_markers <- markers %>%
    filter(grouping_var == "PT") %>%
    mutate(p.adj = p.adjust(p.val, method = "BH")) %>%
    mutate(rank = -log10(p.adj) * abs(logFC) / logFC) %>%
    filter(rank != Inf & rank != -Inf) %>%
    arrange(-rank)

ranked_gene_list <- curr_markers$rank
names(ranked_gene_list) <- curr_markers$gene

# load msigdb geneset
c2 <- gmtPathways(params$c2_gmt_path)

gsea_results <- fgsea::fgseaMultilevel(
    pathways = c2,
    stats = ranked_gene_list,
    minSize = 10,
    maxSize = 5000
)

df <- as.data.frame(gsea_results)
df <- df[!grepl("leadingEdge", colnames(gsea_results))] %>%
    filter(padj < 0.5)
write.csv(
    df,
    file.path(
        params$output_dir,
        paste0("fgseaRes_", pathway_set, "_gene_promoter_cg_pt_vs_te.csv")
    )
)

## plot volcano and gsea results
markers <- markers %>%
    filter(grouping_var == "PT") %>%
    mutate(fdr = p.adjust(p.val, "BH"))

# TODO replace with SuppTable use
gene_lists <- readxl::read_excel(params$genelists_path, skip = 1)


pathway_gobert <- gene_lists %>%
    filter(!is.na(Gobert_Oligodendrocyte_Differentiation_Up)) %>%
    pull(Gobert_Oligodendrocyte_Differentiation_Up)
pathway_kim <- (gene_lists %>%
    filter(!is.na(Kim_All_Disorders_Oligodendrocyte_Abundance_Up)) %>%
    pull(Kim_All_Disorders_Oligodendrocyte_Abundance_Up))
pathway_neuro <- gene_lists %>%
    filter(!is.na(KEGG_Neuroactive_Ligands_and_Receptors)) %>%
    pull(KEGG_Neuroactive_Ligands_and_Receptors)

sig_dn <- cluster_markers %>% filter(Delta < -2.5 & p.val < 0.05)
sig_up <- cluster_markers %>% filter(Delta > 2.5 & p.val < 0.05)

sig_genes_dn <- intersect(
    c(
        pathway_gobert,
        pathway_kim
    ),
    sig_dn$gene
)
sig_genes_up <- intersect(
    pathway_neuro,
    sig_up$gene
)
sig_genes <- c(sig_genes_dn, sig_genes_up)

genes_to_label <- c(
    "CREB5",
    "CCDC61",
    "PEG3",
    "MARCKSL1",
    "FAM171A1",
    "MAST3",
    "RHOA",
    "GPR156",
    "GRM6",
    "CHRNA2",
    "GRIN1",
    "P2RX2",
    "HRH1"
)
cluster_markers$label <- ifelse(
    cluster_markers$gene %in% genes_to_label,
    cluster_markers$gene,
    NA
)

cluster_markers$colour <- case_when(
    cluster_markers$gene %in%
        sig_genes &
        cluster_markers$gene %in% pathway_gobert ~
        "gobert",
    cluster_markers$gene %in%
        sig_genes &
        cluster_markers$gene %in%
            pathway_kim ~
        "kim",
    cluster_markers$gene %in%
        sig_genes &
        cluster_markers$gene %in%
            pathway_neuro ~
        "neuroactive",
    TRUE ~ "black"
)
cluster_markers$alpha_val <- ifelse(cluster_markers$colour == "black", 0.5, 1)
cluster_markers <- cluster_markers %>%
    mutate(label_1 = ifelse(is.na(label), 0, 1)) %>%
    arrange(label_1)

p <- ggplot(cluster_markers, aes(x = Delta, y = -log10(p.val))) +
    geom_vline(xintercept = 2.5, linetype = "dashed") +
    geom_vline(xintercept = -2.5, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_point(aes(colour = colour, alpha = alpha_val), stroke = NA, size = 2) +
    geom_label_repel(
        aes(label = label, colour = colour),
        na.rm = TRUE,
        size = 4,
        force = 1,
        nudge_y = 0.2,
        nudge_x = 0.2,
        max.overlaps = 15
    ) +
    scale_colour_manual(values = c("black", "#FFC759", "#BB4430", "#7EBDC2")) +
    theme_classic() +
    xlab("DNA methylation difference (%)") +
    ylab("-log10(P-Value)") +
    ylim(c(0, 6)) +
    xlim(c(-10, 10))
ggsave(
    plot = p,
    filename = "Fig4h_pt_vs_te_volcano.pdf",
    path = params$plot_dir,
    width = 6,
    height = 8
)

# ---------------------------------------------------------------------------- #
#                                   Figure 4i                                  #
# ---------------------------------------------------------------------------- #
pd <- fgsea::plotEnrichmentData(
    pathway = pathway_kim,
    stats = ranked_gene_list
)
with(
    pd,
    ggplot(data = curve) +
        geom_line(aes(x = rank, y = ES), color = "red4", linewidth = 4) +
        geom_ribbon(
            aes(x = rank, ymin = 0, ymax = ES),
            fill = "red4",
            alpha = 0.5
        ) +
        geom_ribbon(
            data = stats,
            mapping = aes(
                x = rank,
                ymin = 0,
                ymax = stat / maxAbsStat * (spreadES / 4)
            ),
            fill = "grey"
        ) +
        geom_segment(
            data = ticks,
            mapping = aes(
                x = rank,
                y = -spreadES / 16,
                xend = rank,
                yend = spreadES / 16
            ),
            size = 0.3
        ) +
        geom_hline(yintercept = posES, colour = "red4", linetype = "dashed") +
        geom_hline(yintercept = negES, colour = "red4", linetype = "dashed") +
        geom_hline(yintercept = 0, colour = "black") +
        ylim(c(-0.2, 0.4)) +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_line(color = "grey92"),
            plot.title = element_text(hjust = 0.5, size = 30),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            panel.border = element_rect(colour = "black", fill = NA, size = 2)
        ) +
        labs(x = "Rank in Ordered Dataset", y = "Enrichment score (ES)")
) +
    ggtitle(paste("PT vs Tumour GBM cells \n Oligodendrocyte Abundance Up"))
ggsave(filename = "Fig4i_oligo_abundance_gsea.pdf", path = params$plot_dir)

# ---------------------------------------------------------------------------- #
#                                   Figure 4j                                  #
# ---------------------------------------------------------------------------- #
pd <- fgsea::plotEnrichmentData(
    pathway = pathway_neuro,
    stats = ranked_gene_list
)
with(
    pd,
    ggplot(data = curve) +
        geom_line(aes(x = rank, y = ES), color = "#7EBDC2", linewidth = 4) +
        geom_ribbon(
            aes(x = rank, ymin = 0, ymax = ES),
            fill = "#7EBDC2",
            alpha = 0.5
        ) +
        geom_ribbon(
            data = stats,
            mapping = aes(
                x = rank,
                ymin = 0,
                ymax = stat / maxAbsStat * (spreadES / 4)
            ),
            fill = "grey"
        ) +
        geom_segment(
            data = ticks,
            mapping = aes(
                x = rank,
                y = -spreadES / 16,
                xend = rank,
                yend = spreadES / 16
            ),
            size = 0.3
        ) +
        geom_hline(
            yintercept = posES,
            colour = "#7EBDC2",
            linetype = "dashed"
        ) +
        geom_hline(
            yintercept = negES,
            colour = "#7EBDC2",
            linetype = "dashed"
        ) +
        geom_hline(yintercept = 0, colour = "black") +
        ylim(c(-0.2, 0.4)) +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_line(color = "grey92"),
            plot.title = element_text(hjust = 0.5, size = 30),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            panel.border = element_rect(colour = "black", fill = NA, size = 2)
        ) +
        labs(x = "Rank in Ordered Dataset", y = "Enrichment score (ES)")
) +
    ggtitle(paste(
        "PT vs Tumour GBM cells \n Neuroactive Ligand and Receptor Interaction"
    ))
ggsave(filename = "Fig4j_neuroactive_gsea.pdf", path = params$plot_dir)
