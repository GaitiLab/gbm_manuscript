# ---- Code to reproduce Figure S7d ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(
    argparse,
    varhandle,
    log4r,
    ArchR,
    presto,
    ggplot2,
    ggrepel,
    mitch,
    fgsea,
    readr,
    msigdbr,
    escape,
    dittoSeq,
    tidyr,
    dplyr,
    stringr,
    data.table,
    ComplexHeatmap,
    colorRamp2,
)

# Required inputs
params <- list(
    input = "misc/marker_list.csv",
    plot_dir = "output/figures"
)

# Load data
marker_list <- readr::read_csv(params$input)

# Specify the log2FC and FDR cutoffs
curr.log2FC <- 0.1
curr.fdr <- 0.05

# Align significant and insignificant genes
marker_list <- marker_list %>%
    mutate(
        diffaccessible = case_when(
            log2FC > curr.log2FC & FDR <= curr.fdr ~
                "Accessible in Invasive-high OPC/NPC1",
            log2FC < -curr.log2FC & FDR <= curr.fdr ~
                "Accessible in Progenitor-like (NPC + OPC)",
            TRUE ~ "Not significant"
        )
    )

# change the order of the factor
marker_list$diffaccessible <- factor(
    marker_list$diffaccessible,
    levels = c(
        "Accessible in Invasive-high OPC/NPC1",
        "Not significant",
        "Accessible in Progenitor-like (NPC + OPC)"
    )
)

# Generate label for plot annotation
marker_list$dalabel <- NA
marker_list$dalabel[
    marker_list$diffaccessible != "Not significant"
] <- marker_list$name[marker_list$diffaccessible != "Not significant"]

# Overlay the volcano plot with the NOTCH signaling pathway genes
all_gene_sets <- msigdbr(species = "Homo sapiens")
NOTCH <- all_gene_sets %>%
    filter(gs_name %in% c("WP_NOTCH_SIGNALING")) %>%
    pull(gene_symbol)
Oligo <- all_gene_sets %>%
    filter(gs_name %in% c("GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_DN")) %>%
    pull(gene_symbol)

marker_list <- marker_list %>%
    mutate(
        highlight = case_when(
            name %in%
                NOTCH &
                diffaccessible == "Accessible in invasive-high OPC/NPC1" ~
                "NOTCH",
            name %in%
                Oligo &
                diffaccessible == "Accessible in invasive-high OPC/NPC1" ~
                "Oligo",
            TRUE ~ "no"
        )
    )
color_palette <- c(
    "yes" = "#d61f26",
    "no" = "black",
    "specific_gene" = "#00798CFF",
    "NOTCH" = "#89181A",
    "Oligo" = "#F58B1F"
)

# Create the volcano plot with highlighted genes
volc <- ggplot(
    data = marker_list,
    aes(x = log2FC, y = -log10(FDR), label = name)
) +
    geom_point(aes(color = highlight), size = 1, alpha = 0.25, stroke = NA) +
    geom_point(
        data = marker_list[marker_list$highlight != "no", ],
        aes(color = highlight),
        size = 1.5
    ) +
    geom_vline(xintercept = curr.log2FC, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = -curr.log2FC, linetype = "dashed", alpha = 0.5) +
    geom_hline(
        yintercept = -log10(curr.fdr),
        linetype = "dashed",
        alpha = 0.5
    ) +
    expand_limits(
        x = c(
            -ceiling(max(abs(marker_list$log2FC)) * 10) / 10,
            ceiling(max(abs(marker_list$log2FC)) * 10) / 10
        )
    ) +
    geom_label_repel(
        data = marker_list[marker_list$highlight %in% c("NOTCH"), ],
        fill = "white",
        label.padding = 0.1,
        box.padding = 0.5,
        size = 3,
        min.segment.length = 0,
        aes(color = highlight),
        show.legend = FALSE,
        max.overlaps = Inf,
        force = 2,
        nudge_x = ifelse(
            ceiling(max(abs(marker_list$log2FC)) * 10) / 10 > 1,
            0.3,
            0.15
        ),
        nudge_y = 0.3
    ) +
    geom_label_repel(
        data = marker_list[marker_list$name %in% c("TNR", "GPR17", "DSCAM"), ],
        fill = "white",
        label.padding = 0.1,
        box.padding = 0.5,
        size = 3,
        min.segment.length = 0,
        aes(color = highlight),
        show.legend = FALSE,
        max.overlaps = Inf,
        force = 2,
        nudge_x = ifelse(
            ceiling(max(abs(marker_list$log2FC)) * 10) / 10 > 1,
            0.3,
            0.15
        ),
        nudge_y = 0.3
    ) +
    scale_color_manual(
        values = color_palette,
        guide = guide_legend(title = "")
    ) +
    guides(color = guide_legend(title = "")) +
    labs(
        x = expression(Log[2] * " Fold Change (Gene accessibility socre)"),
        y = expression(-Log[10] * " FDR")
    ) +
    theme_minimal() +
    theme(
        legend.key = element_blank(),
        legend.position = "bottom",
        plot.title.position = "plot",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
ggsave(
    filename = paste0(
        "FigS7d_DiffAccess_",
        "all",
        "_low_vs_high_highlight_specific_gene.pdf"
    ),
    path = params$plot_dir,
    units = "in",
    width = 5,
    height = 6
)
