# Code to reproduce Fig 4c,d,e

if(!("pacman" %in% rownames(installed.packages()))){
  install.packages("pacman")
}

pacman::p_load(data.table,
               patchwork,
               stringr,
               ggplot2,
               dplyr,
               tidyr,
               ggrepel,
               ggtext,
               ggpubr,
               fgsea,
               mitch
               )

region_cols <- c(PT = "#0173b2", TE = "#de8f05", TC = "#029e73")

# Enter paths here
output_dir <- ""

# Creating directories needed for outputs
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_dir <- paste0(output_dir, "/path_to_output_directory")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---------------------------------------------------------------------------- #
#                                   Figure 4c                                  #
# ---------------------------------------------------------------------------- #

dmr <- read.table("/GLM_TSSallcells_DATA_TSS_inv_high_Progenitor_like_GLM_MINCPGS_5_GLOBAL_inv_high_vs_pl_like.txt", header = TRUE)

pathway <- gmtPathways("/GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_DN.v2023.2.Hs.gmt")

sig_dn <- dmr %>% filter(DELTA < -0.05 & PVAL < 0.05)
sig_genes <- intersect(pathway[["GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_DN"]], sig_dn$GENE)
dmr$label <- ifelse(dmr$GENE %in% sig_genes, dmr$GENE, NA)
dmr$colour <- ifelse(is.na(dmr$label), "black", "darkorange")
dmr$alpha_val <- ifelse(dmr$colour == "black", 0.8, 1)

p <- ggplot(dmr, aes(x = DELTA*100, y = -log10(PVAL))) + 
    geom_point(aes(colour = colour, alpha = alpha_val), size = 3, stroke = NA) + 
    geom_label_repel(aes(label = label, colour = colour), na.rm = TRUE, size = 4, force = 1.5, nudge_y = 0.2, nudge_x = 0.3) +
    geom_vline(xintercept = 5, linetype = "dashed") +
    geom_vline(xintercept = -5, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_colour_manual(values = c("black", "darkorange")) +
    scale_alpha_continuous(c(0.9,1)) +
    theme_classic() +
    xlab("DNA methylation difference (%)") +
    ylab("-log10(P)") +
    xlim(c(-50, 50))
ggsave(filename = "oligo_dn_volcano_inv_high_vs_pl_like.pdf", path = plot_dir, height = 7, width = 6)

# ---------------------------------------------------------------------------- #
#                                   Figure 4d                                  #
# ---------------------------------------------------------------------------- #

gene_set <- gmtPathways("/c2.cgp.v2023.2.Hs.symbols.gmt")
ranked_genes_response <- dmr %>%
    mutate(rank = -log10(PVAL) * DELTA/abs(DELTA)) %>%
    filter(rank != Inf) %>% filter(rank != -Inf) %>% arrange(-rank)     

ranked_genes_response$ranking <- rank(-ranked_genes_response$rank, ties.method = "first")

response_level_stats <- ranked_genes_response$rank
names(response_level_stats) <- ranked_genes_response$GENE

fgseaRes_ctrl_multilevel <- fgseaMultilevel(pathways = gene_set,
                                            stats = response_level_stats,
                                            minSize=10,
                                            maxSize=2000,
                                            nPermSimple = 10000) %>% 
filter(padj < 0.5) %>% arrange(-abs(NES))

df <- as.data.frame(fgseaRes_ctrl_multilevel)
df <- df[!grepl("leadingEdge", colnames(fgseaRes_ctrl_multilevel))]
write.csv(df, file.path(plot_dir, "fgseaRes.csv"))

pathway <- gmt_import("GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_DN.v2023.2.Hs.gmt")
pd <- plotEnrichmentData(
  pathway = pathway[["GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_DN"]],
  stats = response_level_stats
)
with(pd,
    ggplot(data=curve) +
        geom_line(aes(x=rank, y=ES), color="darkorange", linewidth = 4) +
        geom_ribbon(aes(x=rank, ymin=0, ymax=ES), fill="darkorange", alpha=0.2) + # Added this line
        geom_ribbon(data=stats,
                    mapping=aes(x=rank, ymin=0,
                                ymax=stat/maxAbsStat*(spreadES/4)),
                    fill="grey") +
        geom_segment(data=ticks,
                    mapping=aes(x=rank, y=-spreadES/16,
                                xend=rank, yend=spreadES/16),
                    size=0.3) +
        geom_hline(yintercept=posES, colour="red", linetype="dashed") +
        geom_hline(yintercept=negES, colour="red", linetype="dashed") +
        geom_hline(yintercept=0, colour="black") +
        ylim(c(-0.35, 0.1))+
        theme(
            panel.background = element_blank(),
            panel.grid.major=element_line(color="grey92"),
            plot.title = element_text(hjust = 0.5, size = 30),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            panel.border = element_rect(colour = "black", fill = NA, size = 2)
        ) +
        labs(x="Rank in Ordered Dataset", y="Enrichment score (ES)") +
        ggtitle(paste("Invasive-high vs invasive-low OPC/NPC1 \n Gobert et al. Oligodendrocyte Differentiation Down"))
)

ggsave(filename = "gobert_oligo_diff.pdf", path = "/Users/bensonwu/Downloads", width = 11)

# ---------------------------------------------------------------------------- #
#                                   Figure 4e                                  #
# ---------------------------------------------------------------------------- #
curr_gene_list <- read.csv("de_results.csv") # import de results from fig 1i

rownames(curr_gene_list) <- curr_gene_list$gene

gene_set <- gmtPathways("/c5.go.v2023.2.Hs.symbols.gmt")

ranked_genes_response <- curr_gene_list %>%
    mutate(rank = log2FoldChange) %>%
    filter(rank != Inf) %>% filter(rank != -Inf) %>% arrange(-rank)

ranked_genes_response$ranking <- rank(-ranked_genes_response$rank, ties.method = "first")
response_level_stats <- ranked_genes_response$rank
names(response_level_stats) <- rownames(ranked_genes_response)

pathway <- gmt_import("/GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP.v2023.2.Hs.gmt")
pd <- plotEnrichmentData(
  pathway = pathway[["GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP"]],
  stats = response_level_stats
)
with(pd,
    ggplot(data=curve) +
        geom_line(aes(x=rank, y=ES), color="darkred", linewidth = 4) +
        geom_ribbon(aes(x=rank, ymin=0, ymax=ES), fill="darkred", alpha = 0.5) +
        geom_ribbon(data=stats,
                    mapping=aes(x=rank, ymin=0,
                                ymax=stat/maxAbsStat*(spreadES/4)),
                    fill="grey") +
        geom_segment(data=ticks,
                    mapping=aes(x=rank, y=-spreadES/16,
                                xend=rank, yend=spreadES/16),
                    size=0.3) +
        geom_hline(yintercept=posES, colour="red", linetype="dashed") +
        geom_hline(yintercept=negES, colour="red", linetype="dashed") +
        geom_hline(yintercept=0, colour="black") +
        ylim(c(-0.6, 0.2))+
        theme(
            panel.background = element_blank(),
            panel.grid.major=element_line(color="grey92"),
            plot.title = element_text(hjust = 0.5, size = 30),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            panel.border = element_rect(colour = "black", fill = NA, size = 2)
        ) +
        labs(x="Ranked gene list", y="Enrichment score (ES)")) +
        ggtitle(paste("PT vs Tumour OPC/NPC1-like cells \n Oligodendrocyte Differentiation Up"))
ggsave(filename = "oligo_diff_up_RNA.pdf", path = plot_dir)
