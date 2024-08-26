# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Generate labels for cell type x interaction (binarization)",
    )
    parser$add_argument("--interactions", type = "character", help = "Excel sheet with interactions")
    parser$add_argument("--sample_id", type = "character", help = "Sample ID")
    parser$add_argument("--meta_with_exp", type = "character", help = "Metadata with expression counts (RNA)")
    args <- parser$parse_args()
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Load additional libraries
pacman::p_load(readxl)
meta_with_exp <- readRDS(args$meta_with_exp) %>% column_to_rownames("cell_id")

interactions <- read_excel(args$interactions) %>%
    select(source, target, ligand_complex, receptor_complex) %>%
    distinct() %>%
    mutate(
        source = str_replace_all(source, "Invasive-high OPC/NPC1", "Progenitor_like"),
        target = str_replace_all(target, "Invasive-high OPC/NPC1", "Progenitor_like"),
        cellpair_x_interaction = paste0(source, "__", target, "_x_", ligand_complex, "__", receptor_complex)
    )

labels <- apply(interactions, 1, function(interaction_oi, meta_df) {
    ligand_oi <- interaction_oi[["ligand_complex"]]
    receptor_oi <- interaction_oi[["receptor_complex"]]
    source_oi <- interaction_oi[["source"]]
    target_oi <- interaction_oi[["target"]]
    cellpair_x_interaction_oi <- interaction_oi[["cellpair_x_interaction"]]
    meta_df %>%
        mutate(dummy = case_when(
            # Handle source
            (cell_type == source_oi) & (!!sym(ligand_oi) > 0) ~ paste0(source_oi, ":", ligand_oi),
            (cell_type == source_oi) & (!!sym(ligand_oi) == 0) ~ source_oi,
            # Handle target
            (cell_type == target_oi) & (!!sym(receptor_oi) > 0) ~ paste0(target_oi, ":", receptor_oi),
            (cell_type == target_oi) & (!!sym(receptor_oi) == 0) ~ target_oi,
            .default = "Other"
        )) %>%
        pull(dummy)
}, meta_df = meta_with_exp) %>% data.frame()
colnames(labels) <- interactions %>% pull(cellpair_x_interaction)

log_info("Save list of cell pair x interaction as txt...")
write.table(data.frame(interactions %>% pull(cellpair_x_interaction)),
    col.names = FALSE, quote = FALSE, row.names = FALSE, file = glue("{args$output_dir}/cellpair_x_interaction.txt")
)

log_info("Combine cell position info with labels...")
out <- cbind(meta_with_exp %>% rownames_to_column("cell_id") %>% select(cell_id, x_centroid, y_centroid, cell_type), labels)

log_info("Save as RDS...")
saveRDS(out, file = glue("{args$output_dir}/{args$sample_id}__BANKSY__meta_w_celltype_x_gene_pair_labels__ROI.rds"))

log_info("Finished!")
