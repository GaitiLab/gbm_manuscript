# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
require(GBMutils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Generate table with overlapping interactions with xenium panel",
    )
    parser$add_argument("--xenium_panel", type = "character", help = "Path to Xenium panel Exccel file")
    parser$add_argument("--interactions", type = "character", help = "Path to interactions of interest (Excel file)")
    parser$add_argument("--detected_interactions", type = "character", help = "Path to detected interactions (RDS file)")

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
pacman::p_load(xlsx, readxl)

contains_all_units <- function(row, xenium_genes) {
    length(intersect(
        row[row != ""],
        xenium_genes
    )) == length(row[row != ""])
}

xenium_genes <- read_excel(
    args$xenium_panel,
    skip = 1
) %>% pull(`Gene(10x)`)

interactions <- read_excel(args$interactions) %>%
    separate(complex_interaction,
        into = c("ligand_complex", "receptor_complex"),
        sep = "__"
    ) %>%
    distinct(ligand_complex, receptor_complex)

detected_interactions <- read_excel(args$detected_interactions) %>%
    separate(complex_interaction, into = c("ligand_complex", "receptor_complex"), sep = "__", remove = FALSE) %>%
    separate(source_target, into = c("source", "target"), sep = "__")

# ---- Detected interactions  ----
cell_types_oi <- c("Neuron", "Invasive-high OPC/NPC1", "Progenitor_like")
detected_interactions_oi <- detected_interactions %>%
    filter(
        source %in% cell_types_oi, target %in% cell_types_oi,
        Region == "PT",
        source != target,
        pval_adj < 0.05
    ) %>%
    distinct(ligand_complex, receptor_complex, .keep_all = TRUE)


ligand_units_detected <- detected_interactions_oi %>%
    pull(ligand_complex) %>%
    str_split(., "\\:", simplify = TRUE)
receptor_units_detected <- detected_interactions_oi %>%
    pull(receptor_complex) %>%
    str_split(., "\\:", simplify = TRUE)

ligand_units_detected_in_xenium <- apply(ligand_units_detected, 1, contains_all_units, xenium_genes = xenium_genes)
receptor_units_detected_in_xenium <- apply(receptor_units_detected, 1, contains_all_units, xenium_genes = xenium_genes)
detected_interactions_oi_in_xenium <- detected_interactions_oi[ligand_units_detected_in_xenium & receptor_units_detected_in_xenium, ]


# Handle unique interactions
ligand_units_oi <- interactions %>%
    pull(ligand_complex) %>%
    str_split(., "\\:", simplify = TRUE)
receptor_units_oi <- interactions %>%
    pull(receptor_complex) %>%
    str_split(., "\\:", simplify = TRUE)

ligand_units_oi_in_xenium <- apply(ligand_units_oi, 1, contains_all_units, xenium_genes = xenium_genes)
receptor_units_oi_in_xenium <- apply(receptor_units_oi, 1, contains_all_units, xenium_genes = xenium_genes)
unique_overlap <- interactions[ligand_units_oi_in_xenium & receptor_units_oi_in_xenium, ] %>% mutate(is_unique = TRUE)

common_overlap <- detected_interactions_oi %>%
    filter(
        ligand_complex %in% xenium_genes, receptor_complex %in% xenium_genes,
        !(source == "Progenitor_like" & target == "Invasive-high OPC/NPC1"),
        !(target == "Progenitor_like" & source == "Invasive-high OPC/NPC1"),
    )
# Mark if interaction is 'unique'
combined_overlap <- common_overlap %>% left_join(unique_overlap)

out_filename <- glue("{args$output_dir}/overlap_w_xenium.xlsx")
if (file.exists(out_filename)) {
    file.remove(out_filename)
}
write.xlsx(combined_overlap %>% as.data.frame(), out_filename,
    col.names = TRUE, row.names = FALSE, append = TRUE
)
log_info("Finished!")
