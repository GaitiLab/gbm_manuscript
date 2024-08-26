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
        description = "Annotate clusters/cells and save in metadata, also add gene-counts",
    )
    parser$add_argument("--sample_id", type = "character", help = "Sample ID")
    parser$add_argument("--sheet_name", type = "character", help = "Sheet to use from Excel file ('annot_clusters') for annotation")
    parser$add_argument("--meta", type = "character", help = "path to metadata file")
    parser$add_argument("--cluster_varname", type = "character", help = "Variable in Xenium metadata that contains the clusters")
    parser$add_argument("--exp_counts", type = "character", help = "path to file with sparse gene-count matrix")
    parser$add_argument("--annot_clusters", type = "character", help = "Excel file with annotations")
    parser$add_argument("--cells_oi", type = "character", help = "Path to file with cells of interest")

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

log_info("Load gene counts matrix...")
# Rows = cells, columns = genes
mat <- data.frame(t(data.matrix(readRDS(args$exp_counts)))) %>% rownames_to_column("cell_id")

log_info("Load Excel file with annotated clusters...")
annot_clusters <- read_excel(args$annot_clusters, sheet = args$sheet_name) %>% filter(Sample == args$sample_id)

log_info("Load metadata...")
meta <- readRDS(args$meta)

log_info("Create lookup table...")
lookup_df <- do.call(
    rbind,
    lapply(
        annot_clusters %>%
            dplyr::pull(Label),
        GaitiLabUtils::create_lookup_table,
        annot = annot_clusters,
        label_name = "cell_type"
    )
) %>% dplyr::filter(!is.na(cluster))

log_info("Annotate clusters based on manual annotation...")
meta <- meta %>%
    # Add expression (RNA - counts)
    left_join(mat, by = "cell_id") %>%
    rename(cluster = !!sym(args$cluster_varname)) %>%
    left_join(lookup_df, by = "cluster") %>%
    # If cluster wasn't annotated set to 'Undetermined'
    mutate(cell_type = ifelse(is.na(cell_type), "Undetermined", cell_type)) %>%
    mutate(
        cell_type = case_when(
            # Check for DLL3+ expression, if so, then Progenitor_like
            cluster == 14 & DLL3 > 0 ~
                "Progenitor_like",
            cluster == 14 & ((DLL3 == 0)) ~ "OPC",
            .default = cell_type
        )
    )

log_info("Save metadata + annot + expr...")
saveRDS(
    obj = meta,
    file = glue("{args$output_dir}/{args$sample_id}__BANKSY__meta_annot_w_expr_counts.rds")
)

log_info("Obtain cell IDs of interest from CSV file...")
cells_oi <- data.frame(fread(args$cells_oi, sep = ",", skip = 2)) %>% pull(Cell.ID)

log_info("Only keep cells in region of interest...")
meta <- meta %>% filter(cell_id %in% cells_oi)

log_info("Save metadata + annot + expr...")
saveRDS(
    obj = meta,
    file = glue("{args$output_dir}/{args$sample_id}__BANKSY__meta_annot_w_expr_counts__ROI.rds")
)

log_info("Finished!")
