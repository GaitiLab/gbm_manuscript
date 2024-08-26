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
        description = "Extract (unique) interactions for a specific cell type pair",
    )
    parser$add_argument("--input_file", type = "character", help = "File with filtered aggregated results (402c)")
    parser$add_argument("--condition_varname", type = "character", help = "Variable indicating the 'condition' or 'group'")
    parser$add_argument("--pval_type", type = "character", help = "use 'pval' or 'pval_adj' for filtering", default = "pval_adj")
    parser$add_argument("--condition_oi", type = "character", help = "Condition or group of interest")
    parser$add_argument("--alpha", type = "numeric", help = "Level to use for p-val filtering", default = 0.05)
    parser$add_argument("--output_name", type = "character", help = "Name for output")
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
pacman::p_load(xlsx)

celltype_pairs_oi <- c("Neuron__Invasive-high OPC/NPC1", "Neuron__Progenitor_like", "Invasive-high OPC/NPC1__Neuron", "Progenitor_like__Neuron")

output_filename <- glue("{args$output_dir}/{args$output_name}.xlsx")
if (file.exists(output_filename)) {
    file.remove(output_filename)
}

log_info("Load data + formatting...")
obj <- readRDS(args$input_file) %>%
    # Additional filtering on p-value (aggregate rank)
    filter(!!sym(args$pval_type) < args$alpha) %>%
    separate(source_target, c("source", "target"), sep = "__", remove = FALSE)


obj_undirected <- obj %>%
    # Remove direction by sorting source-target alphabetically
    rowwise() %>%
    mutate(source_target_undirected = paste0(sort(c(source, target)), collapse = "__")) %>%
    # Remove duplicate interactions, when they are found in both directions only keep 1
    distinct(!!sym(args$condition_varname), source_target_undirected, complex_interaction, .keep_all = FALSE)

log_info("Filtering + conversion to wide-format...")
obj_filtered <- obj_undirected %>%
    filter(!!sym(args$condition_varname) == args$condition_oi, source_target_undirected %in% celltype_pairs_oi) %>%
    mutate(dummy = 1) %>%
    pivot_wider(names_from = source_target_undirected, values_from = dummy, values_fill = 0) %>%
    # Mark unique interactions
    mutate(is_unique = ifelse((`Invasive-high OPC/NPC1__Neuron` == 1) & (Neuron__Progenitor_like == 0), 1, 0))

log_info("Get 'unique' interactions...")
unique_interactions <- obj_filtered %>%
    filter(is_unique == 1) %>%
    pull(complex_interaction)

log_info("Formatting...")
unique_interactions_df <- obj %>%
    filter(
        complex_interaction %in% unique_interactions, source_target %in% celltype_pairs_oi,
        !!sym(args$condition_varname) == args$condition_oi
    ) %>%
    data.frame()

log_info("Write data to Excel...")
write.xlsx(unique_interactions_df,
    file = output_filename,
    col.names = TRUE,
    row.names = FALSE,
    append = FALSE
)
log_info("Finished!")
