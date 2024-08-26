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
        description = "Compute average expression",
    )
    parser$add_argument("--input_file", help = "Seurat input file", type = "character")
    parser$add_argument("--assay", type = "character", help = "Assay to use", default = "RNA")
    parser$add_argument("--slot", type = "character", help = "Layer (v5)/slot (v4) to use", default = "counts")
    parser$add_argument("--groupbyvar", type = "character", help = "group_by for AverageExpression")
    parser$add_argument("--is_seurat_v4", type = "numeric", help = "Is Seurat v4 version (default = 1)", default = TRUE)
    parser$add_argument("--prefix", type = "character", help = "Prefix for saving file", default = "")
    parser$add_argument("--apply_is_confident_filter", type = "numeric", help = "Filter for `Confident_Annotation` with gbm_regional_study.rds (only scRNAseq)", default = FALSE)

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
pacman::p_load(Seurat)
if (args$is_seurat_v4) {
    options(Seurat.object.assay.version = "v4")
}
log_info("Load seurat object...")
seurat_obj <- readRDS(args$input_file)
DefaultAssay(seurat_obj) <- args$assay

groupbyvar <- str_split(args$groupbyvar, ",") %>% unlist()

if (args$apply_is_confident_filter) {
    # Only keep confidently annotated cells
    seurat_obj <- subset(seurat_obj, subset = Confident_Annotation)
}

if (args$is_seurat_v4) {
    log_info("Compute Average Expression based on Seurat v4...")

    avg_expr <- AverageExpression(seurat_obj,
        assays = args$assay,
        group.by = groupbyvar,
        slot = args$slot
    )[[args$assay]]
} else {
    log_info("Compute Average Expression based on Seurat v5...")
    if (args$assay == "RNA") {
        seurat_obj <- NormalizeData(seurat_obj)
        seurat_obj <- ScaleData(seurat_obj)
    }

    avg_expr <- AverageExpression(seurat_obj, assays = args$assay, layer = args$slot, group.by = args$groupbyvar)[[args$assay]]
    # %>%
    # data.frame() %>%
    # rownames_to_column("gene")
}

log_info("Save results...")
groupbyvar_str <- paste0(groupbyvar, collapse = "__")

prefix <- ifelse(args$prefix == "", "", paste0(args$prefix, "__"))

saveRDS(avg_expr,
    file = glue("{args$output_dir}/{prefix}mean_exp_by_{groupbyvar_str}_{args$assay}_{args$slot}.rds")
)

log_info("Finished!")
