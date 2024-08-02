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
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/")
    args$input_file <- "output/CCI_CellClass_L2_2_rerun/100_preprocessing/seurat/6234_2895153_A.rds"
    args$groupbyvar <- "Sample,CCI_CellClass_L2_2"
    args$is_seurat_v4 <- TRUE
    args$slot <- "data"
    args$assay <- "RNA"
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
obj <- readRDS(args$input_file)

groupbyvar <- str_split(args$groupbyvar, ",") %>% unlist()

log_info("Compute Average Expression...")
# GBM - For scatterplot
# Only keep `Confident_Annotation` cells
avg_expr <- AverageExpression(subset(obj, subset = Confident_Annotation),
    assays = args$assay,
    group.by = groupbyvar,
    slot = args$slot
)[[args$assay]]

log_info("Save results...")
groupbyvar_str <- paste0(groupbyvar, collapse = "__")

saveRDS(avg_expr,
    file = glue("{args$output_dir}/mean_exp_by_{groupbyvar_str}_{args$assay}_{args$slot}.rds")
)

log_info("COMPLETED!")
