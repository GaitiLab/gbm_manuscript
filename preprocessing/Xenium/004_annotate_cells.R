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
        description = "Annotate clusters/cells and save in metadata, also add gene-count matrix (for gating)",
    )
    parser$add_argument("--sample_id",
        type = "character", help = "Sample ID"
    )
    parser$add_argument("--seurat_obj_path",
        type = "character", help = "Path to Seurat object"
    )
    parser$add_argument("--cluster_varname",
        type = "character",
        help = "Variable in Xenium metadata that contains the clusters"
    )
    parser$add_argument("--annotated_clusters_path",
        type = "character", help = "Excel file with annotations"
    )
    parser$add_argument("--cells_oi_path", type = "character", help = "Path to file with cells of interest")

    params <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    params <- list()
    params$log_level <- 5
    params$output_dir <- "data/Xenium/processed"
    params$annotated_clusters_path <- "misc/xenium_annot.xlsx"
    params$seurat_obj_path <- "data/Xenium/processed/6509_A__BANKSY.rds"
    params$sample_id <- "6509_A"
    params$cells_oi_path <- "misc/6509_A_roi_cell_ids.csv"
}

# Set up logging
logr <- init_logging(log_level = params$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(params$output_dir)

# Load additional libraries
pacman::p_load(Seurat, readxl)


log_info("Load seurat object...")
seurat_obj <- readRDS(params$seurat_obj_path)

log_info("Load Excel file with annotated clusters...")
annot_clusters <- read_excel(params$annotated_clusters_path) %>% mutate(BANKSY_snn_res.0.5 = factor(BANKSY_snn_res.0.5))

log_info("Extract gene-count matrix...")
mat <- data.frame(t(data.matrix(LayerData(seurat_obj, assay = "RNA", layer = "counts"))))

log_info("Annotate cells...")
meta_df <- seurat_obj@meta.data %>%
    # Annotate cells
    left_join(annot_clusters) %>%
    left_join(
        mat %>% rownames_to_column("cell_id"),
        by = "cell_id"
    ) %>%
    # If cluster wasn't annotated set to 'Undetermined'
    mutate(cell_type = ifelse(is.na(cell_type), "Undetermined", cell_type)) %>%
    mutate(
        cell_type = case_when(
            # Check for DLL3+ expression, if so, then Progenitor_like
            BANKSY_snn_res.0.5 == 14 & DLL3 > 0 ~
                "Invasive-high OPC/NPC1",
            BANKSY_snn_res.0.5 == 14 & ((DLL3 == 0)) ~ "OPC",
            .default = cell_type
        )
    )

log_info("Update metadata in seurat object...")
rownames(meta_df) <- meta_df$cell_id
seurat_obj@meta.data <- meta_df

log_info("Save metadata + annot + expr...")
saveRDS(
    obj = seurat_obj@meta.data,
    file = file.path(
        params$output_dir,
        glue("{params$sample_id}__BANKSY__meta_annot_w_expr_counts.rds")
    )
)

log_info("Save updated seurat object...")
saveRDS(
    obj = seurat_obj,
    file = file.path(
        params$output_dir,
        glue("{params$sample_id}__BANKSY__annot_w_expr_counts.rds")
    )
)

log_info("Obtain cell IDs of interest from CSV file...")
cells_oi <- data.frame(
    fread(params$cells_oi_path,
        sep = ",", skip = 2
    )
) %>% pull(Cell.ID)

log_info("Save metadata + annot + expr for ROI...")
saveRDS(
    obj = (subset(seurat_obj, cells = cells_oi))@meta.data,
    file = file.path(
        params$output_dir,
        glue("{params$sample_id}__BANKSY__meta_annot_w_expr_counts__ROI.rds")
    )
)

log_info("Save seurat object for ROI...")
saveRDS(
    obj = subset(seurat_obj, cells = cells_oi),
    file = file.path(
        params$output_dir,
        glue("{params$sample_id}__BANKSY__annot_w_expr_counts__ROI.rds")
    )
)

log_info("Finished!")
