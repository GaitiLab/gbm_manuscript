# ---- Code to reproduce Figure S5e ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(
    GaitiLabUtils,
    GBMutils,
    glue,
    data.table,
    tidyverse,
    stringr,
    readxl,
    ggplot2,
    ggrastr,
    ggtext
)

logr <- GaitiLabUtils::init_logging()

# Required inputs
params <- list(
    # path pointing to the parent directory with sample directories (masks + seurat obj)
    input_dir = "data/visiumhd/processed",
    # path pointing to the spatial_outs directory with the two 6425 samples (A + B) can be downloaded, see publication
    spatial_outs_parent_dir = "data/visiumhd/spatial_outs/",
    plot_dir = "output/figures"
)
GaitiLabUtils::create_dir(params$plot_dir)

sample_ids <- c(
    "6425_A",
    "6425_B"
)

palette_L1 <- GBMutils::load_color_palette("Spatial_CellClass_L1")

current_sample_params <- list()
for (current_sample_id in sample_ids) {
    log_info(glue("Current sample: {current_sample_id}..."))
    current_sample_params$segmentation_masks_path <- file.path(
        params$input_dir,
        current_sample_id,
        paste0(current_sample_id, "__masks.json")
    )
    current_sample_params$seurat_obj_path <- file.path(
        params$input_dir,
        current_sample_id,
        paste0(current_sample_id, ".rds")
    )

    current_sample_params$scale_factors_path <- file.path(
        params$spatial_outs_parent_dir,
        current_sample_id,
        "outs/binned_outputs/square_002um/spatial/scalefactors_json.json"
    )

    log_info("Load scalefactors...")
    scalefactors <-
        fromJSON(file = current_sample_params$scale_factors_path)

    log_info("Load segmentation masks...")
    gdf <- read_sf(current_sample_params$segmentation_masks_path) %>%
        mutate(geometry_px = geometry)
    # Remove coordinate systems, just want to use cartesian system
    st_crs(gdf) <- NA

    log_info("Load Seurat object...")
    seurat_obj <- readRDS(current_sample_params$seurat_obj_path)

    log_info("Combine metadata and segmentation masks...")
    meta_sdf <- sf::st_as_sf(
        (seurat_obj@meta.data) %>%
            left_join(gdf, by = c("cell_id" = "id"))
    ) %>%
        # Convert to px to microns
        mutate(geometry = geometry_px * scalefactors$microns_per_pixel)

    log_info("Plot spatial distribution...")
    p <- ggplot(meta_sdf) +
        geom_sf(aes(fill = Spatial_CellClass_L1), color = NA) +
        GBMutils::GBM_theme() +
        scale_fill_manual(
            values = palette_L1,
            guide = guide_legend(
                override.aes = list(size = 5),
                title = "Cell type"
            )
        ) +
        coord_sf() +
        scale_y_reverse() +
        labs(
            x = "X (µm)",
            y = "Y (µm)",
            subtitle = paste0(
                "No. cells=",
                prettyNum(nrow(meta_sdf), big.mark = ",")
            )
        ) +
        theme(
            aspect.ratio = 1,
            axis.text = element_blank(),
            axis.ticks = element_blank()
        )

    log_info("Save figure as PDF...")
    ggsave(
        plot = p,
        filename = glue("FigS5e_{current_sample_id}_spatial.pdf"),
        path = params$plot_dir,
        width = 8,
        height = 8
    )
}
