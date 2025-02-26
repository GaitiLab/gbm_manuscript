# ---- Code to reproduce Figure 2de ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# Load need packages
pacman::p_load(
    GaitiLabUtils,
    GBMutils,
    glue,
    data.table,
    tidyverse,
    stringr,
    sf,
    rjson,
    Seurat
)
logr <- GaitiLabUtils::init_logging()

params <- list(
    segmentation_masks_path = "data/visiumhd/processed/6425_A/6425_A__masks.json",
    seurat_obj_path = "data/visiumhd/processed/6425_A/6425_A.rds",
    plot_dir = "output/submission/figures",
    k_neighbors = 30,
    lower_limit = 0.10,
    upper_limit = 0.90,
    scale_factors_path = glue(
        "data/visiumhd/spatial_outs/6425_A/outs/binned_outputs/square_002um/spatial/scalefactors_json.json"
    )
)

GaitiLabUtils::create_dir(params$plot_dir)

# ---- Setup ---- #
# Setup color palette
log_info("Load color palette...")
palette_L1 <- GBMutils::load_color_palette("Spatial_CellClass_L1")

cell_type_pair_oi <- c("Neuron", "Malignant")

myfuns_nam <- list(Close = min, Far = max)
signif_level <- 0.05

# Create color palette for distance labels
distance_labels <- setNames(
    c(
        paste0("<", params$lower_limit * 100, "%"),
        paste0(
            params$lower_limit * 100,
            "% <= mean distance <= ",
            params$upper_limit * 100,
            "%"
        ),
        paste0(">", params$upper_limit * 100, "%")
    ),
    c("bottom", "between", "top")
)

palette_dist_gradient <- setNames(
    RColorBrewer::brewer.pal(n = 8, name = "Blues")[c(3, 5, 8)],
    distance_labels
)
palette_malign_dist_grad_x_neuron <- c(
    palette_L1["Neuron"],
    palette_dist_gradient
)


# ---- Create Figure 2d ---- #
log_info("Load scalefactors...")
scalefactors <-
    fromJSON(file = params$scale_factors_path)

log_info("Load segmentation masks...")
gdf <- read_sf(params$segmentation_masks_path) %>%
    mutate(geometry_px = geometry)
# Remove coordinate systems, just want to use cartesian system
st_crs(gdf) <- NA

log_info("Load Seurat object...")
seurat_obj <- readRDS(params$seurat_obj_path)

log_info("Combine metadata and segmentation masks...")
meta_sdf <- sf::st_as_sf(
    (seurat_obj@meta.data) %>%
        left_join(gdf, by = c("cell_id" = "id"))
) %>%
    # Convert to px to microns
    mutate(geometry = geometry_px * scalefactors$microns_per_pixel) %>%
    filter(Spatial_CellClass_L1 %in% cell_type_pair_oi) %>%
    # Convert to px to microns
    mutate(
        vis_label = factor(
            ifelse(Spatial_CellClass_L1 == "Neuron", "Neuron", cat_dist),
            levels = names(colors)
        )
    )
list_val <- lapply(
    myfuns_nam,
    function(f)
        f(
            meta_sdf %>%
                filter(Spatial_CellClass_L1 == "Malignant") %>%
                pull(mean_dist)
        )
)

p <- ggplot() +
    geom_sf(
        data = meta_sdf %>% filter(Spatial_CellClass_L1 == "Malignant"),
        aes(fill = mean_dist),
        color = NA
    ) +
    GBMutils::GBM_theme() +
    scale_fill_distiller(
        palette = "Blues",
        direction = -1,
        breaks = unlist(list_val),
        guide = guide_colorbar(
            override.aes = list(size = 5),
            title.position = "top",
            title = "Malignant cell by its proximity to 30 nearest Neurons",
            barwidth = unit(8, "lines"),
            barheight = unit(0.8, "lines"),
        )
    ) +
    ggnewscale::new_scale_fill() +
    geom_sf(
        data = meta_sdf %>% filter(Spatial_CellClass_L1 == "Neuron"),
        fill = palette_L1[["Neuron"]],
        color = NA
    ) +
    coord_sf() +
    scale_y_reverse() +
    labs(
        x = "X (µm)",
        y = "Y (µm)",
        subtitle = paste0(
            "No. cells=",
            prettyNum(
                meta_sdf %>%
                    filter(Spatial_CellClass_L1 %in% cell_type_pair_oi) %>%
                    nrow(),
                big.mark = ","
            )
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
    filename = glue(
        "Fig2d.pdf"
    ),
    path = params$plot_dir
)
