# ---- Code to reproduce Figure 3a ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# Load required packages
pacman::p_load(
    varhandle,
    ggplot2,
    ggrepel,
    mitch,
    fgsea,
    readr,
    msigdbr,
    escape,
    dittoSeq,
    tidyr,
    dplyr,
    stringr,
    data.table,
    ComplexHeatmap,
    colorRamp2,
    foreach,
    doParallel,
    tidyverse,
    viridis,
    circlize,
    RColorBrewer
)

# Required inputs
params <- list(
    plot_dir = "output/figures",
    # TF activity
    gene_auc_activator_mtx_path = "multiome_results/12_SCENIC_plus/outs/Plots/all_regions/Filtered_TFs/gene_auc_activator_mtx.csv", # SCENIC+ results
    region_auc_activator_mtx_path = "multiome_results/12_SCENIC_plus/outs/Plots/all_regions/Filtered_TFs/region_auc_activator_mtx.csv", # SCENIC+ resultso_label_path = "multiome_results/12_SCENIC_plus/outs/Plots/all_regions/Filtered_TFs/region_TFs_to_label.csv", # Information in Table S2

    # Heatmap data rna and atac
    heatmap_data_rna_path = "multiome_results/12_SCENIC_plus/outs/RNA_heatmap.csv",
    heatmap_data_atac_path = "multiome_results/12_SCENIC_plus/outs/ATAC_heatmap.csv",
    TFs_to_plot_path = "misc/Table S2.xlsx", # can be downloaded online, see publication
)

GaitiLabUtils::create_dir(params$plot_dir)

# ---- Load data  ---- #
OPC_NPC1_markers <- readxl::read_excel(params$genelists_csv_path, skip = 1) %>%
    filter(!is.na(Neftel_OPC), !is.na(Neftel_NPC1)) %>%
    dplyr::select(Neftel_OPC, Neftel_NPC1) %>%
    unlist() %>%
    unname() %>%
    unique()

heatmap_data_rna <- fread(params$heatmap_data_rna_path) |>
    column_to_rownames("TF")
heatmap_data_atac <- fread(params$heatmap_data_atac_path) |>
    column_to_rownames("TF")

# Load TFs to plot
TFs_to_plot <- readxl::read_excel(
    params$TFs_to_plot_path,
    sheet = "GRN",
    skip = 1
)

# ---- Data wrangling ---- #
TF_names_to_plot <- TFs_to_plot %>%
    filter(Cell_type == "Invasive-high OPC/NPC1") %>%
    pull(TF)

# Subset RNA and ATAC heatmap data
heatmap_data_rna_inv <- heatmap_data_rna[TF_names_to_plot, , drop = FALSE]
heatmap_data_atac_inv <- heatmap_data_atac[TF_names_to_plot, , drop = FALSE]

# Define discrete values and manual markers
discrete_values <- 0:5

# Extract the blue half and red half of the RdBu palette and create discrete color functions
full_rd_bu <- rev(brewer.pal(11, "RdBu"))
blue_half <- full_rd_bu[6:1] # "#F7F7F7" "#D1E5F0" "#92C5DE" "#4393C3" "#2166AC" "#053061"
red_half <- full_rd_bu[6:11] # "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B" "#67001F"
rna_discrete_col_fun <- structure(
    blue_half,
    names = as.character(discrete_values)
)
atac_discrete_col_fun <- structure(
    red_half,
    names = as.character(discrete_values)
)

# Label rows with manual markers
manual_inv_markers <- c(
    "OLIG2",
    "MEOX2",
    "NKX6-2",
    "ASCL1",
    "SOX4",
    "TCF4",
    "ZEB1",
    "ZEB2",
    "E2F1",
    "ETV1",
    "HOXD3",
    "MEIS1",
    "SREBF2",
    "CPEB1",
    "E2F2",
    "HOXB3"
)
search_terms <- paste(c(OPC_NPC1_markers, manual_inv_markers), collapse = "|")
label_positions <- grep(search_terms, rownames(heatmap_data_rna_inv))
label_names <- grep(search_terms, rownames(heatmap_data_rna_inv), value = TRUE)

# Row annotation with lines pointing to labels
ha <- rowAnnotation(
    mark = anno_mark(
        at = label_positions,
        labels = label_names,
        side = "right",
        labels_gp = gpar(fontsize = 8),
        link_gp = gpar(col = "black")
    )
)

# Fixed width for the heatmap
heatmap_width <- unit(1, "in")

# Save RNA Heatmap (Blue Half)
pdf(
    file.path(params$plot_dir, "Fig3a_RNA_heatmap_inv.pdf"),
    width = 4,
    height = 8
)
Heatmap(
    heatmap_data_rna_inv,
    name = "RNA",
    col = rna_discrete_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    right_annotation = ha,
    row_title = "TF",
    column_title = "CellType",
    width = heatmap_width,
    heatmap_legend_param = list(
        title = "Number of patients",
        at = discrete_values,
        labels = as.character(discrete_values)
    )
)
dev.off()

# Save ATAC Heatmap (Red Half)
pdf(
    file.path(params$plot_dir, "Fig3a_ATAC_heatmap_inv.pdf"),
    width = 4,
    height = 8
)
Heatmap(
    heatmap_data_atac_inv,
    name = "ATAC",
    col = atac_discrete_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    right_annotation = ha,
    row_title = "TF",
    column_title = "CellType",
    width = heatmap_width,
    heatmap_legend_param = list(
        title = "Number of patients",
        at = discrete_values,
        labels = as.character(discrete_values)
    )
)
dev.off()
