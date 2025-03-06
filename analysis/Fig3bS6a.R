# ---- Code to reproduce Figure 3b and S6a ---- #

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
GaitiLabUtils::set_wd()

# ---- Setup script ---- #

# Load required packages
pacman::p_load(
    argparse,
    varhandle,
    log4r,
    ArchR,
    BSgenome.Hsapiens.UCSC.hg38,
    cowplot,
    stringr,
    dplyr,
    data.table,
    reshape,
    ggpubr,
    rstatix,
    tibble,
    ggrepel,
    ggplot2,
    circlize,
    ComplexHeatmap,
    viridis,
    GaitiLabUtils,
    GBMutils,
    ggh4x
)

# Required inputs
params <- list(
    input = "multiome_results/10_ArchR", # Please use the data provided in the publication and follow ArchR workflow to generate the data
    background_peaks = "per_patient", # Calculate with per patient or all cells sum up - to reproduce the results in the paper, use "per_patient"
    plot_dir = "10_ArchR/Plots",
    expression_matrix = "TF_exp_mtx.csv" # Subset expression matrix for TFs = generated from data from publication and subsetted to TFs of interest
)

params$output_dir <- params$input
plot_dir <- params$plot_dir

log_info("Loading TF deviations...")
all_patient_obj <- readRDS(paste0(plot_dir, "all_patient_obj.rds"))
merged_matrix <- do.call(rbind, all_patient_obj)
head(merged_matrix)

# Load archR project
archr_proj <- loadArchRProject(curr_proj_dir)
cell_type_column <- paste0("Seurat_", params$celltype_column)
patient_column <- paste0("Seurat_", params$patient_column)
region_column <- paste0("Seurat_", params$region_column)
confidence_column <- paste0("Seurat_", params$confidence_column)

# Method: TF with high deviation from ChromVAR in inv high - prog and high correlation to gene expression
log_info("Step 1: Identifying correlated TF motifs and TF expression...")
corGEM_MM <- correlateMatrices(
    ArchRProj = archr_proj,
    useMatrix1 = "GeneExpressionMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "LSI_Combined"
)
fwrite(as.data.frame(corGEM_MM), paste0(plot_dir, "/corGEM_MM.csv"))
corGEM_MM <- fread(paste0(plot_dir, "/corGEM_MM.csv"))

log_info("Step 2: Identifying Deviant TF Motifs...")
if (params$background_peaks == "all_cells") {
    log_info("Using all_cells background peaks...")
    seGroupMotif <- getGroupSE(
        ArchRProj = archr_proj,
        useMatrix = "MotifMatrix",
        groupBy = "Seurat_CCI_CellClass_L2_2"
    )
    seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames == "z", ]

    # Identify the maximum delta in z-score between Invasive-high OPC/NPC1 and Progenitor_like
    head(assay(seZ))
    rowData(seZ)$Inv_delta <- assay(seZ)[, "Invasive-high OPC/NPC1"] -
        assay(seZ)[, "Progenitor-like"]

    log_info(
        "Step 3: Add maximum delta deviation to the correlation data frame..."
    )
    corGEM_MM$Inv_delta <- as.numeric(rowData(seZ)[
        match(corGEM_MM$MotifMatrix_name, rowData(seZ)$name),
        "Inv_delta"
    ])
} else if (params$background_peaks == "per_patient") {
    log_info("Using per_patient background peaks...")
    all_patient_obj <- readRDS(paste0(plot_dir, "/all_patient_obj.rds"))
    merged_matrix <- do.call(rbind, all_patient_obj)
    head(merged_matrix)

    # cell type annotation
    cell_annotation <- getCellColData(
        archr_proj,
        select = c("Seurat_CCI_CellClass_L2_2")
    ) |>
        as.data.frame() |>
        rownames_to_column(var = "CellID")
    head(cell_annotation)

    Inv_delta_mean <- merged_matrix |>
        as.data.frame() |>
        rownames_to_column(var = "CellID") |>
        mutate(
            sample = sub("#.*", "", CellID),
            patient = sub("_.*", "", sample)
        ) |>
        left_join(cell_annotation, by = "CellID") |>
        filter(
            Seurat_CCI_CellClass_L2_2 %in%
                c("Invasive-high OPC/NPC1", "Progenitor-like")
        ) |>
        group_by(patient, Seurat_CCI_CellClass_L2_2) |>
        dplyr::select(-CellID) |>
        summarise_all(mean) |>
        ungroup() |>
        group_by(Seurat_CCI_CellClass_L2_2) |>
        dplyr::select(-patient) |>
        summarise_all(mean) |>
        ungroup() |>
        tidyr::pivot_longer(
            cols = -Seurat_CCI_CellClass_L2_2,
            names_to = "TF",
            values_to = "Mean"
        ) |>
        tidyr::pivot_wider(
            names_from = Seurat_CCI_CellClass_L2_2,
            values_from = Mean
        ) |>
        mutate(Inv_delta = `Invasive-high OPC/NPC1` - `Progenitor-like`)

    log_info(
        "Step 3: Add maximum delta deviation to the correlation data frame..."
    )
    corGEM_MM$Inv_delta <- as.numeric(Inv_delta_mean$Inv_delta[match(
        corGEM_MM$MotifMatrix_name,
        Inv_delta_mean$TF
    )])
}

# Filter out duplicated TFs
corGEM_MM <- corGEM_MM[order(abs(corGEM_MM$cor), decreasing = TRUE), ]
corGEM_MM <- corGEM_MM[
    which(!duplicated(gsub("\\-.*", "", corGEM_MM$MotifMatrix_name))),
]

# TF to label in the plot
high_SCENIC_TF <- motifs_df$TF[motifs_df$Subgroup == "Invasive-high OPC/NPC1"]
high_SCENIC_TF <- paste0(high_SCENIC_TF, "_")

low_SCENIC_TF <- motifs_df$TF[motifs_df$Subgroup == "Differentiated-like"]
low_SCENIC_TF <- paste0(low_SCENIC_TF, "_")
low_SCENIC_TF <- sort(low_SCENIC_TF)

# NPC/OPC markers
NPC_OPC_markers <- c("SOX4", "SOX11", "OLIG1", "ETV1")
NPC_OPC_markers <- paste0(NPC_OPC_markers, "_")

# Manual invasive markers
manual_inv_makers <- c(
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
manual_inv_makers <- paste0(manual_inv_makers, "_")

# AP-1 family TFs
AP1_TF <- c("FOS", "FOSB", "FOSL1", "FOSL2", "FOSL2", "JUN", "JUNB", "JUND")

# Motifs to highlight in ChromVAR
Positive_TF <- corGEM_MM$MotifMatrix_name[
    corGEM_MM$cor > 0 &
        corGEM_MM$Inv_delta > quantile(corGEM_MM$Inv_delta, 0.95)
]
Positive_TF <- Positive_TF[!is.na(Positive_TF)] # remove NA
Positive_TF_name <- substr(Positive_TF, 1, regexpr("_", Positive_TF) - 1) # remove number after _ in the name
Positive_TF_name <- paste0(Positive_TF_name, "_")

# Find high_SCENIC_TF intersect with Positive_TF
motifs_to_plot <- high_SCENIC_TF[high_SCENIC_TF %in% Positive_TF_name]

# Find TFs in each category
high_SCENIC_TF_to_plot <- high_SCENIC_TF[
    !sapply(
        high_SCENIC_TF,
        function(x) any(sapply(motifs_to_plot, function(y) grepl(y, x)))
    )
] # Identified in SCENIC but not in ChromVAR
high_SCENIC_TF_to_plot <- corGEM_MM$MotifMatrix_name[sapply(
    corGEM_MM$MotifMatrix_name,
    function(x) any(sapply(high_SCENIC_TF_to_plot, function(y) grepl(y, x)))
)]
manual_inv_makers_to_plot <- manual_inv_makers[
    !sapply(
        manual_inv_makers,
        function(x) any(sapply(motifs_to_plot, function(y) grepl(y, x)))
    )
] # Identified in SCENIC but not in ChromVAR
manual_inv_makers_to_plot <- corGEM_MM$MotifMatrix_name[sapply(
    corGEM_MM$MotifMatrix_name,
    function(x) any(sapply(manual_inv_makers_to_plot, function(y) grepl(y, x)))
)]
motifs_to_plot <- Positive_TF[sapply(
    Positive_TF,
    function(x) any(sapply(motifs_to_plot, function(y) grepl(y, x)))
)] # Identified in both SCENIC and ChromVAR

# Label the TFs
corGEM_MM$TFRegulator <- "Not significant"
corGEM_MM$TFRegulator[
    corGEM_MM$MotifMatrix_name %in% high_SCENIC_TF_to_plot
] <- "Candidate identified in SCENIC+"
corGEM_MM$TFRegulator[
    corGEM_MM$MotifMatrix_name %in% motifs_to_plot
] <- "Putative regulator"
corGEM_MM$TFRegulator[sapply(
    corGEM_MM$MotifMatrix_name,
    function(x) any(sapply(NPC_OPC_markers, function(y) grepl(y, x)))
)] <- "NPC1/OPC markers"
corGEM_MM$TFRegulator[sapply(
    corGEM_MM$MotifMatrix_name,
    function(x) any(sapply(AP1_TF, function(y) grepl(y, x)))
)] <- "AP-1 family TFs"

# Label only manual_inv_makers_to_plot within "Candidate identified in SCENIC+"
corGEM_MM$Label <- ifelse(
    !corGEM_MM$TFRegulator %in%
        c("Not significant", "Candidate identified in SCENIC+"),
    corGEM_MM$MotifMatrix_name,
    NA
)
manual_inv_makers_indices <- which(
    corGEM_MM$MotifMatrix_name %in% manual_inv_makers_to_plot
)
corGEM_MM$Label[manual_inv_makers_indices] <- corGEM_MM$MotifMatrix_name[
    manual_inv_makers_indices
]

# Fix MotifMatrix_name - delete everything after MotifMatrix_name
corGEM_MM$MotifMatrix_name <- substr(
    corGEM_MM$MotifMatrix_name,
    1,
    regexpr("_", corGEM_MM$MotifMatrix_name) - 1
)

# Visualize the correlation between TF motif and gene expression
ggplot(data.frame(corGEM_MM), aes(cor, Inv_delta, color = TFRegulator)) +
    geom_point(
        data = data.frame(corGEM_MM[
            corGEM_MM$TFRegulator == "Not significant",
        ]),
        aes(cor, Inv_delta, color = TFRegulator),
        size = 1,
        alpha = 0.75
    ) +
    geom_point(
        data = data.frame(corGEM_MM[
            corGEM_MM$TFRegulator != "Not significant",
        ]),
        aes(cor, Inv_delta, color = TFRegulator),
        size = 2
    ) +
    GBM_theme() +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_hline(
        yintercept = quantile(corGEM_MM$Inv_delta, 0.95),
        lty = "dashed",
        color = "darkgrey"
    ) +
    scale_color_manual(
        values = c(
            "Not significant" = "darkgrey",
            "Putative regulator" = "#2B7095",
            "Candidate identified in SCENIC+" = "#EDAE49FF",
            "NPC1/OPC markers" = "#7C9EB5",
            "AP-1 family TFs" = "#C05E00"
        )
    ) +
    geom_label_repel(
        data = data.frame(corGEM_MM[!is.na(corGEM_MM$Label), ]),
        aes(label = MotifMatrix_name),
        size = 3,
        nudge_x = 0.15,
        nudge_y = 0.15,
        max.overlaps = 10
    ) +
    labs(
        y = "TF motif accessibility difference between \ninvasive-high OPC/NPC1 and progenitor-like (Δz-score)",
        x = "Correlation of TF motif accessibility and TF expression"
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(
            min(corGEM_MM$Inv_delta) * 1.05,
            max(corGEM_MM$Inv_delta) * 1.05
        )
    ) +
    scale_x_continuous(
        expand = c(0, 0),
        limits = c(-1, 1)
    ) +
    theme(legend.position = "bottom", legend.title = element_blank())
ggsave(
    paste0(
        plot_dir,
        "/corGEM_MM_TF_Regulator_95_new_SCENIC+_highlight_",
        params$background_peaks,
        ".pdf"
    ),
    width = 7,
    height = 7.5
)


# cell type annotation
cell_annotation <- getCellColData(archr_proj, select = c(cell_type_column)) |>
    as.data.frame() |>
    rownames_to_column(var = "CellID")
head(cell_annotation)

merged_matrix <- merged_matrix |>
    as.data.frame() |>
    rownames_to_column(var = "CellID") |>
    mutate(
        sample = sub("#.*", "", CellID),
        patient = sub("_.*", "", sample)
    ) |>
    left_join(cell_annotation, by = "CellID") |>
    mutate(celltype = get(cell_type_column)) |>
    dplyr::select(-cell_type_column)
colnames(merged_matrix) <- gsub("_[0-9]*", "", colnames(merged_matrix))
merged_matrix <- merged_matrix |>
    tidyr::pivot_longer(
        cols = -c(CellID, sample, patient, celltype),
        names_to = "TF",
        values_to = "ChromVar"
    )

# motif to plot for ChromVAR
motifs_to_plot <- c("TCF12", "TCF4", "TCF3", "ASCL1", "ZEB1", "FOSL1", "FOSL2")

df_neuronal <- merged_matrix |>
    filter(TF %in% motifs_to_plot)

# Load the expression data
exp_mat <- fread(params$expression_matrix)
exp_mat$archr_cellnames <- sub("-1.*", "", exp_mat$archr_cellnames)

for (TF_name in unique(df_neuronal$TF)) {
    log_info("Currently analyzing: ", TF_name)
    curr_TF_chromVAR <- df_neuronal[df_neuronal$TF == TF_name, ] %>%
        mutate(archr_cellnames = CellID)
    curr_TF_chromVAR$archr_cellnames <- sub(
        "-1.*",
        "",
        curr_TF_chromVAR$archr_cellnames
    )
    curr_TF_exp <- exp_mat %>%
        filter(archr_cellnames %in% curr_TF_chromVAR$archr_cellnames) %>%
        dplyr::select(archr_cellnames, TF_name)
    curr_df <- merge(curr_TF_chromVAR, curr_TF_exp, by = "archr_cellnames") %>%
        mutate(
            Sample = sub("#.*", "", archr_cellnames),
            Barcode = sub(".*#", "", archr_cellnames),
            Patient = sub("_.*", "", Sample),
            x = celltype
        )

    curr_df <- curr_df %>%
        group_by(Patient, x) %>%
        dplyr::summarise(
            mean_gex = mean(get(TF_name)),
            mean_acc = mean(ChromVar),
            n_cell = n()
        ) %>%
        ungroup() %>%
        filter(n_cell > 10) # filter out low patient

    curr_df <- curr_df %>%
        group_by(x) %>%
        dplyr::summarise(
            mean_cell_type_gex = mean(mean_gex),
            mean_cell_type_acc = mean(mean_acc),
            total_n = sum(n_cell)
        )

    # Rename cell type
    curr_df <- curr_df %>%
        mutate(
            x = case_when(
                x == "Malignant_NPC1" ~ "NPC1-like",
                x == "Malignant_NPC2" ~ "NPC2-like",
                x == "Malignant_OPC" ~ "OPC-like",
                x == "Invasive-high_OPC_NPC1" ~ "Invasive-high OPC/NPC1",
                x == "Malignant_AC" ~ "AC-like",
                x == "Malignant_MES_AST" ~ "MES-AC-like",
                x == "Malignant_MES_INT" ~ "MES-INT-like",
                x == "Malignant_MES_HYP" ~ "MES-HYP-like"
            )
        )

    color_palette <- c(
        "NPC1-like" = "#86CCE8",
        "NPC2-like" = "#ADD7E5",
        "OPC-like" = "#AECFA5",
        "AC-like" = "#DFB63B",
        "MES-AC-like" = "#F37F72",
        "MES-INT-like" = "#B02425",
        "MES-HYP-like" = "#F1634B",
        "Invasive-high OPC/NPC1" = "#2B7095"
    )

    # Plot with legend
    ggplot(curr_df, aes(x = mean_cell_type_gex, y = mean_cell_type_acc)) +
        geom_point(
            aes(size = total_n, color = "#000000", fill = x),
            shape = 21,
            stroke = 1
        ) + # Set stroke and shape here
        geom_text_repel(aes(label = x), size = 5) +
        scale_color_manual(values = color_palette) + # Use custom color palette
        scale_fill_manual(values = color_palette) + # Use custom color palette for fill
        labs(
            x = "Normalized mean gene expression",
            y = "Motif accessibility (mean z-score)"
        ) +
        GBM_theme() +
        theme(
            legend.position = "bottom",
            legend.title = element_blank()
        ) +
        ggtitle(TF_name)
    ggsave(
        paste0(
            plot_dir,
            "/Neuronal_motif_deviation_score_scatterplot_",
            TF_name,
            method,
            ".pdf"
        ),
        width = 4.5,
        height = 5.5
    )
}
