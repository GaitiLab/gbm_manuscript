# ---------------------------------------------------------------------------- #
#                                Figure 3c S6c                                 #
# ---------------------------------------------------------------------------- #
pacman::p_load(
    argparse,
    varhandle,
    log4r,
    ggplot2,
    ggrepel,
    tidyr,
    dplyr,
    stringr,
    data.table,
    tibble,
    Seurat,
    ggpubr,
    GaitiLabUtils,
    GBMutils
)
# Required inputs

params <- list(
    # Follow SCENIC+ in silico perturbation workflow to generate the input
    input = "", # SCENIC+ in silico perturbation results
    output_dir = "",

    # Path to Seurat object generated using this manuscript's data, raw data and final metadata can be downloaded online, see publication
    merged_obj = "", # Seurat object with original data

    include_regions = "all_regions",
    DE_genes = "inv_sig.csv", # invasive signature genes from Table S2
    TF_to_perturb = "ZEB1"
)

perturbation_dir <- paste0(
    params$output_dir,
    "/perturbation_simulation/",
    params$include_regions
)
curr_seurat_obj <- readRDS(paste0(
    perturbation_dir,
    "/curr_seurat_obj_",
    params$TF_to_perturb,
    ".rds"
))

update_annotation_labels <- function(data, annotation_col) {
    data <- data %>%
        mutate(
            !!annotation_col := case_when(
                !!annotation_col == "Invasive-high OPC/NPC1" ~
                    "Invasive-high OPC/NPC1",
                !!annotation_col == "Malignant_OPC" ~ "Malignant OPC",
                !!annotation_col == "Malignant_NPC1" ~ "Malignant NPC1",
                !!annotation_col == "Malignant_NPC2" ~ "Malignant NPC2",
                !!annotation_col == "Malignant_AC" ~ "Malignant AC",
                !!annotation_col == "Malignant_MES_AST" ~ "Malignant MES-AST",
                !!annotation_col == "Malignant_MES_HYP" ~ "Malignant MES-HYP",
                !!annotation_col == "Malignant_MES_INT" ~ "Malignant MES-INT",
                !!annotation_col == "Progenitor_like" ~ "Progenitor-like",
                !!annotation_col == "Differentiated_like" ~
                    "Differentiated-like",
                TRUE ~ as.character(!!annotation_col) # Retain original value if no match
            )
        )
    return(data)
}

# ---------------------------------------------------------------------------- #
#                                Figure 3c S6c                                 #
# ---------------------------------------------------------------------------- #
plot_invasion_probabilities_per_sample <- function(
    curr_seurat_obj,
    TF_to_perturb,
    perturbation_dir,
    annotation_column) {
    # Ensure the annotation column is a valid column in the Seurat object
    if (!annotation_column %in% colnames(curr_seurat_obj[[]])) {
        stop(
            "The provided 'annotation_column' does not exist in the Seurat object."
        )
    } else {
        annotation_col <- sym(annotation_column)
    }

    # Prepare the data with transition information
    all_cell <- curr_seurat_obj[[]] %>%
        dplyr::select(original_inv, perturbed_inv, !!annotation_col, Sample) %>%
        mutate(
            transition = case_when(
                original_inv < perturbed_inv ~ "Increased invasiveness",
                original_inv > perturbed_inv ~ "Decreased invasiveness",
            )
        )

    # Calculate the transition probability for each cell type
    transition_prob <- all_cell %>%
        group_by(!!annotation_col, transition, Sample) %>%
        summarize(
            count = n(),
            .groups = "drop"
        )

    # Pivot the data to have one row per cell type with columns for increased and decreased invasiveness counts
    transition_prob_per_cell_type <- transition_prob %>%
        pivot_wider(
            names_from = transition,
            values_from = count,
            values_fill = list(count = 0)
        ) %>%
        mutate(
            increased_prob = `Increased invasiveness` /
                (`Increased invasiveness` + `Decreased invasiveness`),
            decreased_prob = `Decreased invasiveness` /
                (`Increased invasiveness` + `Decreased invasiveness`)
        )

    # Calculate the p-value for each cell type
    binomial_results <- transition_prob_per_cell_type %>%
        mutate(
            success = ifelse(
                `Increased invasiveness` > `Decreased invasiveness`,
                `Increased invasiveness`,
                `Decreased invasiveness`
            )
        ) %>%
        group_by(!!annotation_col) %>%
        summarize(
            observed_successes = sum(success == `Increased invasiveness`),
            total_trials = n(),
            p_value = binom.test(
                observed_successes,
                total_trials,
                p = 0.5,
                alternative = "two.sided"
            )$p.value,
            .groups = "drop"
        ) %>%
        mutate(p_value = formatC(p_value, format = "e", digits = 2))
    binomial_results <- update_annotation_labels(
        binomial_results,
        annotation_col
    )

    # Transform the data to long format for plotting
    long_data <- transition_prob_per_cell_type %>%
        dplyr::select(
            !!annotation_col,
            Sample,
            increased_prob,
            decreased_prob
        ) %>%
        pivot_longer(
            cols = c(increased_prob, decreased_prob),
            names_to = "Transition_Type",
            values_to = "Probability"
        ) %>%
        mutate(
            Transition_Type = case_when(
                Transition_Type == "increased_prob" ~ "Increased invasiveness",
                Transition_Type == "decreased_prob" ~ "Decreased invasiveness"
            )
        )
    long_data <- update_annotation_labels(long_data, annotation_col)

    # Calculate mean and standard error of the mean (SEM) for each cell type and transition type
    summary_data <- long_data %>%
        group_by(!!annotation_col, Transition_Type) %>%
        summarize(
            mean_prob = mean(Probability),
            sem_prob = sd(Probability) / sqrt(n()),
            .groups = "drop"
        )
    write.csv(
        summary_data,
        file.path(
            perturbation_dir,
            paste0(
                "change_in_inv_per_sample_mean_",
                TF_to_perturb,
                "_",
                annotation_column,
                ".csv"
            )
        )
    )

    # Rank the cell types based on the decreased invasiveness probability
    rank_by_diff_prob <- summary_data %>%
        filter(Transition_Type == "Decreased invasiveness") %>%
        arrange(desc(mean_prob)) %>%
        pull(!!annotation_col)
    summary_data[[annotation_column]] <- factor(
        summary_data[[annotation_column]],
        levels = rank_by_diff_prob
    )

    # invasiveness color palette
    color_palette <- c(
        "Increased invasiveness" = "#2B7095",
        "Decreased invasiveness" = "#7C9EB5"
    )

    # Plotting
    p <- ggplot(
        summary_data,
        aes(x = !!annotation_col, y = mean_prob, fill = Transition_Type)
    ) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        geom_errorbar(
            aes(ymin = mean_prob - sem_prob, ymax = mean_prob + sem_prob),
            position = position_dodge(width = 0.7),
            width = 0.1,
            color = "black",
            size = 0.5
        ) +
        geom_text(
            data = binomial_results,
            aes(
                x = !!annotation_col,
                y = max(summary_data$mean_prob + summary_data$sem_prob) + 0.05,
                label = p_value
            ),
            inherit.aes = FALSE,
            size = 4,
            color = "black"
        ) +
        scale_fill_manual(values = color_palette) +
        scale_y_continuous(labels = scales::percent_format(scale = 100)) +
        labs(
            title = paste0(
                "Change in invasiveness across cell states upon ",
                TF_to_perturb,
                " KO"
            ),
            x = "Cell state",
            y = "Change in invasiveness"
        ) +
        GBM_theme()

    # Save the plot
    ggsave(
        filename = paste0(
            "change_in_inv_per_sample_",
            TF_to_perturb,
            "_",
            annotation_column,
            ".pdf"
        ),
        plot = p,
        path = perturbation_dir,
        width = 1.5 * length(unique(summary_data[[annotation_column]])) + 3,
        height = 4
    )
}


plot_invasion_difference <- function(
    curr_seurat_obj,
    TF_to_perturb,
    perturbation_dir,
    annotation_column) {
    # Ensure the annotation column is a valid column in the Seurat object
    if (!annotation_column %in% colnames(curr_seurat_obj[[]])) {
        stop(
            "The provided 'annotation_column' does not exist in the Seurat object."
        )
    } else {
        annotation_col <- sym(annotation_column)
    }

    # Prepare the data with transition information
    all_cell <- curr_seurat_obj[[]] %>%
        dplyr::select(original_inv, perturbed_inv, !!annotation_col, Sample) %>%
        mutate(
            transition = case_when(
                original_inv < perturbed_inv ~ "Increased invasiveness",
                original_inv > perturbed_inv ~ "Decreased invasiveness",
            )
        )

    # Calculate the transition probability for each cell type
    transition_prob <- all_cell %>%
        group_by(!!annotation_col, transition, Sample) %>%
        summarize(
            count = n(),
            .groups = "drop"
        )

    # Calculate the transition ratio for each cell type
    transition_counts <- transition_prob %>%
        pivot_wider(
            names_from = transition,
            values_from = count,
            values_fill = 0
        ) %>%
        dplyr::rename(
            Increased = `Increased invasiveness`,
            Decreased = `Decreased invasiveness`
        )
    transition_diffs <- transition_counts %>%
        mutate(Total_n = Increased + Decreased) %>%
        mutate(prop_differece = (Decreased - Increased) / Total_n) %>%
        filter(Total_n >= 20)

    # Calculate the mean and standard error of the mean (SEM) for each cell type
    summary_data <- transition_diffs %>%
        group_by(!!annotation_col) %>%
        summarize(
            mean_difference = mean(prop_differece),
            sem_difference = sd(prop_differece) / sqrt(n()),
            .groups = "drop"
        ) %>%
        arrange(mean_difference)
    summary_data <- update_annotation_labels(summary_data, annotation_col)
    summary_data[[annotation_column]] <- factor(
        summary_data[[annotation_column]],
        levels = summary_data[[annotation_column]]
    )

    ggplot(summary_data, aes(x = !!sym(annotation_col), y = mean_difference)) +
        geom_point(size = 3) + # Points for mean ratios
        geom_errorbar(
            aes(
                ymin = mean_difference - sem_difference,
                ymax = mean_difference + sem_difference
            ),
            width = 0.2
        ) +
        geom_hline(yintercept = 0, linetype = "dashed") + # Reference line difference = 0
        coord_flip() + # Flip axes for better readability
        labs(
            title = "Forest Plot of invasiveness difference",
            x = "Cell state",
            y = "Mean proportion difference between cells with decreased and increased invasiveness"
        ) +
        scale_y_continuous(
            # Adjust y-axis limits and make it symmetric
            limits = c(-1, 1),
            expand = c(0, 0),
            labels = scales::percent_format(accuracy = 1)
        ) +
        GBM_theme()
    ggsave(
        paste0(
            "invasion_difference_",
            TF_to_perturb,
            "_",
            annotation_column,
            ".pdf"
        ),
        path = perturbation_dir,
        width = 8,
        height = 0.5 * length(unique(summary_data[[annotation_column]])) + 2
    )
}

# Plot invasion probabilities per sample
plot_invasion_probabilities_per_sample(
    curr_seurat_obj = curr_seurat_obj,
    TF_to_perturb = params$TF_to_perturb,
    perturbation_dir = perturbation_dir,
    annotation_column = "CCI_CellClass_L2_2"
)

plot_invasion_difference(
    curr_seurat_obj = curr_seurat_obj,
    TF_to_perturb = params$TF_to_perturb,
    perturbation_dir = perturbation_dir,
    annotation_column = "CCI_CellClass_L2_2"
)

# ---------------------------------------------------------------------------- #
#                               Figure 4e 6h                                  #
# ---------------------------------------------------------------------------- #
plot_transition_probabilities_per_sample <- function(
    curr_seurat_obj,
    TF_to_perturb,
    perturbation_dir,
    annotation_column) {
    # Ensure the annotation column is a valid column in the Seurat object
    if (!annotation_column %in% colnames(curr_seurat_obj[[]])) {
        stop(
            "The provided 'annotation_column' does not exist in the Seurat object."
        )
    } else {
        annotation_col <- sym(annotation_column)
    }

    # Prepare the data with transition information
    all_cell <- curr_seurat_obj[[]] %>%
        dplyr::select(original_opc, perturbed_opc, !!annotation_col, Sample) %>%
        mutate(
            transition = case_when(
                original_opc > perturbed_opc ~ "Differentiation",
                original_opc < perturbed_opc ~ "De-differentiation"
            )
        )

    # Calculate the transition probability for each cell type
    transition_prob <- all_cell %>%
        group_by(!!annotation_col, transition, Sample) %>%
        summarize(
            count = n(),
            .groups = "drop"
        )

    # Pivot the data to have one row per cell type with columns for Differentiation and De-differentiation counts
    transition_prob_per_cell_type <- transition_prob %>%
        pivot_wider(
            names_from = transition,
            values_from = count,
            values_fill = list(count = 0)
        ) %>%
        mutate(
            diff_prob = Differentiation /
                (Differentiation + `De-differentiation`),
            dediff_prob = `De-differentiation` /
                (Differentiation + `De-differentiation`)
        ) %>%
        filter(Differentiation + `De-differentiation` >= 20) # Exclude samples with less than 20 cells

    # Calculate the p-value for each cell type
    binomial_results <- transition_prob_per_cell_type %>%
        mutate(
            success = ifelse(
                Differentiation > `De-differentiation`,
                Differentiation,
                `De-differentiation`
            )
        ) %>%
        group_by(!!annotation_col) %>%
        summarize(
            observed_successes = sum(success == Differentiation),
            total_trials = n(),
            p_value = binom.test(
                observed_successes,
                total_trials,
                p = 0.5,
                alternative = "two.sided"
            )$p.value,
            .groups = "drop"
        ) %>%
        mutate(p_value = formatC(p_value, format = "e", digits = 2))
    binomial_results <- update_annotation_labels(
        binomial_results,
        annotation_col
    )

    # Transform the data to long format for plotting
    long_data <- transition_prob_per_cell_type %>%
        dplyr::select(!!annotation_col, Sample, diff_prob, dediff_prob) %>%
        pivot_longer(
            cols = c(diff_prob, dediff_prob),
            names_to = "Transition_Type",
            values_to = "Probability"
        ) %>%
        mutate(
            Transition_Type = case_when(
                Transition_Type == "diff_prob" ~ "Differentiation",
                Transition_Type == "dediff_prob" ~ "De-differentiation"
            )
        )
    long_data <- update_annotation_labels(long_data, annotation_col)

    # Calculate mean and standard error of the mean (SEM) for each cell type and transition type
    summary_data <- long_data %>%
        group_by(!!annotation_col, Transition_Type) %>%
        summarize(
            mean_prob = mean(Probability),
            sem_prob = sd(Probability) / sqrt(n()),
            .groups = "drop"
        )
    write.csv(
        summary_data,
        file.path(
            perturbation_dir,
            paste0(
                "transition_probability_per_sample_mean_",
                TF_to_perturb,
                "_",
                annotation_column,
                ".csv"
            )
        )
    )

    # Rank the cell types based on the differentiation probability
    rank_by_diff_prob <- summary_data %>%
        filter(Transition_Type == "Differentiation") %>%
        arrange(desc(mean_prob)) %>%
        pull(!!annotation_col)
    summary_data[[annotation_column]] <- factor(
        summary_data[[annotation_column]],
        levels = rank_by_diff_prob
    )

    # OPC_dev color palette
    color_palette <- c(
        "Differentiation" = "#F0E9BA",
        "De-differentiation" = "#ED8B22"
    )

    # Plotting
    p <- ggplot(
        summary_data,
        aes(x = !!annotation_col, y = mean_prob, fill = Transition_Type)
    ) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        geom_errorbar(
            aes(ymin = mean_prob - sem_prob, ymax = mean_prob + sem_prob),
            position = position_dodge(width = 0.7),
            width = 0.1,
            color = "black",
            size = 0.5
        ) +
        geom_text(
            data = binomial_results,
            aes(
                x = !!annotation_col,
                y = max(summary_data$mean_prob + summary_data$sem_prob) + 0.05,
                label = p_value
            ),
            inherit.aes = FALSE,
            size = 4,
            color = "black"
        ) +
        scale_fill_manual(values = color_palette) +
        scale_y_continuous(labels = scales::percent_format(scale = 100)) +
        labs(
            title = paste0(
                "Differentiation probability across cell states upon ",
                TF_to_perturb,
                " KO"
            ),
            x = "Cell state",
            y = "Transition probability"
        ) +
        GBM_theme()

    # Save the plot
    ggsave(
        filename = paste0(
            "transition_probability_per_sample_",
            TF_to_perturb,
            "_",
            annotation_column,
            ".pdf"
        ),
        plot = p,
        width = 1.5 * length(unique(summary_data[[annotation_column]])) + 3,
        height = 4,
        path = perturbation_dir
    )
}

plot_transition_difference <- function(
    curr_seurat_obj,
    TF_to_perturb,
    perturbation_dir,
    annotation_column) {
    # Ensure the annotation column is a valid column in the Seurat object
    if (!annotation_column %in% colnames(curr_seurat_obj[[]])) {
        stop(
            "The provided 'annotation_column' does not exist in the Seurat object."
        )
    } else {
        annotation_col <- sym(annotation_column)
    }

    # Prepare the data with transition information
    all_cell <- curr_seurat_obj[[]] %>%
        dplyr::select(original_opc, perturbed_opc, !!annotation_col, Sample) %>%
        mutate(
            transition = case_when(
                original_opc > perturbed_opc ~ "Differentiation",
                original_opc < perturbed_opc ~ "De-differentiation"
            )
        )

    # Calculate the transition probability for each cell type
    transition_prob <- all_cell %>%
        group_by(!!annotation_col, transition, Sample) %>%
        summarize(
            count = n(),
            .groups = "drop"
        )

    # Calculate the transition ratio for each cell type
    transition_counts <- transition_prob %>%
        pivot_wider(
            names_from = transition,
            values_from = count,
            values_fill = 0
        ) %>%
        dplyr::rename(
            Differentiation = Differentiation,
            De_differentiation = `De-differentiation`
        )
    transition_diffs <- transition_counts %>%
        mutate(Total_n = Differentiation + De_differentiation) %>%
        mutate(
            prop_differece = (Differentiation - De_differentiation) / Total_n
        ) %>%
        filter(Total_n >= 20)

    # Calculate the mean and standard error of the mean (SEM) for each cell type
    summary_data <- transition_diffs %>%
        group_by(!!annotation_col) %>%
        summarize(
            mean_difference = mean(prop_differece),
            sem_difference = sd(prop_differece) / sqrt(n()),
            .groups = "drop"
        ) %>%
        arrange(mean_difference)
    summary_data <- update_annotation_labels(summary_data, annotation_col)
    summary_data[[annotation_column]] <- factor(
        summary_data[[annotation_column]],
        levels = summary_data[[annotation_column]]
    )

    ggplot(summary_data, aes(x = !!sym(annotation_col), y = mean_difference)) +
        geom_point(size = 3) + # Points for mean ratios
        geom_errorbar(
            aes(
                ymin = mean_difference - sem_difference,
                ymax = mean_difference + sem_difference
            ),
            width = 0.2
        ) +
        geom_hline(yintercept = 0, linetype = "dashed") + # Reference line difference = 0
        coord_flip() + # Flip axes for better readability
        labs(
            title = "Forest Plot of transition difference",
            x = "Cell state",
            y = "Mean proportion difference between cells that differentiate and de-differentiate"
        ) +
        scale_y_continuous(
            # Adjust y-axis limits and make it symmetric
            limits = c(-1, 1),
            expand = c(0, 0),
            labels = scales::percent_format(accuracy = 1)
        ) +
        GBM_theme()
    ggsave(
        filename = paste0(
            "transition_difference_",
            TF_to_perturb,
            "_",
            annotation_column,
            ".pdf"
        ),
        width = 8,
        path = perturbation_dir,
        height = 0.5 * length(unique(summary_data[[annotation_column]])) + 2
    )
}

# Plot transition probabilities per sample
plot_transition_probabilities_per_sample(
    curr_seurat_obj = curr_seurat_obj,
    TF_to_perturb = params$TF_to_perturb,
    perturbation_dir = perturbation_dir,
    annotation_column = "CCI_CellClass_L2_2"
)

plot_transition_difference(
    curr_seurat_obj = curr_seurat_obj,
    TF_to_perturb = params$TF_to_perturb,
    perturbation_dir = perturbation_dir,
    annotation_column = "CCI_CellClass_L2_2"
)

# ---------------------------------------------------------------------------- #
#                                Figure S6kl                                   #
# ---------------------------------------------------------------------------- #
plot_signature_change_probabilities_per_sample <- function(
    curr_seurat_obj,
    TF_to_perturb,
    signature,
    perturbation_dir,
    annotation_column) {
    # Ensure the annotation column is a valid column in the Seurat object
    if (!annotation_column %in% colnames(curr_seurat_obj[[]])) {
        stop(
            "The provided 'annotation_column' does not exist in the Seurat object."
        )
    } else {
        annotation_col <- sym(annotation_column)
    }

    # Get the original and perturbed signature scores
    original_sig <- sym(paste0(signature, "1"))
    perturbed_sig <- sym(paste0(signature, "_perturbed1"))

    # get title
    title_dict <- c(
        "synaptic_signaling" = "GO:BP Synaptic signaling",
        "synapse" = "GO:CC Synapse",
        "postsynapse" = "GO:CC Postsynapse"
    )
    title <- title_dict[signature]

    # Prepare the data with transition information
    all_cell <- curr_seurat_obj[[]] %>%
        dplyr::select(
            !!original_sig,
            !!perturbed_sig,
            !!annotation_col,
            Sample
        ) %>%
        mutate(
            transition = case_when(
                !!original_sig < !!perturbed_sig ~ "Increased",
                !!original_sig > !!perturbed_sig ~ "Decreased",
            )
        )

    # Calculate the transition probability for each cell type
    transition_prob <- all_cell %>%
        group_by(!!annotation_col, transition, Sample) %>%
        summarize(
            count = n(),
            .groups = "drop"
        )

    # Pivot the data to have one row per cell type with columns for Increased and Decreased counts
    transition_prob_per_cell_type <- transition_prob %>%
        pivot_wider(
            names_from = transition,
            values_from = count,
            values_fill = list(count = 0)
        ) %>%
        mutate(
            increased_prob = Increased / (Increased + Decreased),
            decreased_prob = Decreased / (Increased + Decreased)
        ) %>%
        filter(Increased + Decreased >= 20) # Exclude samples with less than 20 cells

    # Calculate the p-value for each cell type
    binomial_results <- transition_prob_per_cell_type %>%
        mutate(
            success = ifelse(Increased > Decreased, Increased, Decreased)
        ) %>%
        group_by(!!annotation_col) %>%
        summarize(
            observed_successes = sum(success == Increased),
            total_trials = n(),
            p_value = binom.test(
                observed_successes,
                total_trials,
                p = 0.5,
                alternative = "two.sided"
            )$p.value,
            .groups = "drop"
        ) %>%
        mutate(p_value = formatC(p_value, format = "e", digits = 2))
    binomial_results <- update_annotation_labels(
        binomial_results,
        annotation_col
    )

    # Transform the data to long format for plotting
    long_data <- transition_prob_per_cell_type %>%
        dplyr::select(
            !!annotation_col,
            Sample,
            increased_prob,
            decreased_prob
        ) %>%
        pivot_longer(
            cols = c(increased_prob, decreased_prob),
            names_to = "Transition_Type",
            values_to = "Probability"
        ) %>%
        mutate(
            Transition_Type = case_when(
                Transition_Type == "increased_prob" ~ "Increased",
                Transition_Type == "decreased_prob" ~ "Decreased"
            )
        )
    long_data <- update_annotation_labels(long_data, annotation_col)

    # Calculate mean and standard error of the mean (SEM) for each cell type and transition type
    summary_data <- long_data %>%
        group_by(!!annotation_col, Transition_Type) %>%
        summarize(
            mean_prob = mean(Probability),
            sem_prob = sd(Probability) / sqrt(n()),
            .groups = "drop"
        )
    write.csv(
        summary_data,
        file.path(
            perturbation_dir,
            paste0(
                "signature_transition_probability_per_sample_mean_",
                TF_to_perturb,
                "_",
                annotation_column,
                "_",
                signature,
                ".csv"
            )
        )
    )

    # Rank the cell types based on the differentiation probability
    rank_by_diff_prob <- summary_data %>%
        filter(Transition_Type == "Increased") %>%
        arrange(mean_prob) %>%
        pull(!!annotation_col)
    summary_data[[annotation_column]] <- factor(
        summary_data[[annotation_column]],
        levels = rank_by_diff_prob
    )

    # Differentiation color palette
    color_palette <- c("Increased" = "#B2182B", "Decreased" = "#2166AC")

    # Plotting - stacked bar plot
    summary_data <- summary_data %>%
        mutate(
            sign_mean_prob = case_when(
                Transition_Type == "Increased" ~ mean_prob,
                Transition_Type == "Decreased" ~ -mean_prob
            )
        )
    p <- ggplot(
        summary_data,
        aes(x = !!annotation_col, y = sign_mean_prob, fill = Transition_Type)
    ) +
        geom_bar(stat = "identity", position = "stack", width = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_errorbar(
            aes(
                ymin = sign_mean_prob - sem_prob,
                ymax = sign_mean_prob + sem_prob
            ),
            width = 0.2,
            color = "black",
            size = 0.5
        ) +
        geom_text(
            data = binomial_results,
            aes(
                x = !!annotation_col,
                y = max(
                    max(summary_data$sign_mean_prob + summary_data$sem_prob) +
                        0.05,
                    0.9
                ),
                label = p_value
            ),
            inherit.aes = FALSE,
            size = 3,
            color = "black"
        ) +
        scale_fill_manual(values = color_palette) +
        labs(title = title, x = "Cell state", y = "Change in signature") +
        scale_y_continuous(
            labels = scales::percent_format(scale = 100),
            limits = c(-1, 1),
            expand = c(0, 0)
        ) +
        GBM_theme() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()
        )
    ggsave(
        path = perturbation_dir,
        filename = paste0(
            "signature_transition_probability_per_sample_stacked_",
            TF_to_perturb,
            "_",
            annotation_column,
            "_",
            signature,
            ".pdf"
        ),
        plot = p,
        width = 1 * length(unique(summary_data[[annotation_column]])),
        height = 6
    )
}

# Plot signature change probabilities per sample
plot_signature_change_probabilities_per_sample(
    curr_seurat_obj = curr_seurat_obj,
    TF_to_perturb = params$TF_to_perturb,
    perturbation_dir = perturbation_dir,
    signature = "synaptic_signaling",
    annotation_column = "CCI_CellClass_L2_2"
)

plot_signature_change_probabilities_per_sample(
    curr_seurat_obj = curr_seurat_obj,
    TF_to_perturb = params$TF_to_perturb,
    perturbation_dir = perturbation_dir,
    signature = "synapse",
    annotation_column = "CCI_CellClass_L2_2"
)
