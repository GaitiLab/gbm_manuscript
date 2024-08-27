# Code to reproduce Fig 2d

if(!("pacman" %in% rownames(installed.packages()))){
  install.packages("pacman")
}

pacman::p_load(data.table,
               ggplot2,
               ggrepel,
               Seurat,
               stringr,
               dplyr,
               ggpubr,
               foreach,
               doParallel,
               parallel)

# Enter paths here
output_dir <- ""

# Creating directories needed for outputs
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_dir <- paste0(output_dir, "/path_to_output_directory")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# --------------------------------- functions -------------------------------- #

find_neighbours_in_window <- function(adj_table, radius) {
  return (sapply(rownames(adj_table), find_neighbours_in_window_around_spot, adj_table = adj_table, radius = radius))
}

find_neighbours_in_window_around_spot <- function(adj_table, radius, spot_barcode) {

  neighbours <- c(find_neighbouring_spots(adj_table, spot_barcode))
  if (radius == 1) {
    return (neighbours)
  } else {
    neighbours <- c(neighbours, unlist(lapply(neighbours, find_neighbours_in_window_around_spot, adj_table = adj_table, radius = radius - 1)))
    return (unique(neighbours))
  }
}

find_neighbouring_spots <- function(adj_table, spot_barcode) {
  return (na.omit(adj_table[spot_barcode, ]))
}

# Adapted from Greenwald 2024
build_adjacency_table <- function(spot_pos_type) {
  adj_table <- sapply(spot_pos_type$barcodes, function(spot){
    spots_row = spot_pos_type[spot_pos_type$barcodes == spot, "array_row"]
    spots_col = spot_pos_type[spot_pos_type$barcodes == spot, "array_col"]

    # bottom left
    if (spots_col == 0 | spots_row == 0) {
      c1 = NA
    } else {
      curr_spot <- spot_pos_type[spot_pos_type$array_row == (spots_row - 1) & spot_pos_type$array_col == (spots_col - 1), ]
      spot_exists <- length(rownames(curr_spot)) == 1
      if (spot_exists) {
        n1 = spot_pos_type$barcodes[spot_pos_type$array_row == spots_row - 1 & spot_pos_type$array_col == spots_col - 1]
        c1 = n1
      } else {
        c1 = NA
      }
    }

    # bottom right
    if (spots_col == 127 | spots_row == 0) {
      c2 = NA
    } else {
      curr_spot <- spot_pos_type[spot_pos_type$array_row == (spots_row - 1) & spot_pos_type$array_col == (spots_col + 1), ]
      spot_exists <- length(rownames(curr_spot)) == 1
      if (spot_exists) {
        n2 = spot_pos_type$barcodes[spot_pos_type$array_row == spots_row - 1 & spot_pos_type$array_col == spots_col + 1]
        c2 = n2
      } else {
        c2 = NA
      }
    }
    
    # left
    if (spots_col == 0 | spots_col == 1) {
      c3 = NA
    } else {
      curr_spot <- spot_pos_type[spot_pos_type$array_row == (spots_row) & spot_pos_type$array_col == (spots_col - 2), ]
      spot_exists <- length(rownames(curr_spot)) == 1
      if (spot_exists) {
        n3 = spot_pos_type$barcodes[spot_pos_type$array_row == spots_row & spot_pos_type$array_col == spots_col - 2]
        c3 = n3
      } else {
        c3 = NA
      }
    }
    
    # right
    if (spots_col == 126 | spots_col == 127) {
      c4 = NA
    } else {
      curr_spot <- spot_pos_type[spot_pos_type$array_row == (spots_row) & spot_pos_type$array_col == (spots_col + 2), ]
      spot_exists <- length(rownames(curr_spot)) == 1
      if (spot_exists) {
        n4 = spot_pos_type$barcodes[spot_pos_type$array_row == spots_row & spot_pos_type$array_col == spots_col + 2]
        c4 = n4
      } else {
        c4 = NA
      }
    }
    
    # top left
    if (spots_col == 0 | spots_row == 77) {
      c5 = NA
    } else {
      curr_spot <- spot_pos_type[spot_pos_type$array_row == (spots_row + 1) & spot_pos_type$array_col == (spots_col - 1), ]
      spot_exists <- length(rownames(curr_spot)) == 1
      if (spot_exists) {
        n5 = spot_pos_type$barcodes[spot_pos_type$array_row == spots_row + 1 & spot_pos_type$array_col == spots_col - 1]
        c5 = n5
      } else {
        c5 = NA
      }
    }
    
    # top right
    if (spots_col == 127 | spots_row == 77) {
      c6 = NA
    } else {
      curr_spot <- spot_pos_type[spot_pos_type$array_row == (spots_row + 1) & spot_pos_type$array_col == (spots_col + 1), ]
      spot_exists <- length(rownames(curr_spot)) == 1
      if (spot_exists) {
        n6 = spot_pos_type$barcodes[spot_pos_type$array_row == spots_row + 1 & spot_pos_type$array_col == spots_col + 1]
        c6 = n6
      } else {
        c6 = NA
      }
    }
    
    
    return(c(c1,c2,c3,c4,c5,c6))
    
  })
  
  adj_table = t(adj_table)
  row.names(adj_table) = spot_pos_type$barcodes
  
  return(adj_table)
}

# Load DEGs

degs <- read.csv("/de_results.csv")
degs_up <- degs %>% 
  filter(log2FoldChange > 1 & padj < 0.05)
degs_dn <- degs %>% 
  filter(log2FoldChange < -1 & padj < 0.05)
degs_list <- list(degs_up$gene, degs_dn$gene)

# ------------------------- Load Greenwald 2024 data ------------------------- #
# Load Greenwald et al 2024 data
path_to_greenwald <- "/GBM_data/"
data_dirs <- list.dirs(path_to_greenwald, recursive = FALSE)

metadata <- read.csv("/visium_metadata.csv")

regions <- c("inf", "bulk", "T1", "nec")
greenwald <- sapply(data_dirs, function(dir) {
  outs_dir <- file.path(dir, "outs")

  curr_sample <- basename(dir)
  curr_patient <- substr(curr_sample, 1, 6)
  if (grepl("ZH881", curr_sample)) {
    curr_patient <- "ZH881"
  } else if (grepl("ZH916", curr_sample)) {
    curr_patient <- "ZH916"
  }

  so <- Load10X_Spatial(outs_dir)

  sample_metadata <- metadata %>% 
    filter(sample == curr_sample) 
  rownames(sample_metadata) <- sample_metadata$spot_id

  so <- AddMetaData(so, metadata = sample_metadata)

  # filter out cells with NA values as these were removed by preprocessing from paper
  cells_to_keep <- names(which(FALSE == is.na(so$sample)))
  so <- subset(so, cells = cells_to_keep)

  so$dataset <- "Greenwald"
  so$patient <- curr_patient
  so$region <- "NA"
  for (region in regions) {
    if (grepl(region, curr_sample)) {
      so$region <- region
    }
  }

  # Add spatial meta data
  pos_metadata <- read.csv(file.path(outs_dir, "spatial/tissue_positions_list.csv"), header = FALSE)
  pos_metadata <- pos_metadata %>% 
    select(all_of(c("V1", "V2", "V3", "V4"))) %>% 
    filter(V1 %in% colnames(so))
  colnames(pos_metadata) <- c("barcode", "in_tissue", "array_row", "array_col")
  rownames(pos_metadata) <- pos_metadata$barcode

  so <- AddMetaData(so, metadata = pos_metadata)

  so <- SCTransform(so, assay = "Spatial", verbose = FALSE)

  DefaultAssay(so) <- "SCT"
  
  so <- AddModuleScore(so, features = c(degs_list, hypoxia_markers), name = "gs", assay = "SCT")
  so$neuronal_sig <- so$gs1 - so$gs2
  so$hypoxia <- so$gs3

  print(head(so[[]]))

  # Dataset provides 4 classes of cna bins. Take the malignant and mixed-high as the ones to confidently label malignant spots
  malignant_cna_bin <- c("mix_high", "malignant")

  sample_metadata <- so@meta.data %>%
  mutate(mp_l2 = case_when(
    mp %in% c("OPC", "NPC") & neuronal_sig > 0 ~ "inv_up_stem",
    mp %in% c("OPC", "NPC") & neuronal_sig <= 0 ~ "inv_dn_stem",
    TRUE ~ mp
  )) %>%
  mutate(mp_l2 = case_when(
    mp_l2 == "inv_up_stem" & cna_bin %in% malignant_cna_bin ~ "Malignant_inv_up_stem",
    mp_l2 == "inv_up_stem" & !(cna_bin %in% malignant_cna_bin) ~ "NonMalignant_inv_up_stem",
    mp_l2 == "inv_dn_stem" & cna_bin %in% malignant_cna_bin ~ "Malignant_inv_dn_stem",
    mp_l2 == "inv_dn_stem" & !(cna_bin %in% malignant_cna_bin) ~ "NonMalignant_inv_dn_stem",
    mp_l2 == "AC" & cna_bin %in% malignant_cna_bin ~ "Malignant_AC",
    mp_l2 == "AC" & !(cna_bin %in% malignant_cna_bin) ~ "NonMalignant_AC",
    mp_l2 == "MES" & cna_bin %in% malignant_cna_bin ~ "Malignant_MES",
    mp_l2 == "MES" & !(cna_bin %in% malignant_cna_bin) ~ "NonMalignant_MES",
    mp_l2 == "MES_Hyp" & cna_bin %in% malignant_cna_bin ~ "Malignant_MES_Hyp",
    mp_l2 == "MES_Hyp" & !(cna_bin %in% malignant_cna_bin) ~ "NonMalignant_MES_Hyp",
    mp_l2 == "MES_Ast" & cna_bin %in% malignant_cna_bin ~ "Malignant_MES_Ast",
    mp_l2 == "MES_Ast" & !(cna_bin %in% malignant_cna_bin) ~ "NonMalignant_MES_Ast",
    mp_l2 == "chromatin.reg" & cna_bin %in% malignant_cna_bin ~ "Malignant_chromatin.reg",
    mp_l2 == "chromatin.reg" & !(cna_bin %in% malignant_cna_bin) ~ "NonMalignant_chromatin.reg",
    mp_l2 == "metabolism" & cna_bin %in% malignant_cna_bin ~ "Malignant_metabolism",
    mp_l2 == "metabolism" & !(cna_bin %in% malignant_cna_bin) ~ "NonMalignant_metabolism",
    TRUE ~ mp_l2
  ))
  so@meta.data <- metadata

  saveRDS(so, paste0("/", curr_sample, ".rds"))

  write.csv(so[[]], paste0("/", curr_sample, ".csv"))

  return (so)

})

print(greenwald)

# ---------------------------------------------------------------------------- #
#                     Examine co-localization with neurons                     #
# ---------------------------------------------------------------------------- #
# Here, we first define regions around neurons for a given radius. These regions will be used to bin spot locations into two categories 
# which we will then perform over-representation testing on. 

# Keep only samples with neurons for this analysis
greenwald_neuron <- na.omit(sapply(greenwald, function(sample) {
  if ("neuron" %in% unique(sample$mp)) {
    return (TRUE)
  } else {
    return (FALSE)
  }
}))

greenwald_neuron <- greenwald[greenwald_neuron]

print(greenwald_neuron)

# Find neighbours within radius

max_rad <- 5

numCores <- detectCores() - 1
registerDoParallel(cores = numCores)

greenwald_neuron <- foreach(sample = greenwald_neuron, .packages = c("dplyr")) %dopar% {
  spot_pos_type <- sample[[]] %>% select(all_of(c("array_row", "array_col", "mp")))
  spot_pos_type$barcodes <- rownames(spot_pos_type)
  
  adj_table <- build_adjacency_table(spot_pos_type)
  neuron_barcodes <- rownames(spot_pos_type)[spot_pos_type$mp == "neuron"]
  
  for (radius in seq(1, max_rad)) {
    neighbours_table <- find_neighbours_in_window(adj_table, radius)
    neuron_neighbours <- unique(unlist(neighbours_table[names(neighbours_table) %in% neuron_barcodes]))
    sample[[paste0("neuron_neighbour_", radius)]] <- ifelse(rownames(sample[[]]) %in% neuron_neighbours, TRUE, FALSE)
  }
  return(sample)
}

print(greenwald_neuron)

for (sample in greenwald_neuron) {
  p <- SpatialDimPlot(sample, group.by = "neuron_neighbour_3", image.alpha = 0)
  ggsave(filename = paste0(unique(sample$sample), "neuron_neighbour_3_dimplot.pdf"), path = plot_dir)
}

merged_greenwald_neuron <- merge(greenwald_neuron[[1]], greenwald_neuron[-1])

# Testing if relative abundance of inv-high opc/npc1 is greater within defined radii
plot_df <- merged_greenwald_neuron[[]] %>% 
  filter(mp %in% c("OPC", "NPC") & cna_bin %in% c("malignant", "mix_high")) 
write.csv(plot_df, file.path(plot_dir, "plot_df.csv"))

plot_df <- plot_df %>% filter(org1 == "dis")

dfs <- lapply(seq(1, max_rad), function(rad) {
  n_rad <- paste0("neuron_neighbour_", rad)
  df <- plot_df %>%
    group_by_at(c("sample", n_rad)) %>%
    mutate(total = n()) %>%
    ungroup() %>% 
    group_by_at(c("sample", n_rad, "mp_l2")) %>%
    summarise(frac = n() / dplyr::first(total), .groups = "drop") %>%
    select(sample, mp_l2, n_rad, frac)
  df$radius <- rad
  return (df)
})

plot_df <- rbindlist(dfs)
colnames(plot_df) <- c("sample","mp_l2", "in_window", "fraction", "radius")
plot_df$in_window <- factor(plot_df$in_window, levels = c(TRUE, FALSE))
plot_df$radius <- factor(plot_df$radius, levels = seq(1, max_rad))

plot_df <- plot_df %>% 
  filter(in_window == TRUE)

p <- ggplot(plot_df, aes(x = mp_l2, y = fraction, fill = mp_l2)) +
  geom_boxplot() +
  geom_point(position = "identity") +
  theme_classic() +
  facet_wrap(.~radius) +
  stat_compare_means(label = "p.format", method = "wilcox.test", ref.group = "inv_dn_stem")
ggsave(filename = "invasive_sig_neuron_colocal_rad_categorical.pdf", path = plot_dir)

# ---------------------------------------------------------------------------- #
#                          Assess over-representation                          #
# ---------------------------------------------------------------------------- #
n_perm <- 1000

calculate_proportions <- function(data, grouping_col, annotation_col) {
  df <- stats::aggregate(get(grouping_col) ~ get(annotation_col), data = data, FUN = function(x) mean(x))
  colnames(df) <- c(annotation_col, grouping_col)
  return(df)
}

ora <- function(samples, neighbour_res, annotation_col) {
  i <- 1
  sample_p_vals <- foreach(sample = samples, .packages = c("dplyr")) %dopar% {
    set.seed(i)
    i <- i+1
    metadata <- sample[[]]

    observed_proportions <- calculate_proportions(metadata, grouping_col = neighbour_res, annotation_col = "mp_l2")
    mps <- unique(metadata[[annotation_col]])

    permutation_results <- vector("list", length(mps))
    names(permutation_results) <- mps

    # Perform permutations
    for (i in seq(1, n_perm)) {
      shuffled_df <- metadata
      shuffled_df[[annotation_col]] <- sample(shuffled_df[[annotation_col]]) # Shuffle spot labels
      shuffled_proportions <- calculate_proportions(shuffled_df, grouping_col = neighbour_res, annotation_col = "mp_l2")
      
      for (mp in mps) {
        observed <- observed_proportions[[neighbour_res]][observed_proportions[[annotation_col]] == mp]
        shuffled <- shuffled_proportions[[neighbour_res]][shuffled_proportions[[annotation_col]] == mp]
        permutation_results[[mp]] <- c(permutation_results[[mp]], shuffled >= observed)
      }
    }

    # Calculate p-values
    p_values <- sapply(permutation_results, function(x) mean(x))

    sample_name <- unique(metadata$sample)
    result_df <- data.frame(p_val = c(p_values))
    rownames(result_df) <- names(p_values)
    result_df[[annotation_col]] <- names(p_values)

    # Calculate relative abundance of mps
    relative_abundance <- metadata %>% 
      group_by_at(c(annotation_col, neighbour_res)) %>% 
      summarise(mp_neighbour_total = n(), .groups = "drop") %>% 
      group_by(!!sym(annotation_col)) %>% 
      mutate(mp_total = sum(mp_neighbour_total)) %>% 
      filter(!!sym(neighbour_res) == TRUE) %>% 
      ungroup() %>% 
      mutate(total_neuron_prox = sum(mp_neighbour_total), relative_abundance = mp_neighbour_total / total_neuron_prox, 
        norm_relative_abundance = relative_abundance / mp_total) %>% 
      mutate(ordering = match(!!sym(annotation_col), rownames(result_df))) %>% 
      arrange(ordering) %>%
      select(!!sym(annotation_col), norm_relative_abundance)

    relative_abundance$norm_relative_abundance <- relative_abundance$norm_relative_abundance 
    relative_abundance$norm_relative_abundance <- scale(relative_abundance$norm_relative_abundance)
    
    result_df <- full_join(result_df, relative_abundance, by = annotation_col)
    
    write.csv(result_df, file.path(plot_dir, paste0(sample_name, "_ora_colocalization_", neighbour_res, ".csv")))

    return(result_df)
  }

  return(sample_p_vals)
}

r <- 5
results <- ora(greenwald_neuron, neighbour_res = paste0("neuron_neighbour_", r), annotation_col = "mp_l2")
results <- rbindlist(results)
results$padj <- p.adjust(results$p_val, method = "BH")
plot_df <- results %>% 
group_by(mp_l2) %>% 
mutate(signif = padj < 0.1) %>% 
summarise(num_signif = sum(signif), avg_relative_abund = mean(norm_relative_abundance, na.rm = TRUE))

print(plot_df)
write.csv(plot_df, file.path(plot_dir, paste0("neuron_neighbour_", r, "_results.csv")))

# ------------------------ Aggregate results and plot ------------------------ #

res <- list.files(plot_dir, full.names = TRUE)
res <- res[grep("neuron_neighbour_5.csv", res)]

sample_p_vals <- lapply(res, function(x) {
  result <- read.csv(x)
  result$norm_relative_abundance <- scale(result$norm_relative_abundance)                                 
  return (result)
})

results <- rbindlist(sample_p_vals)
results$padj <- p.adjust(results$p_val, method = "BH")
results[is.na(norm_relative_abundance), "norm_relative_abundance"] <- 0
plot_df <- results %>% 
  group_by(mp_l2) %>% 
  mutate(signif = padj < 0.1) %>% 
  summarise(num_signif = sum(signif), avg_relative_abund = mean(norm_relative_abundance, na.rm = TRUE))

print(plot_df)

pct_signif <- sapply(unique(results$mp_l2), function(x) {
  curr_mp_res <- results %>% 
    filter(mp_l2 == x)
  return(sum(curr_mp_res$padj < 0.1) / length(rownames(curr_mp_res)))
})
pct_signif <- as.data.frame(pct_signif)
pct_signif$mp_l2 <- rownames(pct_signif)
plot_df <- full_join(plot_df, pct_signif, by = "mp_l2")

plot_df <- plot_df %>% arrange(-avg_relative_abund) %>% mutate(rank = order(-avg_relative_abund))
plot_df$mp_l2 <- ifelse(plot_df$mp_l2 == "Malignant_inv_up_stem", "Invasive-high OPC/NPC", plot_df$mp_l2)
plot_df$mp_l2 <- ifelse(plot_df$mp_l2 == "Malignant_inv_dn_stem", "Invasive-low OPC/NPC", plot_df$mp_l2)
plot_df$mp_l2 <- ifelse(plot_df$mp_l2 == "neuron", "Neuron", plot_df$mp_l2)
plot_df$mp_l2 <- ifelse(plot_df$mp_l2 == "Malignant_MES_Hyp", "MES-Hypoxia", plot_df$mp_l2)
labels <- plot_df %>% filter(mp_l2 %in% c("Neuron", "Invasive-high OPC/NPC", "MES-Hypoxia", "Invasive-low OPC/NPC"))
plot_df$colour <- ifelse(grepl("^Malignant", plot_df$mp_l2), "Malignant", "Non-malignant")
plot_df$colour <- ifelse(plot_df$mp_l2 == "Invasive-high OPC/NPC", "Invasive-high OPC/NPC", plot_df$colour)
plot_df$colour <- ifelse(plot_df$mp_l2 == "Neuron", "Neuron", plot_df$colour)
plot_df$colour <- factor(plot_df$colour, levels = c("Malignant", "Non-malignant", "Neuron", "Invasive-high OPC/NPC"))
p <- ggplot(plot_df, aes(x = rank, y = avg_relative_abund)) +
  geom_point(aes(colour = colour, size = pct_signif)) +
  scale_size(name = "Significant in fraction of samples", breaks = seq(0, 1, 0.25), range = c(4,10)) +
  geom_point(aes(size = pct_signif, stroke = 2), shape = 21, colour = "black") +
  scale_shape_identity() +
  geom_label_repel(data = labels, aes(label = mp_l2), size = 6, force = 3, max.overlaps = 12, label.size = 0.75, nudge_x = 2, nudge_y = 0.25) + 
  scale_colour_manual(values = c("#C6D8FF","#ffe096", "green4", "orange")) +
  theme_classic() +
  ylab("Scaled Relative Abundance") +
  xlab("Rank") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size = 10)), size = guide_legend())
ggsave(filename = "colocal_rank_plot.pdf", path = plot_dir, width = 13, height = 8)