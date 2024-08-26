# Overview

Assumed working directory is `gbm_manuscript` . In addition, the SLURM workload manager is used for which the configuration has to be checked accordingly.

For each provided bash script check apart from the SLURM configuration, also the path to the conda environments.

Required files:

  + `misc/6509_A_roi_cell_ids.csv`: cells belonging to region of interest manually selected with the [Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer/latest) from 10x Genomics.
  + `misc/manual_annot.xlsx`: table with labels for clusters.

## Preprocessing

1. Run `sbatch preprocessing/Xenium/000_create_seurat_obj_raw.sh`
2. Run `sbatch preprocessing/Xenium/001_qc.sh`
3. Run `sbatch preprocessing/Xenium/002_preprocessing.sh`
4. Run `sbatch preprocessing/Xenium/003_banksy.sh`
5. Run `sbatch preprocessing/Xenium/004_DE_get_cluster_markers.sh`
6. Run `sbatch preprocessing/Xenium/005_annotate_cells.sh`, for this script, an additional file (the gene-count matrix) is required, you can generate this file with: `sbatch preprocessing/Xenium/UTILS-save_expression.sh`.

## Preprocessing for co-localization analyses

Co-localization assessed in two ways:

1. Cell type level
2. Cell type + interaction level

First, run `sbatch preprocessing/Xenium/110_determine_neighbors.sh` .

### Cell type co-localization

1. Run `sbatch preprocessing/Xenium/120_compute_neighbors_metrics_celltype.sh`
2. Run `sbatch preprocessing/Xenium/preprocessing/Xenium/121_combine_neighbors_metrics.R`

### Cell type and interaction co-localization

1. Run `sbatch preprocessing/Xenium/UTILS-prep_cell_type_x_interactions.sh`
2. Run `sbatch preprocessing/Xenium/120_compute_neighbors_metrics__celltype_x_gene_pair.sh`
3. Run `sbatch preprocessing/Xenium/121_combine_neighbors_metrics__celltype_x_gene_pair.sh`

## Preprocessing for mean expression heatmap

Run `sbatch preprocessing/Xenium/UTILS-compute_avg_expr.sh`
