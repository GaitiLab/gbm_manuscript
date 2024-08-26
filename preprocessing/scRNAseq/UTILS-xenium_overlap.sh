#!/usr/bin/env bash
#SBATCH -J UTILS-xenium_overlap.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# ---- Figure 2b-volcano.Rmd ---- #
# Compute average expression data using 'gbm_regional_study' per 'CCI_CellClass_L2_2' for each Sample
# Activating 'cci' environment from https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline
# TODO change this to your directory with the scripts (UTILS-compute_avg_expr.{R, sh})

# TODO to adapted by user
input_file="${PWD}/data/cci_scRNAseq/402c_filtering_aggregated_res.rds"
output_dir="${PWD}/data/cci_scRNAseq/processed"

source "${CONDA_PREFIX}/bin/activate" "cci"

# DO NOT CHANGE when trying to reproducing figure
interactions="${PWD}/data/cci_scRNAseq/processed/interactions_neuron_x_invasive_high_or_PL.xlsx"
xenium_panel="${PWD}/misc/gbm_xenium_414g_version1_aug2023.xlsx"
detected_interactions="${PWD}/data/cci_scRNAseq/interactions_summary.xlsx"

Rscript "${PWD}/preprocessing/scRNAseq/UTILS-xenium_overlap.R" \
    --output_dir ${output_dir} \
    --xenium_panel ${xenium_panel} \
    --interactions ${interactions} \
    --detected_interactions ${detected_interactions}
