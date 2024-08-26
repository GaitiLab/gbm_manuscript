#!/usr/bin/env bash
#SBATCH -J UTILS-extract_interactions.R
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

source "${HOME}/miniforge3/bin/activate" "cci"

# DO NOT CHANGE when trying to reproducing figure
condition_varname="Region"
pval_type="pval_adj"
condition_oi="PT"
alpha=0.05
output_name="interactions_neuron_x_invasive_high_or_PL"

Rscript "${PWD}/preprocessing/scRNAseq/UTILS-extract_interactions.R" \
    --input_file $input_file \
    --output_dir ${output_dir} \
    --condition_varname ${condition_varname} \
    --pval_type ${pval_type} \
    --condition_oi ${condition_oi} \
    --alpha ${alpha} \
    --output_name ${output_name}

