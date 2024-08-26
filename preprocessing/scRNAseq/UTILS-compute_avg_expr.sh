#!/usr/bin/env bash
#SBATCH -J UTILS-compute_avg_expr.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=veryhimem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out


# ---- Figure 2b-volcano.Rmd ---- #
# Compute average expression data using 'gbm_regional_study' per 'CCI_CellClass_L2_2' for each Sample
# Activating 'cci' environment from https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline
# TODO change this to your directory with the scripts (UTILS-compute_avg_expr.{R, sh})

# TODO to adapted by user
input_file="${PWD}/data/cci_scRNAseq/gbm_regional_study.rds"

# TODO: remove later only for testing
input_file="/cluster/projects/gaitigroup/Data/GBM/processed_data/gbm_regional_study.rds"

output_dir="${PWD}/output/cci_scRNAseq/processed"

source "${CONDA_PREFIX}/bin/activate" "cci"

# DO NOT CHANGE when trying to reproducing figure
slot="counts"
assay="RNA"
groupbyvar="Sample,CCI_CellClass_L2_2"
apply_is_confident_filter=1

Rscript "${PWD}/preprocessing/scRNAseq/UTILS-compute_avg_expr.R" \
    --input_file $input_file \
    --output_dir ${output_dir} \
    --slot ${slot} \
    --assay ${assay} \
    --groupbyvar ${groupbyvar} \
    --apply_is_confident_filter ${apply_is_confident_filter}

