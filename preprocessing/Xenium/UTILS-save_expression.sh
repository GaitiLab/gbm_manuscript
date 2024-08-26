#!/usr/bin/env bash
#SBATCH -J UTILS-save_expression.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# TODO to adapted by user
input_file="${PWD}/data/Xenium/processed/6509_A__BANKSY.rds"
output_dir="${PWD}/data/Xenium/processed"

source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

# DO NOT CHANGE when trying to reproducing figure
assay="RNA"
layer="counts"
Rscript "${PWD}/preprocessing/Xenium/UTILS-save_expression.R" \
    --input_file $input_file \
    --output_dir ${output_dir} \
    --layer ${layer} \
    --assay ${assay}