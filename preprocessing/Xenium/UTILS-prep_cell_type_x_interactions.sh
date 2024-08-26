#!/usr/bin/env bash
#SBATCH -J UTILS-prep_cell_type_x_interaction.R
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

sample_id="6509_A"
input_file="${PWD}/data/Xenium/processed/${sample_id}__BANKSY__meta_annot_w_expr_counts__ROI.rds"
output_dir="${PWD}/data/Xenium/processed"
interactions="${PWD}/output/cci_scRNAseq/processed/overlap_w_xenium.xlsx"

source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

Rscript "${PWD}/preprocessing/Xenium/UTILS-prep_cell_type_x_interactions.R" \
    --interactions ${interactions} \
    --sample_id ${sample_id} \
    --meta_with_exp ${input_file} \
    --output_dir ${output_dir}

