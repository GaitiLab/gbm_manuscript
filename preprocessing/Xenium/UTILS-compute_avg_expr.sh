#!/usr/bin/env bash
#SBATCH -J UTILS-compute_avg_expr.R
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

source "${CONDA_PREFIX}/bin/activate" "cci"

# DO NOT CHANGE when trying to reproducing figure
prefix="6509_A__BANKSY"
assay="SCT"
groupbyvar="BANKSY_snn_res.0.5"
apply_is_confident_filter=0
is_seurat_v4=0

# Required for mean expression heatmap Figure
slot="scale.data"
Rscript "${PWD}/preprocessing/scRNAseq/UTILS-compute_avg_expr.R" \
    --input_file $input_file \
    --output_dir ${output_dir} \
    --slot ${slot} \
    --assay ${assay} \
    --groupbyvar ${groupbyvar} \
    --prefix ${prefix} \
    --apply_is_confident_filter ${apply_is_confident_filter} \
    --is_seurat_v4 ${is_seurat_v4}



