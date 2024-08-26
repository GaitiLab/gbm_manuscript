#!/usr/bin/env bash
#SBATCH -J 004_DE_get_cluster_markers.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

sample_id="6509_A"
input_file="${PWD}/data/Xenium/processed/6509_A__BANKSY.rds"
output_dir="${PWD}/data/Xenium/processed"

# Please do not change when trying to regenerate the figures
cluster_label="BANKSY_snn_res.0.5"
assay="SCT"

echo "Activating conda environment..."
source "$CONDA_PREFIX/bin/activate" "r-4.3.3"

Rscript ${PWD}/preprocessing/Xenium/004_DE_get_cluster_markers.R \
    --input_file ${input_file} \
    --output_dir ${output_dir} \
    --cluster_label $cluster_label \
    --sample_id ${sample_id} \
    --assay ${assay}