#!/usr/bin/env bash
#SBATCH -J 003_banksy.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=veryhimem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

sample_id="6509_A"
input_file="${PWD}/data/Xenium/processed/${sample_id}__preproc.rds"
output_dir="${PWD}/data/Xenium/processed"

# Parameters
lambda=0.3
k_geom=10
cluster_res=0.5
features_to_use="all"

echo "Activating conda environment..."
source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

Rscript ${PWD}/preprocessing/Xenium/003_banksy.R \
    --input_file ${input_file} \
    --output_dir ${output_dir} \
    --sample_id ${sample_id} \
    --lambda $lambda \
    --k_geom $k_geom \
    --cluster_res $cluster_res \
    --features_to_use $features_to_use
