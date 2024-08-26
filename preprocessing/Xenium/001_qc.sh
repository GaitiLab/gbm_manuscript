#!/usr/bin/env bash
#SBATCH -J 001_qc.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=00:20:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

sample_id="6509_A"
input_file="${PWD}/data/Xenium/processed/${sample_id}__raw.rds"
meta="${PWD}/data/Xenium/processed/${sample_id}__raw_metadata.rds"
output_dir="${PWD}/data/Xenium/processed"

echo "Activating conda environment..."
source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

Rscript ${PWD}/preprocessing/Xenium/001_qc.R \
    --input_file ${input_file} \
    --meta ${meta} \
    --output_dir ${output_dir} \
    --sample_id ${sample_id}
