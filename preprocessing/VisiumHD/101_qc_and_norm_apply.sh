#!/usr/bin/env bash
#SBATCH -J launch-103_qc_apply.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

job_min=1

# Input/output
input_dir="${PWD}/data/VisiumHD/processed/seurat_objects/raw"
output_dir="${PWD}/data/VisiumHD/processed/seurat_objects/prepped"

# Params
assay="RNA"

# QC filters
min_counts=1
min_features=50
max_ratio_mt=.2

# For H4H array
sample_ids="${PWD}/misc/VisiumHD_sample_info.txt"
partition="himem"
cpus_per_task=1
mem="40G"
time="00:20:00"

job_max=$(wc -l < "${sample_ids}")

echo $job_max

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J 103_qc_apply.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${cpus_per_task}
#SBATCH --mem=${mem}
#SBATCH --partition=${partition}
#SBATCH --time=${time}
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "r-4.3.3"

sample_id=\$(cut -d$'\t' -f1 ${sample_ids} | awk "NR==\${SLURM_ARRAY_TASK_ID}")

input_file="${input_dir}/\${sample_id}.rds"

Rscript "$PWD/preprocessing/VisiumHD/101_qc_and_norm_apply.R" \
    --input_file \${input_file} \
    --output_dir ${output_dir} \
    --min_counts ${min_counts} \
    --min_features ${min_features} \
    --max_ratio_mt ${max_ratio_mt} \
    --assay ${assay} \
    --sample_id \${sample_id}
EOF
