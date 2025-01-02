#!/usr/bin/env bash
#SBATCH -J launch-R102_banksy.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

input_dir="${PWD}/data/VisiumHD/processed/seurat_objects/prepped"
output_dir="${PWD}/data/VisiumHD/processed/seurat_objects/prepped"

# Parameters
lambda=0.2
k_geom=15
cluster_res=1.0
features_to_use="variable"
assay="RNA"

# For H4H array
sample_ids="${PWD}/misc/VisiumHD_sample_info.txt"
partition="himem"
cpus_per_task=1
mem="40G"
time="08:00:00"

# Determine job array limits
# A. Determine number of files
# job_max=$(ls -d -- $sample_dir_all/* | wc -l) 2>/dev/null
# B. Number of lines in a file
job_max=$(wc -l < "${sample_ids}")
# job_max=1
job_min=1
echo $job_max


sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J R102_banksy.R
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

input_file=${input_dir}/\${sample_id}.rds

Rscript ${PWD}/preprocessing/VisiumHD/102_banksy.R \
    --input_file \${input_file} \
    --output_dir ${output_dir} \
    --sample_id \${sample_id} \
    --lambda $lambda \
    --k_geom $k_geom \
    --cluster_res $cluster_res \
    --features_to_use $features_to_use \
    --assay ${assay}
EOF