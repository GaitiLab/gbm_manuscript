#!/usr/bin/env bash
#SBATCH -J launch-003_bin_to_poly.py
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

spatial_outs=/cluster/projects/gaitigroup/Data/GBM_VisiumHD/20240708_LH00244_0121_B22MN7GLT3_Gaiti_Yiyan_Visium
output_dir=${PWD}/data/VisiumHD/processed/003_bin_to_poly

# Expand centroids of bins - to create square bins with size radius_bin
radius_bin=1

# For H4H array
sample_ids="${PWD}/misc/VisiumHD_sample_info.txt"
partition="himem"
cpus_per_task=8
mem="40G"
time="02:00:00"

job_min=1
job_max=$(wc -l < "${sample_ids}")


sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J 003_bin_to_poly.py
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=${partition}
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${cpus_per_task}
#SBATCH --mem=${mem}
#SBATCH --time=${time}
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "stardist-env"

sample_id=\$(cut -d$'\t' -f1 ${sample_ids} | awk "NR==\${SLURM_ARRAY_TASK_ID}")
image_name=\$(cut -d$'\t' -f2 ${sample_ids} | awk "NR==\${SLURM_ARRAY_TASK_ID}")

sample_dir="${spatial_outs}/\${sample_id}/outs/binned_outputs/square_002um"
tissue_positions_path="\${sample_dir}/spatial/tissue_positions.parquet"
scale_factors_path="${spatial_outs}/\${sample_id}/outs/binned_outputs/square_002um/spatial/scalefactors_json.json"

python3 "${PWD}/preprocessing/VisiumHD/003_bin_to_poly.py" \
    --output_dir ${output_dir} \
    --sample_id \${sample_id} \
    --tissue_positions_path \${tissue_positions_path} \
    --radius_bin ${radius_bin} \
    --scale_factors_path \${scale_factors_path}

EOF