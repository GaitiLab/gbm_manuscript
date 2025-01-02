#!/usr/bin/env bash
#SBATCH -J launch-P016_convert_tif_to_poly.py
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=all
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir="$base_dir/Joan/gbm_visiumHD"

segmentation_masks_dir="${work_dir}/output/stardist_segmentation/segmentation_masks/"

# For H4H array
sample_ids="${work_dir}/misc/sample_info.txt"

partition="himem"
cpus_per_task=16
mem="40G"
time="00:30:00"

job_min=1
job_max=$(wc -l < "${sample_ids}")

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J P016_convert_tif_to_poly.py
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

img_path="${segmentation_masks_dir}/\${sample_id}/\${sample_id}__nuclei_segmentation_stardist_masks__expanded.tif"
output_dir="${segmentation_masks_dir}/\${sample_id}"

python3 "Python/002b_convert_tif_to_poly.py" \
    --output_dir \${output_dir}\
    --img_path \${img_path} \
    --n_cores \${SLURM_CPUS_PER_TASK}
EOF
