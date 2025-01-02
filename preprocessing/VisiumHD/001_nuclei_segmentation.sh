#!/usr/bin/env bash
#SBATCH -J launch-001_nuclei_segmentation.py
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# Models 
model_dir="${PWD}/misc/models"
model_name="python_2D_versatile_he"

# Directories
spatial_outs="/cluster/projects/gaitigroup/Data/GBM_VisiumHD/20240708_LH00244_0121_B22MN7GLT3_Gaiti_Yiyan_Visium"
image_dir="/cluster/projects/gaitigroup/Data/GBM_VisiumHD/Visium_Images"
output_dir="${PWD}/VisiumHD/processed/001_stardist_outs"

# Parameters for nuclei segmentation (defaults)
min_p=5
max_p=95

prob_thresh=0.01
nms_thresh=0.001

# For H4H array
sample_ids="${PWD}/misc/VisiumHD_sample_info.txt"
partition="himem"
cpus_per_task=1
mem="60G"
time="08:00:00"


job_min=1
job_max=$(wc -l < "${sample_ids}")

echo $job_max

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J 001_nuclei_segmentation.py
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
sample_id=\$(cut -d$'\t' -f2 ${sample_ids} | awk "NR==\${SLURM_ARRAY_TASK_ID}")

sample_dir="${spatial_outs}/\${sample_id}"
image_path="${image_dir}/\${image_id}.tif"

output_filename="\${sample_id}__nuclei_segmentation_stardist"

echo "Sample ID: \${sample_id}"
echo "Sample dir: \${sample_dir}"
echo "Image path: \${image_path}"

echo "prob_thresh: ${prob_thresh}"
echo "nms_thresh: ${nms_thresh}"

python3 "${PWD}/preprocessing/VisiumHD/001_nuclei_segmentation.py" \
    --output_dir ${output_dir} \
    --model_dir ${model_dir} \
    --model_name ${model_name} \
    --image_path \${image_path} \
    --output_filename \${output_filename} \
    --min_p ${min_p} \
    --max_p ${max_p} \
    --prob_thresh ${prob_thresh} \
    --nms_thresh ${nms_thresh}
EOF