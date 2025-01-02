#!/usr/bin/env bash
#SBATCH -J launch-004_map_barcodes.py
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --partition=himem
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

spatial_outs="/cluster/projects/gaitigroup/Data/GBM_VisiumHD/20240708_LH00244_0121_B22MN7GLT3_Gaiti_Yiyan_Visium"
radius_bin=1

segm_dir="${PWD}/data/VisiumHD/processed/002_nuclei_expansion"
bin_dir="${PWD}/data/VisiumHD/processed/003_bin_to_poly"

output_dir="${PWD}/data/VisiumHD/processed/004_map_barcodes"

# For H4H array
sample_ids="${PWD}/misc/VisiumHD_sample_info.txt"
partition="himem"
cpus_per_task=1
mem="40G"
time="12:00:00"

job_min=1
job_max=$(wc -l < "${sample_ids}")

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J 004_map_barcodes.py
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

sample_dir=${spatial_outs}/\${sample_id}/outs/binned_outputs/square_002um

# Polygons from bins
gdf_barcodes_path="${bin_dir}/\${sample_id}__tissue_positions__poly.parquet"

# Polygons from expanded nuclei (cell segmentation) segmentation
gdf_segmentation_path="${segm_dir}/\${sample_id}__nuclei_segmentation_stardist_masks__expanded__poly.parquet"

python3 "${PWD}/preprocessing/VisiumHD/004_map_barcodes.py" \
    --output_dir ${output_dir}/\${sample_id} \
    --sample_dir \${sample_dir} \
    --gdf_barcodes_path \${gdf_barcodes_path} \
    --gdf_segmentation_path \${gdf_segmentation_path} \
    --sample_id \${sample_id}

EOF