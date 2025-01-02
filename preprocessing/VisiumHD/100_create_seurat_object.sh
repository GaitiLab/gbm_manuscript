#!/usr/bin/env bash
#SBATCH -J launch-100_create_seurat_object.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

input_dir_parent="${PWD}/data/VisiumHD/processed/004_map_barcodes"
scalefactors_dir="/cluster/projects/gaitigroup/Data/GBM_VisiumHD/20240708_LH00244_0121_B22MN7GLT3_Gaiti_Yiyan_Visium"

output_dir="${PWD}/data/VisiumHD/processed/seurat_objects/raw"

# About the data
is_spatial=0
assay="RNA"

# For H4H array
sample_ids="${PWD}/misc/VisiumHD_sample_info.txt"
partition="himem"
cpus_per_task=1
mem="40G"
time="00:30:00"

job_min=1
job_max=$(wc -l < "${sample_ids}")

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J 100_create_seurat_object.R
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
source "\$HOME/miniforge3/bin/activate" "r-4.3.3"

sample_id=\$(cut -d$'\t' -f1 ${sample_ids} | awk "NR==\${SLURM_ARRAY_TASK_ID}")
image_name=\$(cut -d$'\t' -f2 ${sample_ids} | awk "NR==\${SLURM_ARRAY_TASK_ID}")

scalefactors_path="${scalefactors_dir}/\${sample_id}/outs/binned_outputs/square_002um/spatial/scalefactors_json.json"


Rscript "${PWD}/preprocessing/VisiumHD/100_create_seurat_object.R" \
    --input_dir ${input_dir} \
    --output_dir ${output_dir} \
    --sample_id \${sample_id} \
    --is_spatial ${is_spatial} \
    --assay ${assay}
EOF