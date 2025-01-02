#!/usr/bin/env bash
#SBATCH -J launch-R104_annot_seurat_obj.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

cluster_name="BANKSY_snn_res.2"

input_dir="${PWD}/data/VisiumHD/processed/seurat_objects/prepped"
lookup_table_dir=$PWD/misc/annot

# For H4H array
sample_ids="${PWD}/misc/VisiumHD_sample_info.txt"
partition="himem"
cpus_per_task=1
mem="40G"
time="00:15:00"

job_max=$(wc -l < "${sample_ids}")
# job_max=1
job_min=1
echo $job_max

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J R104_annot_seurat_obj.R
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

echo "Sample ID: \${sample_id}"

seurat_obj="${input_dir}/\${sample_id}.rds"
look_up_table="${PWD}/misc/VisiumHD_annot/\${sample_id}.xlsx"

Rscript "${PWD}/preprocessing/VisiumHD/104_annot_seurat_obj.R" \
    --filename \${sample_id} \
    --seurat_obj \${seurat_obj} \
    --output_dir ${input_dir} \
    --look_up_table \${look_up_table} \
    --cluster_name ${cluster_name} \
    --sheet_name ${cluster_name}


EOF