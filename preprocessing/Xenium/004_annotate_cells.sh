#!/usr/bin/env bash
#SBATCH -J 005_annotate_cells.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

sample_id="6509_A"
seurat_obj_path="${PWD}/data/Xenium/processed/6509_A__BANKSY.rds"
annotated_clusters_path="${PWD}/misc/xenium_annot.xlsx"
output_dir="${PWD}/data/Xenium/processed"
cells_oi_path="${PWD}/misc/6509_A_roi_cell_ids.csv"

echo "Activating conda environment..."
source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

Rscript ${PWD}/preprocessing/Xenium/004_annotate_cells.R \
    --output_dir ${output_dir} \
    --sample_id ${sample_id} \
    --cluster_varname ${cluster_varname} \
    --seurat_obj_path ${seurat_obj_path} \
    --annotated_clusters_path ${annot_clusters} \
    --cells_oi ${cells_oi}