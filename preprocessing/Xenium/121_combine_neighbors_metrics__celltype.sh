#!/usr/bin/env bash
#SBATCH -J launch_121_combine_neighbors_metrics_knn
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

input_file="${PWD}/data/Xenium/processed/coloc_analysis/1_celltype_level/6509_A__cell_type__k15nn_mean.rds"
perm_dir="${PWD}/data/Xenium/processed/coloc_analysis/1_celltype_level/cell_type_pairs_perm"
output_dir="${PWD}/data/Xenium/processed/coloc_analysis/1_celltype_level"

sample_id="6509_A"
k_neighbors=15
approach="knn"
group_varname="cell_type"
metric="mean_n_cells"

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J launch_121_combine_neighbors_metrics_knn
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

echo "Activating conda environment..."
source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

Rscript ${PWD}/preprocessing/Xenium/121_combine_neighbors_metrics.R \
    --obs_df ${input_file} \
    --output_dir ${output_dir} \
    --group_varname ${group_varname} \
    --sample_id ${sample_id} \
    --k_neighbors ${k_neighbors} \
    --perm_dir ${perm_dir} \
    --n_cores \${SLURM_CPUS_PER_TASK} \
    --metric ${metric} \
    --approach ${approach}

EOF