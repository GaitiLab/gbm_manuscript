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

job_min=1

input_dir="${PWD}/data/Xenium/processed/coloc_analysis/2_celltype_x_gene_pair_level"
output_dir="${input_dir}"

sample_id="6509_A"
k_neighbors=15
approach="knn"
metric="mean_n_cells"

# Number of interactions
interactions="${PWD}/data/Xenium/processed/cellpair_x_interaction.txt"
n_interactions=$(wc -l < "${interactions}")
n_iter=400

# TODO remove later
n_interactions=2
n_iter=3

job_max=${n_interactions}

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
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

group_varname=\$(awk "NR==\${SLURM_ARRAY_TASK_ID}" ${interactions})

perm_dir="${input_dir}/\${group_varname}_k${k_neighbors}_perm"
input_file="${input_dir}/${sample_id}__\${group_varname}__k${k_neighbors}nn_mean.rds"

echo "Group varname: \${group_varname}"
echo "Permutation directory: \${perm_dir}"
echo "Input file: \${input_file}"

Rscript ${PWD}/preprocessing/Xenium/121_combine_neighbors_metrics.R \
    --obs_df \${input_file} \
    --output_dir ${output_dir} \
    --group_varname \${group_varname} \
    --sample_id ${sample_id} \
    --k_neighbors ${k_neighbors} \
    --perm_dir \${perm_dir} \
    --n_cores \${SLURM_CPUS_PER_TASK} \
    --metric ${metric} \
    --approach ${approach} \
    --n_iter ${n_iter}

EOF
