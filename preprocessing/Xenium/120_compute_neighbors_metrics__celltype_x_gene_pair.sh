#!/usr/bin/env bash
#SBATCH -J launch_120_compute_neighbors_metrics__knn__cellpair_x_interaction
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# Inputs
input_file="${PWD}/data/Xenium/processed/6509_A__BANKSY__meta_w_celltype_x_gene_pair_labels__ROI.rds"
dbscan_file="${PWD}/data/Xenium/processed/coloc_analysis/6509_A__k15nn_dbscan.rds"

sample_id="6509_A"
dbscan_dir="${PWD}/data/Xenium/processed/coloc_analysis/2_celltype_x_gene_pair_level"

k_neighbors=15
approach="knn"

interactions="${PWD}/data/Xenium/processed/cellpair_x_interaction.txt"

# ORIGINAL
iter=250
n_cores=8

# Number of permutations/N x shuffling
job_min=1
job_max=400

# Number of interactions
start_interaction=1
n_interactions=$(wc -l < "${interactions}")

for ((x=${start_interaction}; x<=${n_interactions}; x++))
do
# Gene-pairs
group_varname=$(awk "NR==${x}" ${interactions})
output_dir="${dbscan_dir}/${group_varname}_k${k_neighbors}_perm"


# ---- Observed mean ---- #

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 120_compute_neighbors_metrics__knn__cellpair_x_interaction
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${n_cores}
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

echo "Activating conda environment..."
source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

echo "Gene pair: ${group_varname}"
echo "Output directory: ${output_dir}"

Rscript ${PWD}/preprocessing/Xenium/120_compute_neighbors_metrics.R \
    --meta ${input_file} \
    --output_dir ${dbscan_dir} \
    --group_varname ${group_varname} \
    --sample_id ${sample_id} \
    --k_neighbors ${k_neighbors} \
    --dbscan ${dbscan_file} \
    --n_cores \${SLURM_CPUS_PER_TASK} \
    --approach ${approach}

EOF

# ---- Permutation testing/shuffling labels ---- #
sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J 120_compute_neighbors_metrics__cellpair_x_interaction
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${n_cores}
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

Rscript ${PWD}/preprocessing/Xenium/120_compute_neighbors_metrics.R \
    --meta ${input_file} \
    --output_dir ${output_dir} \
    --group_varname ${group_varname} \
    --sample_id ${sample_id} \
    --k_neighbors ${k_neighbors} \
    --ix \${SLURM_ARRAY_TASK_ID} \
    --dbscan ${dbscan_file} \
    --n_cores \${SLURM_CPUS_PER_TASK} \
    --approach ${approach} \
    --n_iter ${iter}

EOF

done 
