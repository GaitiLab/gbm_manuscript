#!/usr/bin/env bash
#SBATCH -J launch_120_compute_neighbors_metrics__knn
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
input_file="${PWD}/data/Xenium/processed/6509_A__BANKSY__meta_annot_w_expr_counts__ROI.rds"
dbscan_file="${PWD}/data/Xenium/processed/coloc_analysis/6509_A__k15nn_dbscan.rds"

sample_id="6509_A"
dbscan_dir="${PWD}/data/Xenium/processed/coloc_analysis/1_celltype_level"
output_dir="${dbscan_dir}/cell_type_pairs_perm"
group_varname="cell_type"

k_neighbors=15
approach="knn"

# ORIGINAL 
iter=250
n_cores=8

# Number of permutations/N x shuffling
job_min=1
job_max=400

echo $job_max


# ---- Observed mean ---- #

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J 120_compute_neighbors_metrics__knn
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
#SBATCH -J 120_compute_neighbors_metrics__knn
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

Rscript ${PWD}/preprocessing/Xenium//120_compute_neighbors_metrics.R \
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
