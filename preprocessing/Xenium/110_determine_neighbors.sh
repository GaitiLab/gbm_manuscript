#!/usr/bin/env bash
#SBATCH -J launch_110_determine_neighbors.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

# General setup
approach="knn"
k_neighbors=15
sample_id="6509_A"

input_file="${PWD}/data/Xenium/processed/6509_A__BANKSY__meta_annot_w_expr_counts__ROI.rds"
output_dir="${PWD}/data/Xenium/processed/coloc_analysis"

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J 110_determine_neighbors.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=JohnDoe@mail.com
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out

echo "Activating conda environment..."
source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

Rscript "${PWD}/preprocessing/Xenium/110_determine_neighbors.R" \
    --annotated_cells_df ${input_file} \
    --output_dir ${output_dir} \
    --sample_id ${sample_id} \
    --k_neighbors ${k_neighbors} \
    --approach ${approach} 
EOF
