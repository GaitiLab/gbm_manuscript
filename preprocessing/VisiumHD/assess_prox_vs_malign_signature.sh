#!/usr/bin/env bash
#SBATCH -J assess_prox_vs_malign_signature.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out


output_dir="${PWD}/data/VisiumHD/processed/10k"
input_dir="${PWD}/data/VisiumHD"

label="BANKSY_snn_res.2_annot"
n_iter=10000
lower_limit=0.10
upper_limit=0.90
markers="/cluster/projects/gaitigroup/Users/Joan/002_Project_GBM/000_misc/gene_lists/neftel_signatures.rds"
nbin=24

# For H4H array
partition="veryhimem"
cpus_per_task=8
mem="100G"
time="1-00:00:00"

for k_neighbors in 25 30
do 

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -J assess_prox_vs_malign_signature.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${cpus_per_task}
#SBATCH --mem=${mem}
#SBATCH --partition=${partition}
#SBATCH --time=${time}
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "r-4.3.3"

Rscript "${PWD}/preprocessing/VisiumHD/assess_prox_vs_malign_signature.R" \
    --input_dir ${input_dir} \
    --output_dir ${output_dir} \
    --n_iter ${n_iter} \
    --k_neighbors ${k_neighbors} \
    --label ${label} \
    --lower_limit ${lower_limit} \
    --upper_limit ${upper_limit} \
    --markers ${markers} \
    --n_cores \${SLURM_CPUS_PER_TASK} \
    --nbin ${nbin}

EOF

done