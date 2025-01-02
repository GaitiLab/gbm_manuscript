#!/usr/bin/env bash
#SBATCH -J launch_pipeline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Joan.Kant@uhn.ca
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

export NXF_OFFLINE='true'
export NXF_DEFAULT_DSL=2

module load java/18

pipeline_dir="/cluster/projects/gaitigroup/Pipelines/scrnaseq-cellcomm-pipeline"

nextflow run ${pipeline_dir} \
    -profile conda,slurm \
    -params-file "nf-params.yml" \
    --outdir "output" \
    -c "gaitilab.config" \
    -resume
