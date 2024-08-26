# Overview

## Inference of cell-cell-interactions

1. Follow steps 1-3 from the **Quick start** section at the GitHub repository [scrnaseq-cellcomm-pipeline](https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline).
2. Run `bash preprocessing/scRNAseq/local_nf_run_cci_pipeline.sh` (local) or `sbatch preprocessing/scRNAseq/slurm_nf_run_cci_pipeline.sh` (HPC with slurm workload manager). You have to change the following inside the bash script:`project_dir`,  `input_file`,  `output_dir` and `nf_exec`, these are also marked with 'TODO'. The data used for `input_file` can be downloaded from TODO XXXXX. 
When running `slurm_nf_run_cci_pipeline.sh` , please check the configuration and change accordingly, depending on your computing environment.

## Additional files required for generating figures

For Figure 2b-volcano, the following script has to be run: `UTILS-compute_avg_expr.sh`

Needed for Xenium analyses:

1. Run `sbatch preprocessing/scRNAseq/UTILS-extract_interactions.sh`
2. Run `sbatch preprocessing/scRNAseq/UTILS-xenium_overlap.sh`
