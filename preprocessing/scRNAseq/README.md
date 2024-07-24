# Overview

## Inference of cell-cell-interactions

1. Follow steps 1-3 from the **Quick start** section at the GitHub repository [scrnaseq-cellcomm-pipeline](https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline).
2. Run `bash local_nf_run_cci_pipeline.sh` (local) or `slurm_nf_run_cci_pipeline.sh` (HPC with slurm workload manager). You have to change the following inside the bash script:`project_dir`, `input_file`, `output_dir` and `nf_exec`, these are also marked with 'TODO'. The data used for `input_file` can be downloaded from TODO XXXXX. 
When running `slurm_nf_run_cci_pipeline.sh`, please check the configuration and change accordingly, depending on your computing environment.
