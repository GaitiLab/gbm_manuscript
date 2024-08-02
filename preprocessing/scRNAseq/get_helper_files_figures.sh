#!/usr/bin/env bash

# ---- Figure XX-volcano.Rmd ---- #
# Compute average expression data using 'gbm_regional_study' per 'CCI_CellClass_L2_2' for each Sample
# Activating 'cci' environment from https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline
source "${HOME}/miniforge3/bin/activate" "cci"
# TODO to adapted by user
input_file="gbm_regional_study.rds"
output_dir="output"

# DO NOT CHANGE when trying to reproducing figure
slot="counts"
assay="RNA"
groupbyvar="Sample,CCI_CellClass_L2_2"

Rscript "${work_dir}/scripts/UTILS-compute_avg_expr.R" \
    --input_file $input_file \
    --output_dir ${output_dir} \
    --slot ${slot} \
    --assay ${assay} \
    --groupbyvar ${groupbyvar}

