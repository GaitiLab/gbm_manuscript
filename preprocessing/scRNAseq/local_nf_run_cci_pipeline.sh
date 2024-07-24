# TODO change to path to directory with the cloned Github repo.  
project_dir="scrnaseq-cellcomm-pipeline"

echo "PIPELINE CONFIGURATION..."
# Input seurat file
# TODO change to path to 'gbm_regional_study.rds'
input_file="gbm_regional_study.rds"

# Output directory
# TODO change to directory for saving outputs
output_dir="${project_dir}/output/"

# ---- PIPELINE CONFIGURATION ---- #
init_step=1

# Pre-processing
sample_var="Sample"
annot="CCI_CellClass_L2_2"
condition_var="Region"
patient_var="Patient"
min_patients=2
min_cells=50
is_confident=1

# Cell-cell interactions
n_perm=1000
min_pct=0.10
alpha=0.05

# ----  NEXTFLOW CONFIGURATION ---- #
# Path to nextflow executable
# TODO change to path to nextflow executable
nf_exec="${HOME}/nextflow-23.04.3-all"

# Work directory - all executed tasks (processes) are stored here
work_dir="${project_dir}/nf-work"
nf_profile="conda"
# Output directory for: trace, report + timeline by NextFlow
outdir="${project_dir}/nf-logs"

# Create directories
mkdir -p "${output_dir}"
mkdir -p "${project_dir}/nf-logs"

echo "Running pipeline..."
# # Start the pipeline
${nf_exec} run ${project_dir} -with-report -with-trace -resume \
    -profile ${nf_profile} \
    -w ${work_dir} \
    --input_file $input_file \
    --sample_var ${sample_var} \
    --annot ${annot} \
    --min_cells ${min_cells} \
    --n_perm ${n_perm} \
    --min_pct ${min_pct} \
    --alpha $alpha \
    --init_step $init_step \
    --condition_var $condition_var \
    --patient_var $patient_var \
    --min_patients $min_patients \
    --is_confident ${is_confident} \
    --outdir ${outdir} \
    --output_dir ${output_dir}


echo "Done!"
