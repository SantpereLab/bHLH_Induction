#!/bin/bash

#############################################################
###### Summarize quality checks after fastp processing ######
#############################################################

# -----------------
# -- Use MultiQC --
# -----------------

# Load Modules
module load Miniconda3/20240927
multiqc --version

# Set project path
root_dir="/path/to/the/project"

# 

# Run multiqc to gather qc metrics
multiqc "$root_dir/data/02_preprocessed_fastqs/00_temp/" -o "$root_dir/data/02_preprocessed_fastqs/00_temp/z_multiqc_report"
