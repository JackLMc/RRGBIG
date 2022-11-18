#!/bin/bash
#SBATCH -J FastqcMultiqc              		# A single job name for the array
#SBATCH -n 8                       		# Number of cores
#SBATCH -N 1                       		# All cores on one machine
#SBATCH --mem 20000                 		# in MB
#SBATCH -t 2-0:00                  		# Maximum execution time (D-HH:MM)
#SBATCH -o job_%A_%a.log        		# Standard output
#SBATCH -e job_%A_%a.log        		# Standard error
#SBATCH --account=turnerjd-rrg-bioinf-training # Bear account


# This is a slurm job submission script for the production of FastQC reports

# Important head of each script
set -e # Fail the script on the first error
module purge; module load bluebear # For reproducibility. Get the default apps

## Load your reuired Apps
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4

## Create temporary directory for temporary files
BB_WORKDIR=$(mktemp -d /scratch/${USER}_${SLURM_JOBID}.XXXXXX)


## Run the actual script ========================================================
mkdir -p ../../2_quality_control

fastqc `find .. -name '*fastq.gz' -print` --outdir=../../2_quality_control

multiqc ../../2_quality_control -o ../../2_quality_control
