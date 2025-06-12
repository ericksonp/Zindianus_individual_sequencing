#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=erickson

# Code from https://github.com/JimWhiting91/guppy_convergence/blob/main/BayPass/scripts
#I did not actually use this because I ran everything manually
# The script will need to be added to the location put after Rscript because I do not know how you organize your code to be run
# I also split their get_xtx_significance.R into multiple files because Spydur will not let me do terminal commands from R

MASTER=/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales

POP=$2

# Run R script to assess significance and output outlier SNPs -- this might need work because I changed enough in the xtx file
Rscript $MASTER/R/4b_get_xtx_significance.R $POP