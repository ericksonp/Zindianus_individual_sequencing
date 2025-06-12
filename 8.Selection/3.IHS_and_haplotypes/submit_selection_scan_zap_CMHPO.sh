#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --partition=basic


Rscript /scratch/perickso/private/ind_seq/popgen/scripts/selection_scan_zap_CMHPO.R
