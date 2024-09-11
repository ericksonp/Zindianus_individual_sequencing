#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --partition=erickson

Rscript /scratch/perickso/private/ind_seq/sv/scripts/all_scaffolds_inversions_combined.R
