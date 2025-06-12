#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --partition=basic

Rscript /scratch/perickso/private/ind_seq/popgen/scripts/baypass_window_plotting.R
