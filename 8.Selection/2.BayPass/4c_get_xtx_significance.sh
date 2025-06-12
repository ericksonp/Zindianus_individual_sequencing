#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=erickson

# Code from https://github.com/JimWhiting91/guppy_convergence/blob/main/BayPass/scripts
# I did not use this but instead wrote the command (/usr/local/sw/baypass -npop 2 -gfile G.zap_ind_2023_sim -outprefix zap_ind_2023_sim -nthreads 16) in the terminal because it only took about 2 minutes

MASTER=/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales

/usr/local/sw/baypass -npop 2 -gfile $MASTER/G.zap_ind_2023_sim -outprefix $MASTER/zap_ind_2023_sim -nthreads 16