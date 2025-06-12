#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=200:00:00
#SBATCH --partition=basic

######################################
# This script filters the whole genome snp set for linkage and then runs the core model over this filtered dataset to output a covariance matrix
######################################

# Code from https://github.com/JimWhiting91/guppy_convergence/blob/main/BayPass/scripts
# To update this code update the master


#3
#module load ifort/2017.4.196-GCC-6.4.0-2.28

MASTER=/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales

# This script runs BayPass as an array job on the HPC, running each of the 16 subsets on their own CPU
# Assumes that the inputs have already been subsetted
POP=$1

# Run XtX
#spits out zaprionus_individual_2023_3a_baypass.log
/usr/local/sw/baypass -npop 2 -seed zap_ind_2023_core -gfile $MASTER/zap_ind_2023.geno -omegafile $MASTER/zap_ind_2023_core_mat_omega.out  -outprefix $MASTER/zap_ind_2023 -nthreads 8
#baypass -npop 2  -seed ${MOAB_JOBARRAYINDEX} -gfile $MASTER/data/subsets/five_aside_STAR_sub${MOAB_JOBARRAYINDEX}_${POP}.geno -omegafile $MASTER/data/five_aside_STAR_${POP}_CovMatrix.txt -outprefix $MASTER/outputs/${POP}_five_aside_STAR_sub${MOAB_JOBARRAYINDEX} -nthreads 8

# Run with Env covariate e.g. HP-LP BUT WITH AUX model that can account for spatial SNPs
#baypass -npop 10 -seed ${MOAB_JOBARRAYINDEX} -auxmodel -isingbeta 1.0 -gfile $MASTER/data/subsets/${POP}_sub${MOAB_JOBARRAYINDEX}.geno -omegafile $MASTER/data/${POP}_CovMatrix_Final -efile $MASTER/data/${POP}.env -outprefix $MASTER/outputs/${POP}_auxmodel_sub${MOAB_JOBARRAYINDEX} -nthreads 8