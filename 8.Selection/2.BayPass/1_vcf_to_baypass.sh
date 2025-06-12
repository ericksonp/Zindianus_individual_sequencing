#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=erickson

# Script makes input files for baypass based on VCFs, sorted by chromosome
# this is where we start to use the code from https://github.com/JimWhiting91/guppy_convergence/blob/main/BayPass/scripts/01_VCF_to_baypass.sh
# to update this code update the master, dataset and vcf as well as the pop_array and last paste (change them to the populations you want to test)

MASTER=/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales
DATASET=zap_ind_2023
#VCF_DIR=/gpfs/ts0/home/jw962/HP_LP_trials/phasing/phased_vcfs/$DATASET
VCF=/scratch/perickso_shared/alexandra/baypass/data/zaprionus.individual.2023.vcf.gz


# Define the population
pop_array=(MIA VA)

# Loop over popmaps
for POP in "${pop_array[@]}"
do

# Outputs go here
#mkdir $MASTER/data/$POP

# Get counts from VCF files, adding the derived flag means that sites with an ancestral allele will be correctly ordered
vcftools --gzvcf $VCF --counts2 --keep /scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/${POP}_sample_ids.txt --out $MASTER/${DATASET}_${POP}_tmp

# Now we need to edit the format of the outputted file to Baypass standard - 2 cols per SNP, 1 row per SNP, Allele counts
tail -n+2 $MASTER/${DATASET}_${POP}_tmp.frq.count | cut -f 5,6 > $MASTER/${DATASET}_${POP}.geno

rm -f $MASTER/${DATASET}_${POP}_tmp*
done

#combining the whole
paste $MASTER/zap_ind_2023_MIA.geno $MASTER/zap_ind_2023_VA.geno > $MASTER/zap_ind_2023.geno