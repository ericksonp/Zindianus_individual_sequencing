#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=basic

######################################
# This script filters the whole genome snp set for linkage and then runs the core model over this filtered dataset to output a covariance matrix
######################################

# Code from https://github.com/JimWhiting91/guppy_convergence/blob/main/BayPass/scripts
# To update this code update the master, dataset and vcf as well as the POP_N 
# To change the linkage coefficient change the 0.2 to another linkage coefficient in line 27


# Environment
MASTER=/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales

DATASET=zap_ind_2023
#VCF_DIR=/gpfs/ts0/home/jw962/HP_LP_trials/phasing/phased_vcfs/$DATASET
VCF=/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/zaprionus.individual.2023.with.ids.vcf.gz
#VCF=/gpfs/ts0/home/jw962/guppy_research/five_aside_STAR_vcfs/five_aside_STAR_3033083_final.vcf.gz
POP_N=2

# Filter for linkage with plink
/usr/local/sw/plink-2.0-alpha/plink2 --vcf ${VCF} \
--out $MASTER/${DATASET}_plink_out_pruned --indep-pairwise 50 5 0.2 --allow-extra-chr #do we need the outputs file
#the 50 takes windows of 50 variants and then moves 5 and takes the next 50 (45 overlap)
#0.2 determines how much linkage is acceptable

# Merge .geno with SNPs
gunzip -c $VCF | grep -v "#" | cut -f3 | paste - $MASTER/${DATASET}.geno > $MASTER/snp_labelled_geno.txt

# Filter .geno with pruned SNPs
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' $MASTER/${DATASET}_plink_out_pruned.prune.in $MASTER/snp_labelled_geno.txt | cut -f2- > $MASTER/${DATASET}_LD_pruned.geno

# Dry run of core model
/usr/local/sw/baypass -npop $POP_N -gfile $MASTER/${DATASET}_LD_pruned.geno -outprefix $MASTER/${DATASET}_core -nthreads 16
