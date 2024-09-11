#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition=basic
#SBATCH --dependency=afterany:79010


WORKING_FOLDER=/scratch/perickso/private/ind_seq/haplotype_calling
cd $WORKING_FOLDER

source /usr/local/sw/urtools/slurmstuff/condafy.sh
conda init
conda activate erickson

# pixy --stats pi fst dxy \
# --vcf  Scaffold_3.updated.genotypedSNPs.allsites.raw.recode.females.vcf.gz \
# --populations /scratch/perickso/private/ind_seq/popgen/populations_for_pixy_females.txt \
# --window_size 5000 \
# --n_cores 10 \
# --chromosomes Scaffold_3 \
# --output_folder /scratch/perickso/private/ind_seq/popgen/pixy \
# --output_prefix pixy_4pops_5kb_Scaffold_3_females



pixy --stats pi \
--vcf  Scaffold_3.updated.genotypedSNPs.allsites.raw.recode.females.vcf.gz \
--populations /scratch/perickso/private/ind_seq/popgen/subpopulations_for_pixy_females.txt \
--window_size 5000 \
--n_cores 10 \
--chromosomes Scaffold_3 \
--output_folder /scratch/perickso/private/ind_seq/popgen/pixy \
--output_prefix pixy_subpops_10kb_Scaffold_3_females
