#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition=basic
#SBATCH --array=1-5

WORKING_FOLDER=/scratch/perickso/private/ind_seq/haplotype_calling
cd $WORKING_FOLDER


intervals=/scratch/perickso/private/ind_seq/haplotype_calling/intervals.txt

i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

source /usr/local/sw/urtools/slurmstuff/condafy.sh
conda init
conda activate erickson
#
# pixy --stats pi fst dxy \
# --vcf  ${i}.updated.genotypedSNPs.allsites.raw.recode.vcf.gz \
# --populations /scratch/perickso/private/ind_seq/popgen/populations_for_pixy.txt \
# --window_size 5000 \
# --n_cores 10 \
# --chromosomes $i \
# --output_folder /scratch/perickso/private/ind_seq/popgen/pixy \
# --output_prefix pixy_4pops_5kb_$i
#
# pixy --stats pi  \
# --vcf  ${i}.updated.genotypedSNPs.allsites.raw.recode.vcf.gz \
# --populations /scratch/perickso/private/ind_seq/popgen/subpopulations_for_pixy.txt \
# --window_size 10000 \
# --n_cores 10 \
# --chromosomes $i \
# --output_folder /scratch/perickso/private/ind_seq/popgen/pixy \
# --output_prefix pixy_subpops_10kb_$i

pixy --stats pi  \
--vcf  ${i}.updated.genotypedSNPs.allsites.raw.recode.vcf.gz \
--populations /scratch/perickso/private/ind_seq/popgen/subpopulations_for_pixy.txt \
--window_size 5000 \
--n_cores 10 \
--chromosomes $i \
--output_folder /scratch/perickso/private/ind_seq/popgen/pixy \
--output_prefix pixy_subpops_5kb_$i
