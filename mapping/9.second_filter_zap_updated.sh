#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --partition=basic

# This script is the second VCF filtering



WORKING_FOLDER=/scratch/perickso/private/ind_seq/haplotype_calling

cd $WORKING_FOLDER

vcftools --gzvcf $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat.vcf.gz \
   --recode \
   --remove /scratch/perickso/private/ind_seq/haplotype_calling/updated_samples_to_drop.txt \
  --out $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind

  vcftools --vcf  $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind.recode.vcf \
  --min-alleles 2 \
  --max-alleles 2 \
  --max-missing 0.9 \
  --minQ 40 \
  --min-meanDP 10 \
  --max-meanDP 45 \
  --recode \
  --recode-INFO-all \
  --out $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind_secondfilter


  bgzip $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind_secondfilter.recode.vcf
  tabix $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind_secondfilter.recode.vcf.gz

  #find singleton SNPs
  vcftools --gzvcf $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind_secondfilter.recode.vcf.gz --singletons --out singletons.txt
  cut -f1,2 singletons.txt.singletons > singleton_positions.txt


  vcftools --gzvcf $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind_secondfilter.recode.vcf.gz \
     --recode \
     --exclude-positions singleton_positions.txt \
    --out $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind_secondfilter_dropsingletons

bgzip $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind_secondfilter_dropsingletons.recode.vcf
tabix $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind_secondfilter_dropsingletons.recode.vcf.gz
