#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --dependency=afterany:39205


# This script is a pipeline which gather VCFs from all chromosomes.

#Load Modules
PICARD=/usr/local/sw/picard/picard.jar


#Name of pipeline

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/perickso/private/ind_seq/haplotype_calling

#Parameters
#Java
JAVAMEM=450G
CPU=2

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status

# Move to working directory
cd $WORKING_FOLDER

###########################################################################
###########################################################################
# Gather VCFs to make a final VCF

###########################################################################
###########################################################################
#NEW VERSION: WITH INDEL REMOVAL


#make input file:
ls /scratch/perickso/private/ind_seq/haplotype_calling/Scaffold*updated.hardfilterSNP.comeault2020.noINDEL.vcf > /scratch/perickso/private/ind_seq/haplotype_calling/updated.filtered_vcfs.list

java -Xmx$JAVAMEM \
 -jar $PICARD GatherVcfs \
  I=$WORKING_FOLDER/updated.filtered_vcfs.list \
  O=$WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL.vcf


#remove comeault2020 filter
 vcftools --vcf $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL.vcf \
   --remove-filtered-all \
   --recode \
   --out /$WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS

#remove SNPs in annotated repetitive regions

#
bedtools intersect -v -header \
  -a /$WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS.recode.vcf \
  -b /scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.gff \
  > $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat.vcf

##bgzip and tabix
	bgzip $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat.vcf
 	tabix $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat.vcf.gz


echo "done" $(date)
