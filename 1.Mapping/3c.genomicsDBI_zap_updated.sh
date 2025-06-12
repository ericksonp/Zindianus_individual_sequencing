#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=50
#SBATCH --mem=700G
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --array=1-5



# This script will merge gVCFs into a unified database for genotype calling.

#Run R script to generate file list once

if [[ -e "/scratch/perickso/private/ind_seq/Samples_to_haplotype.txt" ]]
then
	echo "starter file exists"
	echo "lets move on"
	date
else
	echo "input file being created"
	Rscript /scratch/perickso/private/ind_seq/scripts/2b.vcfnames.R
	date
fi

#Name of pipeline
PIPELINE=MergegVCFs

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/perickso/private/ind_seq/haplotype_calling

# User defined inputs -- this represents the name of the samples

#combine together two previous finals for new sample mapped
#cat /scratch/perickso/private/ind_seq/Samples_to_haplotype_published.txt /scratch/perickso/private/ind_seq/Samples_to_haplotype.txt > final_sample_to_haplotype.txt

sample_map=/scratch/perickso/private/ind_seq/haplotype_calling/final_sample_to_haplotype.txt

#This file looks like this
#  sample1      sample1.vcf.gz
#  sample2      sample2.vcf.gz
#  sample3      sample3.vcf.gz

#Parameters

#Java
JAVAMEM=600G
CPU=50

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
# Generate Folders and files
###########################################################################
###########################################################################
# this part of the script will check and generate, if necessary, all of the output folders used in the script

if [[ -d "TEMP_MERGEVCF" ]]
then
	echo "Working TEMP_MERGEVCF folder exist"
	echo "lets move on"
	date
else
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/TEMP_MERGEVCF
	date
fi

###########################################################################
###########################################################################
# Merge VCFs using GenomicsDBImport
###########################################################################
###########################################################################

#use slurm array batch id to choose thread to work with
#less /scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.fai | cut -f1 > /scratch/perickso/private/ind_seq/haplotype_calling/intervals.txt

intervals=/scratch/perickso/private/ind_seq/haplotype_calling/intervals.txt

i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

  /usr/local/sw/gatk-4.2.0.0/gatk --java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM}" \
       GenomicsDBImport \
       --genomicsdb-workspace-path $WORKING_FOLDER/Zap_DB_updated_${i} \
       --batch-size 50 \
       --sample-name-map $sample_map \
       --tmp-dir $WORKING_FOLDER/TEMP_MERGEVCF \
       --reader-threads $CPU \
			 -L $i
