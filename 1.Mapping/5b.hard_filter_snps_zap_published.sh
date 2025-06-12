#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --array=1-5
#SBATCH --dependency=afterany:39203

# This script is a pipeline which applies GATK variant recalibration to VCFs.

#Load Modules



#Load local software


#Name of pipeline
PIPELINE=FilterSNPs

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/perickso/private/ind_seq/haplotype_calling

#Reference genome
REFERENCE=/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.gz

#Intervals to analyze
intervals=/scratch/perickso/private/ind_seq/haplotype_calling/intervals.txt


#Parameters
#Java
JAVAMEM=180G
CPU=1

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

#Interval to analyze
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

echo ${i} "is being processed" $(date)

###########################################################################
###########################################################################
# Filter VCFs to retain just SNPs
###########################################################################
###########################################################################

	vcftools \
		--gzvcf $WORKING_FOLDER/${i}.updated.genotyped.raw.vcf.gz \
		--recode \
		--recode-INFO-all \
		--remove-indels \
		--out ${i}.updated.genotypedSNPs.raw

#bgzip and tabix
bgzip ${i}.updated.genotypedSNPs.raw.recode.vcf
tabix ${i}.updated.genotypedSNPs.raw.recode.vcf.gz


###########################################################################
###########################################################################
# Run Filtering using Comeault 2020 parameters
###########################################################################
###########################################################################

gatk --java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM}" VariantFiltration \
      -R $REFERENCE \
      -V $WORKING_FOLDER/${i}.updated.genotypedSNPs.raw.recode.vcf.gz \
			-O $WORKING_FOLDER/${i}.updated.hardfilterSNP.comeault2020.vcf \
			--filter-expression "QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
			--filter-name "Comeault2020_filters"


#bgzip and tabix
	bgzip $WORKING_FOLDER/${i}.updated.hardfilterSNP.comeault2020.vcf
	tabix $WORKING_FOLDER/${i}.updated.hardfilterSNP.comeault2020.vcf.gz


echo ${i} "done" $(date)
