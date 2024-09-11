#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=700G
#SBATCH --time=96:00:00
#SBATCH --partition=medium
#SBATCH --array=1-5
#SBATCH --dependency=afterany:39102

# This script will conduct genotype calling on the GenomeDBI object



#Name of pipeline
PIPELINE=GenotypeGVCFs

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/perickso/private/ind_seq/haplotype_calling

#Reference genome
REFERENCE=/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta

#Intervals to analyze
intervals=/scratch/perickso/private/ind_seq/haplotype_calling/intervals.txt

#Parameters

#Java
JAVAMEM=650G
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
# Identify the Genome database to genotyoe
###########################################################################
###########################################################################

GenomeDB_path=`echo $WORKING_FOLDER/Zap_DB_updated_${i}`

###########################################################################
###########################################################################
# Genotype call the samples in the DBI merged set
###########################################################################
###########################################################################

 gatk --java-options "-Xmx${JAVAMEM}" GenotypeGVCFs \
   -R $REFERENCE \
   -V gendb://$GenomeDB_path \
   -O $WORKING_FOLDER/${i}.updated.genotyped.raw.vcf.gz

echo ${i} "done" $(date)
