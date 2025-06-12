#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --partition=erickson
#SBATCH --array=0-4



#Name of pipeline
PIPELINE=updated_Raw_Quality_Statistics

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/perickso/private/ind_seq/haplotype_calling

#Input file
IN_VCF=$WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat.vcf.gz
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
# Survey quality of final VCF
###########################################################################
###########################################################################

analyses=(\
"--depth" \
"--site-mean-depth" \
"--site-quality" \
"--missing-indv" \
"--missing-site" \
)

vcftools \
--gzvcf $IN_VCF \
`echo ${analyses[${SLURM_ARRAY_TASK_ID}]}` \
--out $PIPELINE

echo "pipeline completed" $(date)
