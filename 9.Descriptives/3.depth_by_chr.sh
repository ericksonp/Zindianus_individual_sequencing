#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --partition=basic
#SBATCH --array=1-5


#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/perickso/private/ind_seq/popgen

#Input file
IN_VCF=$WORKING_FOLDER/zaprionus.individual.2023.vcf.gz

cd $WORKING_FOLDER

###########################################################################
###########################################################################
# Survey quality of final VCF
###########################################################################
###########################################################################

chr=${SLURM_ARRAY_TASK_ID}

vcftools \
--gzvcf $IN_VCF \
--chr Scaffold_${chr} \
--depth \
--out $PIPELINE.${chr}

echo "pipeline completed" $(date)
