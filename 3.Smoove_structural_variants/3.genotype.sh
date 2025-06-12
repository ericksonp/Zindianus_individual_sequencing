#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --partition=basic
#SBATCH --array=2-285

REFERENCE=/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.gz

WORKING_FOLDER=/scratch/perickso/private/ind_seq/sv

SAMPLE_FILE=/scratch/perickso/private/ind_seq/sv/all_samples.txt

i=`sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLE_FILE`

singularity run /opt/containers/smoove.sif smoove genotype \
-d \
-x \
-p 1 \
--name $i \
--outdir $WORKING_FOLDER/results-genotyped/  \
--fasta $REFERENCE \
--vcf $WORKING_FOLDER/zap_all.sites.vcf.gz \
/scratch/perickso/private/ind_seq/RGSM_final_bams/${i}.RG.bam
