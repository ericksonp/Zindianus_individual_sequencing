#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200G
#SBATCH --time=1:00:00
#SBATCH --partition=basic



REFERENCE=/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.gz
WORKING_FOLDER=/scratch/perickso/private/ind_seq/sv


singularity run /opt/containers/smoove.sif smoove merge \
--name zap_all \
-f $REFERENCE \
--outdir $WORKING_FOLDER \
$WORKING_FOLDER/raw-genotyped/*.genotyped.vcf.gz
