#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300G
#SBATCH --time=1:00:00
#SBATCH --partition=basic


REFERENCE=/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.gz
WORKING_FOLDER=/scratch/perickso/private/ind_seq/sv

cd $WORKING_FOLDER

singularity run /opt/containers/smoove.sif smoove paste \
--name zap_all_called_sv \
$WORKING_FOLDER/results-genotyped/*.vcf.gz
