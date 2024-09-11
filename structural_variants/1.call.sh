#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --partition=basic
#SBATCH --array=1-285

#make a file with all sample names
#cat /scratch/perickso/private/raw_data/Fall2020_sequencing/file_list.txt /scratch/perickso/private/raw_data/published/zap_to_download.txt /scratch/perickso/private/raw_data/published/SRR_Acc_List.txt  > /scratch/perickso/private/ind_seq/sv/all_samples.txt

REFERENCE=/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.gz

WORKING_FOLDER=/scratch/perickso/private/ind_seq/sv

SAMPLE_FILE=/scratch/perickso/private/ind_seq/sv/all_samples.txt

#if the final file exists, exit
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $SAMPLE_FILE`

if [[ `wc -l  "$WORKING_FOLDER/raw-genotyped/${i}-smoove.genotyped.vcf.gz"` > 1 ]]
then
	echo "calling already completed, exiting"
	exit
else
	echo "beginning calling "
fi

singularity run /opt/containers/smoove.sif smoove call \
--outdir $WORKING_FOLDER/raw-genotyped \
--exclude /scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.bed  \
--name $i \
--fasta $REFERENCE \
-p 1 \
--genotype /scratch/perickso/private/ind_seq/RGSM_final_bams/${i}.RG.bam
