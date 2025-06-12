#! /bin/bash


#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=0:15:00
#SBATCH --partition=basic
#SBATCH --array=1-266
#SBATCH --dependency=afterany:23680:23454


#cat /scratch/perickso/private/raw_data/Fall2020_sequencing/file_list.txt  /scratch/perickso/private/raw_data/published/zap_to_download.txt > /scratch/perickso/private/ind_seq/all_individual_samples.txt
#mkdir /scratch/perickso/private/ind_seq/coverage

BIO_SAMPLES=/scratch/perickso/private/raw_data/published/SRR_Acc_List.txt
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $BIO_SAMPLES`


samtools coverage /scratch/perickso/private/ind_seq/RGSM_final_bams/${i}.RG.bam > /scratch/perickso/private/ind_seq/coverage/${i}.coverage.txt
