#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --time=48:00:00
#SBATCH --partition=basic
#SBATCH --array=1-268

input=/scratch/perickso/private/ind_seq/popgen/samples_in_full_vcf.txt

samp=`sed -n ${SLURM_ARRAY_TASK_ID}p $input | cut -f 1`

echo "working on ${samp} "

if [[ -e "/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.${samp}.recode.vcf.gz" ]]
then
  echo "gzipped vcf exists, let's move on"
  date
else

  vcftools \
    --gzvcf /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.vcf.gz   \
    --recode \
    --indv ${samp} \
    --out /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.${samp}
  #statements

  bgzip /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.${samp}.recode.vcf
  tabix /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.${samp}.recode.vcf.gz

fi

whatshap phase \
-o /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.${samp}.phased.vcf \
 --reference=/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta \
 /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.${samp}.recode.vcf.gz \
 /scratch/perickso/private/ind_seq/RGSM_final_bams/${samp}.RG.bam
