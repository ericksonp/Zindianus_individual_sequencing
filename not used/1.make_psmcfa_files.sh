#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=200:00:00
#SBATCH --partition=erickson
#SBATCH --array=11

#original was array 1-12 but 11 failed

j=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/perickso/private/ind_seq/popgen/phlash/phlash_pops.txt | cut -f 1`

  echo ${j}
  #make a folder for each population
  #mkdir /scratch/perickso/private/ind_seq/popgen/phlash/${j}
  #cycle through files and make psmcfa files
  while read i ; do
    echo ${i}
    bcftools mpileup -f /scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta /scratch/perickso/private/ind_seq/RGSM_final_bams/${i}.RG.bam | bcftools call -c  | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > /scratch/perickso/private/ind_seq/popgen/phlash/${j}/${i}.v2.fq.gz
    /usr/local/sw/psmc/utils/fq2psmcfa -q 20 /scratch/perickso/private/ind_seq/popgen/phlash/${j}/${i}.v2.fq.gz > /scratch/perickso/private/ind_seq/popgen/phlash/${j}/${i}.psmcfa
    done < /scratch/perickso/private/ind_seq/popgen/smcpp/${j}.samps.smc.txt
