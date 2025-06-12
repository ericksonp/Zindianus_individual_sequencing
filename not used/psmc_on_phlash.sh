#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH --time=20:00:00
#SBATCH --partition=basic
#SBATCH --array=1-180


i=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/perickso/private/ind_seq/popgen/phlash/psmc_validation_input.txt | cut -f 1 `
j=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/perickso/private/ind_seq/popgen/phlash/psmc_validation_input.txt | cut -f 2 `


psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /scratch/perickso/private/ind_seq/popgen/phlash/${i}/${j}.autosome4.psmc  /scratch/perickso/private/ind_seq/popgen/phlash/${i}/${j}.autosome.psmcfa
psmc -N25 -t15 -r5 -p "2+2+25*2+4+6" -o /scratch/perickso/private/ind_seq/popgen/phlash/${i}/${j}.autosome22.psmc  /scratch/perickso/private/ind_seq/popgen/phlash/${i}/${j}.autosome.psmcfa
psmc -N25 -t15 -r5 -p "2+2+25*2+4+6" -o /scratch/perickso/private/ind_seq/popgen/phlash/${i}/${j}.autosome1111.psmc  /scratch/perickso/private/ind_seq/popgen/phlash/${i}/${j}.autosome.psmcfa
