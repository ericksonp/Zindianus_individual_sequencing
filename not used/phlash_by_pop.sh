#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=220G
#SBATCH --time=200:00:00
#SBATCH --partition=erickson
#SBATCH --array=1-8

j=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/perickso/private/ind_seq/popgen/phlash/phlash_pops.txt | cut -f 1`

#conda init bash
source activate phlash
export PYTHONHOME=/usr/local/sw/anaconda/anaconda3/envs/phlash
python3.10 /scratch/perickso/private/ind_seq/popgen/scripts/phlash_by_pop.py ${j}
