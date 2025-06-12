#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=200G
#SBATCH --time=20:00:00
#SBATCH --partition=gpunodes
#SBATCH --array=1

j=`sed -n 1p /scratch/perickso/private/ind_seq/popgen/phlash/phlash_pops.txt | cut -f 1`

source activate phlash
export PYTHONHOME=/usr/local/sw/anaconda/anaconda3/envs/phlash
python3.10 /scratch/perickso/private/ind_seq/popgen/scripts/phlash_test_script.py ${j}
