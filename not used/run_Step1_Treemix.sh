#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=50
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --partition=erickson

cd /scratch/perickso/private/indseq/popgen/treemix
/scratch/perickso/private/ind_seq/popgen/scripts/Step1_Treemix.sh  \
  zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink.prune.pruneddata.treemix.frq.gz \
  50 \
  500 \
  MIA_2019 \
  100 \
  /usr/local/sw/phylip.3.696/consense \
  5pops \
  1 \
  4 \
  10
