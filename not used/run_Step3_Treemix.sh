#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=50
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --partition=erickson

cd /scratch/perickso/private/ind_seq/popgen/treemix
source activate treemix

# /scratch/perickso/private/ind_seq/popgen/scripts/Step3_Treemix.sh  \
#   zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink.prune.pruneddata.treemix.frq.gz \
#   50 \
#   5000 \
#   MIA_2019 \
#   100 \
#   3 \
#   5pops5ksnps \
#   30 \
#   5pops5ksnps_constree.newick \
#   /usr/local/sw/phylip.3.696/consense

  /scratch/perickso/private/ind_seq/popgen/scripts/Step3_Treemix.sh  \
    zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink.prune.pruneddata.treemix.frq.gz \
    50 \
    10 \
    MIA_2019 \
    100 \
    1 \
    5pops10snps \
    30 \
    5pops10snps_constree.newick \
    /usr/local/sw/phylip.3.696/consense
