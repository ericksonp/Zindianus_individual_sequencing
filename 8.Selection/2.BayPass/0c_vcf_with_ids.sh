#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=basic

#I would recomend not rerunning this file but instead copying the zaprionus.individual.2023.with.ids.vcf.gz file to the location you need it

bcftools annotate --set-id '%CHROM\_%POS\' /scratch/perickso_shared/alexandra/baypass/data/zaprionus.individual.2023.vcf.gz -o  /scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/zaprionus.individual.2023.with.ids.vcf.gz -Oz