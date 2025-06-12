#This script filters the entirety of the vcf to make new samples with only the file of interest
# change directories of the .txt files and the vcf of all samples according to your test

#MIA 
cat <<EOF | sbatch
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=basic
bcftools view \
  --samples /scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/MIA_sample_ids.txt --force-samples\
  -O z \
  -o /scratch/perickso_shared/alexandra/BayPassAlexEdited/datallVAFemalesvsFLfemalesa/MIA_samples.vcf.gz \
  /scratch/perickso_shared/alexandra/baypass/data/zaprionus.individual.2023.vcf.gz
EOF

#VA
cat <<EOF | sbatch
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=basic
bcftools view \
  --samples /scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/VA_sample_ids.txt --force-samples\
  -O z \
  -o /scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/VA_samples.vcf.gz \
  /scratch/perickso_shared/alexandra/baypass/data/zaprionus.individual.2023.vcf.gz
EOF