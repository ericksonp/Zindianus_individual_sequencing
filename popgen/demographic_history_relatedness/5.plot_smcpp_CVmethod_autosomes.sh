#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=100G
#SBATCH --time=96:00:00
#SBATCH --partition=basic
#SBATCH --array=1-26

pop=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/perickso/private/ind_seq/popgen/smcpp/pops_for_smcpp.txt | cut -f 1`

pop_list=`paste -sd" " /scratch/perickso/private/ind_seq/popgen/smcpp/pops_for_smcpp.txt`

echo "plotting final plot for ${pop}"

  cd /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/${pop}_cv_10folds_16knots_10000gen_autosomes/
  /opt/containers/smc++ plot /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/${pop}_cv_10folds_16knots_10000gen_autosomes/${pop}_CV10fold_final_plot.png \
  model.final.json \
    -g .08 \
    --csv


echo "plotting folds"
for i in {0..9}; do
  echo ${i}
  cd /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/${pop}_cv_10folds_16knots_10000gen_autosomes/fold${i}/
  /opt/containers/smc++ plot /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/${pop}_cv_10folds_16knots_10000gen_autosomes/fold${i}/${pop}_fold${i}_final_plot.png \
  model.final.json \
    -g .08 \
    --csv
done
