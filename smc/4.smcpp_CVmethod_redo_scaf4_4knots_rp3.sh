#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --partition=basic
#SBATCH --array=1-12

#variable for populations--prepared input files with script smcpp_fileprep.R

#make RepeatMasked bed file:
#gff2bed < /scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.gff > /scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.bed


# touch /scratch/perickso/private/ind_seq/popgen/smcpp/pops_to_rerun
#
# pop_list=`paste -sd" " /scratch/perickso/private/ind_seq/popgen/smcpp/pops_for_smcpp.txt`
#
# for s in $pop_list; do
#   echo $s
#   if  [[ -f "/scratch/perickso/private/ind_seq/popgen/smcpp/${s}/${s}_cv_10folds_16knots_10000gen_autosomes/model.final.json" ]]
# then
#   echo "completed"
# else
#   echo "job did not complete"
# fi
# done

#   echo $s >> /scratch/perickso/private/ind_seq/popgen/smcpp/pops_to_rerun
# fi
# done

#get missing data

#comm -23 <(sort /scratch/perickso/private/ind_seq/popgen/smcpp/pops_for_smcpp.txt) <(sort /scratch/perickso/private/ind_seq/popgen/smcpp/pops_to_rerun) > /scratch/perickso/private/ind_seq/popgen/smcpp/pops_to_rerun_autosomes.txt

pop=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/perickso/private/ind_seq/popgen/phlash/phlash_pops.txt | cut -f 1 `

samps_list=`paste -sd, /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}.samps.smc.txt`

nsamps=`wc -l < /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}.samps.smc.txt`

samps=`paste -sd" " /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}.samps.smc.txt`

echo "starting ${pop}"

# #setting up slurm folder to tell where error occurs
# echo "setting up folders"
# #creating new folders and sub folders for this script
# if [[ -d "/scratch/perickso/private/ind_seq/popgen/smcpp/${pop}" ]]
# then
#   date
# else
#   mkdir /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}
# fi
# #loop for chromosomes 1-5, check to see if folder  is there and make it if not
# for ind in $samps
# do
# echo "setting " ${ind} " as distinguished individual"
#   for chr in {1..5}
#   do
#     echo "chromosome" ${chr}
#
#     #making smc file; remove if one already exists from previous erroneous program
#
#     if [[ -f "/scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/zaprionus.individual.recodeafricaasreference.2023.${pop}_Scaffold_${chr}_${ind}.smc.gz" ]]
#     then
#       echo "removing old file"
#       rm -f /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/zaprionus.individual.recodeafricaasreference.2023.${pop}_Scaffold_${chr}_${ind}.smc.gz
#     fi
#
#     echo "making smc input file for " $chr  $pop $ind
#     /opt/containers/smc++ vcf2smc \
#     -d $ind $ind \
#     --cores 10 \
#     -m /scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.bed.gz \
#     /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.recodeafricaasreference.2023.vcf.gz \
#     /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/zaprionus.individual.recodeafricaasreference.2023.${pop}_Scaffold_${chr}_${ind}.smc.gz \
#     Scaffold_${chr} ${pop}:${samps_list}
#
# done
# done

#note that output was specified here twice...
echo "all files prepped, doing smc calculations"

unset PYTHONHOME
unset PYTHONPATH

  /opt/containers/smc++ cv -o /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/ \
  --timepoints 2 10000 \
   --folds $nsamps \
   --cores 20 \
   --knots 4 \
   -rp 3 \
    2.8e-9 \
    /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/*Scaffold_4*.smc.gz \
    -o /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/${pop}_cv_10folds_4knots_10000gen_scaf4_rp3



    echo "plotting final plot for ${pop}"

      cd /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/${pop}_cv_10folds_4knots_10000gen_scaf4_rp3/
      /opt/containers/smc++ plot /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/${pop}_cv_10folds_4knots_10000gen_scaf4_rp3/${pop}_CV10fold_final_plot.png \
      model.final.json \
        -g .08 \
        --csv


    echo "plotting folds"
    for i in $(seq 1 $nsamps); do
      echo ${i}
      cd /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/${pop}_cv_10folds_4knots_10000gen_scaf4_rp3/fold${i}/
      /opt/containers/smc++ plot /scratch/perickso/private/ind_seq/popgen/smcpp/${pop}/${pop}_cv_10folds_4knots_10000gen_scaf4_rp3/fold${i}/${pop}_fold${i}_final_plot.png \
      model.final.json \
        -g .08 \
        --csv
    done
