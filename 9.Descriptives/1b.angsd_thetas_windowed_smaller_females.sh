#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --partition=johnson


cd /scratch/perickso/private/ind_seq/popgen/angsd

#for x in FL Africa VA-CM VA-HPO; do
# echo "working on:"
# echo $x
#Florida should have been FL, rerun that way
# x=FL

#rerun with x=VA to combine CM and HPO
x=VA
/usr/local/sw/angsd/angsd \
  -bam /scratch/perickso/private/ind_seq/popgen/org_files/${x}_bam.txt \
  -doSaf 1 \
  -anc /scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta \
  -GL 1 \
  -P 1 \
  -out /scratch/perickso/private/ind_seq/popgen/angsd/${x}


  /usr/local/sw/angsd/misc/realSFS \
  /scratch/perickso/private/ind_seq/popgen/angsd/${x}.saf.idx \
  -P 1 \
  -fold 1 > /scratch/perickso/private/ind_seq/popgen/angsd/${x}.sfs

  /usr/local/sw/angsd/misc/realSFS \
  saf2theta \
  /scratch/perickso/private/ind_seq/popgen/angsd/${x}.saf.idx  \
  -outname /scratch/perickso/private/ind_seq/popgen/angsd/${x} \
  -sfs /scratch/perickso/private/ind_seq/popgen/angsd/${x}.sfs \
  -fold 1

  /usr/local/sw/angsd/misc/thetaStat do_stat /scratch/perickso/private/ind_seq/popgen/angsd/${x}.thetas.idx


/usr/local/sw/angsd/misc/thetaStat \
do_stat \
/scratch/perickso/private/ind_seq/popgen/angsd/${x}.thetas.idx  \
-win 10000 \
-step 5000 \
-outnames ${x}.thetasWindow10000_5000



/usr/local/sw/angsd/misc/thetaStat \
do_stat \
/scratch/perickso/private/ind_seq/popgen/angsd/${x}.thetas.idx  \
-win 5000 \
-step 5000 \
-outnames ${x}.thetasWindow5000_5000
