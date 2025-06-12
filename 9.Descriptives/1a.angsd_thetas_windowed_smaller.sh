#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --partition=basic


cd /scratch/perickso/private/ind_seq/popgen/angsd

while read x; do
  echo "working on:"
  echo $x

# /usr/local/sw/angsd/misc/thetaStat \
# do_stat \
# /scratch/perickso/private/ind_seq/popgen/angsd/${x}.thetas.idx  \
# -win 10000 \
# -step 5000 \
# -outnames ${x}.thetasWindow10000_5000
# done < /scratch/perickso/private/ind_seq/popgen/org_files/populations.txt
# 
# /usr/local/sw/angsd/misc/thetaStat \
# do_stat \
# /scratch/perickso/private/ind_seq/popgen/angsd/${x}.thetas.idx  \
# -win 10000 \
# -step 10000 \
# -outnames ${x}.thetasWindow10000_10000
# done < /scratch/perickso/private/ind_seq/popgen/org_files/populations.txt

/usr/local/sw/angsd/misc/thetaStat \
do_stat \
/scratch/perickso/private/ind_seq/popgen/angsd/${x}.thetas.idx  \
-win 5000 \
-step 5000 \
-outnames ${x}.thetasWindow5000_5000
done < /scratch/perickso/private/ind_seq/popgen/org_files/populations.txt
