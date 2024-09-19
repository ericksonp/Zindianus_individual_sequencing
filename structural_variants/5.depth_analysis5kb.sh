#! /bin/bash

#SBATCH -N 1
#SBATCH -c 50
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --partition=basic


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -c $SLURM_CPUS_PER_TASK Rscript /scratch/perickso/private/ind_seq/sv/scripts/depth_analysis_5kbwin.R
