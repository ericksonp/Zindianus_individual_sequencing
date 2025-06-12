#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 50
#SBATCH --mem=360G
#SBATCH --time=100:00:00
#SBATCH --partition=johnson

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -c $SLURM_CPUS_PER_TASK Rscript /scratch/perickso/private/ind_seq/popgen/scripts/LD_decay.R
