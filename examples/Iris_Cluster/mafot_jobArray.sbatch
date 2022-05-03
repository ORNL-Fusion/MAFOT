#!/bin/bash
#SBATCH --array=0-99
#SBATCH --job-name=mafotF
#SBATCH -p preemptable
#SBATCH -o batch_mafot.out
#SBATCH -n 32
#SBATCH -t 120
#SBATCH --mem-per-cpu=1G
##SBATCH --mail-user=wingen@fusion.gat.com
##SBATCH --mail-type=FAIL,ARRAY_TASKS  ## or BEGIN,END,FAIL,ALL,ARRAY_TASKS
#SBATCH --export=ALL
##module purge
module load default-paths
module load mafot
module list
mpirun -n 32 dtfoot_mpi _inner_${SLURM_ARRAY_TASK_ID}.dat co_Er_${SLURM_ARRAY_TASK_ID}