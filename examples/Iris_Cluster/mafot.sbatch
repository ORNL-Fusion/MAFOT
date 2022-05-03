#!/bin/bash
#SBATCH --job-name=mafotL
#SBATCH -p preemptable
#SBATCH -o batch_mafot.out
#SBATCH -n 20
#SBATCH -t 120
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL  ## or BEGIN,END,FAIL,ALL
#SBATCH --export=ALL
##module purge
module load default-paths
module load mafot
module list
mpirun -n 20 dtlaminar_mpi _lam.dat any_tag_of_your_choice