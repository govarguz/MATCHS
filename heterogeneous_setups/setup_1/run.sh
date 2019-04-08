#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./log/tjob.out.%j
#SBATCH -e ./log/tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J test_slurm
# Queue (Partition):
#SBATCH --partition=express
# Number of nodes and MPI tasks per node:
#SBATCH --ntasks-per-node=32
##SBATCH --nodes=1
#SBATCH --ntasks=256
#
#SBATCH --mail-type=none
#SBATCH --mail-user=<userid>@rzg.mpg.de
#
# Wall clock limit:
#SBATCH --time=0:50:00


module purge
module load git/2.16
module load intel/18.0
module load impi/2018.3
module load fftw-mpi/3.3.8
module load anaconda/2/5.1

# Run the program:
source ../../buildVp/ESPRC

srun python polymer_melt_from_restart_adapted.py 1 2


