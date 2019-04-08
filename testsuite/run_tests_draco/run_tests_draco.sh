#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./test_result.%j
#SBATCH -e ./test_result.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J test_slurm
# Queue (Partition):
#SBATCH --partition=express
# Number of nodes and MPI tasks per node:
#SBATCH --ntasks-per-node=32
#  xxx  #SBATCH --nodes=16
#SBATCH --ntasks=27
#
#SBATCH --mail-type=none
#SBATCH --mail-user=<userid>@rzg.mpg.de
#
# Wall clock limit:
#SBATCH --time=0:20:00

module purge

module load git/2.16
module load intel/18.0
module load impi/2018.3
module load fftw-mpi/3.3.8
module load anaconda/2/5.1


# Run the program:
cd ../..
HOMEDIR=`pwd`
source $HOMEDIR/build/ESPRC

#cd $HOMEDIR/examples/polymer_melt_halfcell
cd $HOMEDIR/testsuite/polymer_melt
srun python polymer_melt_test.py 1 dd_duplicate 1 2
srun python polymer_melt_test.py 1 dd_duplicate 2 2
srun python polymer_melt_test.py 1 dd_duplicate 1 1
srun python polymer_melt_test.py 1 dd           2 2
srun python polymer_melt_test.py 1 dd           1 1

#cd $HOMEDIR/testsuite/free_movement/python
#srun python free_movement.py 1 1

