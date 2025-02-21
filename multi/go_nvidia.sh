#!/bin/bash
#SBATCH --account="IscrB_ARESS"
#SBATCH --job-name="cudec"
#SBATCH --time=00:5:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --output=test.out
#SBATCH -p boost_usr_prod
#SBATCH --error=test.err

module load nvhpc/24.3
module load cuda/12.3
module load openmpi/4.1.6--nvhpc--24.3
export LD_LIBRARY_PATH=/leonardo_scratch/large/userexternal/aroccon0/MHIT36_cuDecomp/cuDecomp/build/lib:$LD_LIBRARY_PATH



mpirun -n 2 ./mhit36
