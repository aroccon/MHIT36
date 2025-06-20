#!/bin/bash
#SBATCH --account="tra25_openhack"
#SBATCH --job-name="cudec"
#SBATCH --time=00:10:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --output=test.out
#SBATCH -p boost_usr_prod
#SBATCH --qos=boost_qos_dbg
#SBATCH --error=test.err

module load nvhpc/24.3
module load cuda/12.3
module load openmpi/4.1.6--nvhpc--24.3
#export LD_LIBRARY_PATH=/leonardo_scratch/large/userexternal/aroccon0/MHIT36_cuDecomp/cuDecomp/build/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/leonardo_scratch/large/userexternal/lenzenbe/RE95_256_cuDec/cuDecomp/build/lib:$LD_LIBRARY_PATH
CURRENT_DIR="$(pwd)"
ROOT_DIR="$(dirname "$CURRENT_DIR")/cuDecomp/build/lib"
echo "Using directory: $ROOT_DIR"
export LD_LIBRARY_PATH=$ROOT_DIR:$LD_LIBRARY_PATH


#mpirun -n 4 ./mhit36
mpirun -n 4 nsys profile --trace=cuda,nvtx,mpi -o profile_output_%q{SLURM_PROCID} --stats=true ./mhit36

