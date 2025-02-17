#!/bin/bash
#SBATCH --account="IscrB_*****"
#SBATCH --job-name="test"
#SBATCH --time=00:10:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --output=test.out
#SBATCH -p boost_usr_prod
#SBATCH --error=test.err

module load nvhpc/24.3 

#echo $HOSTNAME > hostname.dat
#./application_name  > out.4096.C.$SLURM_JOBID.dat
./mhit36