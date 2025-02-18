module load nvhpc/24.3
module load cuda/12.3
module load openmpi/4.1.6--nvhpc--24.3
cp Makefile_leonardo Makefile
make clean
make
#mpirun -np 2 ./mhit36 
