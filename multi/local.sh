NVARCH=Linux_x86_64; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/24.3/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/24.3/compilers/bin:$PATH; export PATH
export PATH=$NVCOMPILERS/$NVARCH/24.3/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/24.3/comm_libs/mpi/man
LD_LIBRARY_PATH=/home/milton/MHIT36_cuDecomp/cuDecomp/build/lib
#clean folder output
rm -rf output
mkdir output
cp Makefile_local Makefile
#rm *.dat
make clean
make
mpirun -np 2 ./mhit36 
