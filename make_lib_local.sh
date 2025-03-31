git clone https://github.com/NVIDIA/cuDecomp
cd cuDecomp
mkdir build
cd build
NVARCH=Linux_x86_64; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/24.3/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/24.3/compilers/bin:$PATH; export PATH
export PATH=$NVCOMPILERS/$NVARCH/24.3/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/24.3/comm_libs/mpi/man
cmake ..
make -j
