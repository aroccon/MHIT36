#define CHECK_CUDECOMP_EXIT(f) if (f /= CUDECOMP_RESULT_SUCCESS) call exit(1)

program main
  use cudafor
  use cudecomp
  use cufft
  use mpi

  implicit none

  ! grid dimensions
  integer :: nx, ny, nz ! this is the global grid
  integer :: comm_backend
  integer :: pr, pc
  ! MPI variables
  integer :: rank, ranks, ierr
  integer :: localRank, localComm
  ! cudecomp
  type(cudecompHandle) :: handle
  type(cudecompGridDesc) :: grid_desc
  type(cudecompGridDescConfig) :: config
  type(cudecompGridDescAutotuneOptions) :: options
  integer :: pdims(2) ! pr x pc pencils
  integer :: gdims(3) ! global grid dimensions

  ! MPI initialization, put in rank the local MPI rank number and ranks total number
  ! Same procedura defined in the cuDecomp documentation
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, ranks, ierr)

  call mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, localComm, ierr)
  call mpi_comm_rank(localComm, localRank, ierr)
  ierr = cudaSetDevice(localRank) !assign GPU to MPI rank

  ! Define grid and decomposition
  ! hard coded, then from input
  nx = 64
  ny = 64
  nz = 64
  pr = 2
  pc = 1
  comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P

  CHECK_CUDECOMP_EXIT(cudecompInit(handle, MPI_COMM_WORLD))

  ! config is a struct and pr and pc are the number of pencils along the two directions
  ! gdims is the global grid
  ! create an uninitialized configuration struct and initialize it to defaults using cudecompGridDescConfigSetDefaults. 
  ! Initializing to default values is required to ensure no entries are left uninitialized.
  CHECK_CUDECOMP_EXIT(cudecompGridDescConfigSetDefaults(config))
  pdims = [pr, pc] !pr and pc are the number of pencil along the different directions
  config%pdims = pdims
  gdims = [nx, ny, nz]
  config%gdims = gdims
  config%transpose_comm_backend = comm_backend
  config%transpose_axis_contiguous = .true.

  !do some checks 
  !CHECK_CUDECOMP_EXIT(cudecompGridDescAutotuneOptionsSetDefaults(options))
  !options%dtype = CUDECOMP_DOUBLE_COMPLEX
  !if (comm_backend == 0) then
  !  options%autotune_transpose_backend = .true.
  !endif

  ! initialize cuDecomp with the config file 
  CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_desc, config, options))

  if (rank == 0) then
     write(*,"('Running on ', i0, ' x ', i0, ' process grid ...')") config%pdims(1), config%pdims(2)
     write(*,"('Using ', a, ' backend ...')") cudecompTransposeCommBackendToString(config%transpose_comm_backend)
  end if


  call mpi_finalize(ierr)

end program main
