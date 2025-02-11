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
  type(cudecompPencilInfo) :: piX, piY, piZ  ! size of the pencils in x- y- and z-configuration
  integer(8) :: nElemX, nElemY, nElemZ, nElemWork
  ! cuFFT
  integer :: planX, planY, planZ
  integer :: batchsize
  integer :: status



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
  pc = 2
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

  !
  CHECK_CUDECOMP_EXIT(cudecompGridDescAutotuneOptionsSetDefaults(options))
  options%dtype = CUDECOMP_DOUBLE_COMPLEX
  if (comm_backend == 0) then
    options%autotune_transpose_backend = .true.
  endif

  ! initialize cuDecomp with the config file 
  CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_desc, config, options))

  if (rank == 0) then
     write(*,"('Running on ', i0, ' x ', i0, ' process grid ...')") config%pdims(1), config%pdims(2)
     write(*,"('Using ', a, ' backend ...')") cudecompTransposeCommBackendToString(config%transpose_comm_backend)
  end if

  ! get pencil info
  ! This function returns a pencil struct (piX, piY or piZ) that contains the shape, global lower and upper index bounds (lo and hi), 
  ! size of the pencil, and an order array to indicate the memory layout that will be used (to handle permuted, axis-contiguous layouts).
  ! Additionally, there is a halo_extents data member that indicates the depth of halos for the pencil, by axis, if provided..
  ! Side note:  ! cudecompGetPencilInfo(handle, grid_desc, pinfo_x, 1, [1, 1, 1]) <- in this way the x-pencil also have halo elements
  ! If no halo regions are necessary, a NULL pointer can be provided in place of this array (or omitted)
  ! Pencil info in x-configuration present in PiX (shape,lo,hi,halo_extents,size)
  CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piX, 1))
  nElemX = piX%size !<- number of total elments in x-configuratiion (including halo)
  ! Pencil info in Y-configuration present in PiY
  CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piY, 2))
  nElemY = piY%size
  ! Pencil info in Z-configuration present in PiZ
  CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piZ, 3))
  nElemZ = piZ%size

 ! get workspace size
  CHECK_CUDECOMP_EXIT(cudecompGetTransposeWorkspaceSize(handle, grid_desc, nElemWork))

 ! show the order 1=x, 2=y, 3=z
 ! x-pencils are x,y,z
 ! y-pencils are y,z,x
 ! z-pencils are z,y,x
 ! This is to make the setup of cuFFT easier, cuFFT cannot do FFT along inner directions and also for performance (no stride)
 ! if (rank .eq. 0) then
 !   write(*,*) "Order in X-pencil is", piX%order(1),piX%order(2),piX%order(3)
 !   write(*,*) "Order in Y-pencil is", piY%order(1),piY%order(2),piY%order(3)
 !   write(*,*) "Order in Z-pencil is", piZ%order(1),piZ%order(2),piZ%order(3)
 ! endif

  ! CUFFT initialization
  ! Create plans (forward and backward are the same!)
  batchSize = piX%shape(2)*piX%shape(3) !<- number of FFT (from x-pencil dimension)
  status = cufftPlan1D(planX, nx, CUFFT_Z2Z, batchSize)
  if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating X plan'

  ! it's always 2 and 3 because y-pencil have coordinates y,z,x
  batchSize = piY%shape(2)*piY%shape(3)
  status = cufftPlan1D(planY, ny, CUFFT_Z2Z, batchSize)
  if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating Y plan'

  ! it's always 2 and 3 because y-pencil have coordinates z,y,x
  batchSize = piZ%shape(2)*piZ%shape(3)
  status = cufftPlan1D(planZ, nz, CUFFT_Z2Z, batchSize)
  if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating Z plan'








  call mpi_finalize(ierr)

end program main
