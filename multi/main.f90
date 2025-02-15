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
  integer :: npx, npy, npz
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
  integer :: halo(3) ! halo extensions
  integer :: halo_ext ! 0 no halo, 1 means 1 halo
  type(cudecompPencilInfo) :: piX, piY, piZ  ! size of the pencils in x- y- and z-configuration
  integer(8) :: nElemX, nElemY, nElemZ, nElemWork, nElemWork_halo
  logical :: halo_periods(3)
  ! cuFFT
  integer :: planX, planY, planZ
  integer :: batchsize
  integer :: status
  ! other variables (wavenumber, grid location)
  real(8), allocatable :: x(:), kx(:)
  real(8) :: dx,lx
  integer :: i,j,k,il,jl,kl,ig,jg,kg
  integer, parameter :: Mx = 1, My = 0, Mz = 0
  real(8), device, allocatable :: kx_d(:)
  real(8), parameter :: twopi = 8.0_8*atan(1.0_8)
  ! workign arrays
  complex(8), allocatable :: phi(:), ua(:,:,:)
  complex(8), device, allocatable :: phi_d(:)
  complex(8), pointer, device, contiguous :: work_d(:), work_halo_d(:)
  character(len=40) :: namefile
  ! Code variables
  complex(8), allocatable :: rhsp_complex(:,:,:)
  real(8), allocatable :: rhsp(:,:,:), p(:,:,:)




!########################################################################################################################################
! 1. INITIALIZATION OF MPI AND cuDECOMP AUTOTUNING : START
!########################################################################################################################################
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
  lx=twopi
  nx = 64
  ny = nx
  nz = nx
  pr = 1
  pc = 2
  halo_ext=1
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
  halo = [0, halo_ext, halo_ext] ! no halo along x neeed because is periodic and in physical space i have x-pencil
  ! for transpositions
  config%transpose_comm_backend = comm_backend
  config%transpose_axis_contiguous = .true.
  ! for halo exchanges
  config%halo_comm_backend = CUDECOMP_HALO_COMM_MPI
  ! Setting for periodic halos in all directions (non required to be in config)
  halo_periods = [.true., .true., .true.]


  
  CHECK_CUDECOMP_EXIT(cudecompGridDescAutotuneOptionsSetDefaults(options))
  options%dtype = CUDECOMP_DOUBLE_COMPLEX
  if (comm_backend == 0) then
    options%autotune_transpose_backend = .true.
  endif

  ! initialize cuDecomp with the config file 
  CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_desc, config, options))

  ! Print information on configuration
  if (rank == 0) then
      write(*,"('Running on ', i0, ' x ', i0, ' process grid ...')") config%pdims(1), config%pdims(2)
      write(*,"('Using ', a, ' transpose backend ...')") &
             cudecompTransposeCommBackendToString(config%transpose_comm_backend)
      write(*,"('Using ', a, ' halo backend ...')") &
             cudecompHaloCommBackendToString(config%halo_comm_backend)
  endif

  ! get pencil info
  ! This function returns a pencil struct (piX, piY or piZ) that contains the shape, global lower and upper index bounds (lo and hi), 
  ! size of the pencil, and an order array to indicate the memory layout that will be used (to handle permuted, axis-contiguous layouts).
  ! Additionally, there is a halo_extents data member that indicates the depth of halos for the pencil, by axis.
  ! Side note:  ! cudecompGetPencilInfo(handle, grid_desc, pinfo_x, 1, [1, 1, 1]) <- in this way the x-pencil also have halo elements
  ! If no halo regions are necessary, a NULL pointer can be provided in place of this array (or omitted)
  ! Pencil info in x-configuration present in PiX (shape,lo,hi,halo_extents,size)
  CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piX, 1, halo))
  nElemX = piX%size !<- number of total elments in x-configuratiion (including halo)
  ! Pencil info in Y-configuration present in PiY
  CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piY, 2))
  nElemY = piY%size
  ! Pencil info in Z-configuration present in PiZ
  CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piZ, 3))
  nElemZ = piZ%size

  ! Get workspace sizes
  ! For transpostion
  CHECK_CUDECOMP_EXIT(cudecompGetTransposeWorkspaceSize(handle, grid_desc, nElemWork))
  ! For halo
  CHECK_CUDECOMP_EXIT(cudecompGetHaloWorkspaceSize(handle, grid_desc, 1, halo, nElemWork_halo))


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

  ! define grid
  ! to test the code (impose initial condition, can be then removed)
  allocate(x(nx),kx(nx))
  dx = lx/nx
  x(1)= 0
  do i = 2, nx
     x(i) = x(i-1) + dx
  enddo

  do i = 1, nx/2
    kx(i) = (i-1)*(twoPi/lx)
  enddo
  do i = nx/2+1, nx
    kx(i) = (i-1-nx)*(twoPi/LX)
  enddo
  ! allocate k_d on the device (later on remove and use OpenACC + managed memory?)
  allocate(kx_d, source=kx)

  !########################################################################################################################################
  ! 1. INITIALIZATION AND cuDECOMP AUTOTUNING : END
  !########################################################################################################################################





  !########################################################################################################################################
  ! START STEP 2: ALLOCATE ARRAYS
  !########################################################################################################################################
  ! allocate arrays
  allocate(phi(max(nElemX, nElemY, nElemZ))) !largest among the pencil
  allocate(phi_d, mold=phi) ! phi on device
  allocate(ua(nx, piX%shape(2), piX%shape(3)))
  ! Pressure variable
  allocate(rhsp_complex(piX%shape(1), piX%shape(2), piX%shape(3)))
  allocate(rhsp(piX%shape(1), piX%shape(2), piX%shape(3))) 
  allocate(p(piX%shape(1), piX%shape(2), piX%shape(3))) 


  ! work_d is the device array for transposition (CUDA)
  CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_desc, work_d, nElemWork))
  ! work_halo_d is the device array for halo exchanges (CUDA)
  CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_desc, work_halo_d, nElemWork_halo))
  !########################################################################################################################################
  ! END STEP2: ALLOCATE ARRAYS)
  !########################################################################################################################################





  !########################################################################################################################################
  ! START STEP 3: FLOW AND PHASE FIELD INIT
  !########################################################################################################################################
  ! 3.1 Read from data without halo grid points (avoid out-of-bound)
  ! 3.2 Call halo exchnages along Y and Z for u,v,w and phi
  !########################################################################################################################################
  ! END STEP 3: FLOW AND PHASE FIELD INIT
  !########################################################################################################################################


  !########################################################################################################################################
  ! START STEP 4: PHASE-FIELD SOLVER (EXPLICIT)
  !########################################################################################################################################
  ! 4.1 RHS computation, no need of halo update in thoery
  ! 4.2 Get phi at n+1 
  ! 4.3 Call halo exchnages along Y and Z for phi 
  !########################################################################################################################################
  ! END STEP 4: HASE-FIELD SOLVER (EXPLICIT)
  !########################################################################################################################################


  !########################################################################################################################################
  ! START STEP 5: USTAR COMPUTATION
  !########################################################################################################################################
  ! 5.1 compute rhs 
  ! 5.2 obtain ustar
  ! 5.3 Call halo exchnages along Y and Z for u,v,w
  ! 5.4 Compute divergence (i.e. the rhsp)
  !########################################################################################################################################
  ! END STEP 5: USTAR COMPUTATION 
  !########################################################################################################################################


  !########################################################################################################################################
  ! START STEP 6: POISSON SOLVER FOR PRESSURE
  !########################################################################################################################################
  ! initialize phi and analytical solution
  ! for the moment keep it similar to the cuDecomp example, can be optimized in the future 
  block
   complex(8), pointer :: phi3(:,:,:)
   call c_f_pointer(c_loc(phi), phi3, [piX%shape(1), piX%shape(2), piX%shape(3)])

   !rhsp is a standard array (similar to those that are in the code)
   do kl = 1, piX%shape(3)
      kg = piX%lo(3) + kl - 1 
      do jl = 1, piX%shape(2)
         jg = piX%lo(2) + jl - 1 
         do i = 1, piX%shape(1)
           !rhsp_complex(i,jl,kl) = cmplx(rhsp(i,jl,kl),0.0) !<- use this once solver is ok
           ! For Poisson testing (knwon solutions)
           rhsp_complex(i,jl,kl) = cmplx(sin(Mx*x(i)),0.0) 
           !rhsp_complex(i,jl,kl) = cmplx(sin(Mx*x(i))*sin(My*x(jg))*sin(Mz*x(kg)),0.0)  ! RHS of Poisson equation
         enddo
      enddo
   enddo

   !move these arrays into the device pointer to feed the Poisson solver
   do kl = 1, pix%shape(3)
     kg = piX%lo(3) + kl - 1 
     do jl = 1, pix%shape(2)
        jg = piX%lo(2) + jl - 1 
        do i = 1, pix%shape(1)
          phi3(i,jl,kl) = rhsp_complex(i,jl,kl)  ! RHS of Poisson equation is phi3 (pointer to phi and then copied into the GPU)
           !phi3(i,jl,kl) = 
           ua(i,jl,kl) = -phi3(i,jl,kl)/(Mx**2 + My**2 + Mz**2) ! Solution of Poisson equation
        enddo
     enddo
   enddo
  end block

  ! H2D transfer using CUDA
  phi_d = phi

  ! input rhs 
  write(namefile,'(a,i3.3,a)') 'rhs_',rank,'.dat'
  open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
  write(55) real(phi, KIND=8)
  close(55)

  ! do the FFT3D forward 
  ! phi(x,y,z) -> phi(kx,y,z)
  status = cufftExecZ2Z(planX, phi_d, phi_d, CUFFT_FORWARD)
  if (status /= CUFFT_SUCCESS) write(*,*) 'X forward error: ', status
  ! phi(kx,y,z) -> phi(y,z,kx)
  CHECK_CUDECOMP_EXIT(cudecompTransposeXToY(handle, grid_desc, phi_d, phi_d, work_d, CUDECOMP_DOUBLE_COMPLEX,Pix%halo_extents, [0,0,0]))
  ! phi(y,z,kx) -> phi(ky,z,kx)
  status = cufftExecZ2Z(planY, phi_d, phi_d, CUFFT_FORWARD)
  if (status /= CUFFT_SUCCESS) write(*,*) 'Y forward error: ', status
  ! phi(ky,z,kx) -> phi(z,kx,ky)
  CHECK_CUDECOMP_EXIT(cudecompTransposeYToZ(handle, grid_desc, phi_d, phi_d, work_d, CUDECOMP_DOUBLE_COMPLEX))
  ! phi(z,kx,ky) -> phi(kz,kx,ky)
  status = cufftExecZ2Z(planZ, phi_d, phi_d, CUFFT_FORWARD)
  if (status /= CUFFT_SUCCESS) write(*,*) 'Z forward error: ', status
  ! END of FFT3D forward



  block
    complex(8), device, pointer :: phi3d(:,:,:)
    real(8) :: k2
    !integer :: il, jl, ig, jg
    integer :: offsets(3), xoff, yoff
    integer :: np(3)
    np(piZ%order(1)) = piZ%shape(1)
    np(piZ%order(2)) = piZ%shape(2)
    np(piZ%order(3)) = piZ%shape(3)
    call c_f_pointer(c_devloc(phi_d), phi3d, piZ%shape)
 
    ! divide by -K**2, and normalize
    offsets(piZ%order(1)) = piZ%lo(1) - 1
    offsets(piZ%order(2)) = piZ%lo(2) - 1
    offsets(piZ%order(3)) = piZ%lo(3) - 1
 
    xoff = offsets(1)
    yoff = offsets(2)
    npx = np(1)
    npy = np(2)
    !$cuf kernel do (2)
    do jl = 1, npy
       do il = 1, npx
          jg = yoff + jl
          ig = xoff + il
          do k = 1, nz
             k2 = kx_d(ig)**2 + kx_d(jg)**2 + kx_d(k)**2
             phi3d(k,il,jl) = -phi3d(k,il,jl)/k2/(nx*ny*nz)
          enddo
       enddo
    enddo
 
    ! specify mean (corrects division by zero wavenumber above)
    if (xoff == 0 .and. yoff == 0) phi3d(1,1,1) = 0.0
 
   end block

   ! phi(kz,kx,ky) -> phi(z,kx,ky)
   status = cufftExecZ2Z(planZ, phi_d, phi_d, CUFFT_INVERSE)
  if (status /= CUFFT_SUCCESS) write(*,*) 'Z inverse error: ', status
  ! phi(z,kx,ky) -> phi(ky,z,kx)
  CHECK_CUDECOMP_EXIT(cudecompTransposeZToY(handle, grid_desc, phi_d, phi_d, work_d, CUDECOMP_DOUBLE_COMPLEX))
  ! phi(ky,z,kx) -> phi(y,z,kx)
  status = cufftExecZ2Z(planY, phi_d, phi_d, CUFFT_INVERSE)
  if (status /= CUFFT_SUCCESS) write(*,*) 'Y inverse error: ', status
  ! phi(y,z,kx) -> phi(kx,y,z)
  CHECK_CUDECOMP_EXIT(cudecompTransposeYToX(handle, grid_desc, phi_d, phi_d, work_d, CUDECOMP_DOUBLE_COMPLEX, [0,0,0], Pix%halo_extents))
  ! phi(kx,y,z) -> phi(x,y,z)
  status = cufftExecZ2Z(planX, phi_d, phi_d, CUFFT_INVERSE)
  if (status /= CUFFT_SUCCESS) write(*,*) 'X inverse error: ', status

 ! update halo nodes with pressure (needed for the pressure correction step)
 ! Update X-pencil halos in Y direction
 CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi_d, work_halo_d, CUDECOMP_DOUBLE_COMPLEX, piX%halo_extents, halo_periods, 2))

 ! Update X-pencil halos in Z direction
 CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi_d, work_halo_d, CUDECOMP_DOUBLE_COMPLEX, piX%halo_extents, halo_periods, 3))

 !D2H transfer
 phi = phi_d



 ! check against analytical solution
 block
   complex(8), pointer :: phi3(:,:,:)
   real(8) :: err, maxErr = -1.0, err2
   call c_f_pointer(c_loc(phi), phi3, [piX%shape(1), piX%shape(2), piX%shape(3)])

   !Check errro on complex (ua and phi3 are complex)
   do kl = 1, piX%shape(3)
      do jl = 1, piX%shape(2)
         do i = 1, piX%shape(1)
            err = abs(ua(i,jl,kl)-phi3(i,jl,kl))
            if (err > maxErr) maxErr = err
         enddo
      enddo
   enddo

   ! take back p from phi3
   ! i can span also the halo because they have been already updated
   do kl = 1, piX%shape(3)
      do jl = 1, piX%shape(2)
         do i = 1, piX%shape(1)
            p(i,jl,kl) = real(phi3(i,jl,kl))
         enddo
      enddo
   enddo

   write(*,"('[', i0, '] Max Error: ', e12.6)") rank, maxErr
  end block

  ! For debug
  ! output solution in a file
  write(namefile,'(a,i3.3,a)') 'p_',rank,'.dat'
  open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
  write(55) p
  close(55)

  !########################################################################################################################################
  ! END STEP 7: POISSON SOLVER FOR PRESSURE
  !########################################################################################################################################



  !########################################################################################################################################
  ! START STEP 8: VELOCITY CORRECTION
  !########################################################################################################################################
  ! 5.1 correct velocity 
  ! 5.3 Call halo exchnages along Y and Z for u,v,w
  !########################################################################################################################################
  ! END STEP 8: USTAR COMPUTATION 
  !########################################################################################################################################


call mpi_finalize(ierr)

end program main
