module param
    integer, parameter :: nx=128
    integer :: ny=nx,nz=nx
    double precision :: pi,lx,dx,dxi,ddxi,rhoi,twopi
    integer :: restart,tstart,tfin,dump
    double precision :: dt,mu,rho !flow parameters
    integer :: inflow
    double precision :: f1,f2,f3,k0 ! forcing parameters
    double precision :: radius, sigma, epsr, eps, pos ! phase-field parameters
    real(8) :: times,timef
end module param


module mpivar
   ! MPI variables
   integer :: rank, ranks, ierr
   integer :: localRank, localComm
end module mpivar


module cudecompvar
   use cudecomp
   integer :: npx, npy, npz
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
end module cudecompvar


module velocity
   double precision, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
   double precision, allocatable :: ustar(:,:,:), vstar(:,:,:), wstar(:,:,:)
   double precision, allocatable :: rhsu(:,:,:), rhsv(:,:,:), rhsw(:,:,:)
   complex(8), allocatable :: rhsp_complex(:,:,:)
   real(8), allocatable :: rhsp(:,:,:), p(:,:,:)
   real(8) :: uc,vc,wc,umax,cou
   real(8) :: h11,h12,h13,h21,h22,h23,h31,h32,h33
end module velocity


module phase
   double precision, allocatable :: phi(:,:,:), rhsphi(:,:,:)
   !double precision, allocatable :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
   !double precision, allocatable :: curv(:,:,:), gradphix(:,:,:), gradphiy(:,:,:), gradphiz(:,:,:)
   !double precision, allocatable :: fxst(:,:,:), fyst(:,:,:), fzst(:,:,:)
end module phase


