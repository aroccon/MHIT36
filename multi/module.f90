module param
    integer, parameter :: nx=64
    integer :: ny=nx,nz=nx
    double precision :: pi,lx,dx,dxi,ddxi,rhoi,twopi
    integer :: restart,tstart,tfin,dump
    double precision :: dt,mu,rho !flow parameters
    integer :: inflow
    double precision :: f1,f2,f3,k0 ! forcing parameters
    double precision :: radius, sigma, eps ! phase-field parameters
    real(8) :: times,timef
end module param


module mpivar
   ! MPI variables
   integer :: rank, ranks, ierr
   integer :: localRank, localComm
end module mpivar


module velocity
   double precision, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
   complex(8), allocatable :: rhsp_complex(:,:,:)
   real(8), allocatable :: rhsp(:,:,:), p(:,:,:)
end module velocity


module phase
   double precision, allocatable :: phi(:,:,:), rhsphi(:,:,:)
   !double precision, allocatable :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
   !double precision, allocatable :: curv(:,:,:), gradphix(:,:,:), gradphiy(:,:,:), gradphiz(:,:,:)
   !double precision, allocatable :: fxst(:,:,:), fyst(:,:,:), fzst(:,:,:)
end module phase


