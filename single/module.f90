module param
    integer, parameter :: nx=128
    integer, parameter :: np=100
    double precision :: pi,lx,dx,dxi,ddxi,rhoi
    integer :: restart,tstart,tfin,dump
    double precision :: dt,mu,rho !flow parameters
    integer :: inflow
    double precision :: f1,f2,f3,k0 ! forcing parameters
    double precision :: radius, sigma, eps ! phase-field parameters
end module param


module fastp
    use param
    use cufft
    use iso_c_binding
    integer :: cudaplan_fwd,cudaplan_bwd
    double precision, allocatable :: delsq(:,:,:)
    double precision, allocatable :: kk(:)
    double precision, allocatable :: kx(:,:,:), ky(:,:,:), kz(:,:,:)
    real(c_double), pinned, allocatable :: p(:,:,:), rhsp(:,:,:)
    complex(c_double_complex), pinned, allocatable :: pc(:,:,:)
    complex(c_double_complex), pinned, allocatable :: rhspc(:,:,:)
end module fastp


module velocity
   double precision, allocatable :: div(:,:,:)
   double precision, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
   double precision, allocatable :: ustar(:,:,:), vstar(:,:,:), wstar(:,:,:)
   double precision, allocatable :: rhsu(:,:,:), rhsv(:,:,:), rhsw(:,:,:)
end module velocity


module phase
   double precision, allocatable :: phi(:,:,:), rhsphi(:,:,:)
   double precision, allocatable :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
   double precision, allocatable :: curv(:,:,:), gradphix(:,:,:), gradphiy(:,:,:), gradphiz(:,:,:)
   double precision, allocatable :: fxst(:,:,:), fyst(:,:,:), fzst(:,:,:)
end module phase


module particles
   double precision, allocatable :: xp(:,:), vp(:,:), ufp(:,:), fp(:,:)
   integer, parameter :: ptype=1 !(1=tracer, 2=inertial)
end module particles
