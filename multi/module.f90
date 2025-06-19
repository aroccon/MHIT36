module param
    integer, parameter :: nx=128
    integer :: ny=nx,nz=nx
    double precision :: pi,lx,dx,dxi,ddxi,rhoi,twopi
    integer :: restart,tstart,tfin,dump
    double precision :: gamma, normod
    double precision :: dt,mu,rho !flow parameters
    integer :: inflow, inphi
    double precision :: f1,f2,f3,k0 ! forcing parameters
    double precision :: radius, sigma, epsr, eps, pos, val, epsi, enum ! phase-field parameters
    double precision :: times,timef
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
   double precision, allocatable :: rhsu_o(:,:,:), rhsv_o(:,:,:), rhsw_o(:,:,:)
   complex(8), allocatable :: rhsp_complex(:,:,:)
   double precision, allocatable :: rhsp(:,:,:), p(:,:,:)
   double precision, allocatable :: div(:,:,:)
   double precision :: uc, vc, wc, umax, gumax=1.0d0, cou, alpha, beta
   double precision :: h11, h12, h13, h21, h22, h23, h31, h32, h33
   double precision :: umean, vmean, wmean, gumean, gvmean, gwmean
end module velocity


module phase
   double precision, allocatable :: phi(:,:,:), rhsphi(:,:,:), psidi(:,:,:)
   double precision, allocatable :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
   double precision, allocatable :: chempot(:,:,:), gradphix(:,:,:), gradphiy(:,:,:), gradphiz(:,:,:)
   double precision, allocatable :: fxst(:,:,:), fyst(:,:,:), fzst(:,:,:)
end module phase


! added NVTX for profiing from maxcuda/NVTX_example
module nvtx
use iso_c_binding
implicit none
integer,private :: col(7) = [ int(Z'0000ff00'), int(Z'000000ff'), int(Z'00ffff00'), int(Z'00ff00ff'), int(Z'0000ffff'), int(Z'00ff0000'), int(Z'00ffffff')]
character,private,target :: tempName(256)

type, bind(C):: nvtxEventAttributes
  integer(C_INT16_T):: version=1
  integer(C_INT16_T):: size=48 !
  integer(C_INT):: category=0
  integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
  integer(C_INT):: color
  integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
  integer(C_INT):: reserved0
  integer(C_INT64_T):: payload   ! union uint,int,double
  integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1 
  type(C_PTR):: message  ! ascii char
end type

interface nvtxRangePush
  ! push range with custom label and standard color
  subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
  use iso_c_binding
  character(kind=C_CHAR) :: name(256)
  end subroutine

  ! push range with custom label and custom color
  subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
  use iso_c_binding
  import:: nvtxEventAttributes
  type(nvtxEventAttributes):: event
  end subroutine
end interface

interface nvtxRangePop
  subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
  end subroutine
end interface

contains

subroutine nvtxStartRange(name,id)
  character(kind=c_char,len=*) :: name
  integer, optional:: id
  type(nvtxEventAttributes):: event
  character(kind=c_char,len=256) :: trimmed_name
  integer:: i

  trimmed_name=trim(name)//c_null_char

  ! move scalar trimmed_name into character array tempName
  do i=1,LEN(trim(name)) + 1
     tempName(i) = trimmed_name(i:i)
  enddo


  if ( .not. present(id)) then
    call nvtxRangePush(tempName)
  else
    event%color=col(mod(id,7)+1)
    event%message=c_loc(tempName)
    call nvtxRangePushEx(event)
  end if
end subroutine

subroutine nvtxEndRange
  call nvtxRangePop
end subroutine

end module nvtx

