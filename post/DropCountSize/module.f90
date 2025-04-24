module param
 implicit none
 double precision, parameter :: pi=4.0d0*atan(1.0d0)
 double precision, parameter :: twopi=8.0d0*atan(1.0d0)
 double precision, parameter :: n3o4opi=3.0d0/(4.0d0*pi)
 integer, parameter:: clen=300 !char length
 integer :: nx
 integer :: begin, finish, step
 character(len=clen) :: rootpath
 double precision ::dx,lx,dxtom3
end module param
!*********************************************
module flowvars
 implicit none
 double precision, allocatable, dimension(:,:,:) :: phi
end module flowvars
!*********************************************