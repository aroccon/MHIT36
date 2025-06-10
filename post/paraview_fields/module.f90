module commondata
 integer :: nx, ny, nz
 integer :: nstart,nend,dump
 double precision, parameter :: pi=3.14159265358979
 double precision :: dx,dy,dz
 double precision, allocatable, dimension(:) :: x,y,z
 real, allocatable, dimension(:,:,:) :: u,v,w
end module commondata
