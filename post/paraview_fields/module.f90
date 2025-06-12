module commondata
 integer :: nx, ny, nz
 integer :: nstart, nend, dump
 double precision, parameter :: pi=3.14159265358979
 double precision :: dx,dy,dz
 integer :: uflag, vflag, wflag, phiflag, nfields
 double precision, allocatable, dimension(:) :: x, y, z
 real, allocatable, dimension(:,:,:) :: u, v, w, phi
end module commondata
