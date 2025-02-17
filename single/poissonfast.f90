subroutine poissonfast
use openacc
use cuFFT
use iso_c_binding
use fastp
use param
implicit none
integer :: gerr,i,j,k
double precision :: pm

! A. Roccon 28/02/2024
! 3D FFT-based solution of p_xx + p_yy  + p_zz = rhs;
! with periodic boundary conditions along all directions

! Laplacian matrix acting on the wavenumbers (in principle required only once)
!do i=1,nx
!    do j=1,nx
!        do k=1,nx
!            delsq(i,j,k) = -(kk(i)**2d0 + kk(j)**2d0 + kk(k)**2d0)
!        enddo
!    enddo
!enddo
! Laplacian matrix acting on the wavenumbers
! Avoids solving for the zero wavenumber
!delsq(1,1,1) = 1.d0
!write(*,*) "delsq(30,27,29)", delsq(30,27,29)

!Perform FFT3D forward of the rhsp
!$acc data copyin(rhsp) copyout(p) create(rhspc,pc)
!$acc host_data use_device(rhsp,rhspc)
gerr = gerr + cufftExecD2Z(cudaplan_fwd,rhsp,rhspc)
!$acc end host_data
!!$acc end data


!$acc kernels
do i=1,nx/2+1
    do j=1,nx
        do k=1,nx
            pc(i,j,k)=rhspc(i,j,k)/delsq(i,j,k)
        enddo
    enddo
enddo
!$acc end kernels

!!$acc data copyin(pc) copyout(p)
!$acc host_data use_device(pc,p)
gerr = gerr + cufftExecZ2D(cudaplan_bwd,pc,p)
!$acc end host_data
!$acc end data

! scale p
!$acc kernels
p = p / (nx*nx*nx)
pm=p(1,1,1)
p = p - pm
!$acc end kernels

return
end subroutine





subroutine init_cufft
use openacc
use cuFFT
use iso_c_binding
use fastp
use param
implicit none
integer :: gerr,i,j,k
integer(kind=int_ptr_kind()) :: workSize(1)

! create plans
! Creata plans (forth and back)
!Plan forward
gerr=0
gerr=gerr+cufftCreate(cudaplan_fwd)
gerr=gerr+cufftMakePlan3d(cudaplan_fwd,nx,nx,nx,CUFFT_D2Z,workSize)
if (gerr.ne.0) write(*,*) "Error in cuFFT plan FWD:", gerr

!Plan backward
gerr=0
gerr=gerr+cufftCreate(cudaplan_bwd)
gerr=gerr+cufftMakePlan3d(cudaplan_bwd,nx,nx,nx,CUFFT_Z2D,workSize)
if (gerr.ne.0) write(*,*) "Error in cuFFT plan BWD:", gerr


!create wavenumber
! wavenumbers 
do i=1,nx/2
    kk(i)=i-1
enddo
do i=nx/2+1,nx
    kk(i)=-nx + i -1
enddo

!create delsq
do i=1,nx
    do j=1,nx
        do k=1,nx
            delsq(i,j,k) = -(kk(i)**2d0 + kk(j)**2d0 + kk(k)**2d0)
        enddo
    enddo
enddo

delsq(1,1,1) = 1.d0

end subroutine
