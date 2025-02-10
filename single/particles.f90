
subroutine get_velocity
use particles
use velocity
use param
implicit none
integer :: m,i,j,k,one,ip,jp,kp,im,jm,km
double precision :: xt,yt,zt
double precision :: c1,c2,c3,c4,c5,c6,c7,c8

one=1.d0

! A. Roccon 04/03/2024
!trilinear interpolation 
!compute velocity at particle position
!$acc kernels
do m=1,np
  ! getting the cell where the particle is at
  i = floor(xp(m,1))
  j = floor(xp(m,2))
  k = floor(xp(m,3))
  ! getting fractional coordinates of the particle
  xt = xp(m,1) - dble(i)
  yt = xp(m,2) - dble(j)
  zt = xp(m,3) - dble(k)
  ! cpefficients for trilinear interpolation
  c1 = (one-xt) * (one-yt) * (one-zt)
  c2 = (xt)     * (one-yt) * (one-zt)
  c3 = (one-xt) * (yt)     * (one-zt)
  c4 = (xt)     * (yt)     * (one-zt)
  c5 = (one-xt) * (one-yt) * (zt)
  c6 = (xt)     * (one-yt) * (zt)
  c7 = (one-xt) * (yt)     * (zt)
  c8 = (xt)     * (yt)     * (zt)
  ! account for periodicity 
  ip=i+1
  jp=j+1
  kp=k+1
  im=i-1
  jm=j-1
  km=k-1
  if (ip .gt. nx) ip=1
  if (jp .gt. nx) jp=1
  if (kp .gt. nx) kp=1
  if (im .lt. 1) im=nx
  if (jm .lt. 1) jm=nx
  if (km .lt. 1) km=nx
! Trilinear interpolation (double checkm, from chatGPT and HIT36(
  ufp(m,1)= c1*u(i,j,k)  + c2*u(ip,j,k)  + c3*u(i,jp,k)  + c4*u(ip,jp,k) +&
          + c5*u(i,k,jp) + c6*u(ip,j,kp) + c7*u(i,jp,kp) + c8*u(ip,jp,kp)
  ufp(m,2)= c1*v(i,j,k)  + c2*v(ip,j,k)  + c3*v(i,jp,k)  + c4*v(ip,jp,k) +&
          + c5*v(i,k,jp) + c6*v(ip,j,kp) + c7*v(i,jp,kp) + c8*v(ip,jp,kp)
  ufp(m,3)= c1*w(i,j,k)  + c2*w(ip,j,k)  + c3*w(i,jp,k)  + c4*w(ip,jp,k) +&
          + c5*w(i,k,jp) + c6*w(ip,j,kp) + c7*w(i,jp,kp) + c8*w(ip,jp,kp)
enddo
!$acc end kernels
end subroutine










subroutine move_part
use openacc
use param
use particles
implicit none
integer :: m,n

! Tracers, advance according to local fluid velocity
if (ptype .eq. 1) then
  do m=1,np
    xp(m,1)=xp(m,1) + dt*ufp(m,1)
    xp(m,2)=xp(m,2) + dt*ufp(m,2)
    xp(m,3)=xp(m,3) + dt*ufp(m,3)
  enddo
endif

! Inertial particles
if (ptype .eq. 2) then
! to be implemented
  !do m=1,np
    ! compute forces on particles
    !fp(m,1)=0.d0
    !fp(m,2)=0.d0
    !fp(m,3)=0.d0
    ! update velocity
    !vp(m,1)=vp(m,1) + dt*fp(m,1)
    !vp(m,2)=vp(m,2) + dt*fp(m,2)
    !vp(m,3)=vp(m,3) + dt*fp(m,3)
    ! update position
    !xp(m,1)=xp(m,1) + dt*vp(m,1)
    !xp(m,2)=xp(m,2) + dt*vp(m,2)
    !xp(m,3)=xp(m,3) + dt*vp(m,3)
  !enddo
endif

!check for periodicity
!$acc kernels
do m=1,np
    !particle beyond domain extension
    if (xp(m,1) .gt. lx) xp(m,1)=xp(m,1)-lx
    if (xp(m,2) .gt. lx) xp(m,2)=xp(m,2)-lx
    if (xp(m,3) .gt. lx) xp(m,3)=xp(m,3)-lx
    ! opposite case
    if (xp(m,1) .lt. 0.d0) xp(m,1)=xp(m,1)+lx
    if (xp(m,2) .lt. 0.d0) xp(m,2)=xp(m,2)+lx
    if (xp(m,3) .lt. 0.d0) xp(m,3)=xp(m,3)+lx
enddo 
!$acc end kernels

end subroutine

