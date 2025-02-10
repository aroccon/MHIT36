
!##########################################################################
!###########################################################################
subroutine writefield(t,fieldn)
! Output field, file is written in the src/output folder 

use velocity
use phase
use fastp
implicit none
integer :: t,fieldn
character(len=40) :: namefile

! fieldn=1 means u
! fieldn=2 means v
! fieldn=3 means w
! fieldn=4 means p
! fieldn=5 means phi

if (fieldn .eq. 1) then
write(namefile,'(a,i8.8,a)') './output/u_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) u(:,:,:)
close(55)
endif

if (fieldn .eq. 2) then
write(namefile,'(a,i8.8,a)') './output/v_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) v(:,:,:)
close(55)
endif

if (fieldn .eq. 3) then
write(namefile,'(a,i8.8,a)') './output/w_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) w(:,:,:)
close(55)
endif

if (fieldn .eq. 4) then
write(namefile,'(a,i8.8,a)') './output/p_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) p(:,:,:)
close(55)
endif

if (fieldn .eq. 5) then
write(namefile,'(a,i8.8,a)') './output/phi_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) phi(:,:,:)
close(55)
endif
end subroutine








subroutine readfield(t,fieldn)
! Used in case of fresh start, file is read from the src/init folder 
use velocity
use phase
use fastp
implicit none
integer :: t,fieldn
character(len=40) :: namefile


! fieldn=1 means u
! fieldn=2 means v
! fieldn=3 means w
! fieldn=4 means p
! fieldn=5 means phi

if (fieldn .eq. 1) then
write(namefile,'(a,i8.8,a)') './init/u.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) u
close(55)
endif

if (fieldn .eq. 2) then
write(namefile,'(a,i8.8,a)') './init/v.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) v
close(55)
endif

if (fieldn .eq. 3) then
write(namefile,'(a,i8.8,a)') './init/w.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) w
close(55)
endif

if (fieldn .eq. 4) then
write(namefile,'(a,i8.8,a)') './init/p.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) p
close(55)
endif

if (fieldn .eq. 5) then
write(namefile,'(a,i8.8,a)') './init/phi.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) phi
close(55)
endif
end subroutine






subroutine readfield_restart(t,fieldn)
! Used in case of restart, file is read from the src/output folder (iteration tstart must be present!)
use velocity
use phase
use fastp
implicit none
integer :: t,fieldn
character(len=40) :: namefile

! fieldn=1 means u
! fieldn=2 means v
! fieldn=3 means w
! fieldn=4 means p
! fieldn=5 means phi

if (fieldn .eq. 1) then
write(*,*) "Read u"
write(namefile,'(a,i8.8,a)') './output/u_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) u
close(55)
endif

if (fieldn .eq. 2) then
write(*,*) "Read v"
write(namefile,'(a,i8.8,a)') './output/v_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) v
close(55)
endif

if (fieldn .eq. 3) then
write(*,*) "Read w"
write(namefile,'(a,i8.8,a)') './output/w_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) w
close(55)
endif

if (fieldn .eq. 4) then
write(namefile,'(a,i8.8,a)') './output/p_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) p
close(55)
endif

if (fieldn .eq. 5) then
write(namefile,'(a,i8.8,a)') './output/phi_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) phi
close(55)
endif

end subroutine


















subroutine writepart(t)
use particles
implicit none
integer :: t
character(len=40) :: namefile

!Output for tracers
if (ptype .eq. 1) then
  write(namefile,'(a,i8.8,a)') './output/xp_',t,'.dat'
  open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
  write(55) xp(:,:)
  close(55)
endif

!Output for inertial tracers
if (ptype .eq. 2) then
  write(namefile,'(a,i8.8,a)') './output/xp_',t,'.dat'
  open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
  write(55) xp(:,:)
  close(55)

  write(namefile,'(a,i8.8,a)') './output/vp_',t,'.dat'
  open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
  write(55) vp(:,:)
  close(55)

  !write(namefile,'(a,i8.8,a)') './output/ufp_',t,'.dat'
  !open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
  !write(55) ufp(:,:)
  !close(55)
endif

end subroutine
