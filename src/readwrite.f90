subroutine writefield(t,fieldn)

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