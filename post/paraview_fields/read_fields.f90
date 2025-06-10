subroutine read_fields(nstep)

use commondata

integer :: nstep
character(len=40) :: namedir,namefile
character(len=6) :: numfile
character(len=3) :: setnum
logical :: check



namedir='../results/'
write(numfile,'(i6.6)') nstep



! check if u file exists; if u exists we assume that also v and w exist
namefile=trim(namedir)//'u_'//numfile//'.dat'
!write(*,*) namefile
inquire(file=trim(namefile),exist=check)


allocate(u(nx,ny,nz))
allocate(v(nx,ny,nz))
allocate(w(nx,ny,nz))

if(check.eqv..true.)then
write(*,*) 'Reading step ',nstep,' out of ',nend,' , flow'
  !reading u
  namefile=trim(namedir)//'u_'//numfile//'.dat'
  open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(666) u
  close(666,status='keep')
  !reading v
  namefile=trim(namedir)//'v_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) v
  close(667,status='keep')
  !reading w
  namefile=trim(namedir)//'w_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) w
  close(668,status='keep')
  ! generate paraview output file
  call generate_output(nstep)
endif

deallocate(u)
deallocate(v)
deallocate(w)

return
end
