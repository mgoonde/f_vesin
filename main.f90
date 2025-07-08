program main
  use f_vesin
  implicit none
  integer :: nat
  integer, allocatable :: typ(:)
  real, allocatable :: pos(:,:)
  real :: lat(3,3)
  integer :: i
  character(len=256) :: line
  integer :: n_begin, n_end
  integer :: ierr

  ! read(*,*) nat
  ! read(*,'(a256)') line
  ! n_begin = index(line, "Lattice=") + 9
  ! line = line(n_begin:)
  ! n_end = index(line,'"')-1
  ! line = line(:n_end)
  ! read(line, *) lat
  ! allocate( typ(1:nat) )
  ! allocate( pos(1:3,1:nat))
  ! do i = 1, nat
  !    read(*,*) typ(i), pos(:,i)
  ! end do

  ierr = vesin_compute()


end program main
