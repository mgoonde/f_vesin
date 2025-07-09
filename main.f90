program main
  use f_vesin_wrapper
  implicit none
  integer :: nat
  integer, allocatable :: typ(:)
  real(rp), allocatable :: pos(:,:)
  real(rp) :: lat(3,3)
  integer :: i
  character(len=256) :: line
  integer :: n_begin, n_end
  integer :: ierr
  type( vesin_t ), pointer :: neigh

  read(*,*) nat
  read(*,'(a256)') line
  n_begin = index(line, "Lattice=") + 9
  line = line(n_begin:)
  n_end = index(line,'"')-1
  line = line(:n_end)
  read(line, *) lat
  allocate( typ(1:nat) )
  allocate( pos(1:3,1:nat))
  do i = 1, nat
     read(*,*) typ(i), pos(:,i)
  end do

  neigh => vesin_t( cutoff=4.2_rp, full=.true., return_vectors=.true. )
  ierr = neigh% compute( nat, pos, lat )
  if( ierr /= 0 ) write(*,*) neigh% errmsg

  write(*,*) neigh% length
  deallocate( neigh )
  deallocate( typ, pos )
end program main
