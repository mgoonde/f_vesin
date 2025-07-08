program main
  use f_vesin
  use, intrinsic :: iso_c_binding
  implicit none
  integer :: nat
  integer, allocatable :: typ(:)
  real(c_double), allocatable :: pos(:,:)
  real(c_double) :: lat(3,3)
  integer :: i
  character(len=256) :: line
  integer :: n_begin, n_end
  integer :: ierr
  character(:), allocatable :: errmsg
  type( VesinOptions ) :: options
  type( VesinNeighborList ) :: neigh
  type( f_VesinNeighborList ) :: fneigh

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

  options% cutoff = 4.2_c_double
  options% full = .true.
  options% sorted = .true.
  options% return_vectors = .true.
  ierr = vesin_compute( nat, pos, lat, options, neigh, errmsg, fneigh )
  if( ierr /= 0 ) then
     write(*,*) errmsg
  end if
  call vesin_free( neigh )

  write(*,*) "n pairs:", fneigh% length

  do i = 1, fneigh% length
     write(*,*) fneigh% pairs(:,i), fneigh% vectors(:,i)
  end do


  call fneigh% destroy()
  deallocate( typ, pos )
end program main
