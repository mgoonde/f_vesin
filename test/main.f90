program main
  use m_neighbour
  implicit none
  integer :: nat
  integer, allocatable :: typ(:)
  real(rp), allocatable :: pos(:,:)
  real(rp) :: lat(3,3)
  integer :: i, n
  character(len=256) :: line
  integer :: n_begin, n_end
  integer :: ierr
  type( t_neighbour ) :: neigh
  integer, allocatable :: neigh_ityp(:)
  real(rp), allocatable :: neigh_coords(:,:)

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

  ! initialize
  neigh = t_neighbour()

  ! compute neighbor list
  ierr = neigh% compute( nat, typ, pos, lat, 3.3_rp )
  if( ierr /= 0 ) then
     write(*,*) neigh% errmsg
     error stop
  end if


  ! retrieve neighbor data fo some index
  n = neigh% get( 2, nshell=3, ityplist=neigh_ityp, veclist=neigh_coords, include_idx=.true. )
  if( n .lt. 0 ) then
     write(*,*) neigh% errmsg
     error stop
  end if


  write(*,*) n
  write(*,*)
  do i = 1, n
     write(*,"(2x,i3,3x,3(g12.6,:,1x))") neigh_ityp(i), neigh_coords(:,i)
  end do

  call neigh% destroy()
  deallocate( typ, pos )
  deallocate( neigh_ityp, neigh_coords )
end program main
