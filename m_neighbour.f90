module m_neighbour

  use, intrinsic :: iso_c_binding

  ! import the C-interoperable vesin interface
  use vesin, only: &
       c_VesinOptions => VesinOptions, &
       c_VesinNeighborList => VesinNeighborList, &
       c_vesin_neighbors => vesin_neighbors, &
       c_vesin_free => vesin_free, &
       VesinUnknownDevice, VesinCPU

  implicit none

  private
  public :: rp
  public :: t_neighbour

  ! precision, modify as needed
  integer, parameter :: rp = c_double


  !> A slightly more sophisticated fortran wrapper to `vesin`.
  !!
  !! Use as:
  !!~~~~~~~~~~{.f90}
  !!
  !!    use m_neighbour
  !!    implicit none
  !!    type( t_neighbour ) :: neigh
  !!    integer :: ierr, n
  !!    integer, allocatable :: neig_idx(:), neig_ityp(:)
  !!    real(rp), allocatable :: neig_coords(:,:)
  !!
  !!    ! initialize
  !!    neigh = t_neighbour()
  !!
  !!    ! compute neighbor lists with some value of `rcut`
  !!    ierr = neigh% compute( nat, ityp, pos, lat, rcut )
  !!    ! error
  !!    if( ierr /= 0 ) then
  !!       write(*,*) neigh% errmsg
  !!       stop
  !!    end if
  !!
  !!    ! Get neighbors of `idx=5`, expanded to 2 neighbor shells.
  !!    ! Output the list into `neig_idx`,
  !!    ! the atomic types into `neig_ityp`,
  !!    ! the vectors to `neig_coords`, and
  !!    ! include the original `idx` in the lists.
  !!    n = neigh% get( idx = 5, nbond = 2, include_idx = .true., &
  !!                    list     = neig_idx, &
  !!                    ityplist = neig_ityp, &
  !!                    veclist  = neig_coords )
  !!    ! error
  !!    if( n < 0 ) then
  !!       write(*,*) neigh% errmsg
  !!       stop
  !!    end if
  !!
  !!    ! set an arbitrary list of atomic indices
  !!    list = [ 4, 12, 2, 8, 25 ]
  !!    ! expand the list by 2 neighbor shells
  !!    n = neigh% expand( 2, list=list )
  !!    write(*,*) "expanded list:", list
  !!
  !!    ! compute the list with a different `rcut=rcut_2`
  !!    ierr = neigh% compute( nat, ityp, pos, lat, rcut_2 )
  !!    ! error
  !!    if( ierr /= 0 ) then
  !!       write(*,*) neigh% errmsg
  !!       stop
  !!    end if
  !!
  !!    ! destroy neigh
  !!    call neigh% destroy()
  !!
  !!~~~~~~~~~~
  type, public :: t_neighbour

     ! ----- private ------
     private
     ! vesin stuffs
     type( c_vesinOptions ) :: opts       !< Computation options
     type( c_vesinNeighborList ) :: cdata !< Returned C data
     integer :: device = VesinCPU
     logical :: active = .false.          !< .true. when initialized
     ! pointers to C data, in C precision
     integer(c_size_t) :: length = 0_c_size_t               !< size of the list
     integer( c_size_t ), pointer :: pairs(:,:) => null()   !< shape[2, length]
     integer( c_int32_t ), pointer :: shifts(:,:) => null() !< shape[3, length]
     real( c_double ), pointer :: distances(:) => null()    !< shape[length]
     real( c_double ), pointer :: vectors(:,:) => null()    !< shape[3, length]

     ! local
     integer, allocatable :: cumsum(:) !< cumulative sum of nneig
     ! real(rp) :: lat(3,3) !< lattice for knowing the shifts in expand
     integer, allocatable :: ityp(:) !< copy of ityp



     ! ----- public ------
     character(:), allocatable, public :: errmsg !< error message


   contains
     procedure, public :: compute => t_neighbour_compute
     procedure, public :: get_nn  => t_neighbour_get_nn
     procedure, public :: get     => t_neighbour_get
     procedure, public :: expand  => t_neighbour_expand
     procedure, public :: cluster => t_neighbour_cluster
     ! procedure, public :: get_by_rcut
     ! final :: t_neighbour_destroy
     procedure, public :: destroy => t_neighbour_destroy
  end type t_neighbour



  interface t_neighbour
     procedure :: t_neighbour_construct
  end interface t_neighbour


contains


  function t_neighbour_construct( return_distances )result( self )
    !! call as:
    !!
    !!~~~~~~~~{.f90}
    !! type( t_neighbour ) :: neigh
    !!
    !! neigh = t_neighbour()
    !!~~~~~~~~
    implicit none
    ! type( t_neighbour ), pointer :: self
    type( t_neighbour ) :: self
    logical, intent(in), optional :: return_distances

    ! allocate( t_neighbour::self )

    ! some vesin options
    self%opts% full             = .true.
    self%opts% sorted           = .true.
    self%opts% return_shifts    = .true.
    self%opts% return_distances = .false.
    self%opts% return_vectors   = .true.
    if( present(return_distances)) self%opts% return_distances = return_distances

    ! set as active
    self% active = .true.

  end function t_neighbour_construct


  subroutine t_neighbour_destroy( self )
    !! destroy the `t_neighbour` instance
    implicit none
    ! type( t_neighbour ), intent(inout) :: self
    class( t_neighbour ), intent(inout) :: self
    if( associated(self%pairs))    nullify(self%pairs)
    if( associated(self%shifts))   nullify(self%shifts)
    if( associated(self%distances))nullify(self%distances)
    if( associated(self%vectors))  nullify(self%vectors)
    call c_vesin_free(self%cdata)
    self% active = .false.
    if(allocated(self%cumsum))deallocate(self%cumsum)
    if(allocated(self%ityp))  deallocate(self%ityp)
  end subroutine t_neighbour_destroy



  function t_neighbour_compute( self, nat, typ, pos, lat, rcut ) result( ierr )
    !! compute the neighbor list using `vesin` algorithm.
    implicit none
    class( t_neighbour ), intent(inout) :: self

    !> number of atoms
    integer, intent(in) :: nat

    !> integer atomic types
    integer, intent(in) :: typ(nat)

    !> atomic positions, shape(3,nat)
    real(rp), intent(in) :: pos(3,nat)

    !> lattice vectors in rows
    real(rp), intent(in) :: lat(3,3)

    !> distance cutoff in the same units as `pos`
    real(rp), intent(in) :: rcut

    !> nonzero on error
    integer :: ierr

    logical(c_bool) :: periodic
    type( c_ptr ) :: c_errmsg
    integer :: n, i, idx

    ierr = -1
    if( .not. self% active )then
       self% errmsg = "compute:: t_neighbour instance not initialized."
       return
    end if

    ! set cutoff
    self%opts% cutoff = real( rcut, c_double )
    ! periodic
    periodic = logical( .true., c_bool )

    ! save typ
    allocate( self% ityp, source=typ )
    ! save lat
    ! self% lat = lat

    ierr = int( c_vesin_neighbors( real(pos, c_double), int(nat, c_size_t), &
         real(lat, c_double), periodic, self%device, self%opts, self%cdata, c_errmsg ))
    !
    if( ierr /= 0 ) then
       self%errmsg = "compute:: vesin error:: "//c2f_string(c_errmsg)
       return
    end if

    ! cast c data to pointers
    self% length = self%cdata%length
    if(associated(self%pairs))    nullify(self%pairs)
    if(associated(self%shifts))   nullify(self%shifts)
    if(associated(self%distances))nullify(self%distances)
    if(associated(self%vectors))  nullify(self%vectors)
    !
    n = int( self%cdata%length )
    if(c_associated(self%cdata%pairs))     call c_f_pointer(self%cdata%pairs, self%pairs, shape=[2,n])
    if(c_associated(self%cdata%shifts))    call c_f_pointer(self%cdata%shifts, self%shifts, shape=[3,n])
    if(c_associated(self%cdata%distances)) call c_f_pointer(self%cdata%distances, self%distances, shape=[n])
    if(c_associated(self%cdata%vectors))   call c_f_pointer(self%cdata%vectors, self%vectors, shape=[3,n])

    ! compute the cumsums
    allocate( self% cumsum(0:nat))
    self%cumsum(0)=0
    n = 0
    do idx = 1, nat
       ! count how many neighbors i has; self%pairs is C format, i.e. indices start with 0
       n = n + count( self%pairs(1,:) .eq. idx-1 )
       self% cumsum(idx) = n
    end do

  end function t_neighbour_compute



  function t_neighbour_get( self, idx, list, ityplist, veclist, shiftslist, nbond, include_idx )result(n)
    !! Get the neighbor data of `idx` from neighbour list, up to `nbond` neighbor shells.
    !! Return `n` which is the number of neighbors, or negative on error.
    !! The vector of atom `idx` is NOT included in output by default.
    implicit none
    class( t_neighbour ), intent(inout) :: self

    !> input atomic index
    integer, intent(in)                 :: idx

    !> output list of neighbours to atom `idx`
    integer, allocatable, intent(out), optional  :: list(:)

    !> output list of atomic types neighbor to `idx`
    integer, allocatable, intent(out), optional :: ityplist(:)

    !> output list of vectors neighbour to `idx`
    real(rp), allocatable, intent(out), optional :: veclist(:,:)

    !> output list of box shifts for each neighbor
    integer, allocatable, intent(out), optional :: shiftslist(:,:)

    !> how many bond shells to get (default=1). If `nbond>1`, the output list is not sorted
    integer, intent(in), optional :: nbond

    !> flag to include the atom `idx` in output. If .true., it will be on the first
    !! element of output. Default=.false.
    logical, intent(in), optional :: include_idx

    !> `n`, number of neighbors; if `idx` is invalid, or the
    !! neighbour list has not been computed `n=-1`
    integer :: n

    integer :: nb, i, ndim
    integer, allocatable :: inlist(:), s_inlist(:,:)
    real(rp), allocatable :: v_inlist(:,:)
    logical :: inc_idx

    if(  idx .le. 0 .or. &
         idx .gt. size(self% cumsum) .or. &
         .not.allocated(self%cumsum) .or. &
         .not. self%active) then
       self% errmsg = "get:: idx out of range, or neighbor list not computed."
       ! idx out of range, or neiglist not computed
       n = -1
       return
    end if

    inc_idx = .false.
    if( present(include_idx))inc_idx = include_idx

    nb = 1
    if(present(nbond))nb=nbond
    if( nb < 1 ) then
       n = 0
       return
    end if


    ! create array with just idx
    allocate(inlist(1:1), source=idx )
    ! vecs
    ndim = size(self% vectors,1)
    allocate(v_inlist(1:ndim,1:1) )
    v_inlist(:,:) = 0.0_rp
    ! shifts
    allocate( s_inlist(1:ndim,1:1))
    s_inlist(:,:) = 0

    ! expand it by nb
    n = self% expand( nb, inlist, veclist=v_inlist, shiftslist=s_inlist )

    ! check error
    if( n .lt. 0 ) then
       self%errmsg = self%errmsg//"::get"
       return
    end if


    if(present(veclist)) then
       if( inc_idx ) then
          ! include `idx`
          allocate( veclist, source=v_inlist )
       else
          ! first index of `inlist` is idx, remove it for output
          allocate( veclist(1:ndim,1:n-1))
          do i = 1, n-1
             veclist(:,i) = v_inlist(:,i+1)
          end do
       end if
    end if

    if( present(list)) then
       if( inc_idx ) then
          ! include `idx`
          allocate( list, source=inlist )
       else
          ! remove `idx` which is first element from output list
          allocate(list(1:n-1))
          do i = 1, n-1
             list(i) = inlist(i+1)
          end do
       end if
    end if

    if( present(ityplist)) then
       if( inc_idx ) then
          allocate(ityplist(1:n))
          do i = 1, n
             ityplist(i) = self% ityp( inlist(i) )
          end do
       else
          allocate(ityplist(1:n-1))
          do i = 1, n-1
             ityplist(i) = self% ityp( inlist(i+1) )
          end do
       end if
    end if

    if( present(shiftslist)) then
       if( inc_idx) then
          allocate( shiftslist, source=s_inlist )
       else
          ! skip first idx
          allocate( shiftslist(1:3,1:n-1) )
          do i = 1, n-1
             shiftslist(:,i) = s_inlist(:,i+1)
          end do
       end if
    end if

    if( .not. inc_idx ) n = n - 1

  end function t_neighbour_get


  function t_neighbour_get_nn( self, idx, list, ityplist, veclist, shiftslist, include_idx )result(n)
    !! Get the first neighbor shell of `idx` from neighbour list.
    !! Return `n` which is the number of neighbors.
    !! The vector of atom `idx` is NOT included in output by default.
    implicit none
    class( t_neighbour ), intent(inout) :: self

    !> input atomic index
    integer, intent(in)                 :: idx

    !> output list of neighbours to atom `idx`
    integer, allocatable, intent(out), optional   :: list(:)

    !> output list of atomic types neighbours to `idx`
    integer, allocatable, intent(out), optional :: ityplist(:)

    !> output list of vectors neighbour to `idx`, which is assumed at [0.0, 0.0, 0.0]
    real(rp), allocatable, intent(out), optional  :: veclist(:,:)

    !> output list of vector shifts (periodic images), wrt `idx`
    integer, allocatable, intent(out), optional :: shiftslist(:,:)

    !> flag to include the atom `idx` in output. If .true., it will be on the first
    !! element of output. Default=.false.
    logical, intent(in), optional :: include_idx

    !> `n`, number of first neighbors; if `idx` is invalid, or the
    !! neighbour list has not been computed `n=-1`
    integer :: n


    integer :: i_s, i_e
    logical :: inc_idx
    integer :: ndim, i, midx

    if(  idx .le. 0 .or. &
         idx .gt. size(self%cumsum)-1 .or. &
         .not.allocated(self%cumsum) .or. &
         .not.self%active ) then
       self% errmsg="get_nn:: idx out of bounds, or neighbor list not computed."
       n = -1
       return
    end if

    ! do not include central idx by default
    inc_idx = .false.
    if(present(include_idx))inc_idx=include_idx

    ! starting and ending idx
    i_s = self%cumsum(idx-1)+1
    i_e = self%cumsum(idx)

    n = i_e - i_s + 1

    ndim=3
    if( present(list)) then
       if( inc_idx ) then
          ! including idx, add it to start
          allocate( list(1:n+1) )
          list(1) = idx
          list(2:) = int( self% pairs( 2, i_s:i_e) ) + 1
       else
          ! not including idx
          allocate( list(1:n) )
          list(:) = int( self% pairs(2, i_s:i_e) ) + 1
       end if
    end if

    if(present(ityplist)) then
       if( inc_idx )then
          ! include idx at start
          allocate(ityplist(1:n+1))
          ityplist(1) = self% ityp(idx)
          do i = 1, n
             midx = int( self% pairs( 2, i_s + i - 1 ) ) + 1
             ityplist(1+i) = self% ityp(midx)
          end do
       else
          allocate(ityplist(1:n))
          do i = 1, n
             midx = int( self% pairs( 2, i_s + i - 1) ) + 1
             ityplist(i) = self% ityp(midx)
          end do
       end if
    end if

    if(present(veclist)) then
       if( inc_idx ) then
          ! include idx at start
          allocate( veclist(1:ndim, 1:n+1) )
          veclist(:,1) = 0.0
          veclist(:,2:) = real( self% vectors(:, i_s:i_e ), kind(veclist) )
       else
          allocate( veclist(1:ndim,1:n) )
          veclist(:,:) = real( self%vectors(:, i_s:i_e), kind(veclist))
       end if
    end if

    if(present(shiftslist)) then
       if( inc_idx ) then
          allocate( shiftslist(1:ndim,1:n+1) )
          shiftslist(:,1) = 0
          shiftslist(:,2:) = int( self%shifts(:,i_s:i_e), kind(shiftslist) )
       else
          allocate( shiftslist(1:ndim,1:n) )
          shiftslist(:,:) = int( self%shifts(:,i_s:i_e), kind(shiftslist))
       end if
    end if

    if( inc_idx ) n = n + 1

  end function t_neighbour_get_nn


  function t_neighbour_expand( self, nbond, list, ityplist, veclist, shiftslist ) result( n )
    !! expand the list in input by `nbond` number of bond shells.
    !! NOTE: pbc are not taken into account here, the expansion goes beyond the cell,
    !! which means the same atom index can repeat multiple times on different positions,
    !! as copies of itself by periodicity.
    implicit none
    class( t_neighbour ), intent(inout) :: self

    !> number of bons shells to expand
    integer, intent(in) :: nbond

    !> Assumed allocated on input, containing list to be expanded;
    !! on output: expanded list
    integer, allocatable, intent(inout) :: list(:)

    !> Optionally allocated on input, containing atomic types list;
    !! modified to include the atomic types of expansion on output
    integer, allocatable, intent(inout), optional :: ityplist(:)

    !> Assumed allocated on input, containing veclist;
    !! on output modified to include the expansion vectors
    real(rp), allocatable, intent(inout), optional :: veclist(:,:)

    !> array of shifts
    integer, allocatable, intent(inout), optional :: shiftslist(:,:)

    !> number of elements in output list
    integer :: n
    !
    integer :: i, ntot, ncur, idx, iatm, nprev, nat
    integer, allocatable :: work(:), swork(:,:)
    real(rp), allocatable :: vwork(:,:), vtmp(:,:)
    integer, parameter :: batchsize=50
    integer :: ibond, ii, nl_idx, nn, jj, j
    integer, allocatable :: l_idx(:)
    real(rp), allocatable :: vl_idx(:,:)
    integer, allocatable :: mybox(:,:), stmp(:,:), sl_box(:,:)
    integer :: ibox(3), jbox(3), mbox(3)
    real(rp) :: jvec(3), ivec(3), center(3), origin(3)

    integer, allocatable :: dum(:), dum2(:,:)
    real(rp), allocatable :: rdum2(:,:)
    if(.not.allocated(self%cumsum)) then
       ! neiglist not computed
       self% errmsg="expand:: neighbor list not computed."
       n = -1
       return
    end if

    ! check indices in input list
    nat = size(self%cumsum)
    if( any(list .le. 0) .or. any(list .gt. nat) ) then
       self%errmsg="expand:: some index out of bounds in input list."
       n = -1
       return
    end if

    ! size of input list
    n = size(list)

    ntot = n+batchsize
    ! allocate larger list for work
    allocate( work(1:ntot))
    work(1:n) = list
    ! boxes
    allocate( mybox(1:3,1:ntot), source=0)
    if( present(shiftslist)) mybox(1:3,1:n) = shiftslist(:,:)
    ! vecs
    allocate( vwork(1:3,1:ntot) )
    vwork(:,:) = 0.0_rp
    if( present(veclist)) vwork(:,1:n) = veclist


    ! do i = 1, self%length
    !    write(*,"('idx',i4,2x,'pair',i4,2x,'box',3i3)") i, self%pairs(2,i)+1, self%shifts(:,i)
    ! end do
    ! write(*,*)

    ! current size
    ncur = n
    ! first index
    nprev = 1

    do ibond = 1, nbond
       ! write(*,*) ">> ibond",ibond
       ! current size
       nn = n
       do iatm = nprev, nn
         ! this idx
         idx = work(iatm)
         ! mybox
         mbox = mybox(:,iatm)
         ! origin vec
         origin = vwork(:,iatm)
         ! neiglist of this idx
         nl_idx = self% get_nn( idx, list=l_idx, shiftslist=sl_box, veclist=vl_idx )
         if( nl_idx .lt. 0 ) then
            self%errmsg = self%errmsg//"::expand"
            n = -1
            return
         end if
         !
         ! write(*,"('iatm',i4,1x,'idx',i4,2x,'mybox',3i4)") iatm,idx, mbox
         ! do i = 1, nl_idx
         !     write(*,"(i4,2x,3i4)") l_idx(i), sl_box(:,i)+mbox
         ! end do
         !
         ! check if any l_idx already in work, and in same box:
         ! for each found idx
         jatm_: do j = 1, nl_idx
             !
             jj = l_idx(j)
             jbox = mbox + sl_box(:,j)
             jvec = vl_idx(:,j) + origin
             ! check over all work
             do i = 1, n
                !
                ii = work(i)
                ibox = mybox(:,i)
                ivec = vwork(:,i)
                ! equal idx and equal box == equal atom already in work, skip j
                if( ii.eq.jj .and. all(ibox.eq.jbox) ) then
                   ! write(*,"('ii',i2,2x,'jj',i2,4x,3i3,2x,3i3)") ii, jj, ibox, jbox
                   ! write(*,*) "equal, cyclein j"
                   cycle jatm_
                end if
             end do
             !
             ! check alloc size
             if( n+1 .gt. size(work) ) then
                call move_alloc(work, dum)
                allocate( work(1:n+batchsize))
                work(1:n) = dum(:)
                deallocate(dum)
                !
                call move_alloc(mybox, dum2)
                allocate( mybox(1:3,1:n+batchsize))
                mybox(1:3,1:n) = dum2(:,:)
                deallocate(dum2)
                !
                call move_alloc(vwork, rdum2)
                allocate( vwork(1:3,1:n+batchsize))
                vwork(1:3,1:n) = rdum2(:,:)
                deallocate(rdum2)
             end if
             ! add jj
             n = n + 1
             work(n) = jj
             mybox(:,n) = jbox
             vwork(:,n) = jvec
             ! write(*,"('add iatm',1x,i4,1x,i4,2x,3i4,2x,3f7.3)") n, jj, jbox, jvec
          end do jatm_
       end do ! iatm
       !
       ! next first index
       nprev = nn + 1

       ! write(*,*) "work:", n
       ! do i = 1, n
       !    write(*,"(i4,3x,3f9.4,3x,3i3)") work(i), vwork(:,i), mybox(:,i)
       ! end do

    end do ! ibond


    ! write(*,*) "work:"
    ! write(*,*) n
    ! write(*,*) "properties=species:I:1:pos:R:3:id:I:1"
    ! do i = 1, n
    !    write(*,"(i4,3x,3f9.4,1x,i4,3x,3i3)") self%ityp(work(i)), vwork(:,i), work(i), mybox(:,i)
    ! end do


    ! set output
    ncur = n
    deallocate(list)
    allocate(list(1:ncur))
    list(:) = work(1:ncur)

    if( present(veclist)) then
       deallocate(veclist)
       allocate(veclist(1:size(vwork,1),1:ncur))
       veclist(:,:) = vwork(:,1:ncur)
       center = veclist(:,1)
       do i = 1, n
          veclist(:,i) = veclist(:,i) - center
       end do
    end if

    if( present(ityplist))then
       if(allocated(ityplist))deallocate(ityplist)
       allocate(ityplist(1:ncur))
       do i = 1, ncur
          ityplist(i) = self% ityp( work(i) )
       end do
    end if

    if(present(shiftslist))then
       if(allocated(shiftslist))deallocate(shiftslist)
       allocate(shiftslist(1:3,1:ncur))
       shiftslist(:,:) = mybox(:,1:ncur)
    end if

  end function t_neighbour_expand


  function t_neighbour_cluster( self, idx, list ) result( n )
    !! check if the input `list` is a cluster graph, and return only the complete graph
    !! around `idx`. The order of indices might be lost on output.
    !! NOTE: works purely on indices, no regard for pbc or not.
    implicit none

    class( t_neighbour ), intent(inout) :: self

    !> starting atom index
    integer, intent(in) :: idx

    !> on input the possibly cluster-graph to check, on output the complete cluster
    !! graph around `idx`
    integer, allocatable, intent(inout) :: list(:)

    !> Number of atoms in output cluster.
    !! If `idx` is not present in `list`, return zero
    integer :: n

    integer, allocatable :: work(:), nn(:)
    integer :: cur_size, batchsize, num, i, n_old, nat
    logical :: more
    if(.not.allocated(self%cumsum)) then
       ! neiglist not computed
       self% errmsg="cluster:: neighbor list not computed."
       n = -1
       return
    end if

    ! check indices in input list
    nat = size(self%cumsum)
    if( any(list .le. 0) .or. any(list .gt. nat) ) then
       self%errmsg="cluster:: some index out of bounds in input list."
       n = -1
       return
    end if

    batchsize=size(list,1)
    ! cur_size can never be more than size of input list, no need to realloc

    allocate( work(1:batchsize), source=-1 )

    more = .true.
    ! if idx is not in input list, return zero
    n = 0
    cur_size = 0
    if( .not.any(list .eq. idx)) more = .false.

    if( more ) then
       ! start from first idx
       cur_size = 1
       work( cur_size ) = idx
    end if
    do while( more )
       more = .false.
       ! expand by one bond
       allocate(nn, source = work(1:cur_size))
       n_old = cur_size
       num = self% expand( 1, nn )
       if( num .lt. 0 ) then
          self%errmsg=self%errmsg//":: cluster"
          n=-1
          return
       end if
       ! check if expanded nn on input list
       do i = n_old+1, num
          if( any(list .eq. nn(i)) ) then
             ! add that index to work
             cur_size = cur_size + 1
             work(cur_size) = nn(i)
             more = .true.
          end if
       end do
       deallocate(nn)
    end do

    n = cur_size
    deallocate(list)
    allocate(list, source=work(1:cur_size))

    deallocate( work )
  end function t_neighbour_cluster


  ! transform `type(c_ptr)` string to fortran `character(:),allocatable` string
  function c2f_string(ptr) result(f_string)
    implicit none
    interface
       function c_strlen(str) bind(c, name='strlen')
         use iso_c_binding, only: c_ptr, c_size_t
         implicit none
         type(c_ptr), intent(in), value :: str
         integer(c_size_t) :: c_strlen
       end function c_strlen
    end interface
    type(c_ptr), intent(in) :: ptr
    character(len=:), allocatable :: f_string
    character(len=1, kind=c_char), dimension(:), pointer :: c_string
    integer :: n, i

    if (.not. c_associated(ptr)) then
       f_string = ' '
    else
       n = int(c_strlen(ptr), kind=kind(n))
       call c_f_pointer(ptr, c_string, [n+1])
       allocate( character(len=n)::f_string)
       do i = 1, n
          f_string(i:i) = c_string(i)
       end do
    end if
  end function c2f_string

end module m_neighbour
