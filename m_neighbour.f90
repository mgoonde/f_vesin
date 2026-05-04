module m_neighbour

  use, intrinsic :: iso_c_binding

  ! import the C-interoperable vesin interface
  use vesin_c, only: &
       c_VesinOptions => VesinOptions, &
       c_VesinNeighborList => VesinNeighborList, &
       c_vesin_neighbors => vesin_neighbors, &
       c_vesin_free => vesin_free, &
       c_vesinDevice => VesinDevice, &
       VesinUnknownDevice, VesinCPU, VesinAutoAlgorithm

  implicit none

  private
  public :: rp
  public :: t_neighbour
  public :: t_neighbour_copy

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
  !!    n = neigh% get( idx = 5, nshell = 2, include_idx = .true., &
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
  type :: t_neighbour

     ! ----- private ------
     ! private
     ! vesin stuffs
     type( c_vesinOptions ) :: opts       !< Computation options
     type( c_vesinNeighborList ) :: cdata !< Returned C data
     ! integer :: device = VesinCPU
     type( c_vesinDevice ) :: device
     logical :: initialized = .false.          !< .true. when initialized
     logical :: active = .false.          !< .true. when it contains data
     ! pointers to C data, in C precision
     integer(c_size_t) :: length = 0_c_size_t               !< size of the list
     integer( c_size_t ), pointer :: pairs(:,:) => null()   !< shape[2, length]
     integer( c_int32_t ), pointer :: shifts(:,:) => null() !< shape[3, length]
     real( c_double ), pointer :: distances(:) => null()    !< shape[length]
     real( c_double ), pointer :: vectors(:,:) => null()    !< shape[3, length]

     ! local, private
     logical :: is_copy = .false. !< when .true., the instance is a copy of another,
                                  !! which means pointer data is allocated, not associated
     integer, allocatable :: cumsum(:) !< cumulative sum of nneig
     integer, allocatable :: ityp(:) !< copy of ityp


     ! ----- public ------
     character(:), allocatable, public :: errmsg !< error message


   contains
     procedure, public :: compute => t_neighbour_compute
     procedure, public :: get_nn  => t_neighbour_get_nn
     procedure, public :: get     => t_neighbour_get
     procedure, public :: expand  => t_neighbour_expand
     procedure, public :: cluster => t_neighbour_cluster
     procedure, private :: t_neighbour_get_by_rcut_single
     procedure, private :: t_neighbour_get_by_rcut_list
     generic, public :: get_by_rcut => t_neighbour_get_by_rcut_single
     generic, public :: get_by_rcut => t_neighbour_get_by_rcut_list
     procedure, public :: deactivate => t_neighbour_deactivate
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

    self%device%type = vesinCPU
    self%device%device_id = 0

    ! some vesin options
    self%opts% full             = .true.
    self%opts% sorted           = .true.
    self%opts% return_shifts    = .true.
    self%opts% return_distances = .false.
    self%opts% return_vectors   = .true.
    self%opts% algorithm        = VesinAutoAlgorithm
    if( present(return_distances)) self%opts% return_distances = return_distances

    ! set as initialized
    self% initialized = .true.

  end function t_neighbour_construct


  subroutine t_neighbour_destroy( self )
    !! destroy the `t_neighbour` instance:
    !! deactivate and free underlying C data
    implicit none
    ! type( t_neighbour ), intent(inout) :: self
    class( t_neighbour ), intent(inout) :: self
    call self% deactivate()
    self% initialized = .false.
    if( .not. self%is_copy ) call c_vesin_free(self%cdata)
  end subroutine t_neighbour_destroy

  subroutine t_neighbour_deactivate( self )
    !! deactivate the `t_neighbour` instance:
    !! nullify pointers to C data, and deallocate any local array
    implicit none
    ! type( t_neighbour ), intent(inout) :: self
    class( t_neighbour ), intent(inout) :: self
    if( self%is_copy ) then
       if( associated(self%pairs))    deallocate(self%pairs)
       if( associated(self%shifts))   deallocate(self%shifts)
       if( associated(self%distances))deallocate(self%distances)
       if( associated(self%vectors))  deallocate(self%vectors)
    else
       if( associated(self%pairs))    nullify(self%pairs)
       if( associated(self%shifts))   nullify(self%shifts)
       if( associated(self%distances))nullify(self%distances)
       if( associated(self%vectors))  nullify(self%vectors)
    end if

    if(allocated(self%cumsum))deallocate(self%cumsum)
    if(allocated(self%ityp))  deallocate(self%ityp)
    self% active = .false.
  end subroutine t_neighbour_deactivate



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

    logical(c_bool) :: periodic(3)
    type( c_ptr ) :: c_errmsg
    integer :: n, i, idx

    ierr = -1
    if( .not. self% initialized )then
       self% errmsg = "compute:: t_neighbour instance not initialized."
       return
    end if

    ! if we contain data from previous calc, deactivate and let vesin reuse
    ! its allocation if it can
    if( self% active ) call self% deactivate()

    ! set cutoff
    self%opts% cutoff = real( rcut, c_double )
    ! periodic
    periodic(:) = logical( .true., c_bool )

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

    ! set as active
    self% active = .true.
  end function t_neighbour_compute



  function t_neighbour_get( self, idx, list, ityplist, veclist, shiftslist, nshell, include_idx )result(n)
    !! Get the neighbor data of `idx` from neighbour list, up to `nshell` neighbor shells.
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

    !> how many neighbor shells to get (default=1). If `nshell>1`, the output list is not sorted
    integer, intent(in), optional :: nshell

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
    if(present(nshell))nb=nshell
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
    !! Return `n` which is the number of neighbors, or negative on error.
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


  function t_neighbour_expand( self, nshell, list, ityplist, veclist, shiftslist ) result( n )
    !! expand the list in input by `nshell` number of neighbor shells.
    !! NOTE: pbc are not taken into account here, the expansion goes beyond the cell,
    !! which means the same atom index can repeat multiple times on different positions,
    !! as copies of itself by periodicity.
    implicit none
    class( t_neighbour ), intent(inout) :: self

    !> number of bons shells to expand
    integer, intent(in) :: nshell

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

    !> number of elements in output list, negative on error
    integer :: n
    !
    integer :: i, ntot, ncur, idx, iatm, nprev, nat
    integer, allocatable :: work(:), swork(:,:)
    real(rp), allocatable :: vwork(:,:), vtmp(:,:)
    integer, parameter :: batchsize=50
    integer :: ishell, ii, nl_idx, nn, jj, j
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

    ! check if input list is allocated
    if( .not. allocated(list)) then
       self%errmsg="expand:: input list not allocated."
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

    do ishell = 1, nshell
       ! write(*,*) ">> ishell",ishell
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

    end do ! ishell


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
    !! around `idx`. The returned `list` is strictly a subset of the input `list`.
    !! The order of indices might be lost on output.
    !! NOTE: works purely on indices, no regard for pbc or not.
    implicit none

    class( t_neighbour ), intent(inout) :: self

    !> starting atom index
    integer, intent(in) :: idx

    !> on input the possibly cluster-graph to check, on output the complete cluster
    !! graph around `idx`
    integer, allocatable, intent(inout) :: list(:)

    !> Number of atoms in output cluster.
    !! If `idx` is not present in `list`, return zero.
    !! Negative on error
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

    ! check if input list allocated
    if( .not. allocated(list)) then
       self% errmsg="cluster:: input list not allocated."
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
       ! expand by one shell
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


  function t_neighbour_get_by_rcut_single( self, idx, rcut, nat, ityp, pos, lat, &
       list, ityplist, veclist, shiftslist, include_idx ) result(n)
    !! get neighbor list of specific index atom, within specified atomic system,
    !! and with a cutoff different from the one in `self`.
    implicit none
    class( t_neighbour ), intent(inout) :: self

    !> input atomic index
    integer, intent(in)  :: idx

    !> custom rcut
    real(rp), intent(in) :: rcut

    !> number of atoms
    integer, intent(in) :: nat

    !> atomic types
    integer, intent(in) :: ityp(nat)

    !> atomic positions
    real(rp), intent(in) :: pos(3,nat)

    !> lattice: lat(:,1)=v1, lat(:,2)=v2, lat(:,3)=v3
    real(rp), intent(in) :: lat(3,3)

    !> output list of neighbours to atom `idx`
    integer, allocatable, intent(out), optional  :: list(:)

    !> output list of atomic types neighbor to `idx`
    integer, allocatable, intent(out), optional :: ityplist(:)

    !> output list of vectors neighbour to `idx`
    real(rp), allocatable, intent(out), optional :: veclist(:,:)

    !> output list of box shifts for each neighbor
    integer, allocatable, intent(out), optional :: shiftslist(:,:)

    !> flag to include the atom `idx` in output. If .true., it will be on the first
    !! element of output. Default=.false.
    logical, intent(in), optional :: include_idx

    !> output number of neighbours, or negative on error
    integer :: n

    logical :: inc_idx
    real(rp) :: invlat(3,3)
    integer :: nbox(3)
    integer :: i, j, ii, jj, kk, n_cur, nmax
    real(rp) :: rshift(3), rij(3), ri(3), rj(3)
    integer, allocatable :: slist(:,:)
    integer, allocatable :: tmp(:)
    real(rp), allocatable :: vlist(:,:)
    real(rp) :: dist2, rcut2
    ! real(rp) :: pos2(3,nat)

    if( idx < 1 .or. idx > nat ) then
       n = -1
       self%errmsg="get_by_rcut:: `idx` out of bounds in input."
       return
    end if

    inc_idx = .false.
    if( present(include_idx) ) inc_idx=include_idx

    ! inverse lattice
    call inverse( lat, invlat )

    ! most brute-force algo for single atom.

    ! how many boxes to add in each direction
    ! this counts as from -nbox to +nbox
    do i = 1, 3
       nbox(i) = nint( rcut*norm2(invlat(:,i)) )
    end do
    ! write(*,*) "rcut",rcut
    ! write(*,*) "nbox:",nbox


    rcut2 = rcut*rcut

    n_cur = 0
    ! absolute maximal number of atoms: use for alloc size
    nmax = product(2*nbox+1)*nat
    allocate(tmp(1:nmax))
    allocate(vlist(1:3, 1:nmax))
    allocate(slist(1:3, 1:nmax))
    ! include origin idx on first place
    if( inc_idx ) then
       n_cur = n_cur + 1
       tmp( n_cur ) = idx
       vlist(:, n_cur) = [0.0_rp, 0.0_rp, 0.0_rp]
       slist(:, n_cur) = [0, 0, 0]
    end if

    ! my pos
    ri = pos(:,idx)
    call cart_to_crist( ri, lat, invlat )
    call periodic( ri )
    call crist_to_cart( ri, lat, invlat )

    do j = 1, nat
       rj = pos(:,j)
       call cart_to_crist( rj, lat, invlat )
       call periodic( rj )
       call crist_to_cart( rj, lat, invlat )
       ! for all box shifts
       do ii = -nbox(1), nbox(1)
          do jj = -nbox(2), nbox(2)
             do kk = -nbox(3), nbox(3)
                ! skip self in the same box
                if( ii==0 .and. jj==0 .and. kk==0 .and. j==idx ) cycle
                ! shift (in crist coords)
                rshift = [real(ii,rp), real(jj,rp), real(kk,rp)]
                ! shift in cartesian
                rshift = matmul( lat, rshift )
                rij = rj+rshift - ri
                dist2 = dot_product(rij, rij)
                if( dist2 <= rcut2 ) then
                   ! add j
                   n_cur = n_cur + 1
                   if( n_cur .gt. nmax ) error stop "n_cur exceeds nmax..! should not happen"
                   tmp( n_cur ) = j
                   vlist(:, n_cur ) = rij
                   slist(:, n_cur ) = [ii, jj, kk]
                   ! write(*,*) ii, jj, kk, j, dist2
                end if
             end do
          end do
       end do
    end do

    ! set output
    n = n_cur
    if( present(list)) allocate(list, source=tmp(1:n))
    if( present(veclist)) allocate(veclist, source=vlist(1:3,1:n))
    if( present(shiftslist))allocate(shiftslist, source=slist(1:3,1:n))
    if( present(ityplist)) then
       allocate( ityplist(1:n) )
       do i = 1, n
          ityplist(i) = ityp( tmp(i) )
       end do
    end if


    ! write(*,"(10i4)") tmp
  end function t_neighbour_get_by_rcut_single

  function t_neighbour_get_by_rcut_list( self, idx_list, rcut, nat, ityp, pos, lat, &
       list, ityplist, veclist, shiftslist, include_idx ) result(n)
    !! get neighbor list of specific index atom, within specified atomic system,
    !! and with a cutoff different from the one in `self`.
    implicit none
    class( t_neighbour ), intent(inout) :: self

    !> input list of atomic index
    integer, intent(in)  :: idx_list(:)

    !> custom rcut
    real(rp), intent(in) :: rcut

    !> number of atoms
    integer, intent(in) :: nat

    !> atomic types
    integer, intent(in) :: ityp(nat)

    !> atomic positions
    real(rp), intent(in) :: pos(3,nat)

    !> lattice: lat(:,1)=v1, lat(:,2)=v2, lat(:,3)=v3
    real(rp), intent(in) :: lat(3,3)

    !> output list of neighbours to atom `idx`
    integer, allocatable, intent(out), optional  :: list(:)

    !> output list of atomic types neighbor to `idx`
    integer, allocatable, intent(out), optional :: ityplist(:)

    !> output list of vectors neighbour to `idx`
    real(rp), allocatable, intent(out), optional :: veclist(:,:)

    !> output list of box shifts for each neighbor
    integer, allocatable, intent(out), optional :: shiftslist(:,:)

    !> flag to include the atom `idx` in output. If .true., it will be on the first
    !! element of output. Default=.false.
    logical, intent(in), optional :: include_idx

    !> output number of neighbours, or negative on error
    integer :: n

    integer :: i, j, idx, nlist
    integer :: jj
    integer :: tmp_n, tmp_idx
    integer :: cursize, newsize, n_cur
    integer, allocatable :: tmp_list(:), tmp_ityplist(:), tmp_shiftslist(:,:)
    real(rp), allocatable :: tmp_veclist(:,:)
    integer, allocatable :: tmp2_list(:), tmp2_ityplist(:), tmp2_shiftslist(:,:)
    real(rp), allocatable :: tmp2_veclist(:,:)
    integer, allocatable :: agg_list(:), agg_ityplist(:), agg_shiftslist(:,:)
    real(rp), allocatable :: agg_veclist(:,:)
    integer, allocatable :: jcheck(:)


    ! compute for first index
    idx = idx_list(1)
    tmp_n = self% get_by_rcut( idx, rcut, nat, ityp, pos, lat, &
         list=tmp_list, ityplist=tmp_ityplist, veclist=tmp_veclist, shiftslist=tmp_shiftslist,&
         include_idx=include_idx )

    nlist = size(idx_list)
    allocate( agg_list, source=tmp_list )
    allocate( agg_ityplist, source=tmp_ityplist )
    allocate( agg_shiftslist, source=tmp_shiftslist )
    allocate( agg_veclist, source=tmp_veclist )
    cursize=tmp_n
    ! counter
    n_cur = tmp_n
    ! compute for others in idx_list, shift them by pos(idx)
    do i = 2, nlist
       !
       idx = idx_list(i)
       tmp_n = self% get_by_rcut( idx, rcut, nat, ityp, pos, lat, &
            list=tmp_list, ityplist=tmp_ityplist, veclist=tmp_veclist, shiftslist=tmp_shiftslist,&
            include_idx=include_idx )
       write(*,*) "from idx", idx
       write(*,"(*(i5))") tmp_list

       ! check for aggregator size
       if( n_cur+tmp_n > cursize ) then
          ! realloc
          newsize = cursize + tmp_n
          call move_alloc( agg_list, tmp2_list )
          allocate( agg_list(1:newsize) )
          agg_list(1:cursize) = tmp2_list
          deallocate( tmp2_list )
          !
          call move_alloc( agg_ityplist, tmp2_ityplist )
          allocate( agg_ityplist(1:newsize) )
          agg_ityplist(1:cursize) = tmp2_ityplist
          deallocate( tmp2_ityplist )
          !
          call move_alloc( agg_shiftslist, tmp2_shiftslist )
          allocate( agg_shiftslist(1:3,1:newsize))
          agg_shiftslist(1:3,1:cursize) = tmp2_shiftslist
          deallocate( tmp2_shiftslist )
          !
          call move_alloc( agg_veclist, tmp2_veclist )
          allocate( agg_veclist(1:3,1:newsize))
          agg_veclist(1:3,1:cursize) = tmp2_veclist
          deallocate( tmp2_veclist )
          !
          cursize = newsize
       end if


       j_: do j = 1, tmp_n
          tmp_idx = tmp_list(j)
          if( any(agg_list(1:n_cur)==tmp_list(j)) ) then
             ! check for all instances, if same vec, same shift
             ! jcheck contains all instances of index tmp_list(j) in agg lists
             jcheck = pack( [(i,i=1,n_cur)], [agg_list(1:n_cur)-tmp_idx==0])
             ! write(*,*) "tmp_list(j)",tmp_list(j)
             ! write(*,"(*(i4))") agg_list(1:n_cur)
             ! write(*,*) "jcheck:", jcheck
             do jj = 1, size(jcheck)
                if( agg_ityplist(jj) == tmp_ityplist(j) .and. &
                     all(agg_shiftslist(:,jj) == tmp_shiftslist(:,j)) .and. &
                     all(abs(agg_veclist(:,jj))-abs(tmp_veclist(:,j)) < 1e-6) ) cycle j_
             end do
          end if

          ! write(*,*) "adding", tmp_list(j)
          n_cur = n_cur + 1
          agg_list(n_cur) = tmp_list(j)
          agg_ityplist(n_cur) = tmp_ityplist(j)
          agg_shiftslist(:,n_cur) = tmp_shiftslist(:,j)
          agg_veclist(:,n_cur) = tmp_veclist(:,j) - (pos(:,idx) - pos(:,idx_list(1)))
       end do j_

       deallocate( tmp_list, tmp_shiftslist, tmp_ityplist, tmp_veclist )
    end do

    if( present(list)) allocate(list, source=agg_list(1:n_cur))
    if(present(ityplist))allocate(ityplist, source=agg_ityplist(1:n_cur))
    if(present(veclist))allocate(veclist, source=agg_veclist(1:3,1:n_cur))
    if(present(shiftslist))allocate(shiftslist, source=agg_shiftslist(1:3,1:n_cur))
    n = n_cur

  end function t_neighbour_get_by_rcut_list


  subroutine t_neighbour_copy( from, to )
    !! allocate Fortran pointers, not associate to C-data
    implicit none
    type( t_neighbour ), intent(in) :: from
    ! type( t_neighbour ), pointer, intent(out) :: to
    type( t_neighbour ), intent(out) :: to
    to = t_neighbour()
    ! we are making a hard-copy, indicate that data is allocated here
    to% is_copy = .true.
    to% length = from% length
    if( associated(from%pairs) ) allocate( to%pairs, source=from%pairs )
    if( associated(from%shifts) ) allocate( to%shifts, source=from%shifts )
    if( associated(from%distances) ) allocate( to%distances, source=from%distances )
    if( associated(from%vectors) ) allocate( to%vectors, source=from%vectors )
    if( allocated( from%cumsum)) allocate( to%cumsum, source=from%cumsum )
    if( allocated( from%ityp)) allocate( to%ityp, source=from%ityp )
    ! indicate data is present
    to% active = .true.
  end subroutine t_neighbour_copy



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


  pure subroutine inverse( mat, inv )
    !! inverse of 3x3 matrix
    implicit none
    real(rp), intent(in) :: mat(3,3)
    real(rp), intent(out) :: inv(3,3)

    real(rp) :: det, invdet

    ! calculate the determinant
    det =  mat(1,1)*mat(2,2)*mat(3,3) &
         + mat(1,2)*mat(2,3)*mat(3,1) &
         + mat(1,3)*mat(2,1)*mat(3,2) &
         - mat(1,3)*mat(2,2)*mat(3,1) &
         - mat(1,2)*mat(2,1)*mat(3,3) &
         - mat(1,1)*mat(2,3)*mat(3,2)
    ! invert the determinant
    invdet = 1.0_rp/det
    ! calculate the inverse matrix
    inv(1,1) = invdet  * ( mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2) )
    inv(2,1) = -invdet * ( mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1) )
    inv(3,1) = invdet  * ( mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1) )
    inv(1,2) = -invdet * ( mat(1,2)*mat(3,3) - mat(1,3)*mat(3,2) )
    inv(2,2) = invdet  * ( mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1) )
    inv(3,2) = -invdet * ( mat(1,1)*mat(3,2) - mat(1,2)*mat(3,1) )
    inv(1,3) = invdet  * ( mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2) )
    inv(2,3) = -invdet * ( mat(1,1)*mat(2,3) - mat(1,3)*mat(2,1) )
    inv(3,3) = invdet  * ( mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1) )
  end subroutine inverse

  pure elemental subroutine periodic(c)
    ! periodic boundary condition, for any dimensional vector input in crist coords.
    implicit none
    real(RP), intent(inout) :: c

    if( c .lt. -0.5_rp ) c = c + 1.0_rp
    if( c .ge. 0.5_rp  ) c = c - 1.0_rp
  end subroutine periodic

  pure subroutine crist_to_cart( rij, lat, invlat )
    ! lat is lattice vectors as:
    !
    !  lat(:,1) = a1 a2 a3
    !  lat(:,2) = b1 b2 b3
    !  lat(:,3) = c1 c2 c3
    !
    ! invlat is inverse of lat
    implicit none
    real(rp), intent(inout) :: rij(:)
    real(rp), intent(in) :: lat(:,:)
    real(rp), intent(in) :: invlat(:,:)

    rij = matmul( lat, rij )
  end subroutine crist_to_cart

  pure subroutine cart_to_crist( rij, lat, invlat )
    ! lat is lattice vectors as:
    !
    !  lat(:,1) = a1 a2 a3
    !  lat(:,2) = b1 b2 b3
    !  lat(:,3) = c1 c2 c3
    !
    ! invlat is inverse of lat
    implicit none
    real(rp), intent(inout) :: rij(:)
    real(rp), intent(in) :: lat(:,:)
    real(rp), intent(in) :: invlat(:,:)

    rij = matmul( invlat, rij )
  end subroutine cart_to_crist

end module m_neighbour
