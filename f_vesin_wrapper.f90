module f_vesin_wrapper

  use, intrinsic :: iso_c_binding
  implicit none

  private
  public :: vesin_t
  public :: rp

  !! precision used by caller code (change as needed)
  integer, parameter :: rp = kind(1.0_c_double)


  ! /// Device on which the data can be
  ! enum VesinDevice {
  ! /// Unknown device, used for default initialization and to indicate no
  ! /// allocated data.
  ! VesinUnknownDevice = 0,
  integer( c_int ), parameter :: VesinUnknownDevice = 0
  ! /// CPU device
  ! VesinCPU = 1,
  integer( c_int ), parameter :: VesinCPU = 1
  ! };

  ! /// Options for a neighbor list calculation
  ! struct VesinOptions {
  type, bind(c) ::  VesinOptions

     real( c_double ) :: &
          ! /// Spherical cutoff, only pairs below this cutoff will be included
          ! double cutoff;
          cutoff = 0.0_c_double

     logical( c_bool ) :: &

          ! /// Should the returned neighbor list be a full list (include both `i -> j`
          ! /// and `j -> i` pairs) or a half list (include only `i -> j`)?
          ! bool full;
          full = .false., &

          ! /// Should the neighbor list be sorted? If yes, the returned pairs will be
          ! /// sorted using lexicographic order.
          ! bool sorted;
          sorted = .true., &

          ! /// Should the returned `VesinNeighborList` contain `shifts`?
          ! bool return_shifts;
          return_shifts = .false., &

          ! /// Should the returned `VesinNeighborList` contain `distances`?
          ! bool return_distances;
          return_distances = .false., &

          ! /// Should the returned `VesinNeighborList` contain `vector`?
          ! bool return_vectors;
          return_vectors = .false.

  end type VesinOptions

  ! /// The actual neighbor list
  ! ///
  ! /// This is organized as a list of pairs, where each pair can contain the
  ! /// following data:
  ! ///
  ! /// - indices of the points in the pair;
  ! /// - distance between points in the pair, accounting for periodic boundary
  ! ///   conditions;
  ! /// - vector between points in the pair, accounting for periodic boundary
  ! ///   conditions;
  ! /// - periodic shift that created the pair. This is only relevant when using
  ! ///   periodic boundary conditions, and contains the number of bounding box we
  ! ///   need to cross to create the pair. If the positions of the points are `r_i`
  ! ///   and `r_j`, the bounding box is described by a matrix of three vectors `H`,
  ! ///   and the periodic shift is `S`, the distance vector for a given pair will
  ! ///   be given by `r_ij = r_j - r_i + S @ H`.
  ! ///
  ! /// Under periodic boundary conditions, two atoms can be part of multiple pairs,
  ! /// each pair having a different periodic shift.
  ! struct VESIN_API VesinNeighborList {
  type, bind(c) :: VesinNeighborList

     ! /// Number of pairs in this neighbor list
     ! size_t length;
     integer( c_size_t ) :: length = 0_c_size_t
     ! type( c_ptr ) :: length

     ! /// Device used for the data allocations
     ! VesinDevice device;
     integer( c_int ) :: device = VesinCPU

     ! /// Array of pairs (storing the indices of the first and second point in the
     ! /// pair), containing `length` elements.
     ! size_t (*pairs)[2];
     type( c_ptr ) :: pairs = c_null_ptr

     ! /// Array of box shifts, one for each `pair`. This is only set if
     ! /// `options.return_pairs` was `true` during the calculation.
     ! int32_t (*shifts)[3];
     type( c_ptr ) :: shifts = c_null_ptr

     ! /// Array of pair distance (i.e. distance between the two points), one for
     ! /// each pair. This is only set if `options.return_distances` was `true`
     ! /// during the calculation.
     ! double *distances;
     type( c_ptr ) :: distances = c_null_ptr

     ! /// Array of pair vector (i.e. vector between the two points), one for
     ! /// each pair. This is only set if `options.return_vector` was `true`
     ! /// during the calculation.
     ! double (*vectors)[3];
     type( c_ptr ) :: vectors = c_null_ptr

  end type VesinNeighborList


  ! /// Compute a neighbor list.
  ! ///
  ! /// The data is returned in a `VesinNeighborList`. For an initial call, the
  ! /// `VesinNeighborList` should be zero-initialized (or default-initalized in
  ! /// C++). The `VesinNeighborList` can be re-used across calls to this functions
  ! /// to re-use memory allocations, and once it is no longer needed, users should
  ! /// call `vesin_free` to release the corresponding memory.
  ! ///
  ! /// @param points positions of all points in the system;
  ! /// @param n_points number of elements in the `points` array
  ! /// @param box bounding box for the system. If the system is non-periodic,
  ! ///     this is ignored. This should contain the three vectors of the bounding
  ! ///     box, one vector per row of the matrix.
  ! /// @param periodic is the system using periodic boundary conditions?
  ! /// @param device device where the `points` and `box` data is allocated.
  ! /// @param options options for the calculation
  ! /// @param neighbors non-NULL pointer to `VesinNeighborList` that will be used
  ! ///     to store the computed list of neighbors.
  ! /// @param error_message Pointer to a `char*` that wil be set to the error
  ! ///     message if this function fails. This does not need to be freed when no
  ! ///     longer needed.
  !
  ! C-header:
  !
  !~~~~~~~~~~~~~~~~~{.c}
  ! int VESIN_API vesin_neighbors(
  !     const double (*points)[3],
  !     size_t n_points,
  !     const double box[3][3],
  !     bool periodic,
  !     VesinDevice device,
  !     struct VesinOptions options,
  !     struct VesinNeighborList* neighbors,
  !     const char** error_message
  ! );
  !~~~~~~~~~~~~~~~~~
  interface
     function fvesin_neighbors( &
          points,        &
          n_points,      &
          box,           &
          periodic,      &
          device,        &
          options,       &
          neighbors,     &
          error_message  &
       )result(res)bind( c, name="vesin_neighbors" )
       import :: c_double, c_size_t, c_bool, c_ptr, c_int, VesinOptions, VesinNeighborList
       integer( c_size_t ), value :: n_points
       real( c_double ), intent(in) :: points(3, n_points)
       real( c_double ), intent(in) :: box(3,3)
       logical( c_bool ), value :: periodic
       integer(c_int), value :: device
       type( VesinOptions ), value :: options
       type( VesinNeighborList ) :: neighbors
       type( c_ptr ) :: error_message
       integer( c_int ) :: res
     end function fvesin_neighbors
  end interface


  ! /// Free all allocated memory inside a `VesinNeighborList`, according the it's
  ! /// `device`.
  ! void VESIN_API vesin_free(struct VesinNeighborList* neighbors);
  interface
     subroutine fvesin_free( neighbors ) bind(C, name="vesin_free")
       import :: VesinNeighborList
       type( VesinNeighborList ), value :: neighbors
     end subroutine fvesin_free
  end interface



  !> @details
  !! fortran derived type, holding the input options and output data from Vesin.
  !!
  !! Use as:
  !!~~~~~~~~{.f90}
  !!
  !! program main
  !!   use f_vesin_wrapper
  !!   type( vesin_t ), pointer :: me
  !!
  !!   ! create the instance, set options
  !!   me => vesin_t( cutoff=<val>, return_vectors=<val>, etc)
  !!
  !!   ! launch computation of neighbor list
  !!   ierr = me% compute( nat, pos, box )
  !!   if( ierr/= 0 ) then
  !!      write(*,*) me% errmsg
  !!      stop
  !!   end if
  !!
  !!   ! data is inside `me`:
  !!   write(*,*) me% pairs
  !!
  !!   ! destroy data
  !!   deallocate( me )
  !!
  !! end program main
  !!
  !!~~~~~~~~
  type :: vesin_t

     ! options
     type( VesinOptions ), private :: opts

     !! error message
     character(:), allocatable, public :: errmsg

     !! number of elements in the neighbor list
     integer, public :: length

     !! Array of pairs (storing the indices of the first and second point in the
     !! pair), containing `length` elements.
     integer,  allocatable, public :: pairs(:,:)

     !! Array of box shifts, one for each `pair`. This is only set if
     !! `return_pairs` option was `true` during the calculation.
     integer,  allocatable, public :: shifts(:,:)

     !! Array of pair distance (i.e. distance between the two points), one for
     !! each pair. This is only set if `return_distances` option was `true`
     !! during the calculation.
     real(rp), allocatable, public :: distances(:)

     !! Array of pair vector (i.e. vector between the two points), one for
     !! each pair. This is only set if `return_vectors` option was `true`
     !! during the calculation.
     real(rp), allocatable, public :: vectors(:,:)

   contains
     procedure, public :: compute => vesin_t_compute
     final :: vesin_t_destroy
  end type vesin_t

  ! overload name
  interface vesin_t
     procedure vesin_set_options
  end interface vesin_t


contains

  function vesin_set_options( &
       cutoff, &
       full, &
       sorted, &
       return_shifts, &
       return_distances, &
       return_vectors )result(self)
    !> @details
    !! Construct and set the options for `vesin_t`. This directly sets values to the private
    !! member `vesin_t% opts`, which is a `type( VesinOptions )` instance, used in the
    !! calculation of the neighbor list.
    !!
    !! @param `cutoff`, real ::  Spherical cutoff, only pairs below this cutoff will be included. Default=0.0
    !! @param `full`, logical, optional :: Should the returned neighbor list be a full
    !!       list (include both `i -> j` and `j -> i` pairs) or a half
    !!       list (include only `i -> j`)? Default=.false.
    !! @param `sorted`, logical, optional :: Should the neighbor list be sorted? If yes,
    !!       the returned pairs will be sorted using lexicographic order. Default=.true.
    !! @param `return_shifts`, logical, optional :: Should `vesin_t` contain `shifts`? Default=.false.
    !! @param `return_distances`, logical, optional :: Should `vesin_t` contain `distances`? Default=.false.
    !! @param `return_vectors`, logical, optional :: Should `vesin_t` contain `vectors`? Default=.false.
    !!
    implicit none
    real( rp ), intent(in) :: cutoff
    logical, intent(in), optional :: full
    logical, intent(in), optional :: sorted
    logical, intent(in), optional :: return_shifts
    logical, intent(in), optional :: return_distances
    logical, intent(in), optional :: return_vectors
    type( vesin_t ), pointer :: self

    allocate( vesin_t :: self )

    self% opts% cutoff = real( cutoff, c_double )
    if(present(full)            ) self% opts% full = logical( full, c_bool )
    if(present(sorted)          ) self% opts% sorted = logical( sorted, c_bool )
    if(present(return_shifts)   ) self% opts% return_shifts = logical( return_shifts, c_bool )
    if(present(return_distances)) self% opts% return_distances = logical( return_distances, c_bool )
    if(present(return_vectors)  ) self% opts% return_vectors = logical( return_vectors, c_bool )

  end function vesin_set_options


  function vesin_t_compute( self, nat, pos, box, periodic )result(ierr)
    !> @details
    !! Compute the neighbor list with options provided in `vesin_t% opts`.
    !! The data is recorded first in C format in `type(VesinNeighborList)`, then it is
    !! transferred to fortran format into `self`. The C data is destroyed at the end
    !! of this function.
    !!
    !! @param `nat`, integer :: number of atoms
    !! @param `pos`, real(rp), [3,nat] :: atomic positions
    !! @param `box`, real(rp), [3,3] :: periodic box vectors
    !! @param `periodic`, logical, optional :: flag for (non)-periodic calculation. Default=.true.
    implicit none
    class( vesin_t ), intent(inout) :: self
    integer, intent(in)  :: nat
    real(rp), intent(in) :: pos(3,nat)
    real(rp), intent(in) :: box(3,3)
    logical, intent(in), optional :: periodic
    integer :: ierr

    type( VesinNeighborList ) :: cdata
    logical( c_bool ) :: c_periodic
    integer( c_size_t ) :: c_nat
    real( c_double ) :: c_pos(3,nat), c_box(3,3)
    integer( c_int ) :: dev
    type( c_ptr ) :: c_errmsg

    integer :: n
    integer( c_size_t ), pointer :: pairs(:,:) => null()
    integer( c_int32_t ), pointer :: shifts(:,:) => null()
    real(c_double ), pointer :: distances(:) => null()
    real( c_double ), pointer :: vectors(:,:) => null()

    ! perform periodic calc by default
    c_periodic = .true.
    if(present(periodic)) c_periodic=logical(periodic, c_bool)

    ! input data to C
    c_nat = int( nat, c_size_t )
    c_pos = real( pos, c_double )
    c_box = real( box, c_double )
    dev = VesinCPU
    cdata% device = VesinCPU

    ! compute the neighbor list
    ierr = int( fvesin_neighbors(c_pos, c_nat, c_box, c_periodic, dev, self%opts, cdata, c_errmsg ))
    if( ierr/= 0 ) then
       self% errmsg = c2f_string( c_errmsg )
       return
    end if

    ! transform cdata to fdata
    self%length = int(cdata%length, kind(self%length))
    n = int( cdata%length )
    ! pairs
    if( c_associated(cdata%pairs)) then
       call c_f_pointer( cdata%pairs, pairs, shape=[2,n] )
       allocate( self%pairs, source=int(pairs, kind(self%pairs)) )
    end if
    ! shifts
    if( c_associated(cdata%shifts)) then
       call c_f_pointer( cdata%shifts, shifts, shape=[3,n])
       allocate( self%shifts, source=int(shifts, kind(self%shifts)) )
    end if
    ! distances
    if( c_associated(cdata%distances)) then
       call c_f_pointer( cdata%distances, distances, shape=[n] )
       allocate( self%distances, source=real(distances,kind(self%distances)) )
    end if
    ! vectors
    if( c_associated(cdata%vectors))then
       call c_f_pointer( cdata%vectors, vectors, shape=[3,n])
       allocate( self%vectors, source=real(vectors, kind(self%vectors)) )
    end if
    ! now can free cdata
    call vesin_free(cdata)

  end function vesin_t_compute

  subroutine vesin_free( neighbors )
    implicit none
    type( VesinNeighborList ), intent(inout) :: neighbors
    call fvesin_free( neighbors )
  end subroutine vesin_free

  subroutine vesin_t_destroy( self )
    !! destructor
    implicit none
    type( vesin_t ), intent(inout) :: self
    if( allocated( self%pairs)) deallocate( self%pairs )
    if( allocated( self%shifts)) deallocate( self%shifts )
    if( allocated( self%distances)) deallocate( self%distances )
    if( allocated( self%vectors)) deallocate( self%vectors )
  end subroutine vesin_t_destroy


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


end module f_vesin_wrapper
