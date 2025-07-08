module f_vesin

  use, intrinsic :: iso_c_binding


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
          sorted = .false., &

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
     ! integer( c_size_t ), allocatable :: pairs(:,:)
     ! integer( c_size_t ), pointer :: pairs(:,:)
     type( c_ptr ) :: pairs = c_null_ptr

     ! /// Array of box shifts, one for each `pair`. This is only set if
     ! /// `options.return_pairs` was `true` during the calculation.
     ! int32_t (*shifts)[3];
     ! integer( c_int32_t ), allocatable :: shifts(:,:)
     type( c_ptr ) :: shifts = c_null_ptr

     ! /// Array of pair distance (i.e. distance between the two points), one for
     ! /// each pair. This is only set if `options.return_distances` was `true`
     ! /// during the calculation.
     ! double *distances;
     ! real( c_double ), allocatable :: distances(:)
     type( c_ptr ) :: distances = c_null_ptr

     ! /// Array of pair vector (i.e. vector between the two points), one for
     ! /// each pair. This is only set if `options.return_vector` was `true`
     ! /// during the calculation.
     ! double (*vectors)[3];
     ! real( c_double ), allocatable :: vectors(:,:)
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
       real( c_double ), intent(in) :: points(:, :)
       integer( c_size_t ), value :: n_points
       real( c_double ), intent(in) :: box(3,3)
       logical( c_bool ), value :: periodic
       ! type( VesinDevice ), value :: device
       integer(c_int), value :: device
       type( VesinOptions ), value :: options
       ! type( c_ptr ), value :: neighbors
       type( VesinNeighborList ) :: neighbors
       type( c_ptr ) :: error_message
       integer( c_int ) :: res
     end function fvesin_neighbors
  end interface

contains


  function vesin_compute( nat, pos, box, vesin_opts, vesin_neighbors, errmsg )result(ierr)
    !!
    !!
    implicit none
    integer( c_int ), intent(in) :: nat
    real( c_double ), intent(in) :: pos(3,nat)
    real( c_double ), intent(in) :: box(3,3)
    type( VesinOptions ), intent(in) :: vesin_opts
    type( VesinNeighborList ), intent(out) :: vesin_neighbors
    character(:), allocatable, intent(out) :: errmsg
    integer( c_int ) :: ierr

    logical( c_bool ) :: periodic

    integer( c_size_t ) :: c_nat

    logical( c_bool ) :: c_periodic
    integer(c_int) :: dev

    c_nat = int( nat, c_size_t )

    ! c_periodic = .false.
    c_periodic = .true.


    ! set cpu device
    dev = VesinCPU
    vesin_neighbors% device = VesinCPU

    ierr = fvesin_neighbors( &
         pos,            &
         c_nat,          &
         box,            &
         c_periodic,     &
         dev,            &
         vesin_opts,         &
         vesin_neighbors,    &
         c_error_message )

    if( int(ierr) /= 0 ) errmsg = c2f_string(c_error_message)

  end function vesin_compute


  FUNCTION c2f_string(ptr) RESULT(f_string)
    implicit none
    interface
       function c_strlen(str) bind(c, name='strlen')
         use iso_c_binding, only: c_ptr, c_size_t
         implicit none
         type(c_ptr), intent(in), value :: str
         integer(c_size_t) :: c_strlen
       end function c_strlen
    end interface
    TYPE(c_ptr), INTENT(IN) :: ptr
    CHARACTER(LEN=:), ALLOCATABLE :: f_string
    CHARACTER(LEN=1, KIND=c_char), DIMENSION(:), POINTER :: c_string
    INTEGER :: n, i

    IF (.NOT. C_ASSOCIATED(ptr)) THEN
       f_string = ' '
    ELSE
       n = INT(c_strlen(ptr), KIND=KIND(n))
       ! write(*,*) "strlen",n
       CALL C_F_POINTER(ptr, c_string, [n+1])
       allocate( CHARACTER(LEN=n)::f_string)
       do i = 1, n
          f_string(i:i) = c_string(i)
       end do

       ! f_string = array2string(c_string, n)
    END IF
  END FUNCTION c2f_string

end module f_vesin
