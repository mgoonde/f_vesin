# fortran wrapper to vesin

A slightly more sophisticated fortran wrapper to `vesin` from: https://github.com/Luthaf/vesin

# additional functionality

Adds functions to post-process the neighbor list returned by `vesin`. Example:

```f90
type( t_neighbour ) :: neigh
integer, allocatable :: list(:)
real(rp), allocatable :: coords(:,:)

! initialize neighbor list
neigh = neighbour()

! compute the neighbor list
ierr = neigh% compute( nat, ityp, pos, lat, rcut )
if( ierr /= 0 ) then
   write(*,*) neigh% errmsg
   error stop 1
end if

! get the 1st nearest neighbors shell of `idx`
n = neigh% get_nn( idx, list=list )

! get all up to m-th neighbor shell, and the atomic positions, including the one of `idx`
m = 4
n = neigh% get( idx, nshell=m, list=list, veclist=coords, include_idx=.true. )

! expand the `list` by m neighbor shells
n = neigh% expand( m, list )

! input any `list` of indices, and output just the subset of indices
! which are within the connected cluster of input `idx`
list = [ 3, 5, 19, 21, 6, 123 ]
idx = 5
n = neigh% cluster( idx, list)
write(*,*) list ! example output: 5, 3, 6

! get the neighbor list for a single atom, with a rcut value different
! from the one used in `compute`
idx_new = 2
n = neigh% get_by_rcut( idx_new, rcut_new, nat, ityp, pos, lat, list=list )
```

# compile

```bash
make
```

clean
```bash
make clean
```

# run

```bash
cd test
./main.x < xx.xyz
```
