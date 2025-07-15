FC = gfortran
FFLAGS := -g

VESIN_PATH     := $(HOME)/vesin2
VESIN_BUILD    := ${VESIN_PATH}/build
VESIN_LIBPATH  := ${VESIN_BUILD}/lib
VESIN_FINCLUDE := ${VESIN_PATH}/fortran/include

SRC := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
OBJ := $(SRC)/Obj
MOD := $(SRC)/mod


f_main = $(SRC)/test/main.f90
x_main = $(SRC)/test/main.x

f_neigh = $(SRC)/m_neighbour.f90
o_neigh = $(OBJ)/m_neighbour.o


all: obj test

obj: $(OBJ) $(MOD) $(o_neigh)
test: obj $(x_main)

#
# directories
#
$(OBJ):
	@if [ ! -d $(OBJ) ]; then mkdir $(OBJ) ; fi
$(MOD):
	@if [ ! -d $(MOD) ]; then mkdir $(MOD) ; fi

#
# object
#
$(o_neigh): $(f_neigh)
	$(FC) $(FFLAGS) -J$(MOD) -I$(VESIN_FINCLUDE) -c $^ -o $@


#
# main prog
#
$(x_main): $(f_main) $(o_neigh)
	$(FC) $(FFLAGS) -I$(MOD) -o $@ $^ -L$(VESIN_LIBPATH) -lvesin -lstdc++


clean:
	rm -rf $(OBJ) $(MOD) $(x_main)
