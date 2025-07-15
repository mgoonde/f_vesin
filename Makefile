FC = gfortran
FFLAGS := -g

VESIN_PATH := /home/mgunde/vesin2
VESIN_BUILD := ${VESIN_PATH}/build
VESIN_FINCLUDE := ${VESIN_PATH}/fortran/include

SRC := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
OBJ=$(SRC)/Obj
MOD=$(SRC)/mod


f_main = main.f90
x_main = main.x

f_neigh = $(SRC)/m_neighbour.f90
o_neigh = $(OBJ)/m_neighbour.o

all: $(OBJ) $(MOD) $(x_main)

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
$(x_main): $(o_neigh) $(f_main)
	$(FC) $(FFLAGS) -I$(VESIN_FINCLUDE) -I$(MOD) -o $@ $^ -L$(VESIN_BUILD)/lib -lvesin -lstdc++


clean:
	rm -rf $(OBJ) $(MOD) $(x_main)
