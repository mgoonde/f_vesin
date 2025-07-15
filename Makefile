FC = gfortran
FFLAGS := -g

SRC := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
OBJ := $(SRC)/Obj
MOD := $(SRC)/mod


VESIN_PATH     := $(SRC)/vesin
VESIN_BUILD    := ${VESIN_PATH}/build
VESIN_LIBPATH  := ${VESIN_BUILD}/lib
VESIN_FINCLUDE := ${VESIN_PATH}/fortran/include



f_main = $(SRC)/test/main.f90
x_main = $(SRC)/test/main.x

f_neigh = $(SRC)/m_neighbour.f90
o_neigh = $(OBJ)/m_neighbour.o


all: obj test

obj: vesin $(OBJ) $(MOD) $(o_neigh)
test: obj $(x_main)

#
# directories
#
$(OBJ):
	@if [ ! -d $(OBJ) ]; then mkdir $(OBJ) ; fi
$(MOD):
	@if [ ! -d $(MOD) ]; then mkdir $(MOD) ; fi
$(VESIN_BUILD):
	@if [ ! -d $(VESIN_BUILD) ]; then mkdir $(VESIN_BUILD) ; fi


#
# vesin
#
vesin: ${VESIN_BUILD}
	@cd ${VESIN_BUILD} && cmake -DCMAKE_INSTALL_PREFIX=./ .. && cmake --build . && cmake --install .

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
	rm -rf $(OBJ) $(MOD) $(x_main) $(VESIN_BUILD)
