FC = gfortran
FFLAGS := -g

SRC := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
OBJ := $(SRC)/Obj
MOD := $(SRC)/mod


VESIN_PATH     := $(SRC)/vesin
VESIN_BUILD    := ${VESIN_PATH}/build
VESIN_LIBPATH  := ${VESIN_PATH}/lib
VESIN_INCLUDE  := ${VESIN_BUILD}/fortran/include



f_main = $(SRC)/test/main.f90
x_main = $(SRC)/test/main.x

f_neigh = $(SRC)/m_neighbour.f90
o_neigh = $(OBJ)/m_neighbour.o


all: obj test

obj: submod vesin $(OBJ) $(MOD) $(o_neigh)
test: obj $(x_main)

#
# directories
#
$(OBJ):
	@if [ ! -d $(OBJ) ]; then mkdir $(OBJ) ; fi
$(MOD):
	@if [ ! -d $(MOD) ]; then mkdir $(MOD) ; fi


#
# vesin
#
vesin: submod
	@cd ${VESIN_PATH} && cmake -B ${VESIN_BUILD} -DVESIN_FORTRAN=ON -DCMAKE_INSTALL_PREFIX=./
	@cd ${VESIN_PATH} && cmake --build ${VESIN_BUILD}
	@cd ${VESIN_PATH} && cmake --install ${VESIN_BUILD}

submod:
	@if test ! -d $(VESIN_PATH)/vesin; then \
	git submodule update --init --recursive; fi

#
# object
#
$(o_neigh): $(f_neigh)
	$(FC) $(FFLAGS) -J$(MOD) -I$(VESIN_INCLUDE) -c $^ -o $@


#
# main prog
#
$(x_main): $(f_main) $(o_neigh)
	$(FC) $(FFLAGS) -I$(MOD) -o $@ $^ -L$(VESIN_LIBPATH) -lvesin -lstdc++


clean:
	rm -rf $(OBJ) $(MOD) $(x_main) $(VESIN_BUILD) ${VESIN_LIBPATH} ${VESIN_INCLUDE}
