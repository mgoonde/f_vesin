g++ -I./ -g -c vesin-single-build.cpp
gfortran -g -c f_vesin.f90
gfortran -g -o main.x main.f90 vesin-single-build.o f_vesin.o -lstdc++ -lc
