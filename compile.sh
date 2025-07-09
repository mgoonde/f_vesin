g++ -I./ -g -c vesin-single-build.cpp
gfortran -g -c f_vesin_wrapper.f90
gfortran -I./ -g -o main.x main.f90 vesin-single-build.o f_vesin_wrapper.o -lstdc++
