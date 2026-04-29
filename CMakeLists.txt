cmake_minimum_required(VERSION 3.21)

project(f_vesin LANGUAGES C Fortran)


set(VESIN_FORTRAN ON CACHE BOOL "Build fortran interface to vesin")
add_subdirectory(vesin)

add_library( m_neighbour OBJECT m_neighbour.f90 )
target_link_libraries( m_neighbour vesin )
