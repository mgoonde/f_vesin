cmake_minimum_required(VERSION 3.21)

project(f_vesin LANGUAGES C CXX Fortran)

set(VESIN_FORTRAN ON CACHE BOOL "Build fortran interface to vesin")

include(FetchContent)
FetchContent_Declare( vesin_git
  GIT_REPOSITORY https://github.com/mgoonde/vesin.git
  GIT_TAG doc_fortran
  GIT_SHALLOW 1
)
FetchContent_MakeAvailable(vesin_git)

# add_subdirectory(vesin)

add_library( m_neighbour OBJECT m_neighbour.f90 )

# enable_language(CXX)
target_link_libraries( m_neighbour PUBLIC vesin_fortran )
target_include_directories( m_neighbour INTERFACE ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(test)
