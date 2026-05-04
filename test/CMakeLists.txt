add_executable(main main.f90)
target_link_libraries(main PRIVATE m_neighbour)
set_target_properties(main PROPERTIES LINKER_LANGUAGE Fortran)
