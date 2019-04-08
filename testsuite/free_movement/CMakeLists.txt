add_subdirectory(python)

#add_executable(free_movement free_movement.cpp particle_set.cpp)
#target_link_libraries(free_movement ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${FFTW3_LIBRARIES} ${VAMPIRTRACE_LIBRARIES} _espressopp)

#add_executable(simple_force simple_force.cpp particle_set.cpp)
#target_link_libraries(simple_force ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${FFTW3_LIBRARIES} ${VAMPIRTRACE_LIBRARIES} _espressopp)

#add_executable(free_movement_halfcell free_movement_halfcell.cpp particle_set.cpp)
#target_link_libraries(free_movement_halfcell ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${FFTW3_LIBRARIES} ${VAMPIRTRACE_LIBRARIES} _espressopp)

#add_executable(simple_force simple_force.cpp particle_set.cpp)
#target_link_libraries(simple_force ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${FFTW3_LIBRARIES} ${VAMPIRTRACE_LIBRARIES} _espressopp)


#set (test_parameters -np 12 "./free_movement")
#add_test(NAME test_free_movement COMMAND "mpirun" ${test_parameters})

#set (test_parameters -np 5 "./simple_force")
#add_test(NAME test_simple_force COMMAND "mpirun" ${test_parameters})

#set (test_parameters -np 36 "./free_movement_halfcell")
#add_test(NAME test_free_movement_halfcell COMMAND "mpirun" ${test_parameters})

#set (test_parameters -np 5 "./simple_force")
#add_test(NAME test_simple_force COMMAND "mpirun" ${test_parameters})
