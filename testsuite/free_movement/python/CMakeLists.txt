#add_test(FreeMovementPy ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/free_movement.py)
#set_tests_properties(FreeMovementPy PROPERTIES ENVIRONMENT "${TEST_ENV}")

add_test(FreeMovementPyMPI_27_1_1 mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/free_movement.py 1 1)
add_test(FreeMovementPyMPI_27_1_2 mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/free_movement.py 1 2)
add_test(FreeMovementPyMPI_27_1_3 mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/free_movement.py 1 3)
add_test(FreeMovementPyMPI_27_2_2 mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/free_movement.py 2 2)
add_test(FreeMovementPyMPI_27_2_3 mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/free_movement.py 2 3)
add_test(FreeMovementPyMPI_27_3_3 mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/free_movement.py 3 3)
