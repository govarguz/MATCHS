add_test(PolymerMeltMPI_27_dd_original_1_1   mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/polymer_melt_test.py 1 dd           1 1)
add_test(PolymerMeltMPI_27_dd_original_2_2   mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/polymer_melt_test.py 1 dd           2 2)
add_test(PolymerMeltMPI_27_dd_duplicate_1_1  mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/polymer_melt_test.py 1 dd_duplicate 1 1)
add_test(PolymerMeltMPI_27_dd_duplicate_2_2  mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/polymer_melt_test.py 1 dd_duplicate 2 2)
add_test(PolymerMeltMPI_27_dd_duplicate_1_2  mpirun -np 27 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/polymer_melt_test.py 1 dd_duplicate 1 2)

