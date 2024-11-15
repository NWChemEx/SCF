macro(cxx_mpi_test test_name)
    # Need to call with MPI, so we have to deconstruct the cmaize_add_test and
    # nwx_pybind11_tests methods
    cmaize_add_executable("${test_name}" ${ARGN})
    add_test(
        NAME "${test_name}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
                "${CMAKE_BINARY_DIR}/${test_name}"
    )
endmacro()

macro(python_mpi_test test_name test_script)
    add_test(
        NAME "py_${test_name}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
                "${Python_EXECUTABLE}" 
                "${test_script}"
    )
    nwx_python_path(TEST_PYTHONPATH ${ARGN})
    set_tests_properties(
        "py_${test_name}"
        PROPERTIES ENVIRONMENT "${TEST_PYTHONPATH}"
    )
endmacro()