# Copyright 2024 NWChemEx-Project
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
    if("${BUILD_PYBIND11_PYBINDINGS}")
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
    endif()
endmacro()
