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

cmake_minimum_required(VERSION 3.14)

# Downloads common CMake modules used throughout NWChemEx

#Sets the version to whatever git thinks it is
#include(get_version_from_git)
#get_version_from_git(scf_version "${CMAKE_CURRENT_LIST_DIR}")
project(scf VERSION "1.0.0" LANGUAGES CXX)
include(cmake/get_nwx_cmake.cmake)
include(nwx_versions)
include(get_cmaize)
include(nwx_cxx_api_docs)

set(SCF_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(SCF_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SCF_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")

nwx_cxx_api_docs("${SCF_SOURCE_DIR}" "${SCF_INCLUDE_DIR}")

### Options ###
cmaize_option_list(
    BUILD_TESTING OFF "Should we build the tests?"
    BUILD_PYBIND11_PYBINDINGS ON "Build Python bindings with pybind11?"
    BUILD_TAMM_SCF OFF "Should we build modules that rely on TAMM/Exachem?"
    INTEGRATION_TESTING OFF "Should we build the integration tests?"
)
if("${BUILD_TAMM_SCF}")
    set(DEPENDENCIES simde gauxc tamm exachem libint2 chemcache)
else()
    set(DEPENDENCIES simde gauxc)
endif()

foreach(dependency_i ${DEPENDENCIES})
    include(get_${dependency_i})
endforeach()

cmaize_add_library(
    scf
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/src/scf"
    INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include/scf"
    DEPENDS "${DEPENDENCIES}"
)

if("${BUILD_TAMM_SCF}")
    target_compile_definitions(scf PRIVATE BUILD_TAMM_SCF)
    cmaize_add_executable(
        scf_driver
        SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/examples/driver"
        DEPENDS scf
    )
endif()

include(nwx_pybind11)
nwx_add_pybind11_module(
    ${PROJECT_NAME}
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/src/python"
    DEPENDS "${PROJECT_NAME}"
)

if("${BUILD_TESTING}") 
    include(CTest)
    set(PYTHON_TEST_DIR "${SCF_TESTS_DIR}/python")
    set(CXX_TEST_DIR "${SCF_TESTS_DIR}/cxx")

    include(get_catch2)
    include(cmake/mpi_test.cmake)
    cxx_mpi_test(
        unit_test_scf
        SOURCE_DIR "${CXX_TEST_DIR}/unit_tests"
        INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/src/scf"
        DEPENDS Catch2 scf
    )
    python_mpi_test(
        unit_test_scf
        "${PYTHON_TEST_DIR}/unit_tests/run_unit_tests.py"
        SUBMODULES simde chemist pluginplay parallelzone
    )

    if("${INTEGRATION_TESTING}")
        include(get_nwchemex)
        cxx_mpi_test(
            integration_test_scf
            SOURCE_DIR "${CXX_TEST_DIR}/integration_tests"
            INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/src/scf"
            DEPENDS Catch2 nwchemex
        )

        python_mpi_test(
            integration_test_scf
            "${PYTHON_TEST_DIR}/integration_tests/run_integration_tests.py"
            SUBMODULES nwchemex chemcache simde chemist pluginplay parallelzone
        )
    endif()
endif()
