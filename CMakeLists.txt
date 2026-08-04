# Copyright 2026 NWChemEx-Project
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

# Downloads common CMake modules used throughout NWChemEx. Included before
# project() so nwx_set_version() is available to compute the version.
include(cmake/get_nwx_cmake.cmake)

# Version from scikit-build-core (wheels/sdists) or the latest git tag.
include(nwx_set_version)
nwx_set_version(scf_version "${CMAKE_CURRENT_LIST_DIR}")
project(scf VERSION "${scf_version}" LANGUAGES CXX)

include(disable_in_source_builds)
include(set_default_nwx_options)
option(BUILD_LIBXC             "Build libxc-backed XC functionals?"     ON)

# Documentation target (no-op unless BUILD_DOCS is on)
include(nwx_cxx_api_docs)
nwx_cxx_api_docs("cxx/include" "cxx/src")

### Dependencies ###
# simde is the only PUBLIC dep; gauxc/eigen/gau2grid/libxc are implementation
# details of the XC and eigensolver internals (PRIVATE).
set(scf_private_deps gauxc eigen gau2grid)
if(BUILD_LIBXC)
    list(APPEND scf_private_deps libxc)
endif()

include(get_dependencies)
get_dependencies(simde ${scf_private_deps})

### Library ###
include(nwx_library)
nwx_library(${PROJECT_NAME} "cxx/include" "cxx/src"
    simde PRIVATE ${scf_private_deps})
if(BUILD_LIBXC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE BUILD_LIBXC)
endif()

### Python bindings (no-op unless BUILD_PYBIND11_BINDINGS is on) ###
include(nwx_python_module)
nwx_python_module(${PROJECT_NAME} "cxx/src")

### Tests ###
if(BUILD_TESTING)
    include(nwx_find_mpi)
    nwx_find_mpi()

    # The unit tests white-box XC/eigensolver private headers that pull the
    # external headers (libxc's <xc.h>, gauxc, Eigen), so those externals must
    # be linked into the test executable in addition to scf itself.
    set(scf_test_libs
        ${PROJECT_NAME}
        ${NWX_DEP_TARGET_gauxc}
        ${NWX_DEP_TARGET_eigen}
        ${NWX_DEP_TARGET_gau2grid}
    )
    if(BUILD_LIBXC)
        list(APPEND scf_test_libs ${NWX_DEP_TARGET_libxc})
    endif()

    include(catch2_tests_from_dir)
    catch2_tests_from_dir(
        "test_unit_${PROJECT_NAME}"
        "tests/cxx/unit_tests"
        ${scf_test_libs}
        PRIVATE_INCLUDES "cxx/src" "cxx/src/${PROJECT_NAME}"
    )
    if(BUILD_LIBXC)
        target_compile_definitions(test_unit_${PROJECT_NAME} PRIVATE BUILD_LIBXC)
    endif()

    # The C++ unit tests can also be run under MPI with 2 ranks. Left wired but
    # disabled by default (mirrors ParallelZone); enable if a build needs it.
    include(nwx_mpi_test)
    #nwx_mpi_test(test_unit_${PROJECT_NAME}_mpi test_unit_${PROJECT_NAME} NPROC 2)

    include(nwx_python_test)
    nwx_python_test(
        py_unit_test_${PROJECT_NAME}
        "${CMAKE_CURRENT_LIST_DIR}/tests/python/unit_tests/run_unit_tests.py"
    )

    if(INTEGRATION_TESTING)
        get_dependencies(nwchemex)
        catch2_tests_from_dir(
            "test_integration_${PROJECT_NAME}"
            "tests/cxx/integration_tests"
            ${scf_test_libs} ${NWX_DEP_TARGET_nwchemex}
            PRIVATE_INCLUDES "cxx/src" "cxx/src/${PROJECT_NAME}"
        )
        if(BUILD_LIBXC)
            target_compile_definitions(
                test_integration_${PROJECT_NAME} PRIVATE BUILD_LIBXC
            )
        endif()
        nwx_python_test(
            py_integration_test_${PROJECT_NAME}
            "${CMAKE_CURRENT_LIST_DIR}/tests/python/integration_tests/run_integration_tests.py"
        )
        # nwchemex has no C++ of its own (LANGUAGES NONE), so it never
        # builds a "nwchemex_python" target -- nwx_python_test's PYTHONPATH
        # (built from NWX_PYTHON_MODULE_DIRS, which only tracks fetched
        # deps' compiled <dep>_python targets) doesn't include its
        # pure-Python python/nwchemex package. Append it directly so
        # tests/python/integration_tests/driver/test_scf_driver.py's
        # `import nwchemex` succeeds under ctest too (pytest picks it up
        # separately via tests/python/conftest.py).
        get_test_property(
            py_integration_test_${PROJECT_NAME} ENVIRONMENT _scf_pit_env
        )
        set_tests_properties(py_integration_test_${PROJECT_NAME} PROPERTIES
            ENVIRONMENT "${_scf_pit_env}:${nwchemex_SOURCE_DIR}/python"
        )
        unset(_scf_pit_env)
    endif()
endif()
