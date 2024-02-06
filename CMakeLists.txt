# Copyright 2022 NWChemEx-Project
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
set(VERSION 1.0.0) #TODO: get from git
project(fastscf VERSION "${VERSION}" LANGUAGES CXX)

include(FetchContent)
enable_language(C)
FetchContent_Declare(
    nwx_cmake
    GIT_REPOSITORY https://github.com/NWChemEx/NWXCMake
)
FetchContent_MakeAvailable(nwx_cmake)
list(APPEND CMAKE_MODULE_PATH "${nwx_cmake_SOURCE_DIR}/cmake")

set(
    CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake"
    CACHE STRING "" FORCE
)

include(nwx_versions)
include(get_cmaize)
include(nwx_cxx_api_docs)
nwx_cxx_api_docs("${CMAKE_CURRENT_SOURCE_DIR}/fastscf")

set(DEPENDENCIES "")

cmaize_find_or_build_dependency(
    simde
    URL github.com/NWChemEx/SimDE
    PRIVATE TRUE
    VERSION ${NWX_SIMDE_VERSION}
    BUILD_TARGET simde
    FIND_TARGET nwx::simde
    CMAKE_ARGS BUILD_TESTING=OFF
)
list(APPEND DEPENDENCIES simde)

# Put after SimDE to allow TA to discover BLAS/LAPACK
cmaize_find_or_build_dependency(
    gauxc
    URL github.com/wavefunction91/GauXC
    VERSION ${NWX_GAUXC_VERSION}
    BUILD_TARGET gauxc
    FIND_TARGET gauxc::gauxc
    CMAKE_ARGS BUILD_TESTING=OFF GAUXC_ENABLE_HDF5=OFF
)
list(APPEND DEPENDENCIES gauxc)

cmaize_add_library(
    fastscf
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/src/fastscf"
    INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include/fastscf"
    DEPENDS "${DEPENDENCIES}"
)

cmaize_add_executable(
    fastscf_driver
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/examples/driver"
    DEPENDS fastscf
)
