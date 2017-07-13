# -----------------------------------------------------------------------------
# Numenta Platform for Intelligent Computing (NuPIC)
# Copyright (C) 2017, Numenta, Inc.  Unless you have purchased from
# Numenta, Inc. a separate commercial license for this software code, the
# following terms and conditions apply:
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero Public License version 3 as
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero Public License for more details.
#
# You should have received a copy of the GNU Affero Public License
# along with this program.  If not, see http://www.gnu.org/licenses.
#
# http://numenta.org/licenses/
# -----------------------------------------------------------------------------

# Configure nupic_core and its third party dependencies
#
# INPUT VARIABLES: Passed via -DVAR=VALUE
#    LOCAL_NUPIC_CORE_INSTALL_DIR: Use local nupic.core installed at this location
#                                  instead of building
#    LOCAL_NUPIC_CORE_SOURCE: If the LOCAL_NUPIC_CORE_SOURCE variable is set to
#                             your local nupic_core sources, cmake will use your 
#                             for the version specified in NUPIC_CORE_VERSION
# 
#    NUPIC_CORE_VERSION: Set NUPIC_CORE_VERSION to any released nupic_core 
#                        version to override the version cmake will download.
#                        Default value: "1.0.1"
# 
# OUPUT VARIABLES: Available to the parent scope
#
#   NUPIC_CORE_INCLUDE_DIR: nupic_core include directory
#   NUPIC_CORE_LIB_DIR: nupic_core libraries directory
#   NUPIC_CORE_STATIC_LIB_TARGET: nupic_core static library target
#   GTEST_STATIC_LIB_TARGET: gtest static library target installed by nupic_core
#   SWIG_EXECUTABLE: swig executable installed by nupic_core
#   SWIG_DIR: swig shared directory installed by nupic_core 
# 
# See CapnProto.cmake for other variables and functions  
#

# nupic.core version to use
if (NOT NUPIC_CORE_VERSION)
    set(NUPIC_CORE_VERSION "1.0.1")
endif(NOT NUPIC_CORE_VERSION)

# Check if we should use local sources or download from the releases site 
if(LOCAL_NUPIC_CORE_SOURCE)
    set(nupic_core_url)
    set(NUPIC_CORE_SOURCE_DIR "${LOCAL_NUPIC_CORE_SOURCE}")
else(LOCAL_NUPIC_CORE_SOURCE)
    set(nupic_core_url  "http://releases.numenta.org/nupic.core/${NUPIC_CORE_VERSION}/nupic.core-source-${NUPIC_CORE_VERSION}.tar.gz")
    set(NUPIC_CORE_SOURCE_DIR "${EP_BASE}/Source/nupic_core")
endif(LOCAL_NUPIC_CORE_SOURCE)

# Use nupic.core local installation
if(LOCAL_NUPIC_CORE_INSTALL_DIR)
    add_custom_target(nupic_core)
    # Validate Nupic installation
    if(NOT EXISTS "${LOCAL_NUPIC_CORE_INSTALL_DIR}/include/nupic")
        message(FATAL_ERROR "Invalid nupic.core installation. \
        Make sure nupic.core is installed the LOCAL_NUPIC_CORE_INSTALL_DIR variable \
        points to a valid nupic.core installation")
    endif(NOT EXISTS "${LOCAL_NUPIC_CORE_INSTALL_DIR}/include/nupic")
    set(NUPIC_CORE_INSTALL_DIR "${LOCAL_NUPIC_CORE_INSTALL_DIR}")
else(LOCAL_NUPIC_CORE_INSTALL_DIR)
    # Update file locations based on the external project location
    set(NUPIC_CORE_BUILD_DIR "${EP_BASE}/Build/nupic_core")
    set(NUPIC_CORE_INSTALL_DIR "${EP_BASE}/Install/nupic_core")
    set(NUPIC_CORE_THIRDPARTY_DIR "${NUPIC_CORE_BUILD_DIR}/ThirdParty/Install")

    # Download and build nupic.core
    set_directory_properties(PROPERTIES EP_BASE "${EP_BASE}")
    ExternalProject_Add(
        nupic_core
        SOURCE_DIR ${LOCAL_NUPIC_CORE_SOURCE}
        URL ${nupic_core_url}
        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        CMAKE_GENERATOR ${CMAKE_GENERATOR}
        CMAKE_ARGS
            ${CAPNP_CMAKE_DEFINITIONS}
            -DNUPIC_BUILD_PYEXT_MODULES=OFF
            -DCMAKE_INSTALL_PREFIX=${NUPIC_CORE_INSTALL_DIR}
    )
endif(LOCAL_NUPIC_CORE_INSTALL_DIR)

set(NUPIC_CORE_INCLUDE_DIR "${NUPIC_CORE_INSTALL_DIR}/include")
set(NUPIC_CORE_LIB_DIR "${NUPIC_CORE_INSTALL_DIR}/lib")

# Expose nupic.core directories and targets
set(NUPIC_CORE_INCLUDE_DIR "${NUPIC_CORE_INCLUDE_DIR}" PARENT_SCOPE)
set(NUPIC_CORE_LIB_DIR "${NUPIC_CORE_LIB_DIR}" PARENT_SCOPE)
set(NUPIC_CORE_STATIC_LIB_TARGET "${NUPIC_CORE_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}nupic_core${CMAKE_STATIC_LIBRARY_SUFFIX}" PARENT_SCOPE)
set(NUPIC_CORE_INSTALL_DIR "${NUPIC_CORE_INSTALL_DIR}" PARENT_SCOPE)

# Gtest 
set(GTEST_STATIC_LIB_TARGET "${NUPIC_CORE_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}" PARENT_SCOPE)

# Swig
set(SWIG_EXECUTABLE "${NUPIC_CORE_INSTALL_DIR}/bin/swig${CMAKE_EXECUTABLE_SUFFIX}" PARENT_SCOPE)
set(SWIG_DIR "${NUPIC_CORE_INSTALL_DIR}/share/swig/3.0.2" PARENT_SCOPE)

# Capnp
include(CapnProto)
