# -----------------------------------------------------------------------------
# Numenta Platform for Intelligent Computing (NuPIC)
# Copyright (C) 2015, Numenta, Inc.  Unless you have purchased from
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

# Expose Cap'n Proto from nupic_core external project.
#
# OUTPUT VARIABLES:
#
#   CAPNP_STATIC_LIB_TARGET: name of static library target that contains all of
#                            capnproto library objects.
#
#   CAPNP_INCLUDE_DIRS
#   CAPNP_EXECUTABLE
#   CAPNPC_CXX_EXECUTABLE
#   CAPNP_CMAKE_DEFINITIONS: informational; platform-specific cmake defintions
#                            used by capnproto build
#   CAPNP_COMPILER_DEFINITIONS: list of -D compiler defintions needed by apps
#                               that are built against this library (e.g.,
#                               -DCAPNP_LITE)
#   CAPNP_BINARIES:          Binaries location
#
# EXPORTED FUNCTIONS:
#
#   CREATE_CAPNPC_COMMAND: Create a custom command that runs the capnp compiler.

set_directory_properties(PROPERTIES EP_BASE "${EP_BASE}")

# Output static library target for linking and dependencies
set(CAPNP_STATIC_LIB_TARGET capnp-bundle)

set(CAPNP_INCLUDE_DIRS ${NUPIC_CORE_THIRDPARTY_DIR}/include)
set(CAPNP_EXECUTABLE ${NUPIC_CORE_THIRDPARTY_DIR}/bin/capnp${CMAKE_EXECUTABLE_SUFFIX})
set(CAPNPC_CXX_EXECUTABLE ${NUPIC_CORE_THIRDPARTY_DIR}/bin/capnpc-c++${CMAKE_EXECUTABLE_SUFFIX})

set(CAPNP_COMPILER_DEFINITIONS)


if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(CAPNP_CMAKE_DEFINITIONS -DCAPNP_LITE=1 -DEXTERNAL_CAPNP=1 -DBUILD_TOOLS=OFF)
  # NOTE nupic.core's swig wraps depend on the macro CAPNP_LITE to have a value
  set(CAPNP_COMPILER_DEFINITIONS ${CAPNP_COMPILER_DEFINITIONS} -DCAPNP_LITE=1)
else()
  set(CAPNP_CMAKE_DEFINITIONS -DCAPNP_LITE=0)
endif()


function(CREATE_CAPNPC_COMMAND
         SPEC_FILES SRC_PREFIX INCLUDE_DIR TARGET_DIR OUTPUT_FILES)
  # Create a custom command that runs the capnp compiler on ${SPEC_FILES} and
  # generates ${OUTPUT_FILES} in directory ${TARGET_DIR}

  # nupic_core external project must be installed prior to using capnp 
  set(dependencies ${SPEC_FILES} nupic_core)

  add_custom_command(
    OUTPUT ${OUTPUT_FILES}
    COMMAND ${CAPNP_EXECUTABLE}
        compile -o ${CAPNPC_CXX_EXECUTABLE}:${TARGET_DIR}
        --src-prefix ${SRC_PREFIX} 
        -I ${INCLUDE_DIR}
        -I ${NUPIC_CORE_INCLUDE_DIR}
        ${SPEC_FILES}
    DEPENDS ${dependencies}
    COMMENT "Executing Cap'n Proto compiler"
  )
endfunction(CREATE_CAPNPC_COMMAND)

# Set the relevant variables in the parent scope.
set(CAPNP_STATIC_LIB_TARGET ${CAPNP_STATIC_LIB_TARGET} PARENT_SCOPE)
set(CAPNP_INCLUDE_DIRS ${CAPNP_INCLUDE_DIRS} PARENT_SCOPE)
set(CAPNP_EXECUTABLE ${CAPNP_EXECUTABLE} PARENT_SCOPE)
set(CAPNPC_CXX_EXECUTABLE ${CAPNPC_CXX_EXECUTABLE} PARENT_SCOPE)
set(CAPNP_CMAKE_DEFINITIONS ${CAPNP_CMAKE_DEFINITIONS} PARENT_SCOPE)
set(CAPNP_COMPILER_DEFINITIONS ${CAPNP_COMPILER_DEFINITIONS} PARENT_SCOPE)
set(CAPNP_BINARIES ${NUPIC_CORE_THIRDPARTY_DIR}/bin PARENT_SCOPE)
