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

cmake_minimum_required(VERSION 2.8)
project(htmresearch_core CXX)

set(CMAKE_VERBOSE_MAKEFILE OFF)

set_directory_properties(PROPERTIES EP_BASE "${EP_BASE}")

include(ExternalProject)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")

include(nupic_core)

# Convenience variable that wraps all external include directories.
list(APPEND EXTERNAL_INCLUDE_DIRS ${NUPIC_CORE_INCLUDE_DIR})
list(APPEND EXTERNAL_INCLUDE_DIRS ${APRUTIL1_STATIC_LIB_INC_DIR})
list(APPEND EXTERNAL_INCLUDE_DIRS ${APR1_STATIC_LIB_INC_DIR} "${APR1_STATIC_LIB_INC_DIR}/apr-1")

set(EXTERNAL_INCLUDE_DIRS ${EXTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
