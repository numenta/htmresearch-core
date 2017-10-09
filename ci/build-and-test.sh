#!/bin/bash
# ----------------------------------------------------------------------
# Numenta Platform for Intelligent Computing (NuPIC)
# Copyright (C) 2016, Numenta, Inc.  Unless you have purchased from
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
# ----------------------------------------------------------------------

set -o errexit


USAGE="Usage:

[BUILD_TYPE=Release | Debug] [WHEEL_PLAT=platform] $( basename ${0} )

This script builds and tests the nupic.bindings Python extension.

In Debug builds, also
  - Turns on the Include What You Use check (assumes iwyu is installed)

ASUMPTION: Expects a pristine htmresearch-core source tree without any remnant 
   build artifacts from prior build attempts. Otherwise, behavior is undefined.


INPUT ENVIRONMENT VARIABLES:

  BUILD_TYPE : Specifies build type, which may be either Release or Debug;
               defaults to Release. [OPTIONAL]
  WHEEL_PLAT : Wheel platform name; pass manylinux1_x86_64 for manylinux build;
               leave undefined for all other builds.

OUTPUTS:
  htmresearch-core wheel: On success, the resulting wheel will be located in the
                        subdirectory htmresearch_core_wheelhouse of the source
                        tree's root directory.

  test results: htmresearch-core test results will be located in the subdirectory
                test_results of the source tree's root directory with the
                the following content:

                junit-test-results.xml
                htmlcov/

"

if [[ $1 == --help ]]; then
  echo "${USAGE}"
  exit 0
fi

if [[ $# > 0 ]]; then
  echo "ERROR Unexpected arguments: ${@}" >&2
  echo "${USAGE}" >&2
  exit 1
fi


set -o xtrace


# Apply defaults
BUILD_TYPE=${BUILD_TYPE-"Release"}


HTMRESEARCH_CORE_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

DEST_WHEELHOUSE="${HTMRESEARCH_CORE_ROOT}/htmresearch_core_wheelhouse"

TEST_RESULTS_DIR="${HTMRESEARCH_CORE_ROOT}/test_results"

echo "RUNNING HTMRESEARCH-CORE BUILD: BUILD_TYPE=${BUILD_TYPE}, " \
     "DEST_WHEELHOUSE=${DEST_WHEELHOUSE}" >&2

# Install htmresearch-core dependencies; the htmresearch-core cmake build
# depends on some of them (e.g., numpy).
pip install \
    --ignore-installed \
    -r ${HTMRESEARCH_CORE_ROOT}/bindings/py/requirements.txt

#
# Build htmresearch-core
#

# NOTE without -p to force build failure upon pre-existing build side-effects
mkdir ${HTMRESEARCH_CORE_ROOT}/build
mkdir ${HTMRESEARCH_CORE_ROOT}/build/scripts

cd ${HTMRESEARCH_CORE_ROOT}/build/scripts

# Configure htmresearch-core build
if [[ "$BUILD_TYPE" == "Debug" ]]; then
  EXTRA_CMAKE_DEFINITIONS="-DNTA_COV_ENABLED=ON"

  # Only add iwyu for clang builds
  if [[ $CC == *"clang"* ]]; then
    EXTRA_CMAKE_DEFINITIONS="-DNUPIC_IWYU=ON ${EXTRA_CMAKE_DEFINITIONS}"
  fi
fi

cmake ${HTMRESEARCH_CORE_ROOT} \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    ${EXTRA_CMAKE_DEFINITIONS} \
    -DCMAKE_INSTALL_PREFIX=${HTMRESEARCH_CORE_ROOT}/build/release \
    -DPY_EXTENSIONS_DIR=${HTMRESEARCH_CORE_ROOT}/bindings/py/src/htmresearch_core

# Build htmresearch-core
make install

# Build htmresearch-core python extensions from htmresearch-core build artifacts
if [[ $WHEEL_PLAT ]]; then
  EXTRA_WHEEL_OPTIONS="--plat-name ${WHEEL_PLAT}"
fi

cd ${HTMRESEARCH_CORE_ROOT}
python setup.py bdist_wheel --dist-dir ${DEST_WHEELHOUSE} ${EXTRA_WHEEL_OPTIONS}


#
# Test
#

# Install htmresearch-core before running c++ tests; py_region_test depends on it
pip install \
    --ignore-installed \
    ${DEST_WHEELHOUSE}/htmresearch_core-*.whl


mkdir ${TEST_RESULTS_DIR}
cd ${TEST_RESULTS_DIR}    # so that py.test will deposit its artifacts here

# Run the htmresearch-core c++ tests
${HTMRESEARCH_CORE_ROOT}/build/release/bin/unit_tests --gtest_output=xml:unit_tests_report.xml

# Run tests with pytest options per htmresearch-core/setup.cfg
py.test ${HTMRESEARCH_CORE_ROOT}/bindings/py/tests
