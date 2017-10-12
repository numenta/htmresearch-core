# ----------------------------------------------------------------------
# Numenta Platform for Intelligent Computing (NuPIC)
# Copyright (C) 2017, Numenta, Inc.  Unless you have an agreement
# with Numenta, Inc., for a separate license for this software code, the
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

import os
import shutil
import pkg_resources

from setuptools import setup, find_packages
from distutils.core import Extension


PY_BINDINGS = os.path.dirname(os.path.realpath(__file__))
REPO_DIR = os.path.abspath(os.path.join(PY_BINDINGS, os.pardir, os.pardir))
DARWIN_PLATFORM = "darwin"
LINUX_PLATFORM = "linux"
UNIX_PLATFORMS = [LINUX_PLATFORM, DARWIN_PLATFORM]
WINDOWS_PLATFORMS = ["windows"]



def nupicBindingsPrereleaseInstalled():
  """
  Make an attempt to determine if a pre-release version of nupic.bindings 
  is installed already.

  @return: boolean
  """
  try:
    nupicDistribution = pkg_resources.get_distribution("nupic.bindings")
    if pkg_resources.parse_version(nupicDistribution.version).is_prerelease:
      # A pre-release dev version of nupic.bindings is installed.
      return True
  except pkg_resources.DistributionNotFound:
    pass  # Silently ignore.  The absence of nupic.bindings will be handled by
    # setuptools by default

  return False


def getVersion():
  """
  Get version from local file.
  """
  with open(os.path.join(REPO_DIR, "VERSION"), "r") as versionFile:
    return versionFile.read().strip()



def parse_file(requirementFile):
  try:
    return [
      line.strip()
      for line in open(requirementFile).readlines()
      if not line.startswith("#")
    ]
  except IOError:
    return []



def findRequirements():
  """
  Read the requirements.txt file and parse into requirements for setup's
  install_requirements option.
  """
  PY_BINDINGS = os.path.dirname(os.path.realpath(__file__))
  REPO_DIR = os.path.abspath(os.path.join(PY_BINDINGS, os.pardir, os.pardir))
  requirementsPath = os.path.join(REPO_DIR, "requirements.txt")
  requirements = parse_file(requirementsPath)

  if nupicBindingsPrereleaseInstalled():
    # User has a pre-release version of nupic.bindings installed, which is only
    # possible if the user installed and built nupic.bindings from source and
    # it is up to the user to decide when to update nupic.bindings.  We'll
    # quietly remove the entry in requirements.txt so as to not conflate the
    # two.
    requirements = [req for req in requirements 
                    if "nupic.bindings" not in req]

  return requirements



if __name__ == "__main__":

  PY_BINDINGS = os.path.dirname(os.path.realpath(__file__))
  REPO_DIR = os.path.abspath(os.path.join(PY_BINDINGS, os.pardir, os.pardir))

  # Copy the proto files into the proto Python package.
  destDir = os.path.relpath(os.path.join(PY_BINDINGS, "src", "htmresearch_core", "proto"))

  # Copy the proto files into the proto Python package.
  protoSchemasToCopy = [
    "ApicalTiebreakTemporalMemoryProto.capnp"
  ]


  for schemaFilename in protoSchemasToCopy:
    protoPath = os.path.relpath(os.path.join(REPO_DIR, "src", "nupic", "proto", schemaFilename))
    shutil.copy(protoPath, destDir)

  setup(
    # setuptools replaces "_" with "-", so package name is "htmresearch-core",
    # import namespace "htmresearch_core".
    name="htmresearch-core",
    version=getVersion(),
    # This distribution contains platform-specific C++ libraries, but they are not
    # built with distutils. So we must create a dummy Extension object so when we
    # create a binary file it knows to make it platform-specific.
    ext_modules=[Extension("htmresearch_core.dummy", sources=["dummy.c"])],
    install_requires=findRequirements(),
    package_dir = {"": "src"},
    packages=find_packages("src"),
    package_data={
      "htmresearch_core.proto": ["*.capnp"],
      "htmresearch_core": ["*.so", "*.pyd"],
    },
    description="Numenta's experimental C++ research code",
    author="Numenta",
    author_email="help@numenta.org",
    url="https://github.com/numenta/htmresearch-core",
    long_description = "Python bindings for htmresearch-core.",
    classifiers=[
      "Programming Language :: Python",
      "Programming Language :: Python :: 2",
      "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
      "Operating System :: MacOS :: MacOS X",
      "Operating System :: POSIX :: Linux",
      "Operating System :: Microsoft :: Windows",
      # It has to be "5 - Production/Stable" or else pypi rejects it!
      "Development Status :: 5 - Production/Stable",
      "Environment :: Console",
      "Intended Audience :: Science/Research",
      "Topic :: Scientific/Engineering :: Artificial Intelligence"
    ],
  )
