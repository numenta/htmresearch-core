import os
import shutil
from setuptools import setup, find_packages
from distutils.core import Extension
import pkg_resources


PY_BINDINGS = os.path.dirname(os.path.realpath(__file__))
REPO_DIR = os.path.abspath(os.path.join(PY_BINDINGS, os.pardir, os.pardir))



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

  # Copy the proto files into the proto Python package.
  destDir = os.path.relpath(os.path.join("src", "htmresearch_core", "proto"))

  # Copy the proto files into the proto Python package.
  protoSchemasToCopy = [
    "ExtendedTemporalMemoryProto.capnp"
  ]

  for schemaFilename in protoSchemasToCopy:
    protoPath = os.path.relpath(os.path.join("..", "..", "src", "nupic", "proto", schemaFilename))
    shutil.copy(protoPath, destDir)

  setup(
    # setuptools replaces "_" with "-", so package name is "htmresearch-core",
    # import namespace "htmresearch_core".
    name="htmresearch-core",
    version=getVersion(),
    package_dir = {"": "src"},
    packages=find_packages("src"),
    package_data={
      "htmresearch_core.proto": ["*.capnp"],
      "htmresearch_core": ["*.so", "*.pyd"],
    },
    # This distribution contains platform-specific C++ libraries, but they are not
    # built with distutils. So we must create a dummy Extension object so when we
    # create a binary file it knows to make it platform-specific.
    ext_modules=[Extension("htmresearch_core.dummy", sources=["dummy.c"])],
    install_requires=findRequirements()
  )
