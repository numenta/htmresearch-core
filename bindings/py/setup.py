import os
import shutil
from setuptools import setup, find_packages
from distutils.core import Extension



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
  return parse_file(requirementsPath)



if __name__ == "__main__":

  # Copy the proto files into the proto Python package.
  destDir = os.path.relpath(os.path.join("htmresearch_core", "proto"))

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
    packages=find_packages(),
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
