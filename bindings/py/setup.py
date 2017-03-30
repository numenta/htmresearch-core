import os
import shutil
from setuptools import setup, find_packages
from distutils.core import Extension



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
  )
