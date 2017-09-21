# htmresearch-core

This repository contains the C++ source code for Numenta's [htmresearch](https://github.com/numenta/htmresearch) repository. Please read the htmresearch [README.md](https://github.com/numenta/htmresearch/blob/master/README.md) for more context. All of those disclaimers also apply to this repository.

### Build:

Environment:

 * `$NUPIC_CORE` is the current location of the [nupic.core](https://github.com/numenta/nupic.core) repository that you downloaded from GitHub.
 * `$HTMRESEARCH_CORE` is the current location of this repository that you downloaded from GitHub.

First, build [nupic.core](https://github.com/numenta/nupic.core).

Then:

    mkdir -p $HTMRESEARCH_CORE/build/scripts_release
    cd $HTMRESEARCH_CORE/build/scripts_release
    cmake ../.. -DCMAKE_INSTALL_PREFIX=../release -DCMAKE_BUILD_TYPE=Release -DNUPIC_IWYU=OFF -DLOCAL_NUPIC_CORE_INSTALL_DIR=$NUPIC_CORE/build/release -DPY_EXTENSIONS_DIR=$HTMRESEARCH_CORE/bindings/py/src/htmresearch_core
    make -j6
    make install

### Install htmresearch-core Python library:

    cd $HTMRESEARCH_CORE
    ARCHFLAGS="-arch x86_64" pip install --user -e .
