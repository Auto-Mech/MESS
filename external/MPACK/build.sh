#!/usr/bin/env bash

export CXXFLAGS="${CXXFLAGS} -std=gnu++98"
export BLAS_LIBS=$CONDA_PREFIX/lib/libopenblas.so
export LAPACK_LIBS=$CONDA_PREFIX/lib/libopenblas.so

./configure --prefix=$PREFIX --disable-dependency-tracking --disable-testing

make VERBOSE=1
make install
