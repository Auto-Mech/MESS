#!/usr/bin/env bash

export THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export PROJECT_ROOT=$THIS_DIR/..

export CXXFLAGS="-fvisibility-inlines-hidden -std=c++17 -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -I${CONDA_PREFIX}/include"

cd $PROJECT_ROOT
export CXXFLAGS="${CXXFLAGS} -std=gnu++98"
export BLAS_LIBS=$CONDA_PREFIX/lib/libopenblas.so
export LAPACK_LIBS=$CONDA_PREFIX/lib/libopenblas.so

./configure --prefix=$CONDA_PREFIX --disable-dependency-tracking

make VERBOSE=1
make install
