#!/usr/bin/env bash

export THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export PROJECT_ROOT=$THIS_DIR/..

export CXX=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-c++

export CXXFLAGS="-fvisibility-inlines-hidden -std=c++17 -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -I${CONDA_PREFIX}/include"

mkdir -p $THIS_DIR/build
cd $THIS_DIR/build

cmake $PROJECT_ROOT -DCMAKE_INSTALL_PREFIX=$THIS_DIR -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_CXX_FLAGS="${CXX_FLAGS}"
make VERBOSE=1
make install
ln -s $THIS_DIR/bin/messpf $THIS_DIR/bin/partition_function
ln -s $THIS_DIR/bin/messabs $THIS_DIR/bin/abstraction
ln -s $THIS_DIR/bin/messsym $THIS_DIR/bin/symmetry_number
