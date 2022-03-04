#!/usr/bin/env bash

mkdir -p build
cd build

echo $PREFIX

export CFLAGS="${CFLAGS} -I${CONDA_PREFIX}/include"
export CXXFLAGS="${CXXFLAGS} -I${CONDA_PREFIX}/include -DWITH_MPACK"

cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX

make VERBOSE=1
make install
ln -s $PREFIX/bin/messpf $PREFIX/bin/partition_function
ln -s $PREFIX/bin/messabs $PREFIX/bin/abstraction
ln -s $PREFIX/bin/messsym $PREFIX/bin/symmetry_number
