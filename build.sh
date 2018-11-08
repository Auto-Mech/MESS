#!/usr/bin/env bash

mkdir build
cd build

if [ "$(uname)" == "Darwin" ]; then
    cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$GFORTRAN -DBLAS_LIBRARIES=$PREFIX/lib/libopenblas.dylib -DLAPACK_LIBRARIES=$PREFIX/lib/libopenblas.dylib
elif [ "$(uname)" == "Linux" ]; then
    cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$GFORTRAN -DBLAS_LIBRARIES=$PREFIX/lib/libopenblas.so -DLAPACK_LIBRARIES=$PREFIX/lib/libopenblas.so
fi

make
make install
ln -s $PREFIX/bin/messpf $PREFIX/bin/partition_function
ln -s $PREFIX/bin/messabs $PREFIX/bin/abstraction
ln -s $PREFIX/bin/messsym $PREFIX/bin/symmetry_number
