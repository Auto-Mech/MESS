cmake_minimum_required(VERSION 3.22)
project(MESS LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(MESS_STATIC "Build fully static executables" OFF)
set(MESS_BLAS_SUFFIX "" CACHE STRING
    "Symbol suffix of the ILP64 BLAS/LAPACK library (64_ for OpenBLAS built with SYMBOLSUFFIX=64_)")

if(MESS_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  set(BLA_STATIC ON)
  add_link_options(-static)
endif()

# Dependencies
# MESS passes 64-bit integers to BLAS/LAPACK. For a library that FindBLAS cannot
# find by itself (e.g. libopenblas64_), pass -DBLAS_LIBRARIES and -DLAPACK_LIBRARIES.
set(BLA_SIZEOF_INTEGER 8)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

# SLATEC (vendored subset, default 32-bit Fortran integers)
file(GLOB SLATEC_SOURCES CONFIGURE_DEPENDS external/SLATEC/src/*.f)
add_library(slatec STATIC ${SLATEC_SOURCES})

# MESS library
add_library(messlibs STATIC
  src/libmess/atom.cc
  src/libmess/configuration.cc
  src/libmess/crossrate.cc
  src/libmess/d3.cc
  src/libmess/divsur.cc
  src/libmess/dynamic.cc
  src/libmess/dynlib.cc
  src/libmess/graph_common.cc
  src/libmess/graph_omp.cc
  src/libmess/io.cc
  src/libmess/key.cc
  src/libmess/lapack.cc
  src/libmess/limits.cc
  src/libmess/linpack.cc
  src/libmess/logical.cc
  src/libmess/math.cc
  src/libmess/mess.cc
  src/libmess/model.cc
  src/libmess/monom.cc
  src/libmess/mpack.cc
  src/libmess/multindex.cc
  src/libmess/permutation.cc
  src/libmess/potential.cc
  src/libmess/random.cc
  src/libmess/read.cc
  src/libmess/slatec.cc
  src/libmess/structure.cc
  src/libmess/symmetry.cc
  src/libmess/system.cc
  src/libmess/trajectory.cc
  src/libmess/units.cc
)
target_include_directories(messlibs PUBLIC src/libmess)
target_link_libraries(messlibs PUBLIC
  slatec LAPACK::LAPACK BLAS::BLAS GSL::gsl OpenMP::OpenMP_CXX ${CMAKE_DL_LIBS})

# Rename the BLAS/LAPACK routines MESS calls (dgemm_ -> dgemm_64_, ...)
if(MESS_BLAS_SUFFIX)
  foreach(f IN ITEMS
      dgelsd dgemv dgesdd dgesv dgetrf dgetri dgetrs dppsv dpptrf dpptri dpptrs
      dsbevd dspgvd dspmv dspsv dsptrf dsptri dsptrs dsyev dsyevd zgemm zhpev zhpgvd)
    target_compile_definitions(messlibs PUBLIC ${f}_=${f}_${MESS_BLAS_SUFFIX})
  endforeach()
endif()

# Executables
add_executable(mess src/mess_driver.cc)
add_executable(messpf src/partition_function.cc)
add_executable(messabs src/abstraction.cc)
foreach(exe IN ITEMS mess messpf messabs)
  target_link_libraries(${exe} PRIVATE messlibs)
endforeach()

include(GNUInstallDirs)
install(TARGETS mess messpf messabs)
