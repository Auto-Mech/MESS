project(MESS C CXX)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Simple options
option(USE_MKL "Use Intel MKL instead of OpenBLAS" OFF)
option(USE_INT64 "Use 64-bit integers for BLAS/LAPACK (ILP64 interface)" OFF)
option(USE_DIRECT_CBLAS "Use direct CBLAS instead of GSL's wrapper" OFF)
set(CBLAS_LIBRARY_PATH "" CACHE STRING "Path to specific CBLAS library")

# Set up BLAS/LAPACK libraries
if(USE_MKL)
  message(STATUS "Using Intel MKL")
  if(USE_INT64)
    set(BLAS_LIBRARIES 
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_ilp64.so
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_sequential.so
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so
      m dl
    )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64")
  else()
    set(BLAS_LIBRARIES 
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_sequential.so
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so
      m dl
    )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_32BIT_INTEGERS")
  endif()
  include_directories(/opt/intel/oneapi/mkl/latest/include)
else()
  message(STATUS "Using OpenBLAS")
  find_library(OPENBLAS_LIB NAMES openblas PATHS
    /usr/lib/x86_64-linux-gnu
    /usr/local/lib
  )
  if(OPENBLAS_LIB)
    set(BLAS_LIBRARIES ${OPENBLAS_LIB})
    if(USE_INT64)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_32BIT_INTEGERS")
    endif()
    else()
    message(FATAL_ERROR "OpenBLAS not found")
  endif()
endif()

set(LAPACK_LIBRARIES ${BLAS_LIBRARIES})

# CBLAS will be handled by GSL's internal configuration
set(CBLAS_LIBRARIES "")
message(STATUS "CBLAS will be configured through GSL build")

# Build GSL from source with OpenBLAS/MKL as CBLAS backend
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)
  include(ExternalProject)

set(GSL_URL "https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz")

message(STATUS "Building GSL from source with OpenBLAS/MKL as CBLAS backend")

# Set up CBLAS configuration for GSL
if(USE_MKL)
  set(GSL_CBLAS_CONFIG "--with-cblas=/opt/intel/oneapi/mkl/latest/lib/intel64")
  set(GSL_CBLAS_LIBS "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
  if(USE_INT64)
    set(GSL_CBLAS_LIBS "-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core")
endif()
else()
  # Use OpenBLAS as CBLAS backend
  set(GSL_CBLAS_CONFIG "--with-cblas=${OPENBLAS_LIB}")
  set(GSL_CBLAS_LIBS "${OPENBLAS_LIB}")
endif()

  ExternalProject_Add(
    gsl_external
    URL ${GSL_URL}
  CONFIGURE_COMMAND <SOURCE_DIR>/configure 
    --prefix=${EXTERNAL_INSTALL_LOCATION}
    --disable-shared
    --enable-static
    ${GSL_CBLAS_CONFIG}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )

  set(GSL_INCLUDE_DIRS ${EXTERNAL_INSTALL_LOCATION}/include)
set(GSL_LIBRARIES ${EXTERNAL_INSTALL_LOCATION}/lib/libgsl.a)
  include_directories(${GSL_INCLUDE_DIRS})
set(GSL_FOUND TRUE)
  set(GSL_DOWNLOADED TRUE)

# Find other libraries
find_library(QD_LIBRARIES NAMES qd libqd QUIET)
find_library(SLATEC_LIBRARIES NAMES slatec libslatec QUIET)

# Make SLATEC optional - if not found, create empty library
if(NOT SLATEC_LIBRARIES)
  message(STATUS "SLATEC not found, creating stub library")
  add_library(slatec_stub STATIC)
  target_sources(slatec_stub PRIVATE
    ${CMAKE_SOURCE_DIR}/slatec_stub.c
  )
  set(SLATEC_LIBRARIES slatec_stub)
endif()

# Display configuration
message(STATUS "Configuration:")
message(STATUS "  USE_MKL: ${USE_MKL}")
message(STATUS "  USE_INT64: ${USE_INT64}")
message(STATUS "  USE_DIRECT_CBLAS: ${USE_DIRECT_CBLAS}")
message(STATUS "  BLAS Libraries: ${BLAS_LIBRARIES}")
message(STATUS "  LAPACK Libraries: ${LAPACK_LIBRARIES}")
message(STATUS "  CBLAS Libraries: ${CBLAS_LIBRARIES}")
message(STATUS "  GSL Libraries: ${GSL_LIBRARIES}")

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/libmess)
if(GSL_INCLUDE_DIRS)
  include_directories(${GSL_INCLUDE_DIRS})
endif()

# Find MPI early so we can add sources to MESS_SOURCES
find_package(MPI QUIET)
message(STATUS "MPI_FOUND: ${MPI_FOUND}")
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Source files
set(MESS_SOURCES
  src/libmess/atom.cc
  src/libmess/io.cc
  src/libmess/limits.cc
  src/libmess/math.cc
  src/libmess/symmetry.cc
  src/libmess/d3.cc
  src/libmess/key.cc
  src/libmess/mess.cc
  src/libmess/multindex.cc
  src/libmess/units.cc
  src/libmess/graph_common.cc
  src/libmess/lapack.cc
  src/libmess/mpack.cc
  src/libmess/permutation.cc
  src/libmess/graph_omp.cc
  src/libmess/linpack.cc
  src/libmess/model.cc
  src/libmess/slatec.cc
  src/libmess/crossrate.cc
  src/libmess/random.cc
  src/libmess/read.cc
  src/libmess/divsur.cc
  src/libmess/dynamic.cc
  src/libmess/structure.cc
  src/libmess/configuration.cc
  src/libmess/dynlib.cc
  src/libmess/monom.cc
  src/libmess/logical.cc
  src/libmess/potential.cc
  src/libmess/system.cc
  src/libmess/trajectory.cc
)

# Add MPI communication sources if MPI is available
if(MPI_FOUND)
  message(STATUS "MPI found, adding communication sources")
  list(APPEND MESS_SOURCES
    src/libmess/comm.cc
    src/libmess/auto_comm.cc
    src/libmess/new_comm.cc
  )
  list(LENGTH MESS_SOURCES num_sources)
  message(STATUS "Total MESS_SOURCES after MPI: ${num_sources}")
  message(STATUS "Added MPI communication sources to MESS_SOURCES")
endif()

# Create static library
add_library(messlibs STATIC ${MESS_SOURCES})

# Link libraries to messlibs
target_link_libraries(messlibs
    ${BLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES} 
    ${GSL_LIBRARIES}
  gfortran
  m
  dl
)

# Add GSL dependency
if(GSL_DOWNLOADED)
  add_dependencies(messlibs gsl_external)
  endif()

# Add optional libraries
if(QD_LIBRARIES)
  target_link_libraries(messlibs ${QD_LIBRARIES})
  endif()
if(SLATEC_LIBRARIES)
  target_link_libraries(messlibs ${SLATEC_LIBRARIES})
endif()

# Source files
set(MESS_SOURCES
  src/libmess/atom.cc
  src/libmess/io.cc
  src/libmess/limits.cc
  src/libmess/math.cc
  src/libmess/symmetry.cc
  src/libmess/d3.cc
  src/libmess/key.cc
  src/libmess/mess.cc
  src/libmess/multindex.cc
  src/libmess/units.cc
  src/libmess/graph_common.cc
  src/libmess/lapack.cc
  src/libmess/mpack.cc
  src/libmess/permutation.cc
  src/libmess/graph_omp.cc
  src/libmess/linpack.cc
  src/libmess/model.cc
  src/libmess/slatec.cc
  src/libmess/crossrate.cc
  src/libmess/random.cc
  src/libmess/read.cc
  src/libmess/divsur.cc
  src/libmess/dynamic.cc
  src/libmess/structure.cc
  src/libmess/configuration.cc
  src/libmess/dynlib.cc
  src/libmess/monom.cc
  src/libmess/logical.cc
  src/libmess/potential.cc
  src/libmess/system.cc
  src/libmess/trajectory.cc
)

# Create executables
add_executable(mess src/mess_driver.cc)
add_executable(mess-v2 src/mess_test.cc)
add_executable(messpf src/partition_function.cc)
if(MPI_FOUND)
  add_executable(messabs src/mess_mpi.cc)
  target_link_libraries(messabs MPI::MPI_CXX)
  else()
  message(STATUS "MPI not found, skipping messabs executable")
    endif()
add_executable(messsym src/gumbo.cc)

# Link executables
target_link_libraries(mess messlibs ${CBLAS_LIBRARIES})
target_link_libraries(mess-v2 messlibs ${CBLAS_LIBRARIES})
target_link_libraries(messpf messlibs ${CBLAS_LIBRARIES})
if(MPI_FOUND)
  target_link_libraries(messabs messlibs ${CBLAS_LIBRARIES})
endif()
target_link_libraries(messsym messlibs ${CBLAS_LIBRARIES})

# Test program
add_executable(test_integer_interface test_mkl_ilp64.cpp)
target_link_libraries(test_integer_interface messlibs ${CBLAS_LIBRARIES})