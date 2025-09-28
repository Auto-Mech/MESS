# =============================================================================
# Project Definition
# =============================================================================
cmake_minimum_required(VERSION 3.16)
project(MESS LANGUAGES C CXX Fortran)

# =============================================================================
# Modern CMake Practices & Project-Wide Settings
# =============================================================================
# Set C++ standard for all targets in the project
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Use a dedicated folder for external projects to keep the build directory clean
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/_deps)

# =============================================================================
# Project Options
# =============================================================================
option(USE_MKL "Use Intel MKL instead of OpenBLAS" OFF)
option(USE_INT64 "Use 64-bit integers for BLAS/LAPACK (ILP64 interface)" OFF)
# option(USE_DIRECT_CBLAS "Use direct CBLAS instead of GSL's wrapper" OFF) # This option was unused.

# =============================================================================
# Find Dependencies (BLAS, LAPACK, GSL, MPI, etc.)
# =============================================================================
include(FetchContent)

# --- BLAS / LAPACK Configuration ---
message(STATUS "Configuring BLAS and LAPACK...")

if(USE_MKL)
  message(STATUS "Using Intel MKL")
  if(USE_INT64)
    set(BLAS_LIBRARIES 
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_ilp64.so
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_sequential.so
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so
      m
      dl
    )
  set(LAPACK_LIBRARIES ${BLAS_LIBRARIES})
    message(STATUS "  - Using 64-bit integers (ILP64)")
  else()
    set(BLAS_LIBRARIES 
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_sequential.so
      /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so
      m
      dl
    )
  set(LAPACK_LIBRARIES ${BLAS_LIBRARIES})
    message(STATUS "  - Using 32-bit integers (LP64)")
  endif()
else()
  message(STATUS "Using OpenBLAS")
  set(BLAS_LIBRARIES /usr/lib/x86_64-linux-gnu/libopenblas.so)
  set(LAPACK_LIBRARIES /usr/lib/x86_64-linux-gnu/libopenblas.so)
  if(USE_INT64)
    message(STATUS "  - Using 64-bit integers (ILP64)")
  else()
    message(STATUS "  - Using 32-bit integers (LP64)")
  endif()
endif()

message(STATUS "BLAS Libraries: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK Libraries: ${LAPACK_LIBRARIES}")


# --- GSL (GNU Scientific Library) ---
# First, try to find an installed version. If not found, download and build it.
find_package(GSL QUIET)
if(NOT GSL_FOUND)
  message(STATUS "GSL not found. Building from source...")
  FetchContent_Declare(
    gsl_external
    URL "https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz"
    URL_HASH SHA256=6a99eeed15632c6354895b1dd542ed5a855c0f15d9ad1326c6fe2b2c9e423190
  )
  FetchContent_MakeAvailable(gsl_external)

  # Configure GSL to use the BLAS library we found earlier
  set(GSL_CONFIGURE_OPTIONS
    "--prefix=${EXTERNAL_INSTALL_LOCATION}"
    "--disable-shared"
    "--enable-static"
    "LDFLAGS=-L${BLAS_LIBRARY_DIRS}"
    "LIBS=${BLAS_LIBRARIES}"
  )

  ExternalProject_Add(
    gsl_project
    SOURCE_DIR ${gsl_external_SOURCE_DIR}
    BINARY_DIR ${gsl_external_BINARY_DIR}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure ${GSL_CONFIGURE_OPTIONS}
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND $(MAKE) install
    LOG_CONFIGURE TRUE
    LOG_BUILD TRUE
    LOG_INSTALL TRUE
  )

  set(GSL_INCLUDE_DIRS ${EXTERNAL_INSTALL_LOCATION}/include)
  set(GSL_LIBRARIES ${EXTERNAL_INSTALL_LOCATION}/lib/libgsl.a)
  set(GSL_DOWNLOADED TRUE)
else()
  message(STATUS "Found GSL: ${GSL_LIBRARIES}")
endif()


# --- SLATEC (Common Mathematical Library) ---
# First, try to find an installed version. If not found, download and build it.
find_library(SLATEC_LIBRARY NAMES slatec libslatec)
if(NOT SLATEC_LIBRARY)
    message(STATUS "SLATEC not found. Building from source...")
    FetchContent_Declare(
        slatec_external
        URL "http://www.netlib.org/slatec/slatec_src.tgz"
        URL_HASH SHA1=6100817113c115778eddec5811065be818157554 # SHA1 from CMake docs for this file
    )
    FetchContent_Populate(slatec_external)

    # Glob all Fortran sources from the downloaded content
    file(GLOB SLATEC_SOURCES "${slatec_external_SOURCE_DIR}/*.f")

    # Create a static library from the SLATEC sources
    add_library(slatec_lib STATIC ${SLATEC_SOURCES})
    
    # Set Fortran compiler flags for SLATEC (Fortran 77 compatibility)
    target_compile_options(slatec_lib PRIVATE 
        -std=legacy
        -fPIC
        -O2
        -w  # Suppress warnings for old Fortran code
    )
    
    # Link Fortran runtime
    target_link_libraries(slatec_lib PRIVATE gfortran)
    
    set(SLATEC_LIBRARIES slatec_lib)
    set(SLATEC_DOWNLOADED TRUE)
else()
    message(STATUS "Found SLATEC: ${SLATEC_LIBRARY}")
    set(SLATEC_LIBRARIES ${SLATEC_LIBRARY})
endif()


# --- Other Optional Libraries ---
find_library(QD_LIBRARIES NAMES qd libqd)
find_package(MPI QUIET)

# =============================================================================
# Library Target (messlibs)
# =============================================================================
set(MESS_SOURCES
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

add_library(messlibs STATIC ${MESS_SOURCES})

# Add dependencies if they were built from source
if(GSL_DOWNLOADED)
  add_dependencies(messlibs gsl_project)
endif()

# --- Target Properties: Include Directories ---
target_include_directories(messlibs PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/libmess>
  $<INSTALL_INTERFACE:include>
  ${GSL_INCLUDE_DIRS}
)

# --- Target Properties: Compile Definitions ---
if(USE_INT64)
  target_compile_definitions(messlibs PUBLIC -DUSE_INT64)
else()
  target_compile_definitions(messlibs PUBLIC -DUSE_32BIT_INTEGERS)
endif()

# --- Target Properties: Link Libraries ---
target_link_libraries(messlibs PUBLIC
  ${LAPACK_LIBRARIES} # LAPACK depends on BLAS, so this is often sufficient
    ${GSL_LIBRARIES}
  ${SLATEC_LIBRARIES}
  m dl # Common system libraries
)

# --- Conditionally link optional libraries ---
if(QD_LIBRARIES)
    message(STATUS "Found and linking QD library: ${QD_LIBRARIES}")
    target_link_libraries(messlibs PUBLIC ${QD_LIBRARIES})
else()
    message(STATUS "Optional QD library not found, skipping.")
endif()

# --- Handle MPI ---
if(MPI_FOUND)
  message(STATUS "MPI found, adding communication sources and linking.")
  target_sources(messlibs PRIVATE
    src/libmess/comm.cc
    src/libmess/auto_comm.cc
    src/libmess/new_comm.cc
  )
  target_include_directories(messlibs PUBLIC ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(messlibs PUBLIC MPI::MPI_CXX)
endif()


# =============================================================================
# Executable Targets
# =============================================================================
add_executable(mess src/mess_driver.cc)
add_executable(mess-v2 src/mess_test.cc)
add_executable(messpf src/partition_function.cc)
add_executable(messsym src/gumbo.cc)
add_executable(test_integer_interface test_mkl_ilp64.cpp)

# Link all standard executables to the main library
target_link_libraries(mess messlibs)
target_link_libraries(mess-v2 messlibs)
target_link_libraries(messpf messlibs)
target_link_libraries(messsym messlibs)
target_link_libraries(test_integer_interface messlibs)

# Handle MPI executable separately
if(MPI_FOUND)
  add_executable(messabs src/mess_mpi.cc)
  target_link_libraries(messabs messlibs)
endif()