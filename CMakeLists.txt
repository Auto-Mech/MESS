project(MESS C CXX)
cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Options for automatic dependency handling
option(AUTO_DOWNLOAD_DEPENDENCIES "Automatically download and build missing dependencies" ON)
option(USE_SYSTEM_LIBS "Try to use system libraries before downloading" ON)

# Directory for downloaded and built dependencies
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)
include(ExternalProject)

# Check for direct SLATEC build
if(USE_DIRECT_SLATEC)
  message(STATUS "Using direct SLATEC build")
  if(EXISTS "${CMAKE_SOURCE_DIR}/external/lib/libslatec.a")
    set(SLATEC "${CMAKE_SOURCE_DIR}/external/lib/libslatec.a")
    message(STATUS "Found direct SLATEC: ${SLATEC}")
  else()
    message(WARNING "Direct SLATEC build enabled but library not found. Please run ./slatec_direct.sh first.")
    message(STATUS "Falling back to regular SLATEC build")
    set(USE_DIRECT_SLATEC OFF)
  endif()
endif()

# Function to conditionally download and build a dependency
function(handle_dependency name)
  if(NOT ${name}_FOUND AND AUTO_DOWNLOAD_DEPENDENCIES)
    message(STATUS "Will download and build ${name}")
    set(${name}_FOUND TRUE PARENT_SCOPE)
    set(${name}_DOWNLOADED TRUE PARENT_SCOPE)
  endif()
endfunction()

# Find BLAS and LAPACK
if(USE_SYSTEM_LIBS)
  find_package(BLAS QUIET)
  find_package(LAPACK QUIET)
else()
  set(BLAS_FOUND FALSE)
  set(LAPACK_FOUND FALSE)
endif()

# If BLAS/LAPACK not found, use OpenBLAS which provides both
if((NOT BLAS_FOUND OR NOT LAPACK_FOUND) AND AUTO_DOWNLOAD_DEPENDENCIES)
  message(STATUS "BLAS or LAPACK not found. Will download and build OpenBLAS.")
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
  
  ExternalProject_Add(
    openblas
    GIT_REPOSITORY https://github.com/xianyi/OpenBLAS.git
    GIT_TAG v0.3.21
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DCMAKE_C_FLAGS="-Wno-misleading-indentation"
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND $(MAKE) PREFIX=${EXTERNAL_INSTALL_LOCATION} install
  )
  
  set(BLAS_LIBRARIES ${EXTERNAL_INSTALL_LOCATION}/lib/libopenblas.so)
  set(LAPACK_LIBRARIES ${EXTERNAL_INSTALL_LOCATION}/lib/libopenblas.so)
  set(OPENBLAS_DOWNLOADED TRUE)
else()
  message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
  message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
endif()

# Find Fortran compilation capabilities
# This tells CMake to enable Fortran language for our project
enable_language(Fortran OPTIONAL)

# Find SLATEC
if(USE_SYSTEM_LIBS)
  find_library(SLATEC NAMES slatec libslatec QUIET)
endif()

if(NOT SLATEC AND AUTO_DOWNLOAD_DEPENDENCIES)
  message(STATUS "SLATEC not found. Will download and build it.")

  # Create SLATEC wrapper file content for missing C interfaces
  file(WRITE ${CMAKE_BINARY_DIR}/slatec_wrapper.f "
! Explicit wrapper for SLATEC functions that need C-compatible interfaces
      subroutine davint(x, y, n, xlo, xup, ans, ierr)
      implicit none
      integer n, ierr
      double precision x(n), y(n), xlo, xup, ans
      call davint_(x, y, n, xlo, xup, ans, ierr)
      return
      end

      subroutine ddeabm(f, neq, t, y, tout, info, rtol, atol, idid,
     &                 rwork, lrw, iwork, liw, rpar, ipar)
      external f
      integer neq, info(*), idid, lrw, iwork(*), liw
      double precision t, y(*), tout, rtol(*), atol(*), rwork(*)
      integer ipar(*)
      double precision rpar(*)
      call ddeabm_(f, neq, t, y, tout, info, rtol, atol, idid,
     &            rwork, lrw, iwork, liw, rpar, ipar)
      return
      end

      subroutine dbint4(x, y, n, ibcl, ibcr, fbcl, fbcr, kntopt, t, nc, 
     &                 bc, ws, lws, ierr)
      integer n, ibcl, ibcr, kntopt, nc, lws, ierr
      double precision x(*), y(*), fbcl, fbcr, t(*), bc(*), ws(*)
      call dbint4_(x, y, n, ibcl, ibcr, fbcl, fbcr, kntopt, t, nc, 
     &            bc, ws, lws, ierr)
      return
      end

      function dbvalu(t, a, n, k, ideriv, x, inbv, work)
      integer n, k, ideriv, inbv
      double precision t(*), a(*), x, work(*), dbvalu
      dbvalu = dbvalu_(t, a, n, k, ideriv, x, inbv, work)
      return
      end
")

  # Create build script for SLATEC
  file(WRITE ${CMAKE_BINARY_DIR}/build_slatec.sh "#!/bin/sh
set -e
# Use gfortran as default Fortran compiler if not specified
FC=${CMAKE_Fortran_COMPILER}
if [ -z \"$FC\" ]; then
  FC=gfortran
fi
FFLAGS='-O2 -fPIC -fdefault-real-8 -fdefault-double-8 -std=legacy'
AR=ar

echo 'Compiling SLATEC modules...'
for f in *.f; do
  echo \"Compiling $f...\"
  $FC $FFLAGS -c $f || exit 1
done

echo 'Creating library...'
$AR rcs libslatec.a *.o
echo 'Build complete'
")

  # Set up SLATEC build with more detailed configuration for reproducibility
  set(SLATEC_SOURCE_DIR ${CMAKE_BINARY_DIR}/slatec_src)
  set(SLATEC_BUILD_DIR ${CMAKE_BINARY_DIR}/slatec_build)
  
  include(ExternalProject)
  ExternalProject_Add(
    slatec_external
    URL https://www.netlib.org/slatec/slatec_src.tgz
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/downloads
    SOURCE_DIR ${SLATEC_SOURCE_DIR}
    BINARY_DIR ${SLATEC_BUILD_DIR}
    # Configure step - create a simpler setup
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory ${SLATEC_SOURCE_DIR} ${SLATEC_BUILD_DIR}
                      # Copy wrapper file
                      COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                              ${CMAKE_BINARY_DIR}/slatec_wrapper.f ${SLATEC_BUILD_DIR}/
                      # Copy build script
                      COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                              ${CMAKE_BINARY_DIR}/build_slatec.sh ${SLATEC_BUILD_DIR}/
                      # Make build script executable
                      COMMAND chmod +x ${SLATEC_BUILD_DIR}/build_slatec.sh
    # Build step - run the build script
    BUILD_COMMAND ${SLATEC_BUILD_DIR}/build_slatec.sh
    # Install step - copy built library
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                    ${SLATEC_BUILD_DIR}/libslatec.a ${EXTERNAL_INSTALL_LOCATION}/lib/libslatec.a
    # Disable updates
    UPDATE_COMMAND ""
    # Logging
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    LOG_INSTALL ON
  )
  
  set(SLATEC ${EXTERNAL_INSTALL_LOCATION}/lib/libslatec.a)
  set(SLATEC_DOWNLOADED TRUE)
else()
  message(STATUS "Found SLATEC: ${SLATEC}")
endif()

# Find Fortran libraries needed for SLATEC
# Try to find gfortran in multiple locations and with different possible names
find_library(GFORTRAN_LIBRARY NAMES 
  gfortran libgfortran libgfortran.so.5 libgfortran.so.4 libgfortran.so.3
  PATHS
  /usr/lib
  /usr/lib64
  /usr/lib/gcc
  /usr/local/lib
  /usr/local/lib64
  /usr/lib/x86_64-linux-gnu
  /lib/x86_64-linux-gnu
  # Add location hints from gfortran compiler
  PATH_SUFFIXES gcc/x86_64-linux-gnu/10 gcc/x86_64-linux-gnu/9 gcc/x86_64-linux-gnu/8 gcc/x86_64-linux-gnu/7
)

if(GFORTRAN_LIBRARY)
  message(STATUS "Found gfortran: ${GFORTRAN_LIBRARY}")
else()
  # Fall back to trying to get library location from gfortran itself
  find_program(GFORTRAN_EXECUTABLE gfortran)
  if(GFORTRAN_EXECUTABLE)
    message(STATUS "Found gfortran executable: ${GFORTRAN_EXECUTABLE}")
    execute_process(
      COMMAND ${GFORTRAN_EXECUTABLE} -print-file-name=libgfortran.so
      OUTPUT_VARIABLE GFORTRAN_LIB_PATH
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(GFORTRAN_LIB_PATH)
      get_filename_component(GFORTRAN_LIB_DIR ${GFORTRAN_LIB_PATH} DIRECTORY)
      message(STATUS "gfortran lib directory: ${GFORTRAN_LIB_DIR}")
      # Try to find the library in this directory
      find_library(GFORTRAN_LIBRARY NAMES gfortran libgfortran
        PATHS ${GFORTRAN_LIB_DIR} NO_DEFAULT_PATH)
    endif()
  endif()

  if(GFORTRAN_LIBRARY)
    message(STATUS "Found gfortran via fallback: ${GFORTRAN_LIBRARY}")
  else()
    # If we still can't find it, just hardcode as direct link option
    message(WARNING "gfortran library not found, will try direct linking")
    set(GFORTRAN_LIBRARY "-lgfortran")
  endif()
endif()

# Find optional quadmath library (often needed with gfortran)
find_library(QUADMATH_LIBRARY NAMES quadmath)
if(QUADMATH_LIBRARY)
  message(STATUS "Found quadmath: ${QUADMATH_LIBRARY}")
endif()

# Set up link libraries for SLATEC
set(FORTRAN_LIBRARIES ${GFORTRAN_LIBRARY})
if(QUADMATH_LIBRARY)
  list(APPEND FORTRAN_LIBRARIES ${QUADMATH_LIBRARY})
endif()

# Find GSL
if(USE_SYSTEM_LIBS)
  find_package(GSL QUIET)
endif()

if(NOT GSL_FOUND AND AUTO_DOWNLOAD_DEPENDENCIES)
  message(STATUS "GSL not found. Will download and build it.")
  ExternalProject_Add(
    gsl_external
    URL https://ftp.gnu.org/gnu/gsl/gsl-2.7.tar.gz
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${EXTERNAL_INSTALL_LOCATION}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  set(GSL_INCLUDE_DIRS ${EXTERNAL_INSTALL_LOCATION}/include)
  set(GSL_LIBRARIES ${EXTERNAL_INSTALL_LOCATION}/lib/libgsl.so)
  include_directories(${GSL_INCLUDE_DIRS})
  set(GSL_DOWNLOADED TRUE)
  set(GSL_FOUND TRUE)
else()
  message(STATUS "Found GSL: ${GSL_LIBRARIES}")
endif()

# Find QD
if(USE_SYSTEM_LIBS)
  find_library(QD NAMES qd libqd libqd.a QUIET)
endif()

if(NOT QD AND AUTO_DOWNLOAD_DEPENDENCIES)
  message(STATUS "QD not found. Will download and build it.")
  ExternalProject_Add(
    qd_external
    URL https://www.davidhbailey.com/dhbsoftware/qd-2.3.23.tar.gz
    URL_HASH SHA256=b3eaf41ce413ec08f348ee73e606bd3ff9203e411c377c3c0467f89acf69ee26
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${EXTERNAL_INSTALL_LOCATION}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  set(QD ${EXTERNAL_INSTALL_LOCATION}/lib/libqd.a)
  set(QD_DOWNLOADED TRUE)
else()
  message(STATUS "Found QD: ${QD}")
endif()

# Create a list of dependencies we need to add
set(DEPENDENCIES_TARGETS "")
if(DEFINED OPENBLAS_DOWNLOADED)
  list(APPEND DEPENDENCIES_TARGETS openblas)
endif()
if(DEFINED SLATEC_DOWNLOADED)
  list(APPEND DEPENDENCIES_TARGETS slatec_external)
endif()
if(DEFINED GSL_DOWNLOADED)
  list(APPEND DEPENDENCIES_TARGETS gsl_external)
endif()
if(DEFINED QD_DOWNLOADED)
  list(APPEND DEPENDENCIES_TARGETS qd_external)
endif()

# Set up include directories for downloaded libraries
include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

# Define the library
add_library(messlibs
    ${PROJECT_SOURCE_DIR}/src/libmess/atom.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/io.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/math.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/symmetry.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/d3.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/key.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/mess.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/multindex.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/units.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/graph_common.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/lapack.cc
    # ${PROJECT_SOURCE_DIR}/src/libmess/mpack.cc
    # ${PROJECT_SOURCE_DIR}/src/libmess/mpack_dd.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/permutation.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/graph_omp.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/linpack.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/model.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/slatec.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/crossrate.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/random.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/read.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/divsur.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/dynamic.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/structure.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/configuration.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/dynlib.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/monom.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/logical.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/potential.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/system.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/trajectory.cc)

# Add dependency on external libraries if needed
if(DEPENDENCIES_TARGETS)
  add_dependencies(messlibs ${DEPENDENCIES_TARGETS})
endif()

# Set Fortran-compatible flags for C/C++ code
set_target_properties(messlibs PROPERTIES
  COMPILE_FLAGS "-fPIC -fno-strict-aliasing"
)

# Set up executables
add_executable(mess ${PROJECT_SOURCE_DIR}/src/mess_driver.cc)
add_executable(mess-v2 ${PROJECT_SOURCE_DIR}/src/mess_test.cc)
add_executable(messpf ${PROJECT_SOURCE_DIR}/src/partition_function.cc)
add_executable(messabs ${PROJECT_SOURCE_DIR}/src/abstraction.cc)
add_executable(messsym ${PROJECT_SOURCE_DIR}/src/symmetry_number.cc)

# Link libraries
target_link_libraries(mess
    messlibs 
    ${BLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES} 
    ${GSL_LIBRARIES}
    ${QD} 
    -Wl,--start-group  # Use link groups to handle circular dependencies
    ${SLATEC} 
    ${FORTRAN_LIBRARIES} 
    -lm  # Math library
    -Wl,--end-group
    ${CMAKE_DL_LIBS}
)
target_link_libraries(mess-v2
    messlibs 
    ${BLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES} 
    ${GSL_LIBRARIES}
    ${QD} 
    -Wl,--start-group
    ${SLATEC} 
    ${FORTRAN_LIBRARIES} 
    -lm
    -Wl,--end-group
    ${CMAKE_DL_LIBS}
)
target_link_libraries(messpf
    messlibs 
    ${BLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES} 
    ${GSL_LIBRARIES}
    ${QD} 
    -Wl,--start-group
    ${SLATEC} 
    ${FORTRAN_LIBRARIES} 
    -lm
    -Wl,--end-group
    ${CMAKE_DL_LIBS}
)
target_link_libraries(messabs
    messlibs 
    ${BLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES} 
    ${GSL_LIBRARIES}
    ${QD} 
    -Wl,--start-group
    ${SLATEC} 
    ${FORTRAN_LIBRARIES} 
    -lm
    -Wl,--end-group
    ${CMAKE_DL_LIBS}
)
target_link_libraries(messsym
    messlibs 
    ${BLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES} 
    ${GSL_LIBRARIES}
    ${QD} 
    -Wl,--start-group
    ${SLATEC} 
    ${FORTRAN_LIBRARIES} 
    -lm
    -Wl,--end-group
    ${CMAKE_DL_LIBS}
)

install(TARGETS mess DESTINATION bin)
install(TARGETS mess-v2 DESTINATION bin)
install(TARGETS messpf DESTINATION bin)
install(TARGETS messabs DESTINATION bin)
install(TARGETS messsym DESTINATION bin)

# Create a direct wrapper for missing SLATEC functions
file(WRITE ${CMAKE_BINARY_DIR}/slatec_wrapper.f "
      subroutine davint(x, y, n, xlo, xup, ans, ierr)
      implicit none
      integer n, ierr
      double precision x(n), y(n), xlo, xup, ans
      call davint_(x, y, n, xlo, xup, ans, ierr)
      return
      end

      subroutine ddeabm(f, neq, t, y, tout, info, rtol, atol, idid,
     &                 rwork, lrw, iwork, liw, rpar, ipar)
      external f
      integer neq, info(*), idid, lrw, iwork(*), liw
      double precision t, y(*), tout, rtol(*), atol(*), rwork(*)
      integer ipar(*)
      double precision rpar(*)
      call ddeabm_(f, neq, t, y, tout, info, rtol, atol, idid,
     &            rwork, lrw, iwork, liw, rpar, ipar)
      return
      end

      subroutine dbint4(x, y, n, ibcl, ibcr, fbcl, fbcr, kntopt, t, nc, 
     &                 bc, ws, lws, ierr)
      integer n, ibcl, ibcr, kntopt, nc, lws, ierr
      double precision x(*), y(*), fbcl, fbcr, t(*), bc(*), ws(*)
      call dbint4_(x, y, n, ibcl, ibcr, fbcl, fbcr, kntopt, t, nc, 
     &            bc, ws, lws, ierr)
      return
      end

      function dbvalu(t, a, n, k, ideriv, x, inbv, work)
      integer n, k, ideriv, inbv
      double precision t(*), a(*), x, work(*), dbvalu
      dbvalu = dbvalu_(t, a, n, k, ideriv, x, inbv, work)
      return
      end
")

# Compile the SLATEC wrapper
add_custom_command(
  OUTPUT ${CMAKE_BINARY_DIR}/slatec_wrapper.o
  COMMAND gfortran -c -O2 -fPIC ${CMAKE_BINARY_DIR}/slatec_wrapper.f -o ${CMAKE_BINARY_DIR}/slatec_wrapper.o
  DEPENDS ${CMAKE_BINARY_DIR}/slatec_wrapper.f
  COMMENT "Compiling SLATEC wrapper"
)

# Create a custom target for the wrapper
add_custom_target(slatec_wrapper DEPENDS ${CMAKE_BINARY_DIR}/slatec_wrapper.o)
add_dependencies(messlibs slatec_wrapper)
