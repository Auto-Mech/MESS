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
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND $(MAKE) install
  )
  
  set(BLAS_LIBRARIES ${EXTERNAL_INSTALL_LOCATION}/lib/libopenblas.so)
  set(LAPACK_LIBRARIES ${EXTERNAL_INSTALL_LOCATION}/lib/libopenblas.so)
  set(OPENBLAS_DOWNLOADED TRUE)
else()
  message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
  message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
endif()

# Find SLATEC
if(USE_SYSTEM_LIBS)
  find_library(SLATEC NAMES slatec libslatec QUIET)
endif()

if(NOT SLATEC AND AUTO_DOWNLOAD_DEPENDENCIES)
  message(STATUS "SLATEC not found. Will download and build it.")
  
  # Create a custom download method for SLATEC
  # Use a known good direct source instead of GitHub or Netlib
  ExternalProject_Add(
    slatec_external
    URL https://www.netlib.org/slatec/slatec_src.tgz
    # Don't use hash verification since the file seems to change
    DOWNLOAD_NO_PROGRESS TRUE
    # Configure step - create Makefile
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR> <BINARY_DIR>
                      && ${CMAKE_COMMAND} -E echo "FC = gfortran" > <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "FFLAGS = -O2 -fPIC" >> <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "AR = ar" >> <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "SRCS = $(wildcard *.f)" >> <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "OBJS = $(SRCS:.f=.o)" >> <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "all: libslatec.a" >> <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "libslatec.a: $(OBJS)" >> <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "\t$(AR) rcs $@ $^" >> <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "%.o: %.f" >> <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "\t$(FC) $(FFLAGS) -c $< -o $@" >> <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "clean:" >> <BINARY_DIR>/Makefile
                      && ${CMAKE_COMMAND} -E echo "\trm -f $(OBJS) libslatec.a" >> <BINARY_DIR>/Makefile
    # Build step
    BUILD_COMMAND make -j4
    # Install step
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                    <BINARY_DIR>/libslatec.a ${EXTERNAL_INSTALL_LOCATION}/lib/libslatec.a
    # Disable hash verification
    TLS_VERIFY FALSE
    # Disable update step to avoid re-downloading if it succeeded once
    UPDATE_COMMAND ""
  )
  set(SLATEC ${EXTERNAL_INSTALL_LOCATION}/lib/libslatec.a)
  set(SLATEC_DOWNLOADED TRUE)
else()
  message(STATUS "Found SLATEC: ${SLATEC}")
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
    URL_HASH SHA256=efbbf3785da0e53038be7907500628b466152dbc3c173a87de1b5eba2e19a986
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
    ${PROJECT_SOURCE_DIR}/src/libmess/monom.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/logical.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/potential.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/system.cc
    ${PROJECT_SOURCE_DIR}/src/libmess/trajectory.cc)

# Add dependency on external libraries if needed
if(DEPENDENCIES_TARGETS)
  add_dependencies(messlibs ${DEPENDENCIES_TARGETS})
endif()

# Set up executables
add_executable(mess ${PROJECT_SOURCE_DIR}/src/mess_driver.cc)
add_executable(mess-v2 ${PROJECT_SOURCE_DIR}/src/mess_test.cc)
add_executable(messpf ${PROJECT_SOURCE_DIR}/src/partition_function.cc)
add_executable(messabs ${PROJECT_SOURCE_DIR}/src/abstraction.cc)
add_executable(messsym ${PROJECT_SOURCE_DIR}/src/symmetry_number.cc)

# Link libraries
target_link_libraries(mess
    messlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GSL_LIBRARIES}
    # ${MLAPACK_QD} ${MLAPACK_DD} ${MBLAS_QD} ${MBLAS_DD}
    ${QD} ${SLATEC} ${CMAKE_DL_LIBS}
)
target_link_libraries(mess-v2
    messlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GSL_LIBRARIES}
    # ${MLAPACK_QD} ${MLAPACK_DD} ${MBLAS_QD} ${MBLAS_DD}
    ${QD} ${SLATEC} ${CMAKE_DL_LIBS}
)
target_link_libraries(messpf
    messlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GSL_LIBRARIES}
    # ${MLAPACK_QD} ${MLAPACK_DD} ${MBLAS_QD} ${MBLAS_DD}
    ${QD} ${SLATEC} ${CMAKE_DL_LIBS}
)
target_link_libraries(messabs
    messlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GSL_LIBRARIES}
    # ${MLAPACK_QD} ${MLAPACK_DD} ${MBLAS_QD} ${MBLAS_DD}
    ${QD} ${SLATEC} ${CMAKE_DL_LIBS}
)
target_link_libraries(messsym
    messlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GSL_LIBRARIES}
    # ${MLAPACK_QD} ${MLAPACK_DD} ${MBLAS_QD} ${MBLAS_DD}
    ${QD} ${SLATEC} ${CMAKE_DL_LIBS}
)

install(TARGETS mess DESTINATION bin)
install(TARGETS mess-v2 DESTINATION bin)
install(TARGETS messpf DESTINATION bin)
install(TARGETS messabs DESTINATION bin)
install(TARGETS messsym DESTINATION bin)
