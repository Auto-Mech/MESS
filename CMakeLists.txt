project(MESS C CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(SLATEC REQUIRED NAMES slatec libslatec)
message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
message(STATUS "Found SLATEC: ${SLATEC_LIBRARIES}")
find_library(MPACK NAMES mpack libmpack)
if(MPACK)
    message(STATUS "Compiling with MPACK")
    add_definitions(-DWITH_MPACK)
else()
    message(STATUS "Compiling without MPACK. If you do want MPACK, set -DUSE_MPACK=ON")
endif()

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

add_executable(mess ${PROJECT_SOURCE_DIR}/src/mess_driver.cc)
add_executable(messpf ${PROJECT_SOURCE_DIR}/src/partition_function.cc)
add_executable(messabs ${PROJECT_SOURCE_DIR}/src/abstraction.cc)
add_executable(messsym ${PROJECT_SOURCE_DIR}/src/symmetry_number.cc)

if(USE_MPACK)
    add_library(mpack ${PROJECT_SOURCE_DIR}/src/libmess/mpack.cc)
    target_link_libraries(mess
        messlibs mpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} mlapack_qd
        mlapack_dd mblas_qd mblas_dd qd ${SLATEC} dl)
    target_link_libraries(messpf
        messlibs mpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} mlapack_qd
        mlapack_dd mblas_qd mblas_dd qd ${SLATEC} dl)
    target_link_libraries(messabs
        messlibs mpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} mlapack_qd
        mlapack_dd mblas_qd mblas_dd qd ${SLATEC} dl)
    target_link_libraries(messsym
        messlibs mpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} mlapack_qd
        mlapack_dd mblas_qd mblas_dd qd ${SLATEC} dl)
else()
    target_link_libraries(mess
        messlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} dl)
    target_link_libraries(messpf
        messlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} dl)
    target_link_libraries(messabs
        messlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} dl)
    target_link_libraries(messsym
        messlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${SLATEC} dl)
endif()

install(TARGETS mess DESTINATION bin)
install(TARGETS messpf DESTINATION bin)
install(TARGETS messabs DESTINATION bin)
install(TARGETS messsym DESTINATION bin)
