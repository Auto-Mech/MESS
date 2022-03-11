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
find_package(GSL REQUIRED)
# MPACK
find_library(QD REQUIRED NAMES qd libqd libqd.a)
# find_library(MBLAS_QD REQUIRED NAMES mblas_qd libmblas_qd)
# find_library(MBLAS_DD REQUIRED NAMES mblas_dd libmblas_dd)
# find_library(MLAPACK_QD REQUIRED NAMES mlapack_qd libmlapack_qd)
# find_library(MLAPACK_DD REQUIRED NAMES mlapack_dd libmlapack_dd)
message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
message(STATUS "Found SLATEC: ${SLATEC_LIBRARIES}")
message(STATUS "Found GSL: ${GSL_LIBRARIES}")
message(STATUS "Found QD: ${QD}")
# message(STATUS "Found MBLAS_QD: ${MBLAS_QD}")
# message(STATUS "Found MBLAS_DD: ${MBLAS_DD}")
# message(STATUS "Found MLAPACK_QD: ${MLAPACK_QD}")
# message(STATUS "Found MLAPACK_DD: ${MLAPACK_DD}")

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

add_executable(mess ${PROJECT_SOURCE_DIR}/src/mess_driver.cc)
add_executable(mess-v2 ${PROJECT_SOURCE_DIR}/src/mess_test.cc)
add_executable(messpf ${PROJECT_SOURCE_DIR}/src/partition_function.cc)
add_executable(messabs ${PROJECT_SOURCE_DIR}/src/abstraction.cc)
add_executable(messsym ${PROJECT_SOURCE_DIR}/src/symmetry_number.cc)

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
