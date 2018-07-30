project(MESS CXX Fortran)
cmake_minimum_required(VERSION 3.0)
find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

option(USE_MPACK "Use MPACK for multi-precision linear algebra" OFF)

if(USE_MPACK)
    message(STATUS "Compiling with MPACK")
    add_definitions(-DWITH_MPACK)
else()
    message(STATUS "Compiling without MPACK. If you do want MPACK, set -DUSE_MPACK=ON")
endif()


add_library(mess
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

add_library(slatec
    ${PROJECT_SOURCE_DIR}/external/slatec/d1mach.f
    ${PROJECT_SOURCE_DIR}/external/slatec/d9b0mp.f
    ${PROJECT_SOURCE_DIR}/external/slatec/d9b1mp.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dasum.f
    ${PROJECT_SOURCE_DIR}/external/slatec/davint.f
    ${PROJECT_SOURCE_DIR}/external/slatec/daxpy.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbesi0.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbesj0.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbesj1.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbesk0.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbint4.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbnfac.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbnslv.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbsi0e.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbsk0e.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbspvd.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbspvn.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dbvalu.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dcsevl.f
    ${PROJECT_SOURCE_DIR}/external/slatec/ddeabm.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dderkf.f
    ${PROJECT_SOURCE_DIR}/external/slatec/ddes.f
    ${PROJECT_SOURCE_DIR}/external/slatec/ddot.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dfehl.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dfzero.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dgaus8.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dgeco.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dgefa.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dgefs.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dgesl.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dhstrt.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dhvnrm.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dintp.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dintrv.f
    ${PROJECT_SOURCE_DIR}/external/slatec/drf.f
    ${PROJECT_SOURCE_DIR}/external/slatec/drkfs.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dscal.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dsifa.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dsisl.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dsteps.f
    ${PROJECT_SOURCE_DIR}/external/slatec/dswap.f
    ${PROJECT_SOURCE_DIR}/external/slatec/fdump.f
    ${PROJECT_SOURCE_DIR}/external/slatec/i1mach.f
    ${PROJECT_SOURCE_DIR}/external/slatec/idamax.f
    ${PROJECT_SOURCE_DIR}/external/slatec/initds.f
    ${PROJECT_SOURCE_DIR}/external/slatec/j4save.f
    ${PROJECT_SOURCE_DIR}/external/slatec/xercnt.f
    ${PROJECT_SOURCE_DIR}/external/slatec/xerhlt.f
    ${PROJECT_SOURCE_DIR}/external/slatec/xermsg.f
    ${PROJECT_SOURCE_DIR}/external/slatec/xerprn.f
    ${PROJECT_SOURCE_DIR}/external/slatec/xersve.f
    ${PROJECT_SOURCE_DIR}/external/slatec/xgetua.f)

add_executable(messdr ${PROJECT_SOURCE_DIR}/src/mess_driver.cc)
add_executable(messpf ${PROJECT_SOURCE_DIR}/src/partition_function.cc)
add_executable(messabs ${PROJECT_SOURCE_DIR}/src/abstraction.cc)
set_target_properties(messdr PROPERTIES OUTPUT_NAME "mess")

if(USE_MPACK)
    add_library(mpack ${PROJECT_SOURCE_DIR}/src/libmess/mpack.cc)
    target_link_libraries(messdr
        mess mpack blas lapack mlapack_qd mlapack_dd mblas_qd mblas_dd qd slatec dl)
    target_link_libraries(messpf
        mess mpack blas lapack mlapack_qd mlapack_dd mblas_qd mblas_dd qd slatec dl)
    target_link_libraries(messabs
        mess mpack blas lapack mlapack_qd mlapack_dd mblas_qd mblas_dd qd slatec dl)
else()
    target_link_libraries(messdr mess blas lapack slatec dl)
    target_link_libraries(messpf mess blas lapack slatec dl)
    target_link_libraries(messabs mess blas lapack slatec dl)
endif()

install(TARGETS messdr DESTINATION bin)
install(TARGETS messpf DESTINATION bin)
install(TARGETS messabs DESTINATION bin)
