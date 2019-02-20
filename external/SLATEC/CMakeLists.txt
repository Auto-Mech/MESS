project(MESS Fortran)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-implicit-none")

add_library(slatec SHARED
    ${PROJECT_SOURCE_DIR}/src/d1mach.f
    ${PROJECT_SOURCE_DIR}/src/d9b0mp.f
    ${PROJECT_SOURCE_DIR}/src/d9b1mp.f
    ${PROJECT_SOURCE_DIR}/src/dasum.f
    ${PROJECT_SOURCE_DIR}/src/davint.f
    ${PROJECT_SOURCE_DIR}/src/daxpy.f
    ${PROJECT_SOURCE_DIR}/src/dbesi0.f
    ${PROJECT_SOURCE_DIR}/src/dbesj0.f
    ${PROJECT_SOURCE_DIR}/src/dbesj1.f
    ${PROJECT_SOURCE_DIR}/src/dbesk0.f
    ${PROJECT_SOURCE_DIR}/src/dbint4.f
    ${PROJECT_SOURCE_DIR}/src/dbnfac.f
    ${PROJECT_SOURCE_DIR}/src/dbnslv.f
    ${PROJECT_SOURCE_DIR}/src/dbsi0e.f
    ${PROJECT_SOURCE_DIR}/src/dbsk0e.f
    ${PROJECT_SOURCE_DIR}/src/dbspvd.f
    ${PROJECT_SOURCE_DIR}/src/dbspvn.f
    ${PROJECT_SOURCE_DIR}/src/dbvalu.f
    ${PROJECT_SOURCE_DIR}/src/dcsevl.f
    ${PROJECT_SOURCE_DIR}/src/ddeabm.f
    ${PROJECT_SOURCE_DIR}/src/dderkf.f
    ${PROJECT_SOURCE_DIR}/src/ddes.f
    ${PROJECT_SOURCE_DIR}/src/ddot.f
    ${PROJECT_SOURCE_DIR}/src/dfehl.f
    ${PROJECT_SOURCE_DIR}/src/dfzero.f
    ${PROJECT_SOURCE_DIR}/src/dgaus8.f
    ${PROJECT_SOURCE_DIR}/src/dgeco.f
    ${PROJECT_SOURCE_DIR}/src/dgefa.f
    ${PROJECT_SOURCE_DIR}/src/dgefs.f
    ${PROJECT_SOURCE_DIR}/src/dgesl.f
    ${PROJECT_SOURCE_DIR}/src/dhstrt.f
    ${PROJECT_SOURCE_DIR}/src/dhvnrm.f
    ${PROJECT_SOURCE_DIR}/src/dintp.f
    ${PROJECT_SOURCE_DIR}/src/dintrv.f
    ${PROJECT_SOURCE_DIR}/src/drf.f
    ${PROJECT_SOURCE_DIR}/src/drkfs.f
    ${PROJECT_SOURCE_DIR}/src/dscal.f
    ${PROJECT_SOURCE_DIR}/src/dsifa.f
    ${PROJECT_SOURCE_DIR}/src/dsisl.f
    ${PROJECT_SOURCE_DIR}/src/dsteps.f
    ${PROJECT_SOURCE_DIR}/src/dswap.f
    ${PROJECT_SOURCE_DIR}/src/fdump.f
    ${PROJECT_SOURCE_DIR}/src/i1mach.f
    ${PROJECT_SOURCE_DIR}/src/idamax.f
    ${PROJECT_SOURCE_DIR}/src/initds.f
    ${PROJECT_SOURCE_DIR}/src/j4save.f
    ${PROJECT_SOURCE_DIR}/src/xercnt.f
    ${PROJECT_SOURCE_DIR}/src/xerhlt.f
    ${PROJECT_SOURCE_DIR}/src/xermsg.f
    ${PROJECT_SOURCE_DIR}/src/xerprn.f
    ${PROJECT_SOURCE_DIR}/src/xersve.f
    ${PROJECT_SOURCE_DIR}/src/xgetua.f)

install(TARGETS slatec LIBRARY DESTINATION lib)
