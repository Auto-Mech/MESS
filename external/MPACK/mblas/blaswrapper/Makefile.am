if ENABLE_WRAPPER

lib_LTLIBRARIES = libblaswrapper_gmp.la libblaswrapper_mpfr.la libblaswrapper_qd.la libblaswrapper_dd.la
if ENABLE_DOUBLE
lib_LTLIBRARIES += libblaswrapper_double.la
endif

BLASWARPPER_SOURCES = zaxpy.cpp ../reference/mpackinit.cpp
#zcopy.cpp zdotc.cpp zdotu.cpp zgbmv.cpp \
#zgemm.cpp zgemv.cpp zgerc.cpp zgeru.cpp zhbmv.cpp zhemm.cpp zhemv.cpp \
#zher.cpp zher2.cpp zher2k.cpp zherk.cpp zhpmv.cpp zhpr.cpp zhpr2.cpp \
#zrotg.cpp zscal.cpp zdrot.cpp zdscal.cpp zswap.cpp zsymm.cpp zsyr2k.cpp \
#zsyrk.cpp ztbmv.cpp ztbsv.cpp ztpmv.cpp ztpsv.cpp ztrmm.cpp ztrmv.cpp \
#ztrsm.cpp ztrsv.cpp dzabs1.cpp dzasum.cpp dznrm2.cpp dasum.cpp daxpy.cpp\
#dcopy.cpp ddot.cpp dgbmv.cpp dgemm.cpp dgemv.cpp dger.cpp dnrm2.cpp drot.cpp\
#drotg.cpp drotm.cpp drotmg.cpp dsbmv.cpp dscal.cpp dspmv.cpp \
#dspr.cpp dspr2.cpp dswap.cpp dsymm.cpp dsymv.cpp dsyr.cpp dsyr2.cpp \
#dsyr2k.cpp dsyrk.cpp dtbmv.cpp dtbsv.cpp dtpmv.cpp dtpsv.cpp dtrmm.cpp \
#dtrmv.cpp dtrsm.cpp dtrsv.cpp izamax.cpp idamax.cpp \
#dset.cpp lsame.cpp xerbla.cpp

libblaswrapper_gmp_la_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/$(GMP_INCLUDEDIR) -D___MPACK_BUILD_WITH_GMP___
libblaswrapper_gmp_la_SOURCES = $(BLASWARPPER_SOURCES)

libblaswrapper_mpfr_la_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(top_srcdir)/$(MPFR_INCLUDEDIR) -I$(top_srcdir)/$(MPC_INCLUDEDIR) -D___MPACK_BUILD_WITH_MPFR___
libblaswrapper_mpfr_la_SOURCES = $(BLASWARPPER_SOURCES)

libblaswrapper_qd_la_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/$(QD_INCLUDEDIR) -D___MPACK_BUILD_WITH_QD___
libblaswrapper_qd_la_SOURCES = $(BLASWARPPER_SOURCES)

libblaswrapper_dd_la_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/$(QD_INCLUDEDIR) -D___MPACK_BUILD_WITH_DD___
libblaswrapper_dd_la_SOURCES = $(BLASWARPPER_SOURCES)

if ENABLE_DOUBLE
libblaswrapper_double_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPACK_BUILD_WITH_DOUBLE___
libblaswrapper_double_la_SOURCES = $(BLASWRAPPER_SOURCES)
endif

endif
