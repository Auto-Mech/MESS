#CTBSV, CTPSV might have bugs or intrinsically inaccurate
#Rnrm2 should be checked more
check_PROGRAMS = $(mblas_qd_test_PROGRAMS) 
TESTS = $(check_PROGRAMS)
mblas_qd_test_PROGRAMS = \
arithmetic.debug_qd \
complex.debug_qd \
iCamax.debug.qd_ref iCamax.debug.qd   \
iRamax.debug.qd_ref iRamax.debug.qd   \
Mlsame.debug.qd_ref Mlsame.debug.qd   \
Mxerbla.debug.qd_ref Mxerbla.debug.qd \
CRrot.debug.qd_ref  CRrot.debug.qd    \
CRscal.debug.qd_ref CRscal.debug.qd   \
Caxpy.debug.qd_ref  Caxpy.debug.qd    \
Ccopy.debug.qd_ref  Ccopy.debug.qd    \
Cdotc.debug.qd_ref  Cdotc.debug.qd	\
Cdotu.debug.qd_ref  Cdotu.debug.qd	\
Cgbmv.debug.qd_ref  Cgbmv.debug.qd	\
Cgemm.debug.qd_ref  Cgemm.debug.qd	\
Cgemv.debug.qd_ref  Cgemv.debug.qd	\
Cgerc.debug.qd_ref  Cgerc.debug.qd	\
Cgeru.debug.qd_ref  Cgeru.debug.qd    \
Chbmv.debug.qd_ref  Chbmv.debug.qd	\
Chemm.debug.qd_ref  Chemm.debug.qd	\
Chemv.debug.qd_ref  Chemv.debug.qd	\
Cher.debug.qd_ref   Cher.debug.qd 	\
Cher2.debug.qd_ref  Cher2.debug.qd	\
Cher2k.debug.qd_ref Cher2k.debug.qd	\
Cherk.debug.qd_ref  Cherk.debug.qd	\
Chpmv.debug.qd_ref  Chpmv.debug.qd	\
Chpr.debug.qd_ref   Chpr.debug.qd	\
Chpr2.debug.qd_ref  Chpr2.debug.qd	\
Crotg.debug.qd_ref  Crotg.debug.qd	\
Cscal.debug.qd_ref  Cscal.debug.qd	\
Cswap.debug.qd_ref  Cswap.debug.qd	\
Csymm.debug.qd_ref  Csymm.debug.qd    \
Csyr2k.debug.qd_ref Csyr2k.debug.qd	\
Csyrk.debug.qd_ref  Csyrk.debug.qd	\
Ctbmv.debug.qd_ref  Ctbmv.debug.qd	\
Ctbsv.debug.qd_ref  Ctbsv.debug.qd	\
Ctpmv.debug.qd_ref  Ctpmv.debug.qd	\
Ctpsv.debug.qd_ref  Ctpsv.debug.qd	\
Ctrmm.debug.qd_ref  Ctrmm.debug.qd	\
Ctrmv.debug.qd_ref  Ctrmv.debug.qd	\
Ctrsm.debug.qd_ref  Ctrsm.debug.qd	\
Ctrsv.debug.qd_ref  Ctrsv.debug.qd	\
RCabs1.debug.qd_ref RCabs1.debug.qd	\
RCasum.debug.qd_ref RCasum.debug.qd	\
RCnrm2.debug.qd_ref RCnrm2.debug.qd	\
Rasum.debug.qd_ref  Rasum.debug.qd    \
Raxpy.debug.qd_ref  Raxpy.debug.qd	\
Rcopy.debug.qd_ref  Rcopy.debug.qd	\
Rdot.debug.qd_ref   Rdot.debug.qd	\
Rgbmv.debug.qd_ref  Rgbmv.debug.qd	\
Rgemm.debug.qd_ref  Rgemm.debug.qd	\
Rgemv.debug.qd_ref  Rgemv.debug.qd	\
Rger.debug.qd_ref   Rger.debug.qd	\
Rnrm2.debug.qd_ref  Rnrm2.debug.qd	\
Rrot.debug.qd_ref   Rrot.debug.qd	\
Rrotg.debug.qd_ref  Rrotg.debug.qd	\
Rrotm.debug.qd_ref  Rrotm.debug.qd	\
Rsbmv.debug.qd_ref  Rsbmv.debug.qd	\
Rscal.debug.qd_ref  Rscal.debug.qd	\
Rspmv.debug.qd_ref  Rspmv.debug.qd    \
Rspr.debug.qd_ref   Rspr.debug.qd	\
Rspr2.debug.qd_ref  Rspr2.debug.qd	\
Rswap.debug.qd_ref  Rswap.debug.qd	\
Rsymm.debug.qd_ref  Rsymm.debug.qd	\
Rsymv.debug.qd_ref  Rsymv.debug.qd	\
Rsyr.debug.qd_ref   Rsyr.debug.qd	\
Rsyr2.debug.qd_ref  Rsyr2.debug.qd	\
Rsyr2k.debug.qd_ref Rsyr2k.debug.qd	\
Rsyrk.debug.qd_ref  Rsyrk.debug.qd	\
Rtbmv.debug.qd_ref  Rtbmv.debug.qd	\
Rtbsv.debug.qd_ref  Rtbsv.debug.qd	\
Rtpmv.debug.qd_ref  Rtpmv.debug.qd	\
Rtpsv.debug.qd_ref  Rtpsv.debug.qd	\
Rtrmm.debug.qd_ref  Rtrmm.debug.qd    \
Rtrmv.debug.qd_ref  Rtrmv.debug.qd	\
Rtrsm.debug.qd_ref  Rtrsm.debug.qd	\
Rtrsv.debug.qd_ref  Rtrsv.debug.qd

mblas_qd_testdir = $(prefix)/share/mpack/testing/qd
ref_libdepends = -L$(top_builddir)/mblas/reference -L$(top_builddir)/mblas/reference -lmblas_qd_ref -lmblas_mpfr_ref
opt_libdepends = -L$(top_builddir)/mblas/reference -L$(top_builddir)/mblas/optimized/qd -lmblas_qd -lmblas_mpfr_ref
mplibs=-L$(top_srcdir)/$(MPC_LIBDIR) -L$(top_srcdir)/$(MPFR_LIBDIR) -L$(top_srcdir)/$(GMP_LIBDIR) -L$(top_builddir)/mpfrc++ -lmpfrcxx -lmpc -lmpfr -lgmpxx -lgmp -L$(top_srcdir)/$(QD_LIBDIR) -lqd
blas_libs=$(BLAS_LIBS) $(FLIBS)

debug_srcdepends = ../common/xerbla.cpp ../common/mpack.debug.cpp
debug_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(GMP_INCLUDEDIR) -I$(top_srcdir)/$(MPFR_INCLUDEDIR) -I$(top_srcdir)/$(MPC_INCLUDEDIR) -I$(top_srcdir)/$(QD_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPACK_BUILD_WITH_QD___ -DVERBOSE_TEST

##
arithmetic_debug_qd_SOURCES=../common/arithmetic.debug.cpp $(debug_srcdepends)
arithmetic_debug_qd_CXXFLAGS=$(debug_cxxflags)
arithmetic_debug_qd_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

complex_debug_qd_SOURCES=../common/complex.debug.cpp $(debug_srcdepends)
complex_debug_qd_CXXFLAGS=$(debug_cxxflags)
complex_debug_qd_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

iCamax_debug_qd_ref_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
iCamax_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

iCamax_debug_qd_SOURCES=../common/iCamax.debug.cpp $(debug_srcdepends)
iCamax_debug_qd_CXXFLAGS=$(debug_cxxflags)
iCamax_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

iRamax_debug_qd_ref_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
iRamax_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

iRamax_debug_qd_SOURCES=../common/iRamax.debug.cpp $(debug_srcdepends)
iRamax_debug_qd_CXXFLAGS=$(debug_cxxflags)
iRamax_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Mlsame_debug_qd_ref_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Mlsame_debug_qd_SOURCES=../common/Mlsame.debug.cpp $(debug_srcdepends)
Mlsame_debug_qd_CXXFLAGS=$(debug_cxxflags)
Mlsame_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Mxerbla_debug_qd_ref_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Mxerbla_debug_qd_SOURCES=../common/Mxerbla.debug.cpp $(debug_srcdepends)
Mxerbla_debug_qd_CXXFLAGS=$(debug_cxxflags)
Mxerbla_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

CRrot_debug_qd_ref_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
CRrot_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

CRrot_debug_qd_SOURCES=../common/CRrot.debug.cpp $(debug_srcdepends)
CRrot_debug_qd_CXXFLAGS=$(debug_cxxflags)
CRrot_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

CRscal_debug_qd_ref_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
CRscal_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

CRscal_debug_qd_SOURCES=../common/CRscal.debug.cpp $(debug_srcdepends)
CRscal_debug_qd_CXXFLAGS=$(debug_cxxflags)
CRscal_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Caxpy_debug_qd_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug_qd_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Caxpy_debug_qd_ref_SOURCES=../common/Caxpy.debug.cpp $(debug_srcdepends)
Caxpy_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Caxpy_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Ccopy_debug_qd_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug_qd_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Ccopy_debug_qd_ref_SOURCES=../common/Ccopy.debug.cpp $(debug_srcdepends)
Ccopy_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Ccopy_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cdotc_debug_qd_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cdotc_debug_qd_ref_SOURCES=../common/Cdotc.debug.cpp $(debug_srcdepends)
Cdotc_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cdotc_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cdotu_debug_qd_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cdotu_debug_qd_ref_SOURCES=../common/Cdotu.debug.cpp $(debug_srcdepends)
Cdotu_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cdotu_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cgbmv_debug_qd_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cgbmv_debug_qd_ref_SOURCES=../common/Cgbmv.debug.cpp $(debug_srcdepends)
Cgbmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cgbmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cgemm_debug_qd_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cgemm_debug_qd_ref_SOURCES=../common/Cgemm.debug.cpp $(debug_srcdepends)
Cgemm_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cgemm_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cgemv_debug_qd_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cgemv_debug_qd_ref_SOURCES=../common/Cgemv.debug.cpp $(debug_srcdepends)
Cgemv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cgemv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cgerc_debug_qd_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cgerc_debug_qd_ref_SOURCES=../common/Cgerc.debug.cpp $(debug_srcdepends)
Cgerc_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cgerc_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cgeru_debug_qd_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cgeru_debug_qd_ref_SOURCES=../common/Cgeru.debug.cpp $(debug_srcdepends)
Cgeru_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cgeru_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Chbmv_debug_qd_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Chbmv_debug_qd_ref_SOURCES=../common/Chbmv.debug.cpp $(debug_srcdepends)
Chbmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Chbmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Chemm_debug_qd_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug_qd_CXXFLAGS=$(debug_cxxflags)
Chemm_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Chemm_debug_qd_ref_SOURCES=../common/Chemm.debug.cpp $(debug_srcdepends)
Chemm_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Chemm_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Chemv_debug_qd_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Chemv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Chemv_debug_qd_ref_SOURCES=../common/Chemv.debug.cpp $(debug_srcdepends)
Chemv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Chemv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cher_debug_qd_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cher_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cher_debug_qd_ref_SOURCES=../common/Cher.debug.cpp $(debug_srcdepends)
Cher_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cher_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cher2_debug_qd_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cher2_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cher2_debug_qd_ref_SOURCES=../common/Cher2.debug.cpp $(debug_srcdepends)
Cher2_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cher2_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cher2k_debug_qd_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cher2k_debug_qd_ref_SOURCES=../common/Cher2k.debug.cpp $(debug_srcdepends)
Cher2k_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cher2k_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cherk_debug_qd_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cherk_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cherk_debug_qd_ref_SOURCES=../common/Cherk.debug.cpp $(debug_srcdepends)
Cherk_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cherk_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Chpmv_debug_qd_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Chpmv_debug_qd_ref_SOURCES=../common/Chpmv.debug.cpp $(debug_srcdepends)
Chpmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Chpmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Chpr_debug_qd_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug_qd_CXXFLAGS=$(debug_cxxflags)
Chpr_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Chpr_debug_qd_ref_SOURCES=../common/Chpr.debug.cpp $(debug_srcdepends)
Chpr_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Chpr_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Chpr2_debug_qd_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug_qd_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Chpr2_debug_qd_ref_SOURCES=../common/Chpr2.debug.cpp $(debug_srcdepends)
Chpr2_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Chpr2_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Crotg_debug_qd_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug_qd_CXXFLAGS=$(debug_cxxflags)
Crotg_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Crotg_debug_qd_ref_SOURCES=../common/Crotg.debug.cpp $(debug_srcdepends)
Crotg_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Crotg_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cscal_debug_qd_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cscal_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cscal_debug_qd_ref_SOURCES=../common/Cscal.debug.cpp $(debug_srcdepends)
Cscal_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cscal_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Cswap_debug_qd_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug_qd_CXXFLAGS=$(debug_cxxflags)
Cswap_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Cswap_debug_qd_ref_SOURCES=../common/Cswap.debug.cpp $(debug_srcdepends)
Cswap_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Cswap_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Csymm_debug_qd_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug_qd_CXXFLAGS=$(debug_cxxflags)
Csymm_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Csymm_debug_qd_ref_SOURCES=../common/Csymm.debug.cpp $(debug_srcdepends)
Csymm_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Csymm_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Csyr2k_debug_qd_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug_qd_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Csyr2k_debug_qd_ref_SOURCES=../common/Csyr2k.debug.cpp $(debug_srcdepends)
Csyr2k_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Csyr2k_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Csyrk_debug_qd_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug_qd_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Csyrk_debug_qd_ref_SOURCES=../common/Csyrk.debug.cpp $(debug_srcdepends)
Csyrk_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Csyrk_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Ctbmv_debug_qd_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Ctbmv_debug_qd_ref_SOURCES=../common/Ctbmv.debug.cpp $(debug_srcdepends)
Ctbmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Ctbmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Ctbsv_debug_qd_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Ctbsv_debug_qd_ref_SOURCES=../common/Ctbsv.debug.cpp $(debug_srcdepends)
Ctbsv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Ctbsv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Ctpmv_debug_qd_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Ctpmv_debug_qd_ref_SOURCES=../common/Ctpmv.debug.cpp $(debug_srcdepends)
Ctpmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Ctpmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Ctpsv_debug_qd_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Ctpsv_debug_qd_ref_SOURCES=../common/Ctpsv.debug.cpp $(debug_srcdepends)
Ctpsv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Ctpsv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Ctrmm_debug_qd_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug_qd_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Ctrmm_debug_qd_ref_SOURCES=../common/Ctrmm.debug.cpp $(debug_srcdepends)
Ctrmm_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmm_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Ctrmv_debug_qd_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Ctrmv_debug_qd_ref_SOURCES=../common/Ctrmv.debug.cpp $(debug_srcdepends)
Ctrmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Ctrmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Ctrsm_debug_qd_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug_qd_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Ctrsm_debug_qd_ref_SOURCES=../common/Ctrsm.debug.cpp $(debug_srcdepends)
Ctrsm_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsm_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Ctrsv_debug_qd_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Ctrsv_debug_qd_ref_SOURCES=../common/Ctrsv.debug.cpp $(debug_srcdepends)
Ctrsv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Ctrsv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

RCabs1_debug_qd_ref_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

RCabs1_debug_qd_SOURCES=../common/RCabs1.debug.cpp $(debug_srcdepends)
RCabs1_debug_qd_CXXFLAGS=$(debug_cxxflags)
RCabs1_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

RCasum_debug_qd_ref_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
RCasum_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

RCasum_debug_qd_SOURCES=../common/RCasum.debug.cpp $(debug_srcdepends)
RCasum_debug_qd_CXXFLAGS=$(debug_cxxflags)
RCasum_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

RCnrm2_debug_qd_ref_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

RCnrm2_debug_qd_SOURCES=../common/RCnrm2.debug.cpp $(debug_srcdepends)
RCnrm2_debug_qd_CXXFLAGS=$(debug_cxxflags)
RCnrm2_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rasum_debug_qd_ref_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rasum_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rasum_debug_qd_SOURCES=../common/Rasum.debug.cpp $(debug_srcdepends)
Rasum_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rasum_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Raxpy_debug_qd_ref_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Raxpy_debug_qd_SOURCES=../common/Raxpy.debug.cpp $(debug_srcdepends)
Raxpy_debug_qd_CXXFLAGS=$(debug_cxxflags)
Raxpy_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rcopy_debug_qd_ref_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rcopy_debug_qd_SOURCES=../common/Rcopy.debug.cpp $(debug_srcdepends)
Rcopy_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rcopy_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rdot_debug_qd_ref_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rdot_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rdot_debug_qd_SOURCES=../common/Rdot.debug.cpp $(debug_srcdepends)
Rdot_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rdot_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rgbmv_debug_qd_ref_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rgbmv_debug_qd_SOURCES=../common/Rgbmv.debug.cpp $(debug_srcdepends)
Rgbmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rgbmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rgemm_debug_qd_ref_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rgemm_debug_qd_SOURCES=../common/Rgemm.debug.cpp $(debug_srcdepends)
Rgemm_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rgemm_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rgemv_debug_qd_ref_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rgemv_debug_qd_SOURCES=../common/Rgemv.debug.cpp $(debug_srcdepends)
Rgemv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rgemv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rger_debug_qd_ref_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rger_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rger_debug_qd_SOURCES=../common/Rger.debug.cpp $(debug_srcdepends)
Rger_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rger_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rnrm2_debug_qd_ref_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rnrm2_debug_qd_SOURCES=../common/Rnrm2.debug.cpp $(debug_srcdepends)
Rnrm2_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rnrm2_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rrot_debug_qd_ref_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rrot_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rrot_debug_qd_SOURCES=../common/Rrot.debug.cpp $(debug_srcdepends)
Rrot_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rrot_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rrotg_debug_qd_ref_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rrotg_debug_qd_SOURCES=../common/Rrotg.debug.cpp $(debug_srcdepends)
Rrotg_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rrotg_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rrotm_debug_qd_ref_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rrotm_debug_qd_SOURCES=../common/Rrotm.debug.cpp $(debug_srcdepends)
Rrotm_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rrotm_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rsbmv_debug_qd_ref_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rsbmv_debug_qd_SOURCES=../common/Rsbmv.debug.cpp $(debug_srcdepends)
Rsbmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rsbmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rscal_debug_qd_ref_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rscal_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rscal_debug_qd_SOURCES=../common/Rscal.debug.cpp $(debug_srcdepends)
Rscal_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rscal_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rspmv_debug_qd_ref_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rspmv_debug_qd_SOURCES=../common/Rspmv.debug.cpp $(debug_srcdepends)
Rspmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rspmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rspr_debug_qd_ref_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rspr_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rspr_debug_qd_SOURCES=../common/Rspr.debug.cpp $(debug_srcdepends)
Rspr_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rspr_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rspr2_debug_qd_ref_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rspr2_debug_qd_SOURCES=../common/Rspr2.debug.cpp $(debug_srcdepends)
Rspr2_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rspr2_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rswap_debug_qd_ref_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rswap_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rswap_debug_qd_SOURCES=../common/Rswap.debug.cpp $(debug_srcdepends)
Rswap_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rswap_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rsymm_debug_qd_ref_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rsymm_debug_qd_SOURCES=../common/Rsymm.debug.cpp $(debug_srcdepends)
Rsymm_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rsymm_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rsymv_debug_qd_ref_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rsymv_debug_qd_SOURCES=../common/Rsymv.debug.cpp $(debug_srcdepends)
Rsymv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rsymv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rsyr_debug_qd_ref_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rsyr_debug_qd_SOURCES=../common/Rsyr.debug.cpp $(debug_srcdepends)
Rsyr_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rsyr_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rsyr2_debug_qd_ref_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rsyr2_debug_qd_SOURCES=../common/Rsyr2.debug.cpp $(debug_srcdepends)
Rsyr2_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rsyr2_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rsyr2k_debug_qd_ref_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rsyr2k_debug_qd_SOURCES=../common/Rsyr2k.debug.cpp $(debug_srcdepends)
Rsyr2k_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rsyr2k_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rsyrk_debug_qd_ref_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rsyrk_debug_qd_SOURCES=../common/Rsyrk.debug.cpp $(debug_srcdepends)
Rsyrk_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rsyrk_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rtbmv_debug_qd_ref_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rtbmv_debug_qd_SOURCES=../common/Rtbmv.debug.cpp $(debug_srcdepends)
Rtbmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rtbmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rtbsv_debug_qd_ref_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rtbsv_debug_qd_SOURCES=../common/Rtbsv.debug.cpp $(debug_srcdepends)
Rtbsv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rtbsv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rtpmv_debug_qd_ref_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rtpmv_debug_qd_SOURCES=../common/Rtpmv.debug.cpp $(debug_srcdepends)
Rtpmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rtpmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rtpsv_debug_qd_ref_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rtpsv_debug_qd_SOURCES=../common/Rtpsv.debug.cpp $(debug_srcdepends)
Rtpsv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rtpsv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rtrmm_debug_qd_ref_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rtrmm_debug_qd_SOURCES=../common/Rtrmm.debug.cpp $(debug_srcdepends)
Rtrmm_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rtrmm_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rtrmv_debug_qd_ref_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rtrmv_debug_qd_SOURCES=../common/Rtrmv.debug.cpp $(debug_srcdepends)
Rtrmv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rtrmv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rtrsm_debug_qd_ref_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rtrsm_debug_qd_SOURCES=../common/Rtrsm.debug.cpp $(debug_srcdepends)
Rtrsm_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rtrsm_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

Rtrsv_debug_qd_ref_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug_qd_ref_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug_qd_ref_LDFLAGS=$(ref_libdepends) $(mplibs) $(blas_libs)

Rtrsv_debug_qd_SOURCES=../common/Rtrsv.debug.cpp $(debug_srcdepends)
Rtrsv_debug_qd_CXXFLAGS=$(debug_cxxflags)
Rtrsv_debug_qd_LDFLAGS=$(opt_libdepends) $(mplibs) $(blas_libs)

