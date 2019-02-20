bench_PROGRAMS = \
Rgemv.dd       \
Rgemv.dd_ref   \
Rgemv.qd       \
Rgemv.qd_ref   \
Rgemv.gmp      \
Rgemv.gmp_ref  \
Rgemv.mpfr     \
Rgemv.mpfr_ref \
Rgemv.double   \
Rgemv.double_ref

if !IS_MINGW
bench_PROGRAMS+=dgemv
endif

benchdir=$(prefix)/share/mpack/bench
ddcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(QD_INCLUDEDIR) 
ddref_libdepends = -L$(top_builddir)/mblas/reference -lmblas_dd_ref -L$(top_srcdir)/$(QD_LIBDIR) -lqd
ddopt_libdepends = -L$(top_builddir)/mblas/optimized/dd -lmblas_dd -L$(top_srcdir)/$(QD_LIBDIR) -lqd

qdcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(QD_INCLUDEDIR) 
qdref_libdepends = -L$(top_builddir)/mblas/reference -lmblas_qd_ref -L$(top_srcdir)/$(QD_LIBDIR) -lqd
qdopt_libdepends = -L$(top_builddir)/mblas/optimized/qd -lmblas_qd -L$(top_srcdir)/$(QD_LIBDIR) -lqd

gmpcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(GMP_INCLUDEDIR) 
gmpref_libdepends = -L$(top_builddir)/mblas/reference -lmblas_gmp_ref -L$(top_srcdir)/$(GMP_LIBDIR) -lgmpxx -lgmp
gmpopt_libdepends = -L$(top_builddir)/mblas/optimized/gmp -lmblas_gmp -L$(top_srcdir)/$(GMP_LIBDIR) -lgmpxx -lgmp

mpfrcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(top_srcdir)/$(GMP_INCLUDEDIR) -I$(top_srcdir)/$(MPFR_INCLUDEDIR) -I$(top_srcdir)/$(MPC_INCLUDEDIR)
mpfrref_libdepends = -L$(top_builddir)/mblas/reference -lmblas_mpfr_ref -L$(top_srcdir)/$(MPC_LIBDIR) -L$(top_srcdir)/$(MPFR_LIBDIR) -L$(top_srcdir)/$(GMP_LIBDIR) -L$(top_builddir)/mpfrc++ -lmpfrcxx -lmpfr -lmpc -lgmp
mpfropt_libdepends = -L$(top_builddir)/mblas/optimized/mpfr -lmblas_mpfr -L$(top_srcdir)/$(MPC_LIBDIR) -L$(top_srcdir)/$(MPFR_LIBDIR) -L$(top_srcdir)/$(GMP_LIBDIR) -L$(top_builddir)/mpfrc++ -lmpfrcxx -lmpfr -lmpc -lgmp

doublecxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
doubleref_libdepends = -L$(top_builddir)/mblas/reference -lmblas_double_ref
doubleopt_libdepends = -L$(top_builddir)/mblas/optimized/double -lmblas_double

ref_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
ref_libdepends = $(BLAS_LIBS) $(FLIBS)

Rgemv_dd_ref_SOURCES= Rgemv_dd.cpp
Rgemv_dd_ref_CXXFLAGS= $(ddcxxflags)
Rgemv_dd_ref_LDFLAGS= $(ddref_libdepends)

Rgemv_dd_SOURCES= Rgemv_dd.cpp
Rgemv_dd_CXXFLAGS= $(ddcxxflags)
Rgemv_dd_LDFLAGS= $(ddopt_libdepends)

Rgemv_qd_ref_SOURCES= Rgemv_qd.cpp
Rgemv_qd_ref_CXXFLAGS= $(qdcxxflags)
Rgemv_qd_ref_LDFLAGS= $(qdref_libdepends)

Rgemv_qd_SOURCES= Rgemv_qd.cpp
Rgemv_qd_CXXFLAGS= $(qdcxxflags)
Rgemv_qd_LDFLAGS= $(qdopt_libdepends)

Rgemv_gmp_ref_SOURCES= Rgemv_gmp.cpp
Rgemv_gmp_ref_CXXFLAGS= $(gmpcxxflags)
Rgemv_gmp_ref_LDFLAGS= $(gmpref_libdepends)

Rgemv_gmp_SOURCES= Rgemv_gmp.cpp
Rgemv_gmp_CXXFLAGS= $(gmpcxxflags)
Rgemv_gmp_LDFLAGS= $(gmpopt_libdepends)

Rgemv_mpfr_ref_SOURCES= Rgemv_mpfr.cpp
Rgemv_mpfr_ref_CXXFLAGS= $(mpfrcxxflags)
Rgemv_mpfr_ref_LDFLAGS= $(mpfrref_libdepends)

Rgemv_mpfr_SOURCES= Rgemv_mpfr.cpp
Rgemv_mpfr_CXXFLAGS= $(mpfrcxxflags)
Rgemv_mpfr_LDFLAGS= $(mpfropt_libdepends)

Rgemv_double_SOURCES= Rgemv_double.cpp
Rgemv_double_CXXFLAGS= $(doublecxxflags)
Rgemv_double_LDFLAGS= $(doubleopt_libdepends)

Rgemv_double_ref_SOURCES= Rgemv_double.cpp
Rgemv_double_ref_CXXFLAGS= $(doublecxxflags)
Rgemv_double_ref_LDFLAGS= $(doubleref_libdepends)

dgemv_SOURCES= dgemv.cpp
dgemv_CXXFLAGS= $(ref_cxxflags)
dgemv_LDFLAGS= $(ref_libdepends)
