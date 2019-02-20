bench_PROGRAMS = \
Rcopy.dd       \
Rcopy.dd_ref   \
Rcopy.qd       \
Rcopy.qd_ref   \
Rcopy.gmp      \
Rcopy.gmp_ref  \
Rcopy.mpfr     \
Rcopy.mpfr_ref \
Rcopy.double   \
Rcopy.double_ref

if !IS_MINGW
bench_PROGRAMS+=dcopy
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

Rcopy_dd_ref_SOURCES= Rcopy_dd.cpp
Rcopy_dd_ref_CXXFLAGS= $(ddcxxflags)
Rcopy_dd_ref_LDFLAGS= $(ddref_libdepends)

Rcopy_dd_SOURCES= Rcopy_dd.cpp
Rcopy_dd_CXXFLAGS= $(ddcxxflags)
Rcopy_dd_LDFLAGS= $(ddopt_libdepends)

Rcopy_qd_ref_SOURCES= Rcopy_qd.cpp
Rcopy_qd_ref_CXXFLAGS= $(qdcxxflags)
Rcopy_qd_ref_LDFLAGS= $(qdref_libdepends)

Rcopy_qd_SOURCES= Rcopy_qd.cpp
Rcopy_qd_CXXFLAGS= $(qdcxxflags)
Rcopy_qd_LDFLAGS= $(qdopt_libdepends)

Rcopy_gmp_ref_SOURCES= Rcopy_gmp.cpp
Rcopy_gmp_ref_CXXFLAGS= $(gmpcxxflags)
Rcopy_gmp_ref_LDFLAGS= $(gmpref_libdepends)

Rcopy_gmp_SOURCES= Rcopy_gmp.cpp
Rcopy_gmp_CXXFLAGS= $(gmpcxxflags)
Rcopy_gmp_LDFLAGS= $(gmpopt_libdepends)

Rcopy_mpfr_ref_SOURCES= Rcopy_mpfr.cpp
Rcopy_mpfr_ref_CXXFLAGS= $(mpfrcxxflags)
Rcopy_mpfr_ref_LDFLAGS= $(mpfrref_libdepends)

Rcopy_mpfr_SOURCES= Rcopy_mpfr.cpp
Rcopy_mpfr_CXXFLAGS= $(mpfrcxxflags)
Rcopy_mpfr_LDFLAGS= $(mpfropt_libdepends)

Rcopy_double_SOURCES= Rcopy_double.cpp
Rcopy_double_CXXFLAGS= $(doublecxxflags)
Rcopy_double_LDFLAGS= $(doubleopt_libdepends)

Rcopy_double_ref_SOURCES= Rcopy_double.cpp
Rcopy_double_ref_CXXFLAGS= $(doublecxxflags)
Rcopy_double_ref_LDFLAGS= $(doubleref_libdepends)

if !IS_MINGW
dcopy_SOURCES= dcopy.cpp
dcopy_CXXFLAGS= $(ref_cxxflags)
dcopy_LDFLAGS= $(ref_libdepends)
endif
