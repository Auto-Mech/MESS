bench_PROGRAMS = \
Rgemm.dd       \
Rgemm.dd_ref   \
Rgemm.qd       \
Rgemm.qd_ref   \
Rgemm.gmp      \
Rgemm.gmp_ref  \
Rgemm.mpfr     \
Rgemm.mpfr_ref \
Rgemm.double   \
Rgemm.double_ref

if !IS_MINGW
bench_PROGRAMS+=dgemm
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

Rgemm_dd_ref_SOURCES= Rgemm_dd.cpp
Rgemm_dd_ref_CXXFLAGS= $(ddcxxflags)
Rgemm_dd_ref_LDFLAGS= $(ddref_libdepends)

Rgemm_dd_SOURCES= Rgemm_dd.cpp
Rgemm_dd_CXXFLAGS= $(ddcxxflags)
Rgemm_dd_LDFLAGS= $(ddopt_libdepends)

Rgemm_qd_ref_SOURCES= Rgemm_qd.cpp
Rgemm_qd_ref_CXXFLAGS= $(qdcxxflags)
Rgemm_qd_ref_LDFLAGS= $(qdref_libdepends)

Rgemm_qd_SOURCES= Rgemm_qd.cpp
Rgemm_qd_CXXFLAGS= $(qdcxxflags)
Rgemm_qd_LDFLAGS= $(qdopt_libdepends)

Rgemm_gmp_ref_SOURCES= Rgemm_gmp.cpp
Rgemm_gmp_ref_CXXFLAGS= $(gmpcxxflags)
Rgemm_gmp_ref_LDFLAGS= $(gmpref_libdepends)

Rgemm_gmp_SOURCES= Rgemm_gmp.cpp
Rgemm_gmp_CXXFLAGS= $(gmpcxxflags)
Rgemm_gmp_LDFLAGS= $(gmpopt_libdepends)

Rgemm_mpfr_ref_SOURCES= Rgemm_mpfr.cpp
Rgemm_mpfr_ref_CXXFLAGS= $(mpfrcxxflags)
Rgemm_mpfr_ref_LDFLAGS= $(mpfrref_libdepends)

Rgemm_mpfr_SOURCES= Rgemm_mpfr.cpp
Rgemm_mpfr_CXXFLAGS= $(mpfrcxxflags)
Rgemm_mpfr_LDFLAGS= $(mpfropt_libdepends)

Rgemm_double_SOURCES= Rgemm_double.cpp
Rgemm_double_CXXFLAGS= $(doublecxxflags)
Rgemm_double_LDFLAGS= $(doubleopt_libdepends)

Rgemm_double_ref_SOURCES= Rgemm_double.cpp
Rgemm_double_ref_CXXFLAGS= $(doublecxxflags)
Rgemm_double_ref_LDFLAGS= $(doubleref_libdepends)

if !IS_MINGW
dgemm_SOURCES= dgemm.cpp
dgemm_CXXFLAGS= $(ref_cxxflags)
dgemm_LDFLAGS= $(ref_libdepends)
endif
