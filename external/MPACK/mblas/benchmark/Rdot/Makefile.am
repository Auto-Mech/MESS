bench_PROGRAMS = \
Rdot.dd       \
Rdot.dd_ref   \
Rdot.qd       \
Rdot.qd_ref   \
Rdot.gmp      \
Rdot.gmp_ref  \
Rdot.mpfr     \
Rdot.mpfr_ref \
Rdot.double   \
Rdot.double_ref
if !IS_MINGW
bench_PROGRAMS+=ddot
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

Rdot_dd_ref_SOURCES= Rdot_dd.cpp
Rdot_dd_ref_CXXFLAGS= $(ddcxxflags)
Rdot_dd_ref_LDFLAGS= $(ddref_libdepends)

Rdot_dd_SOURCES= Rdot_dd.cpp
Rdot_dd_CXXFLAGS= $(ddcxxflags)
Rdot_dd_LDFLAGS= $(ddopt_libdepends)

Rdot_qd_ref_SOURCES= Rdot_qd.cpp
Rdot_qd_ref_CXXFLAGS= $(qdcxxflags)
Rdot_qd_ref_LDFLAGS= $(qdref_libdepends)

Rdot_qd_SOURCES= Rdot_qd.cpp
Rdot_qd_CXXFLAGS= $(qdcxxflags)
Rdot_qd_LDFLAGS= $(qdopt_libdepends)

Rdot_gmp_ref_SOURCES= Rdot_gmp.cpp
Rdot_gmp_ref_CXXFLAGS= $(gmpcxxflags)
Rdot_gmp_ref_LDFLAGS= $(gmpref_libdepends)

Rdot_gmp_SOURCES= Rdot_gmp.cpp
Rdot_gmp_CXXFLAGS= $(gmpcxxflags)
Rdot_gmp_LDFLAGS= $(gmpopt_libdepends)

Rdot_mpfr_ref_SOURCES= Rdot_mpfr.cpp
Rdot_mpfr_ref_CXXFLAGS= $(mpfrcxxflags)
Rdot_mpfr_ref_LDFLAGS= $(mpfrref_libdepends)

Rdot_mpfr_SOURCES= Rdot_mpfr.cpp
Rdot_mpfr_CXXFLAGS= $(mpfrcxxflags)
Rdot_mpfr_LDFLAGS= $(mpfropt_libdepends)

Rdot_double_SOURCES= Rdot_double.cpp
Rdot_double_CXXFLAGS= $(doublecxxflags)
Rdot_double_LDFLAGS= $(doubleopt_libdepends)

Rdot_double_ref_SOURCES= Rdot_double.cpp
Rdot_double_ref_CXXFLAGS= $(doublecxxflags)
Rdot_double_ref_LDFLAGS= $(doubleref_libdepends)

if !IS_MINGW
ddot_SOURCES= ddot.cpp
ddot_CXXFLAGS= $(ref_cxxflags)
ddot_LDFLAGS= $(ref_libdepends)
endif
