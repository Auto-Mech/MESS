blasexamples_PROGRAMS = \
Rgemm.dd       \
Rgemm.qd       \
Rgemm.gmp      \
Rgemm.mpfr     \
Rgemm.double

blasexamplesdir=$(prefix)/share/mpack/examples
ddcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(QD_INCLUDEDIR) 
dd_libdepends = -L$(top_builddir)/mblas/optimized/dd -lmblas_dd -L$(top_srcdir)/$(QD_LIBDIR) -lqd

qdcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(QD_INCLUDEDIR) 
qd_libdepends = -L$(top_builddir)/mblas/optimized/qd -lmblas_qd -L$(top_srcdir)/$(QD_LIBDIR) -lqd

gmpcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(GMP_INCLUDEDIR) 
gmp_libdepends = -L$(top_builddir)/mblas/optimized/gmp -lmblas_gmp -L$(top_srcdir)/$(GMP_LIBDIR) -lgmpxx -lgmp

mpfrcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(top_srcdir)/$(GMP_INCLUDEDIR) -I$(top_srcdir)/$(MPFR_INCLUDEDIR) -I$(top_srcdir)/$(MPC_INCLUDEDIR)
mpfr_libdepends = -L$(top_builddir)/mblas/optimized/mpfr -lmblas_mpfr -L$(top_builddir)/mpfrc++ -lmpfrcxx -L$(top_srcdir)/$(MPC_LIBDIR) -L$(top_srcdir)/$(MPFR_LIBDIR) -L$(top_srcdir)/$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

doublecxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -L$(top_builddir)/mblas/optimized/double -lmblas_double

Rgemm_dd_SOURCES= Rgemm_dd.cpp
Rgemm_dd_CXXFLAGS= $(ddcxxflags)
Rgemm_dd_LDFLAGS= $(dd_libdepends)

Rgemm_qd_SOURCES= Rgemm_qd.cpp
Rgemm_qd_CXXFLAGS= $(qdcxxflags)
Rgemm_qd_LDFLAGS= $(qd_libdepends)

Rgemm_gmp_SOURCES= Rgemm_gmp.cpp
Rgemm_gmp_CXXFLAGS= $(gmpcxxflags)
Rgemm_gmp_LDFLAGS= $(gmp_libdepends)

Rgemm_mpfr_SOURCES= Rgemm_mpfr.cpp
Rgemm_mpfr_CXXFLAGS= $(mpfrcxxflags)
Rgemm_mpfr_LDFLAGS= $(mpfr_libdepends)

Rgemm_double_SOURCES= Rgemm_double.cpp
Rgemm_double_CXXFLAGS= $(doublecxxflags)
Rgemm_double_LDFLAGS= $(double_libdepends)

