bench_PROGRAMS = \
Raxpy.dd       \
Raxpy.dd_ref   \
Raxpy.qd       \
Raxpy.qd_ref   \
Raxpy.gmp      \
Raxpy.gmp_ref  \
Raxpy.mpfr     \
Raxpy.mpfr_ref \
Raxpy.double   \
Raxpy.double_ref

if !IS_MINGW
bench_PROGRAMS+=daxpy
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

Raxpy_dd_ref_SOURCES= Raxpy_dd.cpp
Raxpy_dd_ref_CXXFLAGS= $(ddcxxflags)
Raxpy_dd_ref_LDFLAGS= $(ddref_libdepends)

Raxpy_dd_SOURCES= Raxpy_dd.cpp
Raxpy_dd_CXXFLAGS= $(ddcxxflags)
Raxpy_dd_LDFLAGS= $(ddopt_libdepends)

Raxpy_qd_ref_SOURCES= Raxpy_qd.cpp
Raxpy_qd_ref_CXXFLAGS= $(qdcxxflags)
Raxpy_qd_ref_LDFLAGS= $(qdref_libdepends)

Raxpy_qd_SOURCES= Raxpy_qd.cpp
Raxpy_qd_CXXFLAGS= $(qdcxxflags)
Raxpy_qd_LDFLAGS= $(qdopt_libdepends)

Raxpy_gmp_ref_SOURCES= Raxpy_gmp.cpp
Raxpy_gmp_ref_CXXFLAGS= $(gmpcxxflags)
Raxpy_gmp_ref_LDFLAGS= $(gmpref_libdepends)

Raxpy_gmp_SOURCES= Raxpy_gmp.cpp
Raxpy_gmp_CXXFLAGS= $(gmpcxxflags)
Raxpy_gmp_LDFLAGS= $(gmpopt_libdepends)

Raxpy_mpfr_ref_SOURCES= Raxpy_mpfr.cpp
Raxpy_mpfr_ref_CXXFLAGS= $(mpfrcxxflags)
Raxpy_mpfr_ref_LDFLAGS= $(mpfrref_libdepends)

Raxpy_mpfr_SOURCES= Raxpy_mpfr.cpp
Raxpy_mpfr_CXXFLAGS= $(mpfrcxxflags)
Raxpy_mpfr_LDFLAGS= $(mpfropt_libdepends)

Raxpy_double_SOURCES= Raxpy_double.cpp
Raxpy_double_CXXFLAGS= $(doublecxxflags)
Raxpy_double_LDFLAGS= $(doubleopt_libdepends)

Raxpy_double_ref_SOURCES= Raxpy_double.cpp
Raxpy_double_ref_CXXFLAGS= $(doublecxxflags)
Raxpy_double_ref_LDFLAGS= $(doubleref_libdepends)

if !IS_MINGW
daxpy_SOURCES= daxpy.cpp
daxpy_CXXFLAGS= $(ref_cxxflags)
daxpy_LDFLAGS= $(ref_libdepends)
endif
