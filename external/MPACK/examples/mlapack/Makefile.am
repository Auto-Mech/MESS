lapackexamples_PROGRAMS = \
eigenvalue_dd       \
eigenvalue_qd       \
eigenvalue_gmp      \
eigenvalue_mpfr     \
eigenvalue_double   \
inv_mat_dd       \
inv_mat_qd       \
inv_mat_gmp      \
inv_mat_mpfr     \
inv_mat_double   \
hilbert_dd       \
hilbert_qd       \
hilbert_gmp      \
hilbert_mpfr     \
hilbert_double   \
condition_num_gmp \
condition_num_dd  \
condition_num_qd  \
condition_num_mpfr\
condition_num_double\
condition_num_gmp_pocon

lapackexamplesdir=$(prefix)/share/mpack/examples
ddcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(QD_INCLUDEDIR) 
dd_libdepends =  -L$(top_builddir)/mlapack/optimized/dd -lmlapack_dd -L$(top_builddir)/mblas/optimized/dd -lmblas_dd -L$(top_srcdir)/$(QD_LIBDIR) -lqd

qdcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(QD_INCLUDEDIR) 
qd_libdepends = -L$(top_builddir)/mlapack/optimized/qd -lmlapack_qd -L$(top_builddir)/mblas/optimized/qd -lmblas_qd -L$(top_srcdir)/$(QD_LIBDIR) -lqd

gmpcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/$(GMP_INCLUDEDIR) 
gmp_libdepends = -L$(top_builddir)/mlapack/optimized/gmp -lmlapack_gmp -L$(top_builddir)/mblas/optimized/gmp -lmblas_gmp -L$(top_srcdir)/$(GMP_LIBDIR) -lgmpxx -lgmp

mpfrcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(top_srcdir)/$(GMP_INCLUDEDIR) -I$(top_srcdir)/$(MPFR_INCLUDEDIR) -I$(top_srcdir)/$(MPC_INCLUDEDIR)
mpfr_libdepends =  -L$(top_builddir)/mlapack/optimized/mpfr -lmlapack_mpfr -L$(top_builddir)/mblas/optimized/mpfr -lmblas_mpfr -L$(top_builddir)/mpfrc++ -lmpfrcxx -L$(top_srcdir)/$(MPC_LIBDIR) -L$(top_srcdir)/$(MPFR_LIBDIR) -L$(top_srcdir)/$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

doublecxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -L$(top_builddir)/mlapack/optimized/double -lmlapack_double -L$(top_builddir)/mblas/optimized/double -lmblas_double

eigenvalue_dd_SOURCES= eigenvalue_dd.cpp
eigenvalue_dd_CXXFLAGS= $(ddcxxflags)
eigenvalue_dd_LDFLAGS= $(dd_libdepends)

eigenvalue_qd_SOURCES= eigenvalue_qd.cpp
eigenvalue_qd_CXXFLAGS= $(qdcxxflags)
eigenvalue_qd_LDFLAGS= $(qd_libdepends)

eigenvalue_gmp_SOURCES= eigenvalue_gmp.cpp
eigenvalue_gmp_CXXFLAGS= $(gmpcxxflags)
eigenvalue_gmp_LDFLAGS= $(gmp_libdepends)

eigenvalue_mpfr_SOURCES= eigenvalue_mpfr.cpp
eigenvalue_mpfr_CXXFLAGS= $(mpfrcxxflags)
eigenvalue_mpfr_LDFLAGS= $(mpfr_libdepends)

eigenvalue_double_SOURCES= eigenvalue_double.cpp
eigenvalue_double_CXXFLAGS= $(doublecxxflags)
eigenvalue_double_LDFLAGS= $(double_libdepends)

inv_mat_dd_SOURCES= inv_mat_dd.cpp
inv_mat_dd_CXXFLAGS= $(ddcxxflags)
inv_mat_dd_LDFLAGS= $(dd_libdepends)

inv_mat_qd_SOURCES= inv_mat_qd.cpp
inv_mat_qd_CXXFLAGS= $(qdcxxflags)
inv_mat_qd_LDFLAGS= $(qd_libdepends)

inv_mat_gmp_SOURCES= inv_mat_gmp.cpp
inv_mat_gmp_CXXFLAGS= $(gmpcxxflags)
inv_mat_gmp_LDFLAGS= $(gmp_libdepends)

inv_mat_mpfr_SOURCES= inv_mat_mpfr.cpp
inv_mat_mpfr_CXXFLAGS= $(mpfrcxxflags)
inv_mat_mpfr_LDFLAGS= $(mpfr_libdepends)

inv_mat_double_SOURCES= inv_mat_double.cpp
inv_mat_double_CXXFLAGS= $(doublecxxflags)
inv_mat_double_LDFLAGS= $(double_libdepends)

hilbert_dd_SOURCES= hilbert_dd.cpp
hilbert_dd_CXXFLAGS= $(ddcxxflags)
hilbert_dd_LDFLAGS= $(dd_libdepends)

hilbert_qd_SOURCES= hilbert_qd.cpp
hilbert_qd_CXXFLAGS= $(qdcxxflags)
hilbert_qd_LDFLAGS= $(qd_libdepends)

hilbert_gmp_SOURCES= hilbert_gmp.cpp
hilbert_gmp_CXXFLAGS= $(gmpcxxflags)
hilbert_gmp_LDFLAGS= $(gmp_libdepends)

hilbert_mpfr_SOURCES= hilbert_mpfr.cpp
hilbert_mpfr_CXXFLAGS= $(mpfrcxxflags)
hilbert_mpfr_LDFLAGS= $(mpfr_libdepends)

hilbert_double_SOURCES= hilbert_double.cpp
hilbert_double_CXXFLAGS= $(doublecxxflags)
hilbert_double_LDFLAGS= $(double_libdepends)

condition_num_dd_SOURCES= condition_num_dd.cpp
condition_num_dd_CXXFLAGS= $(ddcxxflags)
condition_num_dd_LDFLAGS= $(dd_libdepends)

condition_num_qd_SOURCES= condition_num_qd.cpp
condition_num_qd_CXXFLAGS= $(qdcxxflags)
condition_num_qd_LDFLAGS= $(qd_libdepends)

condition_num_gmp_SOURCES= condition_num_gmp.cpp
condition_num_gmp_CXXFLAGS= $(gmpcxxflags)
condition_num_gmp_LDFLAGS= $(gmp_libdepends)

condition_num_mpfr_SOURCES= condition_num_mpfr.cpp
condition_num_mpfr_CXXFLAGS= $(mpfrcxxflags)
condition_num_mpfr_LDFLAGS= $(mpfr_libdepends)

condition_num_double_SOURCES= condition_num_double.cpp
condition_num_double_CXXFLAGS= $(doublecxxflags)
condition_num_double_LDFLAGS= $(double_libdepends)

condition_num_gmp_pocon_SOURCES= condition_num_gmp_pocon.cpp
condition_num_gmp_pocon_CXXFLAGS= $(gmpcxxflags)
condition_num_gmp_pocon_LDFLAGS= $(gmp_libdepends)

