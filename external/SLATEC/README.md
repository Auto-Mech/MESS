# SLATEC (subset)

These are the routines from the SLATEC Common Mathematical Library, Version 4.1 (July 1993), that MESS calls (`davint`, `dbint4`, `dbvalu`, `ddeabm`) and the routines they call, 24 files in all.
The top-level `CMakeLists.txt` builds them into a static library linked into the MESS executables.
They use the default (32-bit) Fortran `INTEGER`, which matches `Slatec::int_t` in `src/libmess/slatec.hh`.

## Notice

The files come from Giuseppe Borzi's SLATEC distribution:

Giuseppe Borzi'  
e-mail: etana@tiscalinet.it  
http://web.tiscalinet.it/gborzi  
Assistant Professor at the Univ. of Messina - Italia

License:  
General Public License. Slatec is 'Freely distributable' as they say.

For further details, see <http://www.netlib.org/slatec/>.

The subset deliberately leaves out SLATEC's copies of BLAS routines (`ddot`, `daxpy`, …), which none of these call.
They would clash with the ILP64 BLAS in a static link.
