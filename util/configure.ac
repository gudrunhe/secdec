AC_INIT([SecDecUtil], [0.3])
AC_CONFIG_SRCDIR([secdecutil/series.hpp])
AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([acinclude.d]) # define AX_CXX_COMPILE_STDCXX_11 macro (sometimes required for mac)

AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext]) # noext: use -std=c++11 rather than -std=gnu++11

AM_CONDITIONAL([HAVE_SECDEC_CONTRIB], [test x$SECDEC_CONTRIB != x])

AC_OUTPUT(Makefile secdecutil/Makefile tests/Makefile)
