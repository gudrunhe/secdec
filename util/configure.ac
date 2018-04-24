AC_INIT([SecDecUtil], [1.3.1])
AC_CONFIG_SRCDIR([secdecutil/series.hpp])
AM_INIT_AUTOMAKE

dnl define AX_CXX_COMPILE_STDCXX_11 macro (sometimes required for mac)
AC_CONFIG_MACRO_DIR([acinclude.d])

AC_PROG_CXX

dnl noext: use -std=c++11 rather than -std=gnu++11
AX_CXX_COMPILE_STDCXX_11([noext])

AM_CONDITIONAL([HAVE_SECDEC_CONTRIB], [test x$SECDEC_CONTRIB != x])

AC_OUTPUT(Makefile secdecutil/Makefile tests/Makefile)
