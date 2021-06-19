AC_INIT([SecDecUtil], [1.4.5])
AC_CONFIG_SRCDIR([secdecutil/series.hpp])
AM_INIT_AUTOMAKE

dnl check whether to enable cuda support
AC_ARG_ENABLE([cuda],
    [AS_HELP_STRING([--enable-cuda], [Compile and run test cases with cuda. Assumes a cuda-capable compiler.])],
    [AM_CONDITIONAL([SECDEC_WITH_CUDA], [test x$enable_cuda != xno])],
    [
        AS_IF(
            [test "x${CXX}" = "xnvcc"],
            [AM_CONDITIONAL([SECDEC_WITH_CUDA], true)],
            [AM_CONDITIONAL([SECDEC_WITH_CUDA], false)]
        )
    ]
)

dnl define AX_CXX_COMPILE_STDCXX_14 macro (sometimes required for mac)
AC_CONFIG_MACRO_DIR([acinclude.d])

AC_PROG_CXX


AM_CONDITIONAL([HAVE_SECDEC_CONTRIB], [test x$SECDEC_CONTRIB != x])
dnl noext: use -std=c++14 rather than -std=gnu++14
AX_CXX_COMPILE_STDCXX_14([noext])

AC_OUTPUT(Makefile secdecutil/Makefile tests/Makefile)
