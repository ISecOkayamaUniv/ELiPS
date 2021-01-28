#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([elips], [1.0.14])
#AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/mpn.c])
#		src/bls12/bls12_generate_curve.c])
# AC_CONFIG_SRCDIR([src/bn12/bn12_generate_curve.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])
#Initialize libtool with min ver 2.4
LT_INIT
LT_PREREQ([2.4])
#Avoid warning of README
AM_INIT_AUTOMAKE([foreign])

# Checks for programs.
#AC_PROG_CXX
AC_PROG_CC
#AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_CHECK_SIZEOF(char)
AC_CHECK_SIZEOF(short int)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(long long)

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_LIB(gmp, __gmpz_add)
AC_CHECK_FUNCS([gettimeofday memmove memset strtol])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/bls12/Makefile
                 src/bls/Makefile])
AC_OUTPUT
