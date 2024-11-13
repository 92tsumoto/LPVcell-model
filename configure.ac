#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT([UHNpvc],[1.0.2])
# AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([src/system.c])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

# Checks for programs.
#AC_PROG_CC([icc])
AC_PROG_CC([icx])
AM_PROG_CC_C_O

# Checks for libraries.
# FIXME: Replace `main' with a function in `-llapack':
#AC_CHECK_LIB([lapack],[main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m],[main])
AC_CHECK_LIB([X11],[main])
#AC_CHECK_LIB([pthread],[main])
#AC_CHECK_LIB([imf],[main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stdlib.h string.h mathimf.h mkl.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([floor exp pow sqrt sin cos tan])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT([
])
