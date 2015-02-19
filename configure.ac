#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([dispersion], [0.01])
AM_INIT_AUTOMAKE([-Wall foreign])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Checks for programs.

AC_PROG_CC_C99
if test -z "$CC" || test "$CC" = "no"; then
  AC_MSG_ERROR([No suitable C99 compiler found])
fi
AM_PROG_CC_C_O

# add -Wall if possible
AX_CHECK_COMPILER_FLAGS([-Wall], [CFLAGS="$CFLAGS -Wall"])

dnl For Fortran support, uncomment (remove "dnl") the following
AC_PROG_FC
AC_PROG_F77
if test -z "$F77" || test "X$F77" = "Xno"; then
  AC_MSG_ERROR([No suitable F77 compiler found])
fi

# for library support
#AC_PROG_RANLIB
#m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_F77_WRAPPERS
AC_FC_WRAPPERS

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
	src/Makefile])
AC_OUTPUT

echo ""
echo "                       Configure Summary"
echo "                       ================="
echo ""
echo "Compiling Options:"
echo "------------------"
echo "                     C Compiler: $CC"
echo "                       CPPFLAGS: $CPPFLAGS"
echo "                         CFLAGS: $CFLAGS"
echo "                        LDFLAGS: $LDFLAGS"
echo "                           LIBS: $LIBS"
