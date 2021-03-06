AC_INIT([sqt], 0.1)
AC_CONFIG_SRCDIR([src/sqt.h])
AC_CONFIG_MACRO_DIR([m4])

SQT_MAJOR_VERSION=0
SQT_MINOR_VERSION=1
SQT_MICRO_VERSION=0
SQT_INTERFACE_AGE=0
SQT_BINARY_AGE=0
SQT_VERSION=$SQT_MAJOR_VERSION.$SQT_MINOR_VERSION.$SQT_MICRO_VERSION

AC_SUBST(SQT_MAJOR_VERSION)
AC_SUBST(SQT_MINOR_VERSION)
AC_SUBST(SQT_MICRO_VERSION)
AC_SUBST(SQT_VERSION)

# libtool versioning
LT_RELEASE=$SQT_MAJOR_VERSION.$SQT_MINOR_VERSION
LT_CURRENT=`expr $SQT_MICRO_VERSION - $SQT_INTERFACE_AGE`
LT_REVISION=$SQT_INTERFACE_AGE
LT_AGE=`expr $SQT_BINARY_AGE - $SQT_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$SQT_VERSION
PACKAGE=sqt

AC_SUBST(PACKAGE)
AC_SUBST(VERSION)

AM_INIT_AUTOMAKE(subdir-objects)

##$PACKAGE, $VERSION, no-define)

# Specify a configuration file
AM_CONFIG_HEADER(config.h)

AC_DEFINE_UNQUOTED(SQT_MAJOR_VERSION, $SQT_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(SQT_MINOR_VERSION, $SQT_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(SQT_MICRO_VERSION, $SQT_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(SQT_INTERFACE_AGE, $SQT_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(SQT_BINARY_AGE, $SQT_BINARY_AGE, [Binary age])

dnl Initialize libtool
AM_PROG_LIBTOOL

dnl Initialize maintainer mode
#AM_MAINTAINER_MODE

AC_PROG_CC
AM_PROG_CC_C_O

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_ISC_POSIX
AC_HEADER_STDC

AM_PATH_GLIB_2_0(2.4.0, has_glib2=yes, has_glib2=no, gthread gmodule)
if test x$has_glib2 = xyes; then
  glib_cflags=`$PKG_CONFIG glib-2.0 --cflags`
  glib_thread_cflags=`$PKG_CONFIG glib-2.0 --cflags gthread-2.0`
  glib_module_cflags=`$PKG_CONFIG glib-2.0 --cflags gmodule-2.0`
  glib_libs=`$PKG_CONFIG glib-2.0 --libs`
  glib_thread_libs=`$PKG_CONFIG glib-2.0 --libs gthread-2.0`
  glib_module_libs=`$PKG_CONFIG glib-2.0 --libs gmodule-2.0`
  GLIB_LIBS="$glib_libs"
  GLIB_DEPLIBS="$glib_libs"
else
  if test x$with_glib = xyes ; then
    AC_MSG_ERROR([
*** Directory must be specified for --with-glib])
  fi

  if test x$with_glib = x ; then 
    # Look for separately installed glib

    GLIB_REQUIRED_VERSION=1.2.8
    AM_PATH_GLIB($GLIB_REQUIRED_VERSION,,
      AC_MSG_ERROR([
*** GLIB $GLIB_REQUIRED_VERSION or better is required. The latest version of GLIB
*** is always available from ftp://ftp.gtk.org/.]),
      gmodule gthread)

    glib_cflags=`$GLIB_CONFIG --cflags`
    glib_thread_cflags=`$GLIB_CONFIG --cflags gthread`
    glib_module_cflags=`$GLIB_CONFIG --cflags gmodule`
    glib_libs=`$GLIB_CONFIG --libs`
    glib_thread_libs=`$GLIB_CONFIG --libs gthread`
    glib_module_libs=`$GLIB_CONFIG --libs gmodule`
    GLIB_LIBS="$glib_libs"
    GLIB_DEPLIBS="$glib_libs"
  else
    # Use uninstalled glib (assume they got the version right)

    GLIB_CONFIG=$with_glib/glib-config
    if test -x $GLIB_CONFIG ; then 
      :
    else
      AC_MSG_ERROR([GLIB directory ($with_glib) not present or not configured])
    fi

    # For use in wmpi-config
    glib_cflags=`$GLIB_CONFIG --cflags`
    glib_thread_cflags=`$GLIB_CONFIG --cflags gthread`
    glib_module_cflags=`$GLIB_CONFIG --cflags gmodule`
    glib_libs=`$GLIB_CONFIG --libs`
    glib_thread_libs=`$GLIB_CONFIG --libs gthread`
    glib_module_libs=`$GLIB_CONFIG --libs gmodule`

    glib_release=`$GLIB_CONFIG --version | sed 's%\\.[[0-9]]*$%%'`

    # canonicalize relative paths
    case $with_glib in 
      /*)
        glib_dir=$with_glib
        ;;
      *)
        glib_dir="\$(top_builddir)/$with_glib"
        ;;
    esac

    GLIB_CFLAGS="-I$glib_dir -I$glib_dir/gmodule"
    GLIB_LIBS="$glib_dir/libglib.la $glib_dir/gmodule.la"
    GLIB_DEPLIBS=

    AC_SUBST(GLIB_CFLAGS)
    AC_SUBST(GLIB_LIBS)
  fi
fi

AX_GCC_X86_CPU_SUPPORTS(avx, hasavx=yes) 
AX_GCC_X86_CPU_SUPPORTS(avx2, hasavx2=yes) 
AX_GCC_X86_CPU_SUPPORTS(fma, hasfma=yes) 

if test x$hasavx = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mavx"
fi
if test x$hasfma = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mfma"
fi

AC_CHECK_LIB(m, cos)
AC_CONFIG_MACRO_DIR([m4])

AC_OPENMP
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS
ACX_BLAS
ACX_LAPACK

blas_libs="$BLAS_LIBS"
lapack_libs="$LAPACK_LIBS"

#DX_DOXYGEN_FEATURE(ON)
##DX_DOT_FEATURE(ON)
#DX_HTML_FEATURE(ON)
##DX_CHM_FEATURE(OFF)
##DX_CHI_FEATURE(OFF)
##DX_MAN_FEATURE(OFF)
##DX_RTF_FEATURE(OFF)
##DX_XML_FEATURE(OFF)
##DX_PDF_FEATURE(OFF)
##DX_PS_FEATURE(OFF)
#DX_INIT_DOXYGEN(sqt, doc/sqt.dxy)

CFLAGS="$CFLAGS $OPENMP_CFLAGS $SIMD_FLAGS $glib_cflags"
LDFLAGS="$LDFLAGS $glib_libs"

AC_SUBST(blas_libs)
AC_SUBST(lapack_libs)
AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
#AC_SUBST(DX_RULES)

LIBS="$glib_libs $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

AC_CONFIG_FILES([
Makefile
src/Makefile
src/sqt-config
tests/Makefile])

#doc/Makefile])

#doc/html/Makefile])

AC_OUTPUT

