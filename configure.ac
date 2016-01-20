# Copyright (C) 2005 Gareth Murphy <garethcmurphy@gmail.com>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

AC_INIT([mhdvanleer],
        [0.0.1],
        [Gareth Murphy garethcmurphy@gmail.com],
        [mhdvanleer])
AM_INIT_AUTOMAKE([dist-bzip2])
#AC_LANG_(C++)

AM_CONFIG_HEADER(config.h)


#LF_CONFIGURE_CC
#LF_CONFIGURE_CXX
##LF_HOST_TYPE
#LF_SET_WARNINGS



AC_PROG_RANLIB
AC_PROG_CXX
AC_PROG_CC
#AM_PROG_LIBTOOL
AC_ISC_POSIX

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_C_CONST


# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_LIB(m,sqrt)
AC_CHECK_LIB(m,pow)

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([string.h unistd.h])

 AC_DEFINE([CXX_HAS_BUGGY_FOR_LOOPS], [], [Description])
 AC_DEFINE([CXX_HAS_NO_BOOL], [], [Description])
AC_DEFINE([NDEBUG], [], [Description])
AC_DEFINE([YOUR_OS], [], [Description])

dnl HDF5
AC_ARG_WITH(hdf5,
   [  --with-hdf5[=DIR]         use HDF5 package (with optional path DIR) ],
   [with_hdf5="$withval" ])

dnl HDF5: default is yes, build with HDF5, using default path
if test "x$with_hdf5" = "x"; then
  with_hdf5=yes
fi

if test "x$with_hdf5" != "xno"; then
  if test "x$with_hdf5" != "xyes"; then
    LIBS="$LIBS -L$with_hdf5/lib -lhdf5"
    INCLUDES="$INCLUDES -I$with_hdf5/include -I$with_hdf5/include/hdf"
    EXT_INCLUDES="-I$with_hdf5/include -I/usr/"
  else
    INCLUDES="$INCLUDES -I/usr/include/hdf"
    LIBS="$LIBS -lhdf5"
  fi

  AC_CHECK_LIB(hdf5, H5Fopen, [AC_DEFINE([USE_HDF5], [1], [Define if you want to use HDF5])], [
    echo ""
    echo "Error! HDF version 5 is required but was not found"
    echo "       Use --with-hdf5=DIR to specify the HDF5 directory tree"
    echo "       Use --with-hdf5=no  to not use it"
    echo "       Check the README or use configure --help for other libraries needed"
    echo "       (--with-xxxdir = obligatory, --with-xxx = optional (--with-xxx=no to disable))"
    exit -1
    ])
fi




AC_CONFIG_FILES([
   Makefile
   README
   doc/Makefile
   m4/Makefile
   src/Makefile
   src/headers/Makefile
   src/tnt/Makefile
])

AC_OUTPUT
