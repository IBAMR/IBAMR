# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_GTEST],[
echo
echo "=================================================="
echo "Configuring optional package Google Test Framework"
echo "=================================================="

PACKAGE_SETUP_ENVIRONMENT

AC_ARG_ENABLE([gtest],
  AS_HELP_STRING(--enable-gtest,enable support for the optional GTEST library @<:@default=no@:>@),
                 [case "$enableval" in
                    yes)  GTEST_ENABLED=yes ;;
                    no)   GTEST_ENABLED=no ;;
                    *)    AC_MSG_ERROR(--enable-gtest=$enableval is invalid; choices are "yes" and "no") ;;
                  esac],[GTEST_ENABLED=no])

AM_CONDITIONAL([GTEST_ENABLED],[test "$GTEST_ENABLED" = yes])

AC_ARG_WITH([gtest],
  AS_HELP_STRING(--with-gtest=PATH,location of optional GTEST installation),
  [if test "$GTEST_ENABLED" = no ; then
     AC_MSG_WARN([--with-gtest is specified, but support for gtest is disabled])
   else
     if test ! -d "$withval" ; then
       AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-gtest=PATH])
     fi
     GTEST_DIR=$withval
   fi])

if test "$GTEST_ENABLED" = yes; then
  if test x$GTEST_DIR != x ; then
    if test -d "${GTEST_DIR}/include" ; then
      GTEST_CPPFLAGS="-I${GTEST_DIR}/include"
    fi
    if test -d "${GTEST_DIR}/lib" ; then
      GTEST_LDFLAGS="-L${GTEST_DIR}/lib"
    fi
  fi
  if test x$GTEST_DIR = x ; then
      AC_MSG_ERROR([

      !!!! Google Test is enabled, but no path specified !!!!
      
      --enable-gtest was passed, but GTEST_DIR was not specified 
      with --with-gtest=GTEST_DIR
  
  Google Test may be built using Cmake in following manner:
  =============================================================
  git clone https://github.com/google/googletest.git googletest
  cd googletest

  cmake CMakeLists.txt  -DCMAKE_INSTALL_PREFIX=GTEST_DIR
    -DCMAKE_CXX_COMPILER=MPI_DIR/bin/mpicxx
    -DCMAKE_CXX_FLAGS=IBAMR_CXX_FLAGS 
    -DCMAKE_C_COMPILER=MPI_DIR/bin/mpicc 
    -DCMAKE_C_FLAGS=IBAMR_C_FLAGS -Dgtest_build_tests=ON 
    -Dgmock_build_tests=ON -Dgtest_build_samples=ON
    -DBUILD_GMOCK=ON -DBUILD_GTEST=OFF

    make 
    make test
    make install
  ===============================================================
  
  (replace GTEST_DIR, IBAMR_CXX_FLAGS, etc. with actual values)

  See IBAMR_TESTING_GUIDE.md in root IBAMR directory for more details.
  ])
 fi

  CPPFLAGS_PREPEND($GTEST_CPPFLAGS)
  AC_CHECK_HEADER([${GTEST_DIR}/include/gtest/gtest.h],,AC_MSG_ERROR([
  
  !!!! Configure could not find header file gtest.h !!!!
  Please provide to --with-gtest=GTEST_DIR the directory that contains:
  
  GTEST_DIR
  |-include/
    |-gtest/
      +-gtest.h 
  |-lib/
    +-libgtest.a
  
  Google Test may be built using Cmake in following manner:
  =============================================================
  git clone https://github.com/google/googletest.git googletest
  cd googletest

  cmake CMakeLists.txt  -DCMAKE_INSTALL_PREFIX=GTEST_DIR
    -DCMAKE_CXX_COMPILER=MPI_DIR/bin/mpicxx
    -DCMAKE_CXX_FLAGS=IBAMR_CXX_FLAGS 
    -DCMAKE_C_COMPILER=MPI_DIR/bin/mpicc 
    -DCMAKE_C_FLAGS=IBAMR_C_FLAGS -Dgtest_build_tests=ON 
    -Dgmock_build_tests=ON -Dgtest_build_samples=ON
    -DBUILD_GMOCK=ON -DBUILD_GTEST=OFF

    make 
    make test
    make install
  ===============================================================
  
  See IBAMR_TESTING_GUIDE.md in root IBAMR directory for more details.
  ]))

  LDFLAGS_PREPEND($GTEST_LDFLAGS)
  AC_LIB_HAVE_LINKFLAGS([gtest])
  if test "$HAVE_LIBGTEST" = no ; then
    AC_MSG_ERROR([
    !!! Could not find working libgtest.a !!! 
  Please provide to --with-gtest=GTEST_DIR the directory that contains:
  
  GTEST_DIR
  |-include/
    |-gtest/
      +-gtest.h 
  |-lib/
    +-libgtest.a
  
  Google Test may be built using Cmake in following manner:
  =============================================================
  git clone https://github.com/google/googletest.git googletest
  cd googletest

  cmake CMakeLists.txt  -DCMAKE_INSTALL_PREFIX=GTEST_DIR
    -DCMAKE_CXX_COMPILER=MPI_DIR/bin/mpicxx
    -DCMAKE_CXX_FLAGS=IBAMR_CXX_FLAGS 
    -DCMAKE_C_COMPILER=MPI_DIR/bin/mpicc 
    -DCMAKE_C_FLAGS=IBAMR_C_FLAGS -Dgtest_build_tests=ON 
    -Dgmock_build_tests=ON -Dgtest_build_samples=ON
    -DBUILD_GMOCK=ON -DBUILD_GTEST=OFF

    make 
    make test
    make install
  ===============================================================
  
  See IBAMR_TESTING_GUIDE.md in root IBAMR directory for more details.
    ])
  fi

  PACKAGE_CPPFLAGS_PREPEND($GTEST_CPPFLAGS)
  PACKAGE_LDFLAGS_PREPEND($GTEST_LDFLAGS)
  PACKAGE_LIBS_PREPEND("$LIBGTEST")
 
else
  AC_MSG_NOTICE([Optional package Google Test Framework is DISABLED])
fi

PACKAGE_RESTORE_ENVIRONMENT

])
