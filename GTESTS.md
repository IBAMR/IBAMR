
#IBAMR tests#
## Dependencies #
To run IBAMR's automated tests using `make gtest` or `make gtest-long` the [Google Test Framework](https://github.com/google/googletest) must be installed using the same compiler and compiler flags as the IBAMR project. 

><font size=3>**Note:** The Google Test Framework is often referred to as **gtest**

Discussion of why this is necessary can be found in [Google Test README](googletest/googletest/README.md).

In short, it is necessary to use the same compiler and same compiler flags as you will be using for your IBAMR build to avoid possibly aberrant runtime behavior.

Below is an example of how you would build gtest as a user on an Linux system with all IBAMR dependencies already in place as outlined in the [IBAMR building wiki](https://github.com/IBAMR/IBAMR/wiki/Building). 

This build technique requires `cmake` which is available through most package managers. Additionally users can find pre-compiled binaries and instructions for building from source at [Cmake project's install page](https://cmake.org/install/). 

*If you would like to build googletest using only the autotools, please refer to the appendix at the end of this document.*

```
cd $HOME/sfw/linux/
git clone https://github.com/google/googletest.git googletest
cd googletest

cmake CMakeLists.txt \
-DBUILD_GMOCK=ON \
-DBUILD_GTEST=OFF \
-DCMAKE_CXX_COMPILER=$HOME/sfw/linux/openmpi/1.10.2/bin/mpicxx \
-DCMAKE_CXX_FLAGS=-Wall \
-DCMAKE_C_COMPILER=$HOME/sfw/linux/openmpi/1.10.2/bin/mpicc \
-DCMAKE_C_FLAGS=-Wall \
-DCMAKE_INSTALL_PREFIX=$HOME/sfw/linux/googletest/ \
-Dgtest_build_tests=ON \
-Dgmock_build_tests=ON \
-Dgtest_build_samples=ON

make 
make test
make install
```

> You should set `DC_MAKE_CXX_COMPILER`, `DCMAKE_C_COMPILER`, `DCMAKE_C_FLAGS` and `DCMAKE_CXX_FLAGS` to the same values as you have build IBAMR with. You can see your most recent options by looking looking at the `config.log` in your IBAMR build directory.

Once the gtest library is installed, IBAMR must be configured with the LDFLAGS variable pointing to it. i.e. add 

>```--enable-gtest --with-gtest=$HOME/sfw/linux/googletest ```

to your `configure` invocation.


## How to run the tests ##

Once the gtest library is installed and IBAMR is configured properly, all the tests can be run by invoking `make gtest` or `make gtest-long` in the root build directory. 

#### `make gtest` vs. `make gtest-all` ####
* `make gtest` is intended as a smoke test/sanity check that can be run in less than ten minutes on an average machine used for development (assuming the tests are already compiled). 

* `make gtest` compiles a selection of 2D and 3D tests, and of these selected only runs the 2D tests.

*  `make gtest-all` compiles and runs every single test (2D and 3D) and takes a significant amount of time (exact runtime depends on your system resources).

## Running Individual Tests ##

Each test currently included in the library is based off of example applications and uses the same source code as the examples, but runs the simulation from within a gtest application and analyzes the results. To run an individual test, navigate to the example of interest and run `make gtest` from the command line to compile the test. Then the application can be run with `./test2d input2d.test` or `./test3d input3d.test` from the command line .

## Interpreting Results ##

Every gtest application included in the IBAMR library returns a result indicating whether or not the example runs.

Below is output from an individual test application:

```
[==========] Running 7 tests from 1 test case.
[----------] Global test environment set-up.
[----------] 7 tests from IBFE_explicit_ex0_2d
[ RUN      ] IBFE_explicit_ex0_2d.example_runs
[       OK ] IBFE_explicit_ex0_2d.example_runs (0 ms)
[ RUN      ] IBFE_explicit_ex0_2d.u_L1Norm
[       OK ] IBFE_explicit_ex0_2d.u_L1Norm (0 ms)
[ RUN      ] IBFE_explicit_ex0_2d.u_L2Norm
[       OK ] IBFE_explicit_ex0_2d.u_L2Norm (0 ms)
[ RUN      ] IBFE_explicit_ex0_2d.u_MaxNorm
[       OK ] IBFE_explicit_ex0_2d.u_MaxNorm (0 ms)
[ RUN      ] IBFE_explicit_ex0_2d.p_L1Norm
[       OK ] IBFE_explicit_ex0_2d.p_L1Norm (0 ms)
[ RUN      ] IBFE_explicit_ex0_2d.p_L2Norm
[       OK ] IBFE_explicit_ex0_2d.p_L2Norm (0 ms)
[ RUN      ] IBFE_explicit_ex0_2d.p_MaxNorm
[       OK ] IBFE_explicit_ex0_2d.p_MaxNorm (0 ms)
[----------] 7 tests from IBFE_ex0_2d (0 ms total)

[----------] Global test environment tear-down
[==========] 7 tests from 1 test case ran. (1 ms total)
[  PASSED  ] 7 tests.
```

Notice that each test name follows the convention:
`NAME_OF_TEST_APPLICATION.NAME_OF_TEST`

**Example interpretation 1:**
```
[ RUN      ] IBFE_explicit_ex0_2d.example_runs
[       OK ] IBFE_explicit_ex0_2d.example_runs (0 ms)
```
The first part of the test name indicates that the test was of the application found in examples/IBFE/explicit/ex0.
The second part of the name, "example_runs", indicates that the test is simply reporting back if the example was able to execute without any runtime errors.

**Example interpretation 2:**
```
[ RUN      ] IBFE_explicit_ex0_2d.u_MaxNorm
[       OK ] IBFE_explicit_ex0_2d.u_MaxNorm (0 ms)
[ RUN      ] IBFE_explicit_ex0_2d.p_L1Norm
[       OK ] IBFE_explicit_ex0_2d.p_L1Norm (0 ms)
```

Again, the first part of the test names indicates that the tests were of the application found in examples/IBFE/explicit/ex0.

Here the test `u_MaxNorm` is observing whether or not the test application's "Max Norm" (also known as the infinity norm) of the velocity is within an acceptable margin of error from an analytical or laboratory result. 

The second test is performing a similar calculation except this time with the L1 norm of the pressure.

#### _Notation_ ####

**L1_Norm**: The one-norm (also known as the L1-norm, `1 norm, or mean norm)  is defined as the sum of the absolute values of its components.

**L2_Norm**:  The two-norm (also known as the L2-norm, 2-norm, mean-square norm, or least-squares norm) is defined as the square root of the sum of the squares of the absolute values of its components.

**MaxNorm**:  The infinity norm (also known as the L∞-norm, `∞-norm, max norm, or uniform norm) is defined as the maximum of the absolute values of its
components.

## Adding tests to existing applications ##

Adding new tests to examples with existing gtest applications, for example to `$IBAMR_DIR/examples/IBFE/explicit/ex1/`, can be done by adding tests to the `test_main.cpp` file in that directory following the pattern:

```
TEST(TEST_CASE_NAME, your_unique_test_name_here) {
	// do stuff
	// make some assertions
}
```

Many types  of assertions are available and are well documented in the [Google Test Primer](https://github.com/google/googletest/blob/master/googletest/docs/Primer.md) and the [Google Test Advanced Guide](https://github.com/google/googletest/blob/master/googletest/docs/AdvancedGuide.md).

If you desire to collect additional data from the example, it is suggested you pass a reference to the object that you'd like to test to the example. An example of this can be found in:
`$IBAMR_DIR/examples/IBFE/explicit/ex0/test_main.cpp `.

Any setup that should be shared between all tests can be performed in the `main( )` method in test_main.cpp before the call to `return RUN_ALL_TESTS(  );`.

```cpp
int main( int argc, char** argv ) {
     /* these lines processes any command line arguments to googletest 
      * and then prepare to hand the remaining arguments that were not
      * valid googletest arguments to the example 
      */
    testing::InitGoogleTest( &argc, argv );
    ex_argc = argc;
    ex_argv = argv; 

	// this line runs the simulation returns true if no run-time errors occur
    ex_runs = run_example(ex_argc, ex_argv);
    
    /* All tests have access to anything declared or manipulated here, 
     * as the tests themselves are actually run within this main method
     * upon the invocation of RUN_ALL_TESTS( )
     */
    
    return RUN_ALL_TESTS( );
}
```
## Creating new test applications ##

The essential elements to setting up your new google test application are:

```
// test_main.cpp
#include <gtest/gtest.h>
//#include "your class/header etc. here"

//function prototypes, variable declarations
// example:
static const double REL_ERROR = 1.0e-8;

// Set names of test based on if compiled with 2D or 3D libraries
#if (NDIM == 2)
#define TEST_CASE_NAME name_of_class_2d
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME name_of_class_3d
#endif


int main( int argc, char** argv ) {
    testing::InitGoogleTest( &argc, argv );

    // any shared objects or data that will be refered to 
    // within individual test cases
    // can be created or calculated here
   
    return RUN_ALL_TESTS( );
}

```

If you would like to contribute your tests to IBAMR, please ensure that functionality of any other make targets do not require the user to have gtest enabled, i.e. all applications that rely on google test should be made using the `make gtest` target.
```
## Makefile.am
## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

GTEST_DRIVER   = test_main.cpp

GTESTS   =
EXTRA_PROGRAMS =
if LIBMESH_ENABLED
if SAMRAI2D_ENABLED
if GTEST_ENABLED
GTESTS   += test2d
endif
EXTRA_PROGRAMS += $(GTESTS)
endif
endif

test2d_CXXFLAGS = $(main2d_CXXFLAGS)
test2d_LDADD = $(main2d_LDADD)
test2d_SOURCES = $(GTEST_DRIVER)

if GTEST_ENABLED
gtest: $(GTESTS)

else
gtest:
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	@echo " To run tests, Google Test Framework must be enabled.                "
	@echo " Configure IBAMR with additional options:                            "
	@echo "                      --enable-gtest --with-gtest=path               "
	@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
endif
```

## Output ##
By default, the gtest applications will only output to stdout.

The **`$GTEST_OUTPUT`** environment variable is available to make the gtest applications generate xml reports of test results. Its behavior is detailed in [the advanced guide to gtest](https://github.com/google/googletest/blob/48ee8e98abc950abd8541e15550b18f8f6cfb3a9/googletest/docs/V1_7_AdvancedGuide.md#generating-an-xml-report). 

In short,

> To generate the XML report, set the `GTEST_OUTPUT` environment variable or the `--gtest_output` flag to the string `"xml:_path_to_output_file_"`, which will create the file at the given location. You can also just use the string `"xml"`, in which case the output can be found in the `test_detail.xml` file in the current directory.

This is especially useful when using the [JUnit plugin for Jenkins](https://wiki.jenkins-ci.org/display/JENKINS/JUnit+Plugin) which knows how to parse these xml files. 

The following is a build script for a job running tests and reporting the results using the [JUnit plugin](https://wiki.jenkins-ci.org/display/JENKINS/JUnit+Plugin):
```bash
export HOME=/srv/sfw
export BOOST_ROOT=$HOME/linux/boost/1.61.0
export PETSC_ARCH=linux-debug
export PETSC_DIR=$HOME/petsc/3.7.2
export GTEST_OUTPUT="xml:$WORKSPACE/"
export SAMRAI_DIR="$HOME/samrai/2.4.4/linux-debug"
rm -rf *.xml
./.jenkins_quicktest

```

With this set, all test results will be amalgamated into one report the the plugin can use to determine the build status and generate graphical reports in the job view. Old test results must be removed at the beginning of each build in order to avoid cluttering up the workspace.

> _**Note:**_`$WORKSPACE` *is one of the [environment variable is available from Jenkins](https://wiki.jenkins-ci.org/display/JENKINS/Building+a+software+project#Buildingasoftwareproject-JenkinsSetEnvironmentVariables) and resolves to the absolute path to the workspace*

## APPENDIX: How to build libgtest without cmake or libtool ##
This method is not recommended, but possible.

Building libgtest with cmake method is preferable, as it is the currently maintained build system for the Google Test Framework. However, there are are autotools scripts available in the interior googletest directory. 

If you do not have `cmake` available and do not wish to install it, here is an example of how you might build libgtest using the autotools:

First, download the googletest repository from git:
```
cd $HOME/sfw/linux/
git clone https://github.com/google/googletest.git googletest
cd googletest
```

IBAMR configure is expecting a static library (a library with ".a" file ending), but the googletest Makefile.am creates ["libtool libraries"](http://stackoverflow.com/questions/1238035/what-are-libtools-la-file-for). 

So you must go into the Makefile.am and change everything that says **something like**:
```
lib_LTLIBRARIES = lib/libgtest.la lib/libgtest_main.la

lib_libgtest_la_SOURCES = src/gtest-all.cc
```
**to** 
```
lib_LIBRARIES = lib/libgtest.a lib/libgtest_main.a

lib_libgtest_a_SOURCES = src/gtest-all.cc

```

This example is representative but not exhaustive. There are several more occurrences where you must replace `LTLIBRARIES` with `LIBRARIES` and switch `la` to `a`.

```
cd $HOME/sfw/linux
cd $HOME/sfw/linux/
git clone https://github.com/google/googletest.git googletest
cd googletest/googletest

autoreconf --force --verbose --install -I config -I m4
./configure CXXFLAGS=-Wall FCFLAGS=-Wall CC="ccache $HOME/linux/openmpi/1.10.2/bin/mpicc" CXX="ccache $HOME/linux/openmpi/1.10.2/bin/mpicxx" FC=$HOME/linux/openmpi/1.10.2/bin/mpif90 CPPFLAGS=-DOMPI_SKIP_MPICXX

make
``` 

Now there should be a header at the path: 

`$HOME/sfw/linux/googletest/googletest/include/gtest/gtest.h`
 
 and a static library at the path:
 
 `$HOME/sfw/linux/googletest/googletest/lib/libgtest.a`.

Now, when configuring with IBAMR you will provide this inner directory to configure:
 `--enable-gtest --with-gtest=$HOME/sfw/linux/googletest/googletest`


