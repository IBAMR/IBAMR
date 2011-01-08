/*
 * tests/qd_timer.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2004
 *
 * Contains code to time basic operations.
 */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <qd/qd_real.h>
#include <qd/fpu.h>
#include "tictoc.h"

using std::cout;
using std::cerr;
using std::endl;
using std::sqrt;
using std::strcmp;
using std::setw;
using std::setprecision;
using std::fixed;

// Global flags passed to the main program.
static bool flag_test_double = false;
static bool flag_test_dd = false;
static bool flag_test_qd = false;
static bool flag_verbose = false;
static int  long_factor = 1;

template <class T>
class TestSuite {
public:
  void test1();
  void test2();
  void test3();
  void test4();
  void test5();
  void test6();
  void test7();
  void test8();
  void test9();
  void testall();
  T pi();
};

template <class T>
T TestSuite<T>::pi() { return T::_pi; }

template <>
double TestSuite<double>::pi() { return 3.141592653589793116; }

void print_timing(double nops, double t) {
  double mops = 1.0e-6 * nops / t;
  cout << fixed;
  cout << setprecision(6) << setw(10) << 1.0 / mops << " us";
  cout << setprecision(4) << setw(10) << mops << " mop/s" << endl;
}

template <class T>
void TestSuite<T>::test1() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing addition..." << endl;
  }

  int n = 100000, i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a1 = 1.0 / T(7.0);
  T a2 = 1.0 / T(11.0);
  T a3 = 1.0 / T(13.0);
  T a4 = 1.0 / T(17.0);
  T b1 = 0.0, b2 = 0.0, b3 = 0.0, b4 = 0.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    b1 += a1;
    b2 += a2;
    b3 += a3;
    b4 += a4;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << b1+b2+b3+b4 << endl;
    cout << 4*n << " operations in " << t << " s." << endl;
  } else {
    cout << "   add: ";
  }

  print_timing(4.0*n, t);
}

template <class T>
void TestSuite<T>::test2() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing multiplication ..." << endl;
  }

  int n = 100000, i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a1 = 1.0 + 1.0 / T(static_cast<double>(n));
  T a2 = 1.0 + 2.0 / T(static_cast<double>(n));
  T a3 = 1.0 + 3.0 / T(static_cast<double>(n));
  T a4 = 1.0 + 4.0 / T(static_cast<double>(n));
  T b1 = 1.0, b2 = 1.0, b3 = 1.0, b4 = 1.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    b1 *= a1;
    b2 *= a2;
    b3 *= a3;
    b4 *= a4;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << b1+b2+b3+b4 << endl;
    cout << 4*n << " operations in " << t << " s." << endl;
  } else {
    cout << "   mul: ";
  }

  print_timing(4.0*n, t);
}

template <class T>
void TestSuite<T>::test3() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing division ..." << endl;
  }

  int n = 100000, i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a1 = 1.0 + 1.0 / T(static_cast<double>(n));
  T a2 = 1.0 + 2.0 / T(static_cast<double>(n));
  T a3 = 1.0 + 3.0 / T(static_cast<double>(n));
  T a4 = 1.0 + 4.0 / T(static_cast<double>(n));
  T b1 = 1.0, b2 = 1.0, b3 = 1.0, b4 = 1.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    b1 /= a1;
    b2 /= a2;
    b3 /= a3;
    b4 /= a4;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << b1+b2+b3+b4 << endl;
    cout << 4*n << " operations in " << t << " s." << endl;
  } else {
    cout << "   div: ";
  }

  print_timing(4.0*n, t);
}

template <class T>
void TestSuite<T>::test4() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing square root ..." << endl;
  }

  int n = 10000, i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
  T b1 = 1.0 + pi();
  T b2 = 2.0 + pi();
  T b3 = 3.0 + pi();
  T b4 = 4.0 + pi();

  tic(&tv);
  for (i = 0; i < n; i++) {
    a1 = sqrt(a1 + b1);
    a2 = sqrt(a2 + b2);
    a3 = sqrt(a3 + b3);
    a4 = sqrt(a4 + b4);
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << a1+a2+a3+a4 << endl;
    cout << 4*n << " operations in " << t << " s." << endl;
  } else {
    cout << "  sqrt: ";
  }

  print_timing(4.0*n, t);
}

template <class T>
void TestSuite<T>::test5() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing sin ..." << endl;
  }

  int n = 4000, i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a = 0.0;
  T b = 3.0 * pi() / static_cast<double>(n);
  T c = 0.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    a += b;
    c += sin(a);
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << c << endl;
    cout << n << " operations in " << t << " s." << endl;
  } else {
    cout << "   sin: ";
  }

  print_timing(n, t);
}

template <class T>
void TestSuite<T>::test6() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing log ..." << endl;
  }

  int n = 1000, i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a = 0.0;
  T c = exp(T(-50.1));
  T d = exp(T(100.2) / double(n));

  tic(&tv);
  for (i = 0; i < n; i++) {
    a = a + log(c);
    c *= d;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "a = " << a << endl;
    cout << n << " operations in " << t << " s." << endl;
  } else {
    cout << "   log: ";
  }

  print_timing(n, t);
}

template <class T>
void TestSuite<T>::test7() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing dot ..." << endl;
  }

  int n = 100000, i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a1 = 1.0 / T(7.0);
  T a2 = 1.0 / T(11.0);
  T a3 = 1.0 / T(13.0);
  T a4 = 1.0 / T(17.0);
  T b1 = 1.0 - T(1.0) / static_cast<double>(n);
  T b2 = 1.0 - T(2.0) / static_cast<double>(n);
  T b3 = 1.0 - T(3.0) / static_cast<double>(n);
  T b4 = 1.0 - T(4.0) / static_cast<double>(n);
  T x1 = 1.0, x2 = 1.0, x3 = 1.0, x4 = 1.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    x1 = a1 + b1 * x1;
    x2 = a2 + b2 * x2;
    x3 = a3 + b3 * x3;
    x4 = a4 + b4 * x4;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << x1+x2+x3+x4 << endl;
    cout << 8*n << " operations in " << t << " s." << endl;
  } else {
    cout << "   dot: ";
  }

  print_timing(8.0*n, t);
}

template <class T>
void TestSuite<T>::test8() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing exp ..." << endl;
  }

  int n = 1000, i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a = 0.0;
  T c = -5.0;
  T d = 10.0 / static_cast<double>(n);

  tic(&tv);
  for (i = 0; i < n; i++) {
    a = a + exp(c);
    c += d;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "a = " << a << endl;
    cout << n << " operations in " << t << " s." << endl;
  } else {
    cout << "   exp: ";
  }

  print_timing(n, t);
}

template <class T>
void TestSuite<T>::test9() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing cos ..." << endl;
  }

  int n = 4000, i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a = 0.0;
  T b = 3.0 * pi() / static_cast<double>(n);
  T c = 0.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    a += b;
    c += cos(a);
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << c << endl;
    cout << n << " operations in " << t << " s." << endl;
  } else {
    cout << "   cos: ";
  }

  print_timing(n, t);
}

template <class T>
void TestSuite<T>::testall() {
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
  test8();
  test9();
}

void print_usage() {
  cout << "qd_test [-h] [-dd] [-qd] [-all]" << endl;
  cout << "  Performs timing tests of the quad-double library." << endl;
  cout << "  By default, double-double and quad-double arithmetics" << endl;
  cout << "  are timed." << endl;
  cout << endl;
  cout << "  -h -help  Prints this usage message." << endl;
  cout << "  -double   Time arithmetic of double." << endl;
  cout << "  -dd       Time arithmetic of double-double." << endl;
  cout << "  -qd       Time arithmetic of quad-double." << endl;
  cout << "  -all      Perform both double-double and quad-double tests." << endl;
  cout << "  -v        Verbose output." << endl;
  cout << "  -long     Perform a longer timing loop." << endl;
}

int main(int argc, char *argv[]) {
  unsigned int old_cw;
  fpu_fix_start(&old_cw);

  /* Parse the arguments. */
  char *arg;
  for (int i = 1; i < argc; i++) {
    arg = argv[i];
    if (strcmp(arg, "-h") == 0 || strcmp(arg, "-help") == 0) {
      print_usage();
      std::exit(0);
    } else if (strcmp(arg, "-double") == 0) {
      flag_test_double = true;
    } else if (strcmp(arg, "-dd") == 0) {
      flag_test_dd = true;
    } else if (strcmp(arg, "-qd") == 0) {
      flag_test_qd = true;
    } else if (strcmp(arg, "-all") == 0) {
      flag_test_double = flag_test_dd = flag_test_qd = true;
    } else if (strcmp(arg, "-v") == 0) {
      flag_verbose = true;
    } else if (strcmp(arg, "-long") == 0) {
      long_factor *= 10;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
    }
  }

  /* If no flag, test both double-double and quad-double. */
  if (!flag_test_double && !flag_test_dd && !flag_test_qd) {
    flag_test_dd = true;
    flag_test_qd = true;
  }

  if (flag_test_double) {
    TestSuite<double> test;

    cout << endl;
    cout << "Timing double" << endl;
    cout << "-------------" << endl;
    test.testall();
  }

  if (flag_test_dd) {
    TestSuite<dd_real> test;

    cout << endl;
    cout << "Timing dd_real" << endl;
    cout << "--------------" << endl;
    test.testall();
  }

  if (flag_test_qd) {
    TestSuite<qd_real> test;

    cout << endl;
    cout << "Timing qd_real" << endl;
    cout << "--------------" << endl;
    test.testall();
  }
  
  fpu_fix_end(&old_cw);
  return 0;
}

