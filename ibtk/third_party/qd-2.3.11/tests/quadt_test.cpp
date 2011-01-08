/*
 * tests/quadt_test.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This class contains a test suite for the quadt integration 
 * code (see quadt.h).
 */

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <qd/fpu.h>
#include <qd/qd_real.h>

#include "tictoc.h"

using std::cout;
using std::cerr;
using std::endl;
using std::abs;
using std::exp;
using std::log;
using std::sqrt;
using std::cos;
using std::atan;
using namespace qd;

#include "quadt.h"

/** Various flags passed to the main program. */
static bool flag_verbose   = false;
static bool flag_test_d    = false;
static bool flag_test_dd   = false;
static bool flag_test_qd   = false;
static bool flag_last_only = false;

template <class T> 
class constants {
public:
  static const T &pi;
  static const T &pi2;
  static const T &pi4;
  static const T &log2;
};

template <class T>
const T &constants<T>::pi   = T::_pi;
template <class T>
const T &constants<T>::pi2  = T::_pi2;
template <class T>
const T &constants<T>::pi4  = T::_pi4;
template <class T>
const T &constants<T>::log2 = T::_log2;

template <>
const double &constants<double>::pi = 3.14159265358979;
template <>
const double &constants<double>::pi2 = 1.57079632679490;
template <>
const double &constants<double>::pi4 = 0.785398163397448;
template <>
const double &constants<double>::log2 = 0.693147180559945;


/* Sample Functions */
template <class T>
class CircleFunction {
  T r, r2;
public:
  CircleFunction(T radius) {
    r = radius;
    r2 = sqr(r);
  }
  T operator() (T x) const {
    if (abs(x) >= r)
      return 0;
    return sqrt(r2 - sqr(x));
  }
};

template <class T>
class SecantFunction {
public:
  T operator() (T x) const {
    return 1.0 / cos(x);
  }
};

template <class T>
class TestFunction1 {
public:
  T operator() (T x) const {
    return x * log(1.0 + x);
  }
};

template <class T>
class TestFunction2 {
public:
  T operator() (T x) const {
    return sqr(x) * atan(x);
  }
};

template <class T>
class TestFunction3 {
public:
  T operator() (T x) const {
    if (x <= 0.0)
      return 0.0;
    return sqr(log(x));
  }
};

template <class T>
class TestFunction4 {
public:
  T operator() (T x) const {
    T tmp;
    if (x >= constants<T>::pi2)
      return 0.0;
    tmp = tan(x);
    if (tmp < 0.0)
      return 0.0;
    return sqrt(tmp);
  }
};

template <class T>
class TestFunction5 {
public:
  T operator() (T x) const {
    T t;

    if (x <= 0.0)
      return 0.0;

    if (x > 0.00146) {
      T rt = 1.0 / x;
      t = 1.0 / (exp(rt) * sqrt(rt) * sqr(x));
    } else
      t = 0.0;

    t = 1.0 / (exp(x) * sqrt(x)) + t;

    return t;
  }
};

template <class T>
void convert(char *s, T *x) {
  *x = s;
}

template <>
void convert(char *s, double *x) {
  *x = atof(s);
}

template <class T>
class CosineProduct {
private:
  T coeff[32];
public:
  CosineProduct() {
    const char *f_name = "coeff.dat";
    char s[100];

    FILE *f = fopen(f_name, "r");
    if (f == NULL) {
      cerr << "Failed to open coefficient file " << f_name << "." << endl;
      exit(-1);
    }

    for (int i = 0; i < 32; i++) {
      fscanf(f, "%s", s);
      convert(s, &coeff[i]);
    }
    
    fclose(f);
  }

  T operator() (T x) const {
    T xx = (x + 1.0) * 0.5;
    T val = cos(2.0 * xx);
    T tmp = 0.0;
    T xp = sqr(xx);
    T x2 = xp;
    for (int i = 1; i < 512; i++) {
      val *= cos(xx / static_cast<double>(i));
    }
    for (int i = 0; i < 32; i++, xp *= x2) {
      tmp += coeff[i] * xp;
    }
    val *= exp(tmp);

    return val;
  }
};

template <class T>
class InvCosineProduct {
private:
  T coeff[32];
public:
  InvCosineProduct() {
    const char *f_name = "coeff.dat";
    char s[100];

    FILE *f = fopen(f_name, "r");
    if (f == NULL) {
      cerr << "Failed to open coefficient file " << f_name << "." << endl;
      exit(-1);
    }

    for (int i = 0; i < 32; i++) {
      fscanf(f, "%s", s);
      convert(s, &coeff[i]);
    }
    
    fclose(f);
  }

  T operator() (T x) const {
    T xx = (x + 1.0) * 0.5;

    if (xx < 0.005)
      return 0.0;

    T inv_x = 1.0 / xx;
    T val = cos(2.0 * inv_x);
    T tmp = 0.0;
    T xp = sqr(inv_x);
    T x2 = xp;
    for (int i = 1; i < 512; i++) {
      val *= cos(inv_x / static_cast<double>(i));
    }
    for (int i = 0; i < 32; i++, xp *= x2) {
      tmp += coeff[i] * xp;
    }
    val *= exp(tmp);
    val *= x2;

    return val;
  }
};

template <class T>
class quadt_tester {
private:
  double eps;
  quadt<T> *q;
public:
  quadt_tester(double eps) { 
    this->eps = eps; 
    q = new quadt<T>(eps);
  }

  ~quadt_tester() {
    delete q;
  }

  template <class F>
  void test_integral(F &f, T a, T b, T truth);

  void test();
};

template <class T> template <class F>
void quadt_tester<T>::test_integral(F &f, T a, T b, T truth) {
  int r;
  T result;
  double err_est, err;
  double tol = eps * 1024;
  
  r = q->integrate(f, a, b, tol, result, err_est);
  err = abs(to_double(result - truth));
  if (flag_verbose) {
    cout << "    Result: " << result << endl;
    cout << "     Truth: " << truth << endl;
    cout << "Est. Error: " << err_est << endl;
    cout << "True Error: " << err << endl;
    cout << endl;
  }
}

template <class T>
void quadt_tester<T>::test() {
  CosineProduct<T> f8;
  InvCosineProduct<T> invf8;

  if (!flag_last_only) {
    CircleFunction<T> f1(1.0);
    SecantFunction<T> f2;
    TestFunction1<T> f3;
    TestFunction2<T> f4;
    TestFunction3<T> f5;
    TestFunction4<T> f6;
    TestFunction5<T> f7;
    
    cout << "Test 1." << endl;
    test_integral(f1, T(-1.0), T(1.0), constants<T>::pi2);
    
    cout << "Test 2." << endl;
    test_integral(f2, T(0.0), constants<T>::pi4, log(1.0 + sqrt(T(2.0))));
    
    cout << "Test 3." << endl;
    test_integral(f3, T(0.0), T(1.0), T(0.25));
    
    cout << "Test 4." << endl;
    test_integral(f4, T(0.0), T(1.0), 
                  constants<T>::pi4/3.0 - T(1.0) / 6.0 + log(T(2.0)) / 6.0);
    
    cout << "Test 5." << endl;
    test_integral(f5, T(0.0), T(1.0), T(2.0));
    
    cout << "Test 6." << endl;
    test_integral(f6, T(0.0), constants<T>::pi2,
                  constants<T>::pi2 * sqrt(T(2.0)));
    
    cout << "Test 7." << endl;
    test_integral(f7, T(0.0), T(1.0), sqrt(constants<T>::pi));
  }

  cout << "Test 8." << endl;
    
  double tol = eps * 1024;
  double err;
  T r, r1, r2;
  T truth = constants<T>::pi4 * 0.5;

  q->integrate_u(f8, tol, r1, err);

  r1 *= 0.5;
  if (flag_verbose) {
    cout << "  Result 1: " << r1 << endl;
    cout << "Est. Error: " << err << endl;
  }

  q->integrate_u(invf8, tol, r2, err);

  r2 *= 0.5;
  if (flag_verbose) {
    cout << "  Result 2: " << r2 << endl;
    cout << "Est. Error: " << err << endl;
  }

  r = r1 + r2;
  err = abs(to_double(r - truth));

  if (flag_verbose) {
    cout << "    Result: " << r << endl;
    cout << "     Truth: " << truth << endl;
    cout << "True Error: " << err << endl;
    cout << endl;  
  }

}

template <class T>
void test_quadt(double eps) {

  tictoc tv;
  double tm1, tm2;

  tic(&tv);
  quadt_tester<T> tester (eps);
  tm1 = toc(&tv);

  tic(&tv);
  tester.test();
  tm2 = toc(&tv);

  cout << "Setup CPU Time = " << tm1 << endl;
  cout << " Test CPU Time = " << tm2 << endl;
  cout << "Total CPU Time = " << tm1 + tm2 << endl;
}

void print_usage() {
  cout << "quadt_test [-dd] [-qd] [-all] [-v] [-x]" << endl;
  cout << "  Performs a selected set of integration using " << endl;
  cout << "  the quad-double library." << endl;
  cout << endl;
  cout << "  -h -help  Print this usage message." << endl;
  cout << "  -d        Perform quadrature test with regular double precision." << endl;
  cout << "  -dd       Perform quadrature test with double-double." << endl;
  cout << "  -qd       Perform quadrature test with quad-double." << endl;
  cout << "            This is the default." << endl;
  cout << "  -all      Perform quadrature test with all three precision types." << endl;
  cout << "  -v" << endl;
  cout << "  -verbose  Prints out detailed test results." << endl;
  cout << "  -x        Perform only the last test, an interesting quadrature" << endl;
  cout << "            whose value is *very* close to pi / 8 (see paper)." << endl;
}
  
int main(int argc, char **argv) {

  char *arg;
  
  /* Parse the command-line flags. */
  for (int i = 1; i < argc; i++) {
    arg = argv[i];
    if (strcmp(arg, "-h") == 0 || strcmp(arg, "-help") == 0) {
      print_usage();
      exit(0);
    } else if (strcmp(arg, "-d") == 0) {
      flag_test_d = true;
    } else if (strcmp(arg, "-dd") == 0) {
      flag_test_dd = true;
    } else if (strcmp(arg, "-qd") == 0) {
      flag_test_qd = true;
    } else if (strcmp(arg, "-all") == 0) {
      flag_test_d = flag_test_dd = flag_test_qd = true;
    } else if (strcmp(arg, "-v") == 0 || strcmp(arg, "-verbose") == 0) {
      flag_verbose = true;
    } else if (strcmp(arg, "-x") == 0) {
      flag_last_only = true;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
    }
  }


  unsigned int old_cw;
  fpu_fix_start(&old_cw);

  if (!flag_test_d && !flag_test_dd && !flag_test_qd)
    flag_test_qd = true;

  double _eps = 1.11022302462516e-16;
  if (flag_test_d)
    test_quadt<double>  (_eps);

  if (flag_test_dd)
    test_quadt<dd_real> (dd_real::_eps);

  if (flag_test_qd) 
    test_quadt<qd_real> (qd_real::_eps);

  fpu_fix_end(&old_cw);
  return 0;
}



