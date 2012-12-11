/*
 * tests/pslq_test.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * A driver for the pslq program which exercises the double-double and 
 * quad-double library.
 */

#include <cmath>
#include <iostream>
#include <ctime>
#include <cstring>
#include <limits>
#include <iomanip>
#include <qd/fpu.h>

#include "tictoc.h"
#include "pslq.h"

using std::cout;
using std::cerr;
using std::endl;
using std::strcmp;

int g_verbose = 0;
bool flag_double_pslq = false;
bool flag_dd_pslq = false;
bool flag_qd_pslq = false;

/* Computes the value of the given n-th degree polynomial at point x
   where the (n+1) coefficients of the polynomial is given in a. */
template <class T>
T polyeval(T *a, int n, T &x, double &err_bnd) {
  /* Use Horner's evaluation scheme. */

  T t = a[n];
  err_bnd = std::abs(to_double(t)) * 0.5;
  for (int i = n-1; i >= 0; i--) {
    t *= x;
    t += a[i];
    err_bnd *= std::abs(to_double(x));
    err_bnd += std::abs(to_double(t));
  }
  err_bnd = (2.0 * err_bnd - to_double(t)) * std::numeric_limits<T>::epsilon();

  return t;
}

double nroot(double x, int n) {
  return std::pow(x, 1.0 / n);
}

bool is_int(double x) {
  return (std::abs(x) <= std::numeric_limits<int>::max() &&
          static_cast<int>(x) == x);
}

/* Sets r = 2^(1/p) + 3^(1/q) and tries to recover the algebraic
 * polynomial of degree pq performing PSLQ on 1, r, r^2, ..., r^n. */
template <class T>
bool pslq_test(int p, int q, double eps, int max_itr = 100000) {
  T *x, *b;
  T r = nroot(T(2.0), p) + nroot(T(3.0), q);
  T t;
  int err;
  tictoc tv;
  double tm;
  int n = p * q + 1;
  std::ios_base::fmtflags fmt = cout.flags();

  b = new T[n];
  x = new T[n];

  /* Fill in vector x with powers of r. */
  x[0] = 1.0;
  x[1] = r;
  t = r*r;
  for (int i = 2; i < n; i++, t *= r) x[i] = t;
  
  cout << "  testing pslq_test(" << p << ", " << q << ") ..." << endl;
  if (g_verbose) cout << std::setprecision(std::numeric_limits<T>::digits10) << "    r = " << r << endl;

  /* Construct algebraic relation */
  tic(&tv);
  err = pslq<T>(x, n, b, eps, max_itr);
  tm = toc(&tv);

  cout << "    elapsed time = " << std::setprecision(4) << tm << " seconds." << endl;
  cout << std::right << std::setprecision(2) << std::fixed;
  if (!err) {
    if (g_verbose) {
      cout << "    polynomial: ";
      for (int i = 0; i < n; i++) { 
        if (i > 0) cout << "                ";
        cout << std::setprecision(0) << std::setw(24) << b[i] << endl;
      }
    }

    /* Check if r satisfies the polynomial. */
    double err_bnd;
    t = abs(polyeval<T>(b, n-1, r, err_bnd));
    err = t > 10.0 * err_bnd;
    cout << std::scientific << std::setprecision(4);
    if (err || g_verbose) {
      cout << "    residual    = " << t << endl;
      cout << "    error bound = " << err_bnd << endl;
    }
  }

  delete [] x;
  delete [] b;

  if (err)
    cout << "  test FAILED." << endl;
  else
    cout << "  test passed." << endl;
  cout << endl;

  cout.flags(fmt);
  return !err;
}

/* We need this since Sun C++ compiler seems to miscompile when
 * eps parameter is given default (templated) argument. */
template <class T>
bool pslq_test(int p, int q) {
  return pslq_test<T>(p, q, std::numeric_limits<T>::epsilon());
}

void print_usage() {
  cout << "pslq_test [-h] [-n N] [-d] [-dd] [-qd] [-all] [-verbose]" << endl;
  cout << "  Performs the PSLQ algorithm on 1, r, r^2, ..., r^{n-1}" << endl;
  cout << "  where r is a root of a constructed integer coefficient" << endl;
  cout << "  polynomial.  PSLQ algorithm should reconstruct the polynomial" << endl;
  cout << "  in most cases where the degree is not too high and the" << endl;
  cout << "  polynomial is irreducible over the rationals." << endl;
  cout << endl;
  cout << "  -h -help  Print this usage message and exit." << endl;
  cout << "  -d        Perform PSLQ with double precision (53 bit mantissa)." << endl;
  cout << "  -dd       Perform PSLQ with double-double precision." << endl;
  cout << "            (about 106 bits of significand)." << endl;
  cout << "  -qd       Perform PSLQ with quad-double precision." << endl;
  cout << "            (about 212 bits of significand).  This is the default." << endl;
  cout << "  -all      Perform PSLQ with all three precisions above." << endl;
  cout << "  -verbose" << endl;
  cout << "  -v        Increase verbosity." << endl;
}

int main(int argc, char **argv) {
  char *arg;

  /* Parse the command-line arguments. */
  for (int i = 1; i < argc; i++) {
    arg = argv[i];
    if (strcmp(arg, "-h") == 0 || strcmp(arg, "-help") == 0) {
      print_usage();
      return 0;
    } else if (strcmp(arg, "-d") == 0) {
      flag_double_pslq = true;
    } else if (strcmp(arg, "-dd") == 0) {
      flag_dd_pslq = true;
    } else if (strcmp(arg, "-qd") == 0) {
      flag_qd_pslq = true;
    } else if (strcmp(arg, "-all") == 0) {
      flag_double_pslq = flag_dd_pslq = flag_qd_pslq = true;
    } else if (strcmp(arg, "-v") == 0 || strcmp(arg, "-verbose") == 0) {
      g_verbose++;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
    }
  }

  if (!flag_double_pslq && !flag_dd_pslq && !flag_qd_pslq) {
    flag_dd_pslq = true;
    flag_qd_pslq = true;
  }

  unsigned int old_cw;
  fpu_fix_start(&old_cw);

  bool pass = true;
  if (flag_double_pslq) {
    cout << "Performing double-precision PSLQ." << endl;
    pass &= pslq_test<double>(2, 2);
    pass &= pslq_test<double>(3, 2);
  }

  if (flag_dd_pslq) {
    cout << "Performing double-double precision PSLQ." << endl;
    pass &= pslq_test<dd_real>(2, 2);
    pass &= pslq_test<dd_real>(2, 3);
    pass &= pslq_test<dd_real>(2, 4);
    pass &= pslq_test<dd_real>(3, 3);
    pass &= pslq_test<dd_real>(2, 5);
  }

  if (flag_qd_pslq) {
    cout << "Performing quad-double precision PSLQ." << endl;
    pass &= pslq_test<qd_real>(3, 3);
    pass &= pslq_test<dd_real>(2, 5);
    pass &= pslq_test<qd_real>(4, 3);
    pass &= pslq_test<qd_real>(2, 6);
    pass &= pslq_test<qd_real>(2, 7);
    pass &= pslq_test<qd_real>(3, 5);
  }

  if (pass)
    cout << "All tests passed." << endl;
  else
    cout << "Some tests FAILED." << endl;

  fpu_fix_end(&old_cw);
  return (pass ? 0 : 1);
}


