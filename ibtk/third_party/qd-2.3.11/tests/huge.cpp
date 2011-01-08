/*
 * tests/huge.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2007
 *
 * This contains tests for check for accuracy when dealing with numbers 
 * near overflow.
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <qd/qd_real.h>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

// Global flags passed to the main program.
static bool flag_test_dd = false;
static bool flag_test_qd = false;
bool flag_verbose = false;

bool print_result(bool result) {
  if (result)
    cout << "Test passed." << endl;
  else
    cout << "Test FAILED." << endl;
  return result;
}

void print_usage() {
  cout << "qd_test [-h] [-dd] [-qd] [-all]" << endl;
  cout << "  Tests output of large numbers." << endl;
  cout << endl;
  cout << "  -h -help  Prints this usage message." << endl;
  cout << "  -dd       Perform tests with double-double types." << endl;
  cout << "  -qd       Perform tests with quad-double types." << endl;
  cout << "            This is the default." << endl;
  cout << "  -all      Perform both double-double and quad-double tests." << endl;
  cout << "  -v" << endl;
  cout << "  -verbose  Print detailed information for each test." << endl;
}

bool check(string str, string true_str) {
  bool pass = (str == true_str);
  if (!pass) {
    cout << "     fail: " << str << endl;
    cout << "should be: " << true_str << endl;
  } else if (flag_verbose) {
    cout << "     pass: " << str << endl;
  }
  return pass;
}

template <class T>
bool test_huge() {
  bool pass = true;
  int digits = T::_ndigits - 1;
  T x = T::_pi * T("1.0e290");

  string pi_str = T::_pi.to_string(digits, 0, std::ios_base::fixed);
  if (flag_verbose) cout << pi_str << endl;
  for (int i = 0; i < 18; i++, x *= 10.0) {
    std::ostringstream os;
    os << pi_str << "e+" << (290 + i);
    pass &= check(x.to_string(digits), os.str());
  }

  x = -T::_pi * T("1.0e290");
  pi_str = "-" + pi_str;
  for (int i = 0; i < 18; i++, x *= 10.0) {
    std::ostringstream os;
    os << pi_str << "e+" << (290 + i);
    pass &= check(x.to_string(digits), os.str());
  }
  
  return pass;
}

template <class T>
bool test_max(string true_str) {
  bool pass = true;
  int digits = T::_ndigits - 1;
  pass &= check(T::_max.to_string(digits), true_str);
  pass &= check((-T::_max).to_string(digits), "-" + true_str);
  return pass;
}

int main(int argc, char *argv[]) {
  
  bool pass = true;
  unsigned int old_cw;
  fpu_fix_start(&old_cw);

  /* Parse the arguments. */
  for (int i = 1; i < argc; i++) {
    string arg(argv[i]);

    if (arg == "-h" || arg == "-help") {
      print_usage();
      exit(0);
    } else if (arg == "-dd") {
      flag_test_dd = true;
    } else if (arg == "-qd") {
      flag_test_qd = true;
    } else if (arg == "-all") {
      flag_test_dd = flag_test_qd = true;
    } else if (arg == "-v" || arg == "-verbose") {
      flag_verbose = true;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
    }
  }

  /* If no flag, test both double-double and quad-double. */
  if (!flag_test_dd && !flag_test_qd) {
    flag_test_dd = true;
    flag_test_qd = true;
  }

  cout << "Testing output of huge numbers..." << endl;

  if (flag_test_dd) {
    cout << endl;
    cout << "Testing dd_real ..." << endl;
    pass &= test_huge<dd_real>();
    pass &= test_max<dd_real>("1.797693134862315807937289714053e+308");
    print_result(pass);
  }

  if (flag_test_qd) {
    cout << endl;
    cout << "Testing qd_real ..." << endl;
    pass &= test_huge<qd_real>();
    pass &= test_max<qd_real>(
        "1.7976931348623158079372897140530286112296785259868571699620069e+308");
    print_result(pass);
  }
  
  fpu_fix_end(&old_cw);
  return (pass ? 0 : 1);
}

