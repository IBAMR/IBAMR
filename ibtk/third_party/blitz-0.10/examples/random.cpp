/*****************************************************************************
 * random.cpp        Blitz++ random numbers example
 */

#include <random/uniform.h>
#include <blitz/numinquire.h>
#include <time.h>
#include <iostream>
#include <iomanip>

using namespace ranlib;
using namespace blitz;

// workaround for broken streams in Compaq cxx, can't handle long double
#if defined(__DECCXX)
#define LD_HACK(x) static_cast<double>(x)
#else
#define LD_HACK(x) x
#endif

template<typename T>
void printRandoms()
{
  Uniform<T> x;
  //x.seed((unsigned int)time(0));
  x.seed(5);
  int N=5;
  for (int i = 0; i < N; ++i) 
    cout << setprecision(digits10(T())) << LD_HACK(x.random()) << endl;

  cout << endl;
}

int main()
{
// test get/set state interface
  Uniform<double> x;
  Uniform<double>::T_state S = x.getState();
  x.setState(S);
  std::string str = x.getStateString();
  x.setState(str);

  cout << "Some random float: " << endl;
  printRandoms<float>();

  cout << "Some random doubles: " << endl;
  printRandoms<double>();

  cout << "Some random long doubles: " << endl;
  printRandoms<long double>();

  return 0;
}

