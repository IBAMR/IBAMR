/*
 * tests/quadt.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This contains a C++ class template for a tanh-sinh quadrature
 * algorithm, which employs the transformation
 *
 *   t  <--  tanh (sinh (x))
 *
 * This quadrature scheme is suitable for any function that is 
 * continuous, infinitely differentiable and integrable on a finite
 * open interval.  It can also be used for certain integrals on
 * infinite intervals by making a suitable change of variable.
 * While this routine is not quite as efficient as Gaussian quadrature, 
 * it can be used for functions with an integrable singularity at one 
 * or both of the endpoints.  Further, this scheme has the advantage
 * that function evaluation at one level are all utilized at the next
 * level, thus saving significant computation.
 * 
 * This program is based on David Bailey's tquadt.f program (written
 * in Fortran 90).  C++ conversion, quad-double precision support, 
 * and few other changes have been added.
 */
#ifndef _QUADT_CC_
#define _QUADT_CC_

/* Suppose we are given the integral
 *
 *        / 1
 *   I =  |     f(x) dx
 *        / -1
 *
 * Then the substitution  t = tanh (sinh (x))  gives
 * 
 *   dt = (tanh (sinh (x))' dx = (sech (sinh (x)))^2 * cosh(x) dx
 *
 * Also the limit point x = 1 corresponds to t = 
 */
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <qd/qd_real.h>

template<class T>
class quadt {
public:
  /* Constructor.  This will create a tanh-sinh quadrature class.
     Parameters
       eps    -- The machine epsilon of the variable type to be used. */
  quadt(double eps);

  /* Destructor.  This will take care of disposing internal tables, etc. */
  ~quadt();

  /* Computes the integral of the function f from -1 to 1.
     The class F is any class with overloaded operator() (T &). */
  template <class F>
  int integrate_u(const F &f, double tol, T &result, double &err);
  
  /* Computes the integral of the function f from a to b.
     The class F is any class with overloaded operator() (T &). */
  template <class F>
  int integrate(const F &f, T a, T b, double tol, T &result, double &err);

private:
  int max_level;
  double initial_width, final_width;
  int table_size;
  double eps;

  /* Pre-computed quadrature points. */
  T *weights;
  T *points;

  /* Scales and translates the given function f from the
     interval  [a, b] to [-1, 1] so it can be evaluated using
     the tanh-sinh substitution.                             */
  template <class F>
  class UnitFunction {
  private:
    F f;
    T offset, h;
  public:
    UnitFunction(const F &f, const T &a, const T &b) : f(f) {
      offset = 0.5 * (a + b);
      h = (b - a) * 0.5;
    }
    T operator()(T x) const {
      return f(offset + h * x) * h;
    }
  };

  /* Initializes the weight and abcissa table. */
  void init_table();
};



/*-**** Class Template Implementations ****-*/
template <class T>
quadt<T>::quadt(double eps) {
  max_level = 11;
  initial_width = 0.5;
  final_width = std::ldexp(initial_width, -max_level+1);
  table_size = static_cast<int>(2.0 * 7.0 / final_width);
  this->eps = eps;

  init_table();
}

template <class T>
quadt<T>::~quadt() {
  delete [] weights;
  delete [] points;
}

template <class T>
void quadt<T>::init_table() {

  weights = new T[table_size];
  points  = new T[table_size];

  double h = initial_width * 2.0;
  double dt;
  double t;
  int i = 0;
  T sinh_t, cosh_t, sinh_s, cosh_s;
  T x, w;
  for (int level = 1; level <= max_level; level++, h *= 0.5) {
    t = h * 0.5;
    dt = (level == 1) ? t : h;
    for (;; t += dt) {
      sincosh(T(t), sinh_t, cosh_t);
      sincosh(sinh_t, sinh_s, cosh_s);
      x = sinh_s / cosh_s;
      //      w = (cosh_t / cosh_s) / cosh_s;
      w = (cosh_t / sqr(cosh_s));

      if (x == 1.0 || w < eps) {
        weights[i++] = 0.0;
        break;
      }

      points[i] = x;
      weights[i] = w;
      i++;
    }
  }

}

template <class T> template <class F>
int quadt<T>::integrate_u(const F &f, double tol, 
                          T &result, double &err) {
  T r1, r2, r3, s;
  T x, w;
  int level;
  double h = initial_width;
  bool conv = false;
  int i = 0;
  
  r1 = r2 = r3 = 0.0;
  s = f(T(0.0));
  for (level = 1; level <= max_level; level++, h *= 0.5) {

    /* Compute the integral */
    for (;;) {
      x = points[i];
      w = weights[i];
      i++;
      if (w == 0.0)
        break;
      s += w * (f(x) + f(-x));
    }

    r1 = s * h;

    /* Check for convergence. */
    if (level > 2) {
      double e1, e2, d1, d2;

      e1 = abs(to_double(r1 - r2));
      if (e1 == 0.0)
        err = eps;
      else {
        e2 = abs(to_double(r1 - r3));
        d1 = log(e1);
        d2 = log(e2);
        
        err = exp(d1 * d1 / d2);
      }

      cout << " level = " << level << endl;
      cout << "     r = " << r1 << endl;
      cout << "   err = " << err << endl;

      if (err < abs(r1) * tol) {
        conv = true;
        break;
      }
    }

    r2 = r1;
    r3 = r2;
  }
  
  if (level > max_level)
    puts("Level exhausted.");
  
  result = r1;
  if (!conv) {
    /* No convergence. */
    return -1;
  }
  
  return 0;
}

template <class T> template <class F>
int quadt<T>::integrate(const F &f, T a, T b, double tol, 
                        T &result, double &err) {
  if (a == -1.0 && b == 1.0)
    return integrate_u(f, tol, result, err);
  else {
    UnitFunction<F> unit_f(f, a, b);
    return integrate_u< UnitFunction<F> >(unit_f, tol, result, err);
  }
}

#endif  /* _QUADT_CC_ */


