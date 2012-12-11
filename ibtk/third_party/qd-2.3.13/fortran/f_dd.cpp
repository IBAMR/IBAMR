/*
 * src/c_dd.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains the C wrapper functions for double-double precision arithmetic.
 * This can be used from Fortran code.
 */
#include "config.h"
#ifdef HAVE_FORTRAN
#include <cstring>

#include "config.h"
#include <qd/dd_real.h>
#include <qd/c_dd.h>

#define f_dd_add          FC_FUNC_(f_dd_add, F_DD_ADD)
#define f_dd_add_dd_d     FC_FUNC_(f_dd_add_dd_d, F_DD_ADD_DD_D)

#define f_dd_sub          FC_FUNC_(f_dd_sub, F_DD_SUB)
#define f_dd_sub_dd_d     FC_FUNC_(f_dd_sub_dd_d, F_DD_SUB_DD_D)
#define f_dd_sub_d_dd     FC_FUNC_(f_dd_sub_d_dd, F_DD_SUB_D_DD)

#define f_dd_mul          FC_FUNC_(f_dd_mul, F_DD_MUL)
#define f_dd_mul_dd_d     FC_FUNC_(f_dd_mul_dd_d, F_DD_MUL_DD_D)

#define f_dd_div          FC_FUNC_(f_dd_div, F_DD_DIV)
#define f_dd_div_dd_d     FC_FUNC_(f_dd_div_dd_d, F_DD_DIV_DD_D)
#define f_dd_div_d_dd     FC_FUNC_(f_dd_div_d_dd, F_DD_DIV_D_DD)

#define f_dd_sqrt         FC_FUNC_(f_dd_sqrt, F_DD_SQRT)
#define f_dd_sqr          FC_FUNC_(f_dd_sqr, F_DD_SQR)

#define f_dd_abs          FC_FUNC_(f_dd_abs, F_DD_ABS)

#define f_dd_npwr         FC_FUNC_(f_dd_npwr, F_DD_NPWR)
#define f_dd_nroot        FC_FUNC_(f_dd_nroot, F_DD_NROOT)

#define f_dd_nint         FC_FUNC_(f_dd_nint, F_DD_NINT)
#define f_dd_aint         FC_FUNC_(f_dd_aint, F_DD_AINT)
#define f_dd_floor        FC_FUNC_(f_dd_floor, F_DD_FLOOR)
#define f_dd_ceil         FC_FUNC_(f_dd_ceil, F_DD_CEIL)

#define f_dd_exp          FC_FUNC_(f_dd_exp, F_DD_EXP)
#define f_dd_log          FC_FUNC_(f_dd_log, F_DD_LOG)
#define f_dd_log10        FC_FUNC_(f_dd_log10, F_DD_LOG10)

#define f_dd_sin          FC_FUNC_(f_dd_sin, F_DD_SIN)
#define f_dd_cos          FC_FUNC_(f_dd_cos, F_DD_COS)
#define f_dd_tan          FC_FUNC_(f_dd_tan, F_DD_TAN)

#define f_dd_asin         FC_FUNC_(f_dd_asin, F_DD_ASIN)
#define f_dd_acos         FC_FUNC_(f_dd_acos, F_DD_ACOS)
#define f_dd_atan         FC_FUNC_(f_dd_atan, F_DD_ATAN)
#define f_dd_atan2        FC_FUNC_(f_dd_atan2, F_DD_ATAN2)

#define f_dd_sinh         FC_FUNC_(f_dd_sinh, F_DD_SINH)
#define f_dd_cosh         FC_FUNC_(f_dd_cosh, F_DD_COSH)
#define f_dd_tanh         FC_FUNC_(f_dd_tanh, F_DD_TANH)

#define f_dd_asinh        FC_FUNC_(f_dd_asinh, F_DD_ASINH)
#define f_dd_acosh        FC_FUNC_(f_dd_acosh, F_DD_ACOSH)
#define f_dd_atanh        FC_FUNC_(f_dd_atanh, F_DD_ATANH)

#define f_dd_sincos       FC_FUNC_(f_dd_sincos, F_DD_SINCOS)
#define f_dd_sincosh      FC_FUNC_(f_dd_sincosh, F_DD_SINCOSH)

#define f_dd_swrite       FC_FUNC_(f_dd_swrite, F_DD_SWRITE)
#define f_dd_write        FC_FUNC_(f_dd_write, F_DD_WRITE)

#define f_dd_neg          FC_FUNC_(f_dd_neg, F_DD_NEG)
#define f_dd_rand         FC_FUNC_(f_dd_rand, F_DD_RAND)
#define f_dd_comp         FC_FUNC_(f_dd_comp, F_DD_COMP)
#define f_dd_comp_dd_d    FC_FUNC_(f_dd_comp_dd_d, F_DD_COMP_DD_D)

#define f_dd_pi           FC_FUNC_(f_dd_pi, F_DD_PI)
#define f_dd_nan          FC_FUNC_(f_dd_nan, F_DD_NAN)

#define TO_DOUBLE_PTR(a, ptr) ptr[0] = a.x[0]; ptr[1] = a.x[1];

extern "C" {

/* add */
void f_dd_add(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = dd_real(a) + dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_dd_add_dd_d(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = dd_real(a) + *b;
  TO_DOUBLE_PTR(cc, c);
}


/* sub */
void f_dd_sub(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = dd_real(a) - dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_dd_sub_dd_d(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = dd_real(a) - *b;
  TO_DOUBLE_PTR(cc, c);
}
void f_dd_sub_d_dd(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = *a - dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}


/* mul */
void f_dd_mul(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = dd_real(a) * dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_dd_mul_dd_d(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = dd_real(a) * *b;
  TO_DOUBLE_PTR(cc, c);
}


/* div */
void f_dd_div(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = dd_real(a) / dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_dd_div_dd_d(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = dd_real(a) / *b;
  TO_DOUBLE_PTR(cc, c);
}
void f_dd_div_d_dd(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = *a / dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}


void f_dd_sqrt(const double *a, double *b) {
  dd_real bb;
  bb = sqrt(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_sqr(const double *a, double *b) {
  dd_real bb;
  bb = sqr(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_dd_abs(const double *a, double *b) {
  dd_real bb;
  bb = abs(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_dd_npwr(const double *a, const int *n, double *b) {
  dd_real bb;
  bb = npwr(dd_real(a), *n);
  TO_DOUBLE_PTR(bb, b);
}

void f_dd_nroot(const double *a, const int *n, double *b) {
  dd_real bb;
  bb = nroot(dd_real(a), *n);
  TO_DOUBLE_PTR(bb, b);
}

void f_dd_nint(const double *a, double *b) {
  dd_real bb;
  bb = nint(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_aint(const double *a, double *b) {
  dd_real bb;
  bb = aint(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_floor(const double *a, double *b) {
  dd_real bb;
  bb = floor(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_ceil(const double *a, double *b) {
  dd_real bb;
  bb = ceil(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_dd_log(const double *a, double *b) {
  dd_real bb;
  bb = log(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_log10(const double *a, double *b) {
  dd_real bb;
  bb = log10(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_exp(const double *a, double *b) {
  dd_real bb;
  bb = exp(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_dd_sin(const double *a, double *b) {
  dd_real bb;
  bb = sin(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_cos(const double *a, double *b) {
  dd_real bb;
  bb = cos(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_tan(const double *a, double *b) {
  dd_real bb;
  bb = tan(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_dd_asin(const double *a, double *b) {
  dd_real bb;
  bb = asin(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_acos(const double *a, double *b) {
  dd_real bb;
  bb = acos(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_atan(const double *a, double *b) {
  dd_real bb;
  bb = atan(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_dd_atan2(const double *a, const double *b, double *c) {
  dd_real cc;
  cc = atan2(dd_real(a), dd_real(b));
  TO_DOUBLE_PTR(cc, c);
}

void f_dd_sinh(const double *a, double *b) {
  dd_real bb;
  bb = sinh(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_cosh(const double *a, double *b) {
  dd_real bb;
  bb = cosh(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_tanh(const double *a, double *b) {
  dd_real bb;
  bb = tanh(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_dd_asinh(const double *a, double *b) {
  dd_real bb;
  bb = asinh(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_acosh(const double *a, double *b) {
  dd_real bb;
  bb = acosh(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_dd_atanh(const double *a, double *b) {
  dd_real bb;
  bb = atanh(dd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_dd_sincos(const double *a, double *s, double *c) {
  dd_real ss, cc;
  sincos(dd_real(a), ss, cc);
  TO_DOUBLE_PTR(ss, s);
  TO_DOUBLE_PTR(cc, c);
}

void f_dd_sincosh(const double *a, double *s, double *c) {
  dd_real ss, cc;
  sincosh(dd_real(a), ss, cc);
  TO_DOUBLE_PTR(ss, s);
  TO_DOUBLE_PTR(cc, c);
}

/* Writes a dd_real into a character array of length maxlen, with
 * the given precision.   The rest of the array will be filled with
 * spaces.  Parameter maxlen should at least be precision + 7 
 * characters.  Prec can be zero to put out the defaut number of 
 * digits. */
void f_dd_swrite(const double *a, int *precision, char *s, int *maxlen) {
  int prec = *precision;
  if (prec <= 0 || prec > dd_real::_ndigits) prec = dd_real::_ndigits;
  std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0);
  std::string str = dd_real(a).to_string(prec, 0, fmt, false, true);

  int len = 0;
  if (a[0] < 0.0) {
    strncpy(&s[len], str.c_str(), *maxlen - len);
  } else {
    s[len++] = ' ';
    strncpy(&s[len], str.c_str(), *maxlen - len);
  }

  len += str.length();
  for (int i = len; i < *maxlen; i++)
    s[i] = ' ';
}

void f_dd_write(const double *a) {
  std::cout << dd_real(a) << std::endl;
}

void f_dd_neg(const double *a, double *b) {
  b[0] = -a[0];
  b[1] = -a[1];
}

void f_dd_rand(double *a) {
  dd_real aa;
  aa = ddrand();
  TO_DOUBLE_PTR(aa, a);
}

void f_dd_comp(const double *a, const double *b, int *result) {
  dd_real aa(a), bb(b);
  if (aa < bb)
    *result = -1;
  else if (aa > bb)
    *result = 1;
  else 
    *result = 0;
}

void f_dd_comp_dd_d(const double *a, const double *b, int *result) {
  dd_real aa(a);
  if (aa < *b)
    *result = -1;
  else if (aa > *b)
    *result = 1;
  else 
    *result = 0;
}

void f_dd_comp_d_dd(const double *a, const double *b, int *result) {
  dd_real bb(b);
  if (*a < bb)
    *result = -1;
  else if (*a > bb)
    *result = 1;
  else 
    *result = 0;
}

void f_dd_pi(double *a) {
  TO_DOUBLE_PTR(dd_real::_pi, a);
}

void f_dd_nan(double *a) {
  TO_DOUBLE_PTR(dd_real::_nan, a);
}

}
#endif /* HAVE_FORTRAN */
