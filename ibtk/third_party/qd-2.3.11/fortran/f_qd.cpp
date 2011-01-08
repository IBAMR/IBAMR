/*
 * src/c_qd.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains C wrapper function for quad-double precision arithmetic.
 * This can be used from fortran code.
 */
#include "config.h"
#ifdef HAVE_FORTRAN

#include <cstring>

#include "config.h"
#include <qd/qd_real.h>
#include <qd/c_qd.h>

#define f_qd_add          FC_FUNC_(f_qd_add, F_QD_ADD)
#define f_qd_add_qd_dd    FC_FUNC_(f_qd_add_qd_dd, F_QD_ADD_QD_DD)
#define f_qd_add_qd_d     FC_FUNC_(f_qd_add_qd_d, F_QD_ADD_QD_D)

#define f_qd_sub          FC_FUNC_(f_qd_sub, F_QD_SUB)
#define f_qd_sub_dd_qd    FC_FUNC_(f_qd_sub_dd_qd, F_QD_SUB_DD_QD)
#define f_qd_sub_qd_dd    FC_FUNC_(f_qd_sub_qd_dd, F_QD_SUB_QD_DD)
#define f_qd_sub_d_qd     FC_FUNC_(f_qd_sub_d_qd, F_QD_SUB_D_QD)
#define f_qd_sub_qd_d     FC_FUNC_(f_qd_sub_qd_d, F_QD_SUB_QD_D)

#define f_qd_mul          FC_FUNC_(f_qd_mul, F_QD_MUL)
#define f_qd_mul_qd_dd    FC_FUNC_(f_qd_mul_qd_dd, F_QD_MUL_QD_DD)
#define f_qd_mul_qd_d     FC_FUNC_(f_qd_mul_qd_d, F_QD_MUL_QD_D)

#define f_qd_div          FC_FUNC_(f_qd_div, F_QD_DIV)
#define f_qd_div_dd_qd    FC_FUNC_(f_qd_div_dd_qd, F_QD_DIV_DD_QD)
#define f_qd_div_qd_dd    FC_FUNC_(f_qd_div_qd_dd, F_QD_DIV_QD_DD)
#define f_qd_div_d_qd     FC_FUNC_(f_qd_div_d_qd, F_QD_DIV_D_QD)
#define f_qd_div_qd_d     FC_FUNC_(f_qd_div_qd_d, F_QD_DIV_QD_D)

#define f_qd_sqrt         FC_FUNC_(f_qd_sqrt, F_QD_SQRT)
#define f_qd_sqr          FC_FUNC_(f_qd_sqr, F_QD_SQR)

#define f_qd_abs          FC_FUNC_(f_qd_abs, F_QD_ABS)

#define f_qd_npwr         FC_FUNC_(f_qd_npwr, F_QD_NPWR)
#define f_qd_nroot        FC_FUNC_(f_qd_nroot, F_QD_NROOT)

#define f_qd_nint         FC_FUNC_(f_qd_nint, F_QD_NINT)
#define f_qd_aint         FC_FUNC_(f_qd_aint, F_QD_AINT)
#define f_qd_floor        FC_FUNC_(f_qd_floor, F_QD_FLOOR)
#define f_qd_ceil         FC_FUNC_(f_qd_ceil, F_QD_CEIL)

#define f_qd_exp          FC_FUNC_(f_qd_exp, F_QD_EXP)
#define f_qd_log          FC_FUNC_(f_qd_log, F_QD_LOG)
#define f_qd_log10        FC_FUNC_(f_qd_log10, F_QD_LOG10)

#define f_qd_sin          FC_FUNC_(f_qd_sin, F_QD_SIN)
#define f_qd_cos          FC_FUNC_(f_qd_cos, F_QD_COS)
#define f_qd_tan          FC_FUNC_(f_qd_tan, F_QD_TAN)

#define f_qd_asin         FC_FUNC_(f_qd_asin, F_QD_ASIN)
#define f_qd_acos         FC_FUNC_(f_qd_acos, F_QD_ACOS)
#define f_qd_atan         FC_FUNC_(f_qd_atan, F_QD_ATAN)
#define f_qd_atan2        FC_FUNC_(f_qd_atan2, F_QD_ATAN2)

#define f_qd_sinh         FC_FUNC_(f_qd_sinh, F_QD_SINH)
#define f_qd_cosh         FC_FUNC_(f_qd_cosh, F_QD_COSH)
#define f_qd_tanh         FC_FUNC_(f_qd_tanh, F_QD_TANH)

#define f_qd_asinh        FC_FUNC_(f_qd_asinh, F_QD_ASINH)
#define f_qd_acosh        FC_FUNC_(f_qd_acosh, F_QD_ACOSH)
#define f_qd_atanh        FC_FUNC_(f_qd_atanh, F_QD_ATANH)

#define f_qd_sincos       FC_FUNC_(f_qd_sincos, F_QD_SINCOS)
#define f_qd_sincosh      FC_FUNC_(f_qd_sincosh, F_QD_SINCOSH)

#define f_qd_swrite       FC_FUNC_(f_qd_swrite, F_QD_SWRITE)
#define f_qd_write        FC_FUNC_(f_qd_write, F_QD_WRITE)

#define f_qd_neg          FC_FUNC_(f_qd_neg, F_QD_NEG)
#define f_qd_rand         FC_FUNC_(f_qd_rand, F_QD_RAND)
#define f_qd_comp         FC_FUNC_(f_qd_comp, F_QD_COMP)
#define f_qd_comp_qd_d    FC_FUNC_(f_qd_comp_qd_d, F_QD_COMP_QD_D)

#define f_qd_pi           FC_FUNC_(f_qd_pi, F_QD_PI)
#define f_qd_nan          FC_FUNC_(f_qd_nan, F_QD_NAN)

#define TO_DOUBLE_PTR(a, ptr) ptr[0] = a.x[0]; ptr[1] = a.x[1]; \
                              ptr[2] = a.x[2]; ptr[3] = a.x[3];

extern "C" {



/* add */
void f_qd_add(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) + qd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_add_qd_dd(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) + dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_add_qd_d(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) + *b;
  TO_DOUBLE_PTR(cc, c);
}



/* sub */
void f_qd_sub(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) - qd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_sub_qd_dd(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) - dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_sub_dd_qd(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = dd_real(a) - qd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_sub_qd_d(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) - *b;
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_sub_d_qd(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = *a - qd_real(b);
  TO_DOUBLE_PTR(cc, c);
}



/* mul */
void f_qd_mul(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) * qd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_mul_qd_dd(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) * dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_mul_qd_d(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) * *b;
  TO_DOUBLE_PTR(cc, c);
}



/* div */
void f_qd_div(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) / qd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_div_qd_dd(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) / dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_div_dd_qd(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = dd_real(a) / qd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_div_qd_d(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = qd_real(a) / *b;
  TO_DOUBLE_PTR(cc, c);
}
void f_qd_div_d_qd(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = *a / qd_real(b);
  TO_DOUBLE_PTR(cc, c);
}




/* selfadd */
void f_qd_selfadd(const double *a, double *b) {
  qd_real bb(b);
  bb += qd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_selfadd_dd(const double *a, double *b) {
  qd_real bb(b);
  bb += dd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_selfadd_d(const double *a, double *b) {
  qd_real bb(b);
  bb += *a;
  TO_DOUBLE_PTR(bb, b);
}



/* selfsub */
void f_qd_selfsub(const double *a, double *b) {
  qd_real bb(b);
  bb -= qd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_selfsub_dd(const double *a, double *b) {
  qd_real bb(b);
  bb -= dd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_selfsub_d(const double *a, double *b) {
  qd_real bb(b);
  bb -= *a;
  TO_DOUBLE_PTR(bb, b);
}



/* selfmul */
void f_qd_selfmul(const double *a, double *b) {
  qd_real bb(b);
  bb *= qd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_selfmul_dd(const double *a, double *b) {
  qd_real bb(b);
  bb *= dd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_selfmul_d(const double *a, double *b) {
  qd_real bb(b);
  bb *= *a;
  TO_DOUBLE_PTR(bb, b);
}



/* selfdiv */
void f_qd_selfdiv(const double *a, double *b) {
  qd_real bb(b);
  bb /= qd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_selfdiv_dd(const double *a, double *b) {
  qd_real bb(b);
  bb /= dd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_selfdiv_d(const double *a, double *b) {
  qd_real bb(b);
  bb /= *a;
  TO_DOUBLE_PTR(bb, b);
}



void f_qd_sqrt(const double *a, double *b) {
  qd_real bb;
  bb = sqrt(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_sqr(const double *a, double *b) {
  qd_real bb;
  bb = sqr(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_qd_abs(const double *a, double *b) {
  qd_real bb;
  bb = abs(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_qd_npwr(double *a, int *n, double *b) {
  qd_real bb;
  bb = npwr(qd_real(a), *n);
  TO_DOUBLE_PTR(bb, b);
}

void f_qd_nroot(double *a, int *n, double *b) {
  qd_real bb;
  bb = nroot(qd_real(a), *n);
  TO_DOUBLE_PTR(bb, b);
}

void f_qd_nint(const double *a, double *b) {
  qd_real bb;
  bb = nint(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_aint(const double *a, double *b) {
  qd_real bb;
  bb = aint(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_floor(const double *a, double *b) {
  qd_real bb;
  bb = floor(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_ceil(const double *a, double *b) {
  qd_real bb;
  bb = ceil(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_qd_log(const double *a, double *b) {
  qd_real bb;
  bb = log(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_log10(const double *a, double *b) {
  qd_real bb;
  bb = log10(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_exp(const double *a, double *b) {
  qd_real bb;
  bb = exp(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_qd_sin(const double *a, double *b) {
  qd_real bb;
  bb = sin(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_cos(const double *a, double *b) {
  qd_real bb;
  bb = cos(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_tan(const double *a, double *b) {
  qd_real bb;
  bb = tan(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_qd_asin(const double *a, double *b) {
  qd_real bb;
  bb = asin(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_acos(const double *a, double *b) {
  qd_real bb;
  bb = acos(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_atan(const double *a, double *b) {
  qd_real bb;
  bb = atan(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_qd_atan2(const double *a, const double *b, double *c) {
  qd_real cc;
  cc = atan2(qd_real(a), qd_real(b));
  TO_DOUBLE_PTR(cc, c);
}

void f_qd_sinh(const double *a, double *b) {
  qd_real bb;
  bb = sinh(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_cosh(const double *a, double *b) {
  qd_real bb;
  bb = cosh(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_tanh(const double *a, double *b) {
  qd_real bb;
  bb = tanh(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_qd_asinh(const double *a, double *b) {
  qd_real bb;
  bb = asinh(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_acosh(const double *a, double *b) {
  qd_real bb;
  bb = acosh(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void f_qd_atanh(const double *a, double *b) {
  qd_real bb;
  bb = atanh(qd_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void f_qd_sincos(const double *a, double *s, double *c) {
  qd_real ss, cc;
  sincos(qd_real(a), ss, cc);
  TO_DOUBLE_PTR(cc, c);
  TO_DOUBLE_PTR(ss, s);
}

void f_qd_sincosh(const double *a, double *s, double *c) {
  qd_real ss, cc;
  sincosh(qd_real(a), ss, cc);
  TO_DOUBLE_PTR(cc, c);
  TO_DOUBLE_PTR(ss, s);
}

/* Writes a dd_real into a character array of length maxlen, with
 * the given precision.   The rest of the array will be filled with
 * spaces.  Parameter maxlen should at least be precision + 7 
 * characters.  Prec can be zero to put out the defaut number of 
 * digits. */
void f_qd_swrite(const double *a, int *precision, char *s, int *maxlen) {
  int prec = *precision;
  if (prec <= 0 || prec > qd_real::_ndigits) prec = qd_real::_ndigits;
  std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0);
  std::string str = qd_real(a).to_string(prec, 0, fmt, false, true);

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

void f_qd_write(const double *a) {
  std::cout << qd_real(a) << std::endl;
}

void f_qd_neg(const double *a, double *b) {
  b[0] = -a[0];
  b[1] = -a[1];
  b[2] = -a[2];
  b[3] = -a[3];
}

void f_qd_rand(double *a) {
  qd_real aa;
  aa = qdrand();
  TO_DOUBLE_PTR(aa, a);
}

void f_qd_comp(const double *a, const double *b, int *result) {
  qd_real aa(a), bb(b);
  if (aa < bb)
    *result = -1;
  else if (aa > bb)
    *result = 1;
  else 
    *result = 0;
}

void f_qd_comp_qd_d(const double *a, const double *b, int *result) {
  qd_real aa(a);
  if (aa < *b)
    *result = -1;
  else if (aa > *b)
    *result = 1;
  else 
    *result = 0;
}

void f_qd_comp_d_qd(const double *a, const double *b, int *result) {
  qd_real bb(b);
  if (*a < bb)
    *result = -1;
  else if (*a > bb)
    *result = 1;
  else 
    *result = 0;
}

void f_qd_pi(double *a) {
  TO_DOUBLE_PTR(qd_real::_pi, a);
}

void f_qd_nan(double *a) {
  TO_DOUBLE_PTR(qd_real::_nan, a);
}

}

#endif /* HAVE_FORTRAN */

