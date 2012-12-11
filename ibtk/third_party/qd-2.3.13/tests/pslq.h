/*
 * tests/pslq.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 *
 * Implementation of PSLQ Integer Relation Detection Algorithm
 * due to H. R. P. Ferguson and D. H. Bailey.  See
 *
 * A new polynomial time algorithm for finding relations among
 * real numbers, Supercomputing Research Center Tech Report 
 * SRC-93-093 (March 1993).                                 
 *
 * This code is based in part on David Bailey's F90 version.
 */

#include <cmath>
#include <cstdio>

#include <qd/dd_real.h>
#include <qd/qd_real.h>

using std::sqrt;
using std::abs;

using namespace qd;

#define MIN(a, b) ( ((a)<(b)) ? (a) : (b) )
#define MAX(a, b) ( ((a)>(b)) ? (a) : (b) )
#define SWAP(a, b)  { t = a; a = b; b = t; }
#define SQR(a) ( (a)*(a) )

/* Create an n-by-m matrix of T.  Diagonal elements
   are initialized to diag, while all other elements are
   initialized to elem.                                      */
template <class T>
T **new_matrix(int nr_rows, int nr_cols, 
                    T diag = 0.0, T elem = 0.0) {
  T **m = new T *[nr_rows];
  int i, j;

  for (i = 0; i < nr_rows; i++)
    m[i] = new T[nr_cols];
  for (i = 0; i < nr_rows; i++)
    for (j = 0; j < nr_cols; j++)
      m[i][j] = (i == j) ? diag : elem;
  return m;
}

/* Create an n-vector of T.  Each elements are initialized to elem. */
template <class T>
T *new_vector(int n, T elem = 0.0) {
  T *v = new T[n];
  for (int i = 0; i < n; i++)
    v[i] = elem;
  return v;
}

/* Destroys a matrix m. */
template <class T>
void delete_matrix(T **m, int nr_rows) {
  for (int i = 0; i < nr_rows; i++)
    delete [] m[i];
  delete [] m;
}

/* Destroys a vector v. */
template <class T>
void delete_vector(T *v) {
  delete [] v;
}

/* The parameter gamma used in the PSLQ algorithm. */
static const double gam = 1.2;

/* Perform PSLQ integer relation algorithm to find any
   linear relation among the n numbers in the vector x.
   It returns the coefficients found in the vector r.
   The parameter eps provides the precision of type T. */   
template <class T>
int pslq(const T *x, int n, T *r, double eps, int max_itr) {
  T *s = new_vector<T>(n);
  T *y = new_vector<T>(n);
  T **a = new_matrix<T>(n, n, 1.0);
  T **b = new_matrix<T>(n, n, 1.0);
  T **h = new_matrix<T>(n, n-1, 0.0);
  T t;
  double teps = 16.0 * eps;
  int i, j, k;
  int err = 0;

  /* Compute the vector s */
  t = x[n-1] * x[n-1];
  s[n-1] = abs(x[n-1]);
  for (i = n-2; i >= 0; i--) {
    t += x[i] * x[i];
    s[i] = sqrt(t);
  }

  /* Normalize the vector x, s. */
  t = s[0];
  for (i = 0; i < n; i++)
    y[i] = x[i] / t;
  s[0] = 1.0;
  for (i = 1; i < n; i++)
    s[i] /= t;

  /* Construct matrix H. */
  for (i = 0; i < n; i++) {
    for (j = 0; j <= MIN(i, n-2); j++) {
      h[i][j] = (i == j) ? s[j+1]/s[j] : - y[i]*y[j] / (s[j] * s[j+1]);
    }
  }
  
  /* Reduce matrix H. */
  for (i = 1; i < n; i++) {
    for (j = i-1; j >= 0; j--) {
      t = nint(h[i][j] / h[j][j]);
      y[j] += t * y[i];
      for (k = 0; k <= j; k++)
        h[i][k] -= t * h[j][k];
      for (k = 0; k < n; k++) {
        a[i][k] -= t * a[j][k];
        b[k][j] += t * b[k][i];
      }
    }
  }
  
  int m;
  int itr = 0;
  bool done = false;

  while(!done) {

    itr++;
    
    /* Select m such that gam^i * |H_ii| is maximal when i = m. */
    T m_val = -1.0;
    T g = gam;
    m = -1;

    for (i = 0; i < n-1; i++, g *= gam) {
      t = abs(g * h[i][i]);
      if (t > m_val) {
        m_val = t;
        m = i;
      }
    }

    if (m < 0) { 
      /* This shouldn't happen. */
      err = 1;
      break;
    }

    /* Exchange entries m and m+1 of y, 
                   rows m and m+1 of A and H, 
                columns m and m+1 of B.        */
    SWAP(y[m], y[m+1])
    for (i = 0; i < n; i++) {
      SWAP(a[m][i], a[m+1][i])
    }
    for (i = 0; i < n-1; i++) {
      SWAP(h[m][i], h[m+1][i])
    }
    for (i = 0; i < n; i++) {
      SWAP(b[i][m], b[i][m+1])
    }


    /* Remove the corner on H diagonal. */
    if (m < n-2) {
      T t0, t1, t2, t3, t4;
      t0 = sqrt(SQR(h[m][m]) + SQR(h[m][m+1]));
      t1 = h[m][m] / t0;
      t2 = h[m][m+1] / t0;
      for (i = m; i < n; i++) {
        t3 = h[i][m];
        t4 = h[i][m+1];
        h[i][m] = t1 * t3 + t2 * t4;
        h[i][m+1] = t1 * t4 - t2 * t3;
      }
    }

    /* Reduce H. */
    for (i = m+1; i < n; i++) {
      for (j = MIN(i-1, m+1); j >= 0; j--) {
        t = nint(h[i][j]/h[j][j]);
        y[j] += t * y[i];
        for (k = 0; k <= j; k++) {
          h[i][k] -= t * h[j][k];
        }
        for (k = 0; k < n; k++) {
          a[i][k] -= t * a[j][k];
          b[k][j] += t * b[k][i];
        }
      }
    }


    /* Norm bound */
    m_val = -1.0e308;
    for (j = 0; j < n-1; j++) {
      t = abs(h[j][j]);
      if (t > m_val)
        m_val = t;
    }

    /* Check the y vector for zeros. */
    for (i = 0; i < n; i++) {
      t = abs(y[i]);
      if (t < teps) {
        m = i;
        done = true;
        break;
      }
    }

    if (itr > max_itr) {
      done = true;
      err = -1;
    }

  }  /* while */

  /* Get the coefficients. */
  if (err == 0) {
    for (i = 0; i < n; i++) {
      r[i] = b[i][m];
    }
  }

  delete_matrix<T>(h, n);
  delete_matrix<T>(a, n);
  delete_matrix<T>(b, n);
  delete_vector<T>(y);
  delete_vector<T>(s);

  return err;
}

