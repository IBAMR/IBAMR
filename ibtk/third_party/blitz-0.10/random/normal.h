// -*- C++ -*-
// $Id$

/*
 * This is a modification of the Kinderman + Monahan algorithm for
 * generating normal random numbers, due to Leva:
 *
 * J.L. Leva, Algorithm 712. A normal random number generator, ACM Trans.
 * Math. Softw.  18 (1992) 454--455.
 *
 * http://www.acm.org/pubs/citations/journals/toms/1992-18-4/p449-leva/
 *
 * Note: Some of the constants used below look like they have dubious
 * precision.  These constants are used for an approximate bounding
 * region test (see the paper).  If the approximate test fails,
 * then an exact region test is performed.
 *
 * Only 0.012 logarithm evaluations are required per random number
 * generated, making this method comparatively fast.
 *
 * Adapted to C++ by T. Veldhuizen.
 */

#ifndef BZ_RANDOM_NORMAL
#define BZ_RANDOM_NORMAL

#ifndef BZ_RANDOM_UNIFORM
 #include <random/uniform.h>
#endif

BZ_NAMESPACE(ranlib)

template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class NormalUnit : public UniformOpen<T,IRNG,stateTag>
{
public:
    typedef T T_numtype;

  NormalUnit() {}

  explicit NormalUnit(unsigned int i) :
    UniformOpen<T,IRNG,stateTag>(i) {};

    T random()
    {
        const T s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472;
        const T r1 = 0.27597, r2 = 0.27846;

        T u, v;

        for (;;) {
          // Generate P = (u,v) uniform in rectangle enclosing
          // acceptance region:
          //   0 < u < 1
          // - sqrt(2/e) < v < sqrt(2/e)
          // The constant below is 2*sqrt(2/e).

          u = this->getUniform();
          v = 1.715527769921413592960379282557544956242L
              * (this->getUniform() - 0.5);

          // Evaluate the quadratic form
          T x = u - s;
          T y = fabs(v) - t;
          T q = x*x + y*(a*y - b*x);

          // Accept P if inside inner ellipse
          if (q < r1)
            break;

          // Reject P if outside outer ellipse
          if (q > r2)
            continue;

          // Between ellipses: perform exact test
          if (v*v <= -4.0 * log(u)*u*u)
            break;
        }

        return v/u;
    }

};


template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class Normal : public NormalUnit<T,IRNG,stateTag> {

public:
    typedef T T_numtype;

    Normal(T mean, T standardDeviation)
    {
        mean_ = mean;
        standardDeviation_ = standardDeviation;
    }

  Normal(T mean, T standardDeviation, unsigned int i) :
    NormalUnit<T,IRNG,stateTag>(i)
  {
        mean_ = mean;
        standardDeviation_ = standardDeviation;
  };

    T random()
    {
        return mean_ + standardDeviation_
           * NormalUnit<T,IRNG,stateTag>::random();
    }

private:
    T mean_;
    T standardDeviation_;
};

BZ_NAMESPACE_END

#endif // BZ_RANDOM_NORMAL
