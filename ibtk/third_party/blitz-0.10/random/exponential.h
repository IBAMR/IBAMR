// -*- C++ -*-
// $Id$

/*
 * This generator uses the straightforward transformation
 *  x = - log(y) * m
 *
 * to turn a uniform (0,1) y into an exponentially distributed
 * variable x.  x has density function
 *
 * f(x) = (1/m) exp(-(1/m)x)  (x > 0)
 *
 * and mean m.
 *
 * NEEDS_WORK: Adapt the method of Ahrens and Dieter.  This will
 * require extending the precision of the constants.
 *
 * Ahrens, J.H. and Dieter, U.  Computer Methods for Sampling From the
 * Exponential and Normal Distributions. Comm. ACM, 15,10 (Oct. 1972), p. 873.
 */

#ifndef BZ_RANDOM_EXPONENTIAL
#define BZ_RANDOM_EXPONENTIAL

#ifndef BZ_RANDOM_UNIFORM
 #include <random/uniform.h>
#endif

BZ_NAMESPACE(ranlib)

template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class ExponentialUnit : public UniformOpen<T,IRNG,stateTag>
{
public:
    typedef T T_numtype;

	ExponentialUnit() {}

  explicit ExponentialUnit(unsigned int i) :
    UniformOpen<T,IRNG,stateTag>(i) {};

    T random()
    {
        return - log(UniformOpen<T,IRNG,stateTag>::random());
    }
};

template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class Exponential : public ExponentialUnit<T,IRNG,stateTag> {

public:
    typedef T T_numtype;

    Exponential(T mean)
    {
        mean_ = mean;
    }

  Exponential(T mean, unsigned int i) :
    UniformOpen<T,IRNG,stateTag>(i)
  {
        mean_ = mean;
  };

    T random()
    {
        return mean_ * ExponentialUnit<T,IRNG,stateTag>::random();
    }

private:
    T mean_;
};

BZ_NAMESPACE_END

#endif // BZ_RANDOM_EXPONENTIAL
