// -*- C++ -*-
// $Id$

/*
 * F distribution
 *
 * This code has been adapted from RANDLIB.C 1.3, by
 * Barry W. Brown, James Lovato, Kathy Russell, and John Venier.
 * Code was originally by Ahrens and Dieter (see above).
 *
 * Adapter's notes:
 * BZ_NEEDS_WORK: how to handle seeding for the two gamma RNGs if
 * independentState is used?
 * BZ_NEEDS_WORK: code for handling possible overflow when xden
 * is tiny seems a bit flaky.
 */

#ifndef BZ_RANDOM_F
#define BZ_RANDOM_F

#ifndef BZ_RANDOM_GAMMA
 #include <random/gamma.h>
#endif

BZ_NAMESPACE(ranlib)

template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class F {
public:
    typedef T T_numtype;

    F(T numeratorDF, T denominatorDF)
    {
        setDF(numeratorDF, denominatorDF);
        mindenom = 0.085 * blitz::tiny(T());
    }

  F(T numeratorDF, T denominatorDF, unsigned int i) :
    ngamma(i), dgamma(i)
    {
        setDF(numeratorDF, denominatorDF);
        mindenom = 0.085 * blitz::tiny(T());
    }

    void setDF(T _dfn, T _dfd)
    {
        BZPRECONDITION(_dfn > 0.0);
        BZPRECONDITION(_dfd > 0.0);
        dfn = _dfn;
        dfd = _dfd;

        ngamma.setMean(dfn/2.0);
        dgamma.setMean(dfd/2.0);
    }

    T random()
    {
        T xnum = 2.0 * ngamma.random() / dfn;
        T xden = 2.0 * ngamma.random() / dfd;

        // Rare event: Will an overflow probably occur?
        if (xden <= mindenom)
        {
            // Yes, just return huge(T())
            return blitz::huge(T());
        }

        return xnum / xden;
    }

  void seed(IRNG_int s, IRNG_int r)
    {
        // This is such a bad idea if independentState is used. Ugh.
        // If sharedState is used, it is merely inefficient (the
        // same RNG is seeded twice).

      // yes it's unacceptable -- changed to using two seeds / Patrik
      // in fact should probably be two uncorrelated IRNGs...

        ngamma.seed(s);
        dgamma.seed(r);
    }

protected:
    Gamma<T,IRNG,stateTag> ngamma, dgamma;
    T dfn, dfd, mindenom;
};

BZ_NAMESPACE_END

#endif // BZ_RANDOM_F
