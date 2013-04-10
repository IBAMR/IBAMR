// -*- C++ -*-
// $Id$

/*
 * Chi^2 distribution
 *
 * This code has been adapted from RANDLIB.C 1.3, by
 * Barry W. Brown, James Lovato, Kathy Russell, and John Venier.
 * Code was originally by Ahrens and Dieter (see above).
 *
 * Adapter's notes:
 */

#ifndef BZ_RANDOM_CHISQUARE
#define BZ_RANDOM_CHISQUARE

#ifndef BZ_RANDOM_GAMMA
 #include <random/gamma.h>
#endif

BZ_NAMESPACE(ranlib)

template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class ChiSquare : public Gamma<T,IRNG,stateTag>
{
public:
    typedef T T_numtype;

    ChiSquare(T df)
      : Gamma<T,IRNG,stateTag>(df/2.0) // isn't this redundant with setDF call?
    {
        setDF(df);
    }

  ChiSquare(T df, unsigned int i)
    : Gamma<T,IRNG,stateTag>(df/2.0, i)
    {
        setDF(df);
    }

    void setDF(T _df)
    {
        BZPRECONDITION(_df > 0.0);
        df = _df;
        Gamma<T,IRNG,stateTag>::setMean(df/2.0);
    }

    T random()
    {
        return 2.0 * sgamma();
    }

protected:
    T sgamma()
    {
        return Gamma<T,IRNG,stateTag>::random();
    }

    T df;
};

BZ_NAMESPACE_END

#endif // BZ_RANDOM_CHISQUARE
