// -*- C++ -*-
/***************************************************************************
 * random/uniform.h       Uniform IRNG wrapper class
 *
 * $Id$
 *
 * Copyright (C) 1997-2011 Todd Veldhuizen <tveldhui@acm.org>
 *
 * This file is a part of Blitz.
 *
 * Blitz is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 *
 * Blitz is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Blitz.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Suggestions:          blitz-devel@lists.sourceforge.net
 * Bugs:                 blitz-support@lists.sourceforge.net
 *
 * For more information, please see the Blitz++ Home Page:
 *    https://sourceforge.net/projects/blitz/
 *
 ***************************************************************************/

#ifndef BZ_RANDOM_UNIFORM_H
#define BZ_RANDOM_UNIFORM_H

#include <random/default.h>

#ifndef FLT_MANT_DIG
 #include <float.h>
#endif

BZ_NAMESPACE(ranlib)

/*****************************************************************************
 * UniformClosedOpen generator: uniform random numbers in [0,1).
 *****************************************************************************/

template<typename T = defaultType, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class UniformClosedOpen { };

// These constants are 1/2^32, 1/2^64, 1/2^96, 1/2^128
const long double norm32open = .2328306436538696289062500000000000000000E-9L;
const long double norm64open = .5421010862427522170037264004349708557129E-19L;
const long double norm96open = .1262177448353618888658765704452457967477E-28L;
const long double norm128open = .2938735877055718769921841343055614194547E-38L;


template<typename IRNG, typename stateTag>
class UniformClosedOpen<float,IRNG,stateTag>
  : public IRNGWrapper<IRNG,stateTag>
{

public:
    typedef float T_numtype;

  UniformClosedOpen() {};
  UniformClosedOpen(unsigned int i) :
    IRNGWrapper<IRNG,stateTag>::IRNGWrapper(i) {};

    float random()
    {
#if FLT_MANT_DIG > 96
 #if (FLT_MANT_DIG > 128) && !defined(BZ_IGNORE_RNG_ERRORS)
  #error RNG code assumes float mantissa is <= 128 bits (not true for your platform).  Use -DBZ_IGNORE_RNG_ERRORS to ignore this warning.
 #endif
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();
        IRNG_int i4 = this->irng_.random();
        return i1 * norm128open + i2 * norm96open + i3 * norm64open
            + i4 * norm32open;
#elif FLT_MANT_DIG > 64
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();
        return i1 * norm96open + i2 * norm64open + i3 * norm32open;
#elif FLT_MANT_DIG > 32
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        return i1 * norm64open + i2 * norm32open;
#else
        IRNG_int i1 = this->irng_.random();
        return i1 * norm32open;
#endif
    }

    float getUniform()
    { return random(); }
};

template<typename IRNG, typename stateTag>
class UniformClosedOpen<double,IRNG,stateTag>
  : public IRNGWrapper<IRNG,stateTag>
{
public:
    typedef double T_numtype;

  UniformClosedOpen() {}
  UniformClosedOpen(unsigned int i) :
    IRNGWrapper<IRNG,stateTag>::IRNGWrapper(i) {};

    double random()
    {
#if DBL_MANT_DIG > 96
 #if (DBL_MANT_DIG > 128) && !defined(BZ_IGNORE_RNG_ERRORS)
  #error RNG code assumes double mantissa is <= 128 bits (not true for your platform).  Use -DBZ_IGNORE_RNG_ERRORS to ignore this warning.
 #endif

        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();
        IRNG_int i4 = this->irng_.random();
        return i1 * norm128open + i2 * norm96open + i3 * norm64open
            + i4 * norm32open;
#elif DBL_MANT_DIG > 64
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();
        return i1 * norm96open + i2 * norm64open + i3 * norm32open;
#elif DBL_MANT_DIG > 32
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        return i1 * norm64open + i2 * norm32open;
#else
        IRNG_int i1 = this->irng_.random();
        return i1 * norm32open;
#endif

    }

    double getUniform() { return random(); }
};

template<typename IRNG, typename stateTag>
class UniformClosedOpen<long double,IRNG,stateTag>
  : public IRNGWrapper<IRNG,stateTag>
{
public:
    typedef long double T_numtype;

  UniformClosedOpen() {};
  UniformClosedOpen(unsigned int i) :
    IRNGWrapper<IRNG,stateTag>::IRNGWrapper(i) {};

    long double random()
    {
#if LDBL_MANT_DIG > 96
 #if (LDBL_MANT_DIG > 128) && !defined(BZ_IGNORE_RNG_ERRORS)
  #error RNG code assumes long double mantissa is <= 128 bits (not true for your platform).  Use -DBZ_IGNORE_RNG_ERRORS to ignore this warning.
#endif

        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();
        IRNG_int i4 = this->irng_.random();
        return i1 * norm128open + i2 * norm96open + i3 * norm64open
            + i4 * norm32open;
#elif LDBL_MANT_DIG > 64
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();
        return i1 * norm96open + i2 * norm64open + i3 * norm32open;
#elif LDBL_MANT_DIG > 32
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        return i1 * norm64open + i2 * norm32open;
#else
        IRNG_int i1 = this->irng_.random();
        return i1 * norm32open;
#endif
    }

    long double getUniform() { return random(); }
};

// For people who don't care about open or closed: just give them
// ClosedOpen (this is the fastest).

template<class T, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class Uniform : public UniformClosedOpen<T,IRNG,stateTag>
{
public:
  Uniform() {};
  Uniform(unsigned int i) :
    UniformClosedOpen<T,IRNG,stateTag>::UniformClosedOpen(i) {}
 };

/*****************************************************************************
 * UniformClosed generator: uniform random numbers in [0,1].
 *****************************************************************************/

// This constant is 1/(2^32-1)
const long double norm32closed = .2328306437080797375431469961868475648078E-9L;

// These constants are 2^32/(2^64-1) and 1/(2^64-1), respectively.

const long double norm64closed1 =
    .23283064365386962891887177448353618888727188481031E-9L;
const long double norm64closed2 =
    .54210108624275221703311375920552804341370213034169E-19L;

// These constants are 2^64/(2^96-1), 2^32/(2^96-1), and 1/(2^96-1)
const long double norm96closed1 = .2328306436538696289062500000029387358771E-9L;
const long double norm96closed2 =
    .5421010862427522170037264004418131333707E-19L;
const long double norm96closed3 =
    .1262177448353618888658765704468388886588E-28L;

// These constants are 2^96/(2^128-1), 2^64/(2^128-1), 2^32/(2^128-1) and
// 1/(2^128-1).
const long double norm128clos1 = .2328306436538696289062500000000000000007E-9L;
const long double norm128clos2 = .5421010862427522170037264004349708557145E-19L;
const long double norm128clos3 = .1262177448353618888658765704452457967481E-28L;
const long double norm128clos4 = .2938735877055718769921841343055614194555E-38L;


template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class UniformClosed { };

template<typename IRNG, typename stateTag>
class UniformClosed<float,IRNG,stateTag> : public IRNGWrapper<IRNG,stateTag> {

public:
    typedef float T_numtype;

  UniformClosed() {};
  UniformClosed(unsigned int i) :
    IRNGWrapper<IRNG,stateTag>::IRNGWrapper(i) {};

    float random()
    {
#if FLTMANT_DIG > 96
 #if (FLT_MANT_DIG > 128) && !defined(BZ_IGNORE_RNG_ERRORS)
  #error RNG code assumes float mantissa is <= 128 bits (not true for your platform).  Use -DBZ_IGNORE_RNG_ERRORS to ignore this warning.
 #endif
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();
        IRNG_int i4 = this->irng_.random();

        return i1 * norm128clos1 + i2 * norm128clos2
          + i3 * norm128clos3 + i4 * norm128clos4;
#elif FLT_MANT_DIG > 64
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();

        return i1 * norm96closed1 + i2 * norm96closed2
          + i3 * norm96closed3;
#elif FLT_MANT_DIG > 32
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();

        return i1 * norm64closed1 + i2 * norm64closed2;
#else
        IRNG_int i = this->irng_.random();
        return i * norm32closed;
#endif

    }

    float getUniform()
    { return random(); }
};

template<typename IRNG, typename stateTag>
class UniformClosed<double,IRNG,stateTag> : public IRNGWrapper<IRNG,stateTag> {

public:
    typedef double T_numtype;

  UniformClosed() {};
  UniformClosed(unsigned int i) :
    IRNGWrapper<IRNG,stateTag>::IRNGWrapper(i) {};

    double random()
    {
#if DBL_MANT_DIG > 96
 #if (DBL_MANT_DIG > 128) && !defined(BZ_IGNORE_RNG_ERRORS)
  #error RNG code assumes double mantissa is <= 128 bits (not true for your platform).  Use -DBZ_IGNORE_RNG_ERRORS to ignore this warning.
 #endif
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();
        IRNG_int i4 = this->irng_.random();

        return i1 * norm128clos1 + i2 * norm128clos2
          + i3 * norm128clos3 + i4 * norm128clos4;
#elif DBL_MANT_DIG > 64
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();

        return i1 * norm96closed1 + i2 * norm96closed2
          + i3 * norm96closed3;
#elif DBL_MANT_DIG > 32
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();

        return i1 * norm64closed1 + i2 * norm64closed2;
#else
        IRNG_int i = this->irng_.random();
        return i * norm32closed;
#endif

    }

    double getUniform()
    { return random(); }
};

template<typename IRNG, typename stateTag>
class UniformClosed<long double,IRNG,stateTag>
  : public IRNGWrapper<IRNG,stateTag> {

public:
    typedef long double T_numtype;

  UniformClosed() {};
  UniformClosed(unsigned int i) :
    IRNGWrapper<IRNG,stateTag>::IRNGWrapper(i) {};

    long double random()
    {
#if LDBL_MANT_DIG > 96
 #if (LDBL_MANT_DIG > 128) && !defined(BZ_IGNORE_RNG_ERRORS)
  #error RNG code assumes long double mantissa is <= 128 bits (not true for your platform).  Use -DBZ_IGNORE_RNG_ERRORS to ignore this warning.
 #endif
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();
        IRNG_int i4 = this->irng_.random();

        return i1 * norm128clos1 + i2 * norm128clos2
          + i3 * norm128clos3 + i4 * norm128clos4;
#elif LDBL_MANT_DIG > 64
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();
        IRNG_int i3 = this->irng_.random();

        return i1 * norm96closed1 + i2 * norm96closed2
          + i3 * norm96closed3;
#elif LDBL_MANT_DIG > 32
        IRNG_int i1 = this->irng_.random();
        IRNG_int i2 = this->irng_.random();

        return i1 * norm64closed1 + i2 * norm64closed2;
#else
        IRNG_int i = this->irng_.random();
        return i * norm32closed;
#endif
    }

    long double getUniform()
    { return random(); }

};

/*****************************************************************************
 * UniformOpen generator: uniform random numbers in (0,1).
 *****************************************************************************/

template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class UniformOpen : public UniformClosedOpen<T,IRNG,stateTag> {
  public:
    typedef T T_numtype;

  UniformOpen() {};
  UniformOpen(unsigned int i) :
    UniformClosedOpen<T,IRNG,stateTag>(i) {};

    T random()
    {
        // Turn a [0,1) into a (0,1) interval by weeding out
        // any zeros.
        T x;
        do {
          x = UniformClosedOpen<T,IRNG,stateTag>::random();
        } while (x == 0.0L);

        return x;
    }

    T getUniform()
    { return random(); }

};

/*****************************************************************************
 * UniformOpenClosed generator: uniform random numbers in (0,1]
 *****************************************************************************/

template<typename T = double, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class UniformOpenClosed : public UniformClosedOpen<T,IRNG,stateTag> {

public:
    typedef T T_numtype;

  UniformOpenClosed() {};
  UniformOpenClosed(unsigned int i) :
    UniformClosedOpen<T,IRNG,stateTag>(i) {};


    T random()
    {
        // Antithetic value: taking 1-X where X is [0,1) results
        // in a (0,1] distribution.
        return 1.0 - UniformClosedOpen<T,IRNG,stateTag>::random();
    }

    T getUniform()
    { return random(); }
};

BZ_NAMESPACE_END

#endif // BZ_RANDOM_UNIFORM_H
