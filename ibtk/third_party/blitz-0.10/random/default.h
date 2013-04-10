// -*- C++ -*-
/***************************************************************************
 * random/default.h       Default IRNG wrapper class
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

#ifndef BZ_RANDOM_DEFAULT_H
#define BZ_RANDOM_DEFAULT_H

#include <random/mt.h>

BZ_NAMESPACE(ranlib)

// Some terminology:
// IRNG = Integer Random Number Generator.  IRNGs generate random
//        integers, which are used to create floating-point random
//        numbers.
// RNG  = Random Number Generator.  RNGs use IRNGs to create floating-
//        point random numbers following desired distributions.

typedef float defaultType;

// These are type tags.  A RNG with sharedState shares an IRNG
// with other RNGs.  An RNG with independentState
// contains its own IRNG.  Generally, sharedState RNGs should be
// used.

struct sharedState { };
struct independentState { };
typedef sharedState defaultState;

typedef unsigned int IRNG_int;


// IRNGWrapper handles shared and independent state IRNGs.
// If a class inherits from IRNGWrapper<IRNG,sharedState>,
// it gets a static IRNG (i.e. the IRNG state is shared among
// all RNGs); if it inherits from IRNGWrapper<IRNG,independentState>,
// it gets an independent IRNG (the IRNG state is encapsulated
// in the RNG, and is not shared among RNGs).

template<typename IRNG, typename state>
class IRNGWrapper {
};

template<typename IRNG>
class IRNGWrapper<IRNG,sharedState> {

public:
    void seed(IRNG_int x)
    { irng_.seed(x); }

  void seed(std::vector<IRNG_int> x)
    { irng_.seed(x); }

    typedef typename IRNG::T_state T_state;
    T_state getState() const { return irng_.getState(); }
    std::string getStateString() const { return irng_.getStateString(); }
    void setState(const T_state& s) { irng_.setState(s); }
    void setState(const std::string& s) { irng_.setState(s); }

protected:
    static IRNG irng_;
};

template<typename IRNG>
IRNG IRNGWrapper<IRNG,sharedState>::irng_;

template<typename IRNG>
class IRNGWrapper<IRNG,independentState> {

public:
  IRNGWrapper() {};
  IRNGWrapper(unsigned int i) : irng_(MersenneTwisterCreator::create(i)) {};

    void seed(IRNG_int x)
    { irng_.seed(x); }

  void seed(std::vector<IRNG_int> x)
    { irng_.seed(x); }

    typedef typename IRNG::T_state T_state;
    T_state getState() const { return irng_.getState(); }
    std::string getStateString() const { return irng_.getStateString(); }
    void setState(const T_state& s) { irng_.setState(s); }
    void setState(const std::string& s) { irng_.setState(s); }

protected:
    IRNG irng_;
};

// defaultIRNG is a type alias for the default Integer Random
// Number Generator (IRNG).

typedef MersenneTwister defaultIRNG;

BZ_NAMESPACE_END

#endif // BZ_RANDOM_DEFAULT_H

