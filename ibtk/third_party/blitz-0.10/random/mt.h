// -*- C++ -*-

/*
 * $Id$
 *
 * A C-program for MT19937: Integer version (1998/4/6)
 *  genrand() generates one pseudorandom unsigned integer (32bit)
 * which is uniformly distributed among 0 to 2^32-1  for each
 * call. sgenrand(seed) set initial values to the working area
 * of 624 words. Before genrand(), sgenrand(seed) must be
 * called once. (seed is any 32-bit integer except for 0).
 *   Coded by Takuji Nishimura, considering the suggestions by
 * Topher Cooper and Marc Rieffel in July-Aug. 1997.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later
 * version.
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Library General Public License for more details.
 * You should have received a copy of the GNU Library General
 * Public License along with this library; if not, write to the
 * Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 * Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
 * When you use this, send an email to: matumoto@math.keio.ac.jp
 * with an appropriate reference to your work.
 *
 * REFERENCE
 * M. Matsumoto and T. Nishimura,
 * "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform
 * Pseudo-Random Number Generator",
 * ACM Transactions on Modeling and Computer Simulation,
 * Vol. 8, No. 1, January 1998, pp 3--30.
 *
 * See
 *     http://www.math.keio.ac.jp/~matumoto/emt.html
 * and
 *     http://www.acm.org/pubs/citations/journals/tomacs/1998-8-1/p3-matsumoto/
 *
 */

#ifndef BZ_RAND_MT
#define BZ_RAND_MT

#include <blitz/blitz.h>

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iterator>

#ifndef UINT_MAX
  #include <limits.h>
#endif

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#endif

BZ_NAMESPACE(ranlib)

#if UINT_MAX < 4294967295U
  typedef unsigned long twist_int;  // must be at least 32 bits
#else
  typedef unsigned int twist_int;
#endif

class MersenneTwister
{
public:

private:

#if defined(BZ_HAVE_NAMESPACES) && defined(BZ_HAVE_STD)
  typedef std::vector<twist_int> State;
#else
  typedef vector<twist_int> State;
#endif
  typedef State::size_type SizeType;
  typedef State::iterator Iter;

  // Implements Step 2 and half of Step 3 in the MN98 description
  struct BitMixer {
    BitMixer() : s0(0), K(0x9908b0df) {};
    BitMixer(twist_int k) : s0(0), K(k) {};

    // Return 0 if lsb of s1=0, a if lsb of s1=1.
    inline twist_int low_mask (twist_int s1) const {
      // This does not actually result in a branch because it's
      // equivalent to ( -(s1 & 1u) ) & K
      return (s1&1u) ? K : 0u;
    }

    // Return y>>1 in MN98 (step 2 + part of 3).
    // y = (x[i] AND u) OR (x[i+1 mod n] AND ll), where s0=x[i] and
    // s1=x[i+1 mod n]
    inline twist_int high_mask (twist_int s1) const {
      return ( (s0&0x80000000) | (s1&0x7fffffff) ) >>1;
    }

    // Calculate (half of) step 3 in MN98.
    inline twist_int operator() (twist_int s1) {
      // (y>>1) XOR (0 if lsb of s1=0, a if lsb of s1=1)
      // x[i+m] is XORed in reload
      // (Note that it is OK to call low_mask with s1 (x[i+1]) and not
      // with y, like MN98 does, because s1&1 == y&1 by construction.
      twist_int r = high_mask(s1) ^ low_mask(s1);
      s0 = s1;
      return r;
    }
    twist_int s0; // this is "x[i]" in the MN98 description
    const twist_int K; // MN98 "a" vector
  };

enum { N = 624,
       PF = 397, // MN98 "m"
       reference_seed = 4357 };

  void initialize()
  {
    S.resize(N);
    I = S.end();
  }

public:
  MersenneTwister() : b_(0x9D2C5680), c_(0xEFC60000)
  {
    initialize();
    seed();

    // There is a problem: static initialization + templates do not
    // mix very well in C++.  If you have a static member
    // of a class template, there is no guarantee on its order iin
    // static initialization.  This MersenneTwister class is used
    // elsewhere as a static member of a template class, and it is
    // possible (in fact, I've done so) to create a static initializer
    // that will invoke the seed() method of this object before its
    // ctor has been called (result: crash).
    // ANSI C++ is stranger than fiction.
    // Currently the documentation forbids using RNGs from
    // static initializers.  There doesn't seem to be a good
    // fix.
  }

  MersenneTwister(twist_int aa, twist_int bb, twist_int cc) :
    twist_(aa), b_(bb), c_(cc)
  {
    initialize();
    seed();
  }

  MersenneTwister(twist_int initial_seed) : b_(0x9D2C5680), c_(0xEFC60000)
  {
    initialize();
    seed(initial_seed);
  }

  MersenneTwister(twist_int aa, twist_int bb, twist_int cc,
		  twist_int initial_seed) : twist_(aa), b_(bb), c_(cc)
  {
    initialize();
    seed(initial_seed);
  }

  // Seed. Uses updated seed algorithm from mt19937ar. The old
  // algorithm would yield sequences that were close from seeds that
  // were close.
  void seed (twist_int seed = reference_seed)
  {
    // seed cannot equal 0
    if (seed == 0)
      seed = reference_seed;

    S[0] = seed & 0xFFFFFFFF;
    for (SizeType mti=1; mti<S.size(); ++mti) {
      S[mti] = (1812433253U * (S[mti-1] ^ (S[mti-1] >> 30)) + mti);
      S[mti] &= 0xffffffffU;
    }

    reload();
  }

  // Seed by array, swiped directly from mt19937ar. Gives a larger
  // initial seed space.
  void seed (State seed_vector)
  {
    SizeType i, j, k;
    seed(19650218U);
    i=1; j=0;
    const SizeType N=S.size();
	const SizeType n=seed_vector.size();
    k = (N>n ? N : n);
    for (; k; k--) {
      S[i] = (S[i] ^ ((S[i-1] ^ (S[i-1] >> 30)) * 1664525U))
	+ seed_vector[j] + j; /* non linear */
      S[i] &= 0xffffffffU; /* for WORDSIZE > 32 machines */
      i++; j++;
      if (i>=N) { S[0] = S[N-1]; i=1; }
      if (j>=n) j=0;
    }
    for (k=N-1; k; k--) {
      S[i] = (S[i] ^ ((S[i-1] ^ (S[i-1] >> 30)) * 1566083941UL))
	- i; /* non linear */
      S[i] &= 0xffffffffU; /* for WORDSIZE > 32 machines */
      i++;
      if (i>=N) { S[0] = S[N-1]; i=1; }
    }

    S[0] = 0x80000000U; /* MSB is 1; assuring non-zero initial array */

    reload();
  }

  // generate a full new x array
  void reload (void)
  {
    // This check is required because it is possible to call random()
    // before the constructor.  See the note above about static
    // initialization.

    Iter p0 = S.begin();
    Iter pM = p0 + PF;
    twist_ (S[0]); // set x[i]=x[0] in the twister (prime the pump)
    for (Iter pf_end = S.begin()+(N-PF); p0 != pf_end; ++p0, ++pM)
      // mt[kk] = mt[kk+m] XOR ((y>>1)XOR(mag01), as calc by BitMixer)
      *p0 = *pM ^ twist_ (p0[1]);

    // This is the "modulo part" where kk+m rolls over
    pM = S.begin();
    for (Iter s_end = S.begin()+(N-1); p0 != s_end; ++p0, ++pM)
      *p0 = *pM ^ twist_ (p0[1]);
    // and final element where kk+1 rolls over
    *p0 = *pM ^ twist_ (S[0]);

    I = S.begin();
  }

  inline twist_int random (void)
  {
    if (I >= S.end()) reload();
    // get next word from array
    twist_int y = *I++;
    // Step 4+5 in MN98, multiply by tempering matrix
    y ^= (y >> 11);
    y ^= (y <<  7) & b_;
    y ^= (y << 15) & c_;
    y ^= (y >> 18);
    return y;
  }

  // functions for getting/setting state
  class mt_state {
    friend class MersenneTwister;
  private:
    State S;
    State::difference_type I;
  public:
    mt_state() { }
	mt_state(State s, State::difference_type i) : S(s), I(i) { }
    mt_state(const std::string& s) {
      std::istringstream is(s);
      is >> I;
      S = State(std::istream_iterator<twist_int>(is),
		std::istream_iterator<twist_int>());
      assert(!S.empty());
    }
    operator bool() const { return !S.empty(); }
    std::string str() const {
      if (S.empty())
	return std::string();
      std::ostringstream os;
      os << I << " ";
      std::copy(S.begin(), S.end(),
		std::ostream_iterator<twist_int>(os," "));
      return os.str();
    }
#ifdef BZ_HAVE_BOOST_SERIALIZATION
    friend class boost::serialization::access;
    /** Serialization operator. Relies on the ability to serialize
	std::vector. */
    template<class T_arch>
    void serialize(T_arch& ar, const unsigned int version) {
      ar & S & I;
  };
#endif

  };

  typedef mt_state T_state;
  T_state getState() const { return T_state(S, I-S.begin()); }
  std::string getStateString() const {
    T_state tmp(S, I-S.begin());
    return tmp.str();
  }
  void setState(const T_state& s) {
    if (!s) {
      std::cerr << "Error: state is empty" << std::endl;
      return;
    }
    S = s.S;
    I = S.begin() + s.I;
  }
  void setState(const std::string& s) {
    T_state tmp(s);
    setState(tmp);
  }

private:
  BitMixer twist_;
  const twist_int b_,c_;

  State   S;
  Iter    I;
};


/** This class creates MersenneTwisters with different parameters
    indexed by and ID number. */
class MersenneTwisterCreator
{
public:
  static MersenneTwister create(unsigned int i) {
    // We only have n different parameter sets
    i = i % n;
    return MersenneTwister(a_[i], b_[i], c_[i]);
  };

private:
  static const unsigned int n=48;
  static const twist_int a_[n];
  static const twist_int b_[n];
  static const twist_int c_[n];
};

BZ_NAMESPACE_END

#endif // BZ_RAND_MT
