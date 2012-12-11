// -*- C++ -*-
/***************************************************************************
 * random/discrete-uniform.h       Discrete uniform IRNG wrapper class
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

#ifndef BZ_RANDOM_DISCRETE_UNIFORM_H
#define BZ_RANDOM_DISCRETE_UNIFORM_H

#include <random/default.h>

BZ_NAMESPACE(ranlib)

template<typename T = unsigned int, typename IRNG = defaultIRNG,
    typename stateTag = defaultState>
class DiscreteUniform : public IRNGWrapper<IRNG,stateTag>
{
public:
    typedef T T_numtype;

    DiscreteUniform(T n)
    {
        BZPRECONDITION(n < 4294967295U);
        n_ = n;
    }

  DiscreteUniform(T n, unsigned int i) :
    IRNGWrapper<IRNG,stateTag>::IRNGWrapper(i)
    {
        BZPRECONDITION(n < 4294967295U);
        n_ = n;
    }

    T random()
    {
        return this->irng_.random() % n_;
    }

private:
    T n_;
};

BZ_NAMESPACE_END

#endif // BZ_RANDOM_DISCRETE_UNIFORM_H
