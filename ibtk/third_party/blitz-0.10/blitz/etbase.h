// -*- C++ -*-
/***************************************************************************
 * blitz/etbase.h    Declaration of the ETBase<T> class
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

#ifndef BZ_ETBASE_H
#define BZ_ETBASE_H

BZ_NAMESPACE(blitz)

template<typename T>
class ETBase {
public:
  typedef T T_unwrapped;

    ETBase()
    { }

    ETBase(const ETBase<T>&)
    { }

    T& unwrap() { return static_cast<T&>(*this); }

    const T& unwrap() const { return static_cast<const T&>(*this); }

    ETBase<T>& wrap() { return static_cast<ETBase<T>&>(*this); }

    const ETBase<T>& wrap() const { return static_cast<const ETBase<T>&>(*this); }
};

BZ_NAMESPACE_END

#endif // BZ_ETBASE_H

