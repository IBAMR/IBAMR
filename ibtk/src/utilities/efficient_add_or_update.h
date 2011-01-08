// Filename: IBTK_Utilities.h
// Created on 06 Mar 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_efficient_add_or_update
#define included_efficient_add_or_update

//////////////////////////////////////////////////////////////////////////////

namespace IBTK
{
/*!
 * \brief A templated function that adds to or updates the contents of a STL map
 * in an efficient manner.
 *
 * \note See Item 24 in Scott Meyers' book "Effective STL."
 */
template<typename MapType, typename KeyArgType, typename ValueArgType>
typename MapType::iterator
efficient_add_or_update(
    MapType& m,
    const KeyArgType& k,
    const ValueArgType& v)
{
    // Find either where k is or where k should be in the map.
    typename MapType::iterator lb = m.lower_bound(k);

    // If k belongs to the map, update the value of the pair and return an
    // iterator to the pair.
    //
    // Otherwise, add pair(k,v) to m and return an iterator to the new map
    // element.
    if (lb != m.end() && !m.key_comp()(k,lb->first) )
    {
        lb->second = v;
        return lb;
    }
    else
    {
        typedef typename MapType::value_type MVT;
        return m.insert(lb,MVT(k,v));
    }
}// efficient_add_or_update
}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_efficient_add_or_update
