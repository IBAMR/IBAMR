// Filename: LNodeIndexIterator.I
// Created on 11 Dec 2009 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
LNodeIndexIterator::LNodeIndexIterator()
    : d_box(),
      d_index_it(),
      d_node_set(),
      d_node_it()
{
    // intentionally blank
    return;
}// LNodeIndexIterator

inline
LNodeIndexIterator::LNodeIndexIterator(
    const LNodeIndexIterator& that)
    : d_box(that.d_box),
      d_index_it(that.d_index_it),
      d_node_set(that.d_node_set),
      d_node_it(that.d_node_it)
{
    // intentionally blank
    return;
}// LNodeIndexIterator

inline
LNodeIndexIterator::~LNodeIndexIterator()
{
    // intentionally blank
    return;
}// ~LNodeIndexIterator

inline LNodeIndexIterator&
LNodeIndexIterator::operator=(
    const LNodeIndexIterator& that)
{
    if (this != &that)
    {
        d_box      = that.d_box;
        d_index_it = that.d_index_it;
        d_node_set = that.d_node_set;
        d_node_it  = that.d_node_it;
    }
    return *this;
}// operator=

inline bool
LNodeIndexIterator::operator==(
    const LNodeIndexIterator& that)
{
    return ((d_node_set == NULL && that.d_node_set == NULL) ||
            (d_box      == that.d_box &&
             d_index_it == that.d_index_it &&
             d_node_set == that.d_node_set &&
             d_node_it  == that.d_node_it));
}// operator==

inline bool
LNodeIndexIterator::operator!=(
    const LNodeIndexIterator& that)
{
    return !(*this == that);
}// operator!=

inline LNodeIndexIterator&
LNodeIndexIterator::operator++()
{
    if (d_node_set == NULL) return *this;
    if (++d_node_it != d_node_set->end()) return *this;
    d_index_it++;
    if (d_index_it)
    {
        d_node_set = &(*d_index_it);
        while (d_index_it && !d_box.contains(d_index_it.getIndex()))
        {
            d_index_it++;
        }

        if (d_index_it && d_box.contains(d_index_it.getIndex()))
        {
            d_node_set = &(*d_index_it);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!d_node_set->empty());
#endif
            d_node_it = d_node_set->begin();
        }
        else
        {
            d_node_set = NULL;
        }
    }
    else
    {
        d_node_set = NULL;
    }
    return *this;
}// operator++

inline LNodeIndexIterator
LNodeIndexIterator::operator++(
    int)
{
    LNodeIndexIterator tmp(*this);
    ++(*this);
    return tmp;
}// operator++

inline LNodeIndex&
LNodeIndexIterator::operator*() const
{
    return getLNodeIndex();
}// operator*

inline LNodeIndex&
LNodeIndexIterator::getLNodeIndex() const
{
    return **d_node_it;
}// getLNodeIndex

inline const SAMRAI::hier::Index<NDIM>&
LNodeIndexIterator::getCellIndex() const
{
    return d_index_it.getIndex();
}// getCellIndex

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
