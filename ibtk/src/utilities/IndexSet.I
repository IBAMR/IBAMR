// Filename: IndexSet.I
// Created on 30 Apr 2010 by Boyce Griffith
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

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IndexSet::IndexSet()
    : index_set()
{
    // intentionally blank
    return;
}// IndexSet

inline
IndexSet::IndexSet(
    const IndexSet& from)
    : index_set(from.index_set)
{
    // intentionally blank
    return;
}// IndexSet

inline
IndexSet::~IndexSet()
{
    // intentionally blank
    return;
}// ~IndexSet

inline IndexSet&
IndexSet::operator=(
    const IndexSet& that)
{
    if (this == &that) return *this;  // check for self-assignment
    index_set = that.index_set;
    return *this;
}// operator=

inline void
IndexSet::copySourceItem(
    const SAMRAI::hier::Index<NDIM>& src_index,
    const SAMRAI::hier::IntVector<NDIM>& src_offset,
    const IndexSet& src_item)
{
    index_set = src_item.index_set;
    return;
}// copySourceItem

inline size_t
IndexSet::getDataStreamSize() const
{
    return SAMRAI::tbox::AbstractStream::sizeofInt()*(index_set.size()+1);
}// getDataStreamSize

inline void
IndexSet::packStream(
    SAMRAI::tbox::AbstractStream& stream) const
{
    stream << int(index_set.size());
    if (!index_set.empty())
    {
        std::vector<int> index_set_vec(index_set.begin(), index_set.end());
        stream.pack(&index_set_vec[0],index_set.size());
    }
    return;
}// packStream

inline void
IndexSet::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int size;
    stream >> size;
    if (size > 0)
    {
        std::vector<int> index_set_vec(size);
        stream.unpack(&index_set_vec[0],size);
        index_set.clear();
        index_set.insert(index_set_vec.begin(), index_set_vec.end());
    }
    return;
}// unpackStream

inline void
IndexSet::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database) const
{
    database->putInteger("index_set.size()",index_set.size());
    if (!index_set.empty())
    {
        std::vector<int> index_set_vec(index_set.begin(), index_set.end());
        database->putIntegerArray("index_set_vec",&index_set_vec[0],index_set.size());
    }
    return;
}// putToDatabase

inline void
IndexSet::getFromDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database)
{
    int size = database->getInteger("index_set.size()");
    if (size > 0)
    {
        std::vector<int> index_set_vec(size);
        database->getIntegerArray("index_set_vec",&index_set_vec[0],size);
        index_set.clear();
        index_set.insert(index_set_vec.begin(), index_set_vec.end());
    }
    return;
}// getFromDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
