// Filename: ParallelEdgeMap.h
// Created on 28 Jun 2010 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_ParallelEdgeMap
#define included_ParallelEdgeMap

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <utility>

#include "tbox/DescribedClass.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class ParallelEdgeMap is a utility class for managing edge maps (i.e.,
 * maps from vertices to links between vertices) in parallel.
 */
class ParallelEdgeMap : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    ParallelEdgeMap();

    /*!
     * \brief Destructor.
     */
    ~ParallelEdgeMap();

    /*!
     * \brief Add an edge to the edge map.
     *
     * \return The master node index to be associated with the edge in the edge
     * map.
     *
     * \note This method is not collective (i.e., it does not have to be called
     * by all MPI tasks); however, it is necessary to call the collective
     * function ParallelEdgeMap::communicateData() to finalize all parallel
     * communication.
     *
     * \note By default, the master index associated with each edge is the
     * vertex with minimum index in the link.
     */
    int addEdge(const std::pair<int, int>& link, int mastr_idx = -1);

    /*!
     * \brief Remove an edge from the edge map.
     *
     * \note This method is not collective (i.e., it does not have to be called
     * by all MPI tasks); however, it is necessary to call the collective
     * function ParallelEdgeMap::communicateData() to finalize all parallel
     * communication.
     *
     * \note The master index argument is optional and is only used as a hint to
     * attempt to find the link in the link table.
     */
    void removeEdge(const std::pair<int, int>& link, int mastr_idx = -1);

    /*!
     * \brief Communicate data to (re-)initialize the edge map.
     */
    void communicateData();

    /*!
     * \brief Return a const reference to the edge map.
     */
    const std::multimap<int, std::pair<int, int> >& getEdgeMap() const;

private:
    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     *
     * \note This constructor is not implemented and should not be used.
     */
    ParallelEdgeMap(const ParallelEdgeMap& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ParallelEdgeMap& operator=(const ParallelEdgeMap& that);

    // Member data.
    std::multimap<int, std::pair<int, int> > d_edge_map;
    std::multimap<int, std::pair<int, int> > d_pending_additions, d_pending_removals;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ParallelEdgeMap
