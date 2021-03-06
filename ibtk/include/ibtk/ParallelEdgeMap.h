// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_ParallelEdgeMap
#define included_IBTK_ParallelEdgeMap

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "tbox/DescribedClass.h"

#include <map>
#include <utility>

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
    ParallelEdgeMap() = default;

    /*!
     * \brief Destructor.
     */
    virtual ~ParallelEdgeMap() = default;

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
    ParallelEdgeMap(const ParallelEdgeMap& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ParallelEdgeMap& operator=(const ParallelEdgeMap& that) = delete;

    // Member data.
    std::multimap<int, std::pair<int, int> > d_edge_map;
    std::multimap<int, std::pair<int, int> > d_pending_additions, d_pending_removals;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_ParallelEdgeMap
