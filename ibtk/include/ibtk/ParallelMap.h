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

#ifndef included_IBTK_ParallelMap
#define included_IBTK_ParallelMap

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <map>
#include <vector>

namespace IBTK
{
class Streamable;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class ParallelMap is a utility class for associating integer keys with
 * arbitrary data items in parallel.
 */
class ParallelMap : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    ParallelMap() = default;

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    ParallelMap(const ParallelMap& from) = default;

    /*!
     * \brief Destructor.
     */
    virtual ~ParallelMap() = default;

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ParallelMap& operator=(const ParallelMap& that);

    /*!
     * \brief Add an item with the specified key to the map.
     *
     * \note This method is not collective (i.e., it does not have to be called
     * by all MPI tasks); however, it is necessary to call the collective
     * function ParallelMap::communicateData() to finalize all parallel
     * communication.
     *
     * \note The underling map data structure is \em not updated until the
     * collective method communicateData() is called, even for \em serial runs.
     */
    void addItem(int key, SAMRAI::tbox::Pointer<Streamable> item);

    /*!
     * \brief Remove an item from the map.
     *
     * \note This method is not collective (i.e., it does not have to be called
     * by all MPI tasks); however, it is necessary to call the collective
     * function ParallelMap::communicateData() to finalize all parallel
     * communication.
     *
     * \note The underling map data structure is \em not updated until the
     * collective method communicateData() is called, even for \em serial runs.
     */
    void removeItem(int key);

    /*!
     * \brief Communicate data to (re-)initialize the map.
     */
    void communicateData();

    /*!
     * \brief Return a const reference to the map.
     */
    const std::map<int, SAMRAI::tbox::Pointer<Streamable> >& getMap() const;

private:
    // Member data.
    std::map<int, SAMRAI::tbox::Pointer<Streamable> > d_map;
    std::map<int, SAMRAI::tbox::Pointer<Streamable> > d_pending_additions;
    std::vector<int> d_pending_removals;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_ParallelMap
