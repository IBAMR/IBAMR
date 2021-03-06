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

#ifndef included_IBTK_ParallelSet
#define included_IBTK_ParallelSet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "tbox/DescribedClass.h"

#include <set>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class ParallelSet is a utility class for storing collections of
 * integer keys in parallel.
 */
class ParallelSet : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    ParallelSet() = default;

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    ParallelSet(const ParallelSet& from) = default;

    /*!
     * \brief Destructor.
     */
    virtual ~ParallelSet() = default;

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ParallelSet& operator=(const ParallelSet& that);

    /*!
     * \brief Add an item with the specified key to the set.
     *
     * \note This method is not collective (i.e., it does not have to be called
     * by all MPI tasks); however, it is necessary to call the collective
     * function ParallelSet::communicateData() to finalize all parallel
     * communication.
     *
     * \note The underling set data structure is \em not updated until the
     * collective method communicateData() is called, even for \em serial runs.
     */
    void addItem(int key);

    /*!
     * \brief Remove an item from the set.
     *
     * \note This method is not collective (i.e., it does not have to be called
     * by all MPI tasks); however, it is necessary to call the collective
     * function ParallelSet::communicateData() to finalize all parallel
     * communication.
     *
     * \note The underling set data structure is \em not updated until the
     * collective method communicateData() is called, even for \em serial runs.
     */
    void removeItem(int key);

    /*!
     * \brief Communicate data to (re-)initialize the set.
     */
    void communicateData();

    /*!
     * \brief Return a const reference to the set.
     */
    const std::set<int>& getSet() const;

private:
    // Member data.
    std::set<int> d_set;
    std::vector<int> d_pending_additions, d_pending_removals;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_ParallelSet
