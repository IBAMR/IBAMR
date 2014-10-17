// Filename: ParallelSet.h
// Created on 03 Mar 2011 by Boyce Griffith
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

#ifndef included_ParallelSet
#define included_ParallelSet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <set>
#include <vector>

#include "tbox/DescribedClass.h"

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
    ParallelSet();

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    ParallelSet(const ParallelSet& from);

    /*!
     * \brief Destructor.
     */
    ~ParallelSet();

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

#endif //#ifndef included_ParallelSet
