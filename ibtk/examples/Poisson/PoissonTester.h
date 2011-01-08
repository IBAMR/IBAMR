// Filename: PoissonTester.h
// Created on 12 Feb 2005 by Boyce Griffith
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

#ifndef included_PoissonTester
#define included_PoissonTester

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <BasePatchHierarchy.h>
#include <BasePatchLevel.h>
#include <StandardTagAndInitStrategy.h>
#include <VariableContext.h>
#include <tbox/Pointer.h>

// NAMESPACE
using namespace SAMRAI;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Tester for Poisson preconditioners.
 */
class PoissonTester
    : public mesh::StandardTagAndInitStrategy<NDIM>
{
public:
    /*!
     * \brief Default constructor.
     */
    PoissonTester();

    /*!
     * \brief Destructor.
     */
    ~PoissonTester();

    int
    getUIndex()
        { return d_U_idx; }

    int
    getVIndex()
        { return d_V_idx; }

    int
    getFIndex()
        { return d_F_idx; }

    ///
    ///  The following routines:
    ///
    ///      initializeLevelData(),
    ///      resetHierarchyConfiguration()
    ///
    ///  are concrete implementations of functions declared in the
    ///  mesh::StandardTagAndInitStrategy abstract base class.
    ///

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     */
    void
    initializeLevelData(
        const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level=tbox::Pointer<hier::PatchLevel<NDIM> >(NULL),
        const bool allocate_data=true);

    /*!
     * Reset data after the hierarchy has changed (e.g., due to regridding) and
     * the data has been initialized on the new levels.
     */
    void
    resetHierarchyConfiguration(
        const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level);

protected:

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PoissonTester(
        const PoissonTester& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PoissonTester&
    operator=(
        const PoissonTester& that);

    /*
     * Patch data descriptor indices and data contexts.
     */
    tbox::Pointer<hier::VariableContext> d_context;
    int d_U_idx, d_V_idx, d_F_idx;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "PoissonTester.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PoissonTester
