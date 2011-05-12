// Filename: IBDataPostProcessor.h
// Created on 24 Sep 2008 by Boyce Griffith
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

#ifndef included_IBDataPostProcessor
#define included_IBDataPostProcessor

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>

// SAMRAI INCLUDES
#include <PatchHierarchy.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBDataPostProcessor provides a generic interface for specifying
 * post-processing code for use in an IB computation.
 */
class IBDataPostProcessor
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IBDataPostProcessor();

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBDataPostProcessor();

    /*!
     * \brief Setup the data needed to post-process the data on the specified
     * level of the patch hierarchy.
     *
     * \note A default empty implementation is provided.
     */
    virtual void
    initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool initial_time,
        IBTK::LDataManager* const lag_manager);

    /*!
     * \brief Post-process data on the patch hierarchy.
     */
    virtual void
    postProcessData(
        const int u_idx,
        const int p_idx,
        const int f_idx,
        std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > F_data,
        std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > X_data,
        std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > U_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level_number,
        const int finest_level_number,
        const double data_time,
        IBTK::LDataManager* const lag_manager) = 0;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBDataPostProcessor(
        const IBDataPostProcessor& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBDataPostProcessor&
    operator=(
        const IBDataPostProcessor& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBDataPostProcessor.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBDataPostProcessor
