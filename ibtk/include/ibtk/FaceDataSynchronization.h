// Filename: FaceDataSynchronization.h
// Created on 03 Feb 2011 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#ifndef included_IBTK_FaceDataSynchronization
#define included_IBTK_FaceDataSynchronization

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <utility>
#include <vector>

#include "CartesianGridGeometry.h"
#include "CoarsenAlgorithm.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "RefineAlgorithm.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FaceDataSynchronization encapsulates the operations required to
 * "synchronize" face-centered values defined at patch boundaries.
 */
class FaceDataSynchronization : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Class FaceDataSynchronization::SynchronizationTransactionComponent
     * encapsulates options for filling ghost cell values via class
     * FaceDataSynchronization.
     */
    class SynchronizationTransactionComponent
    {
    public:
        friend class FaceDataSynchronization;

        /*!
         * \brief Default constructor.
         */
        inline SynchronizationTransactionComponent(int data_idx = -1, std::string coarsen_op_name = "NONE")
            : d_data_idx(data_idx), d_coarsen_op_name(std::move(coarsen_op_name))
        {
            // intentionally blank
            return;
        } // SynchronizationTransactionComponent

        /*!
         * \brief Copy constructor.
         *
         * \param from The value to copy to this object.
         */
        inline SynchronizationTransactionComponent(const SynchronizationTransactionComponent& from)
            : d_data_idx(from.d_data_idx), d_coarsen_op_name(from.d_coarsen_op_name)
        {
            // intentionally blank
            return;
        } // SynchronizationTransactionComponent

        /*!
         * \brief Assignment operator.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        inline SynchronizationTransactionComponent& operator=(const SynchronizationTransactionComponent& that)
        {
            if (this != &that)
            {
                d_data_idx = that.d_data_idx;
                d_coarsen_op_name = that.d_coarsen_op_name;
            }
            return *this;
        } // operator=

        /*!
         * \brief Destructor.
         */
        inline ~SynchronizationTransactionComponent()
        {
            // intentionally blank
            return;
        } // ~SynchronizationTransactionComponent

    private:
        // Data.
        int d_data_idx;
        std::string d_coarsen_op_name;
    };

    /*!
     * \brief Default constructor.
     */
    FaceDataSynchronization();

    /*!
     * \brief Destructor.
     */
    ~FaceDataSynchronization() override;

    /*!
     * \brief Setup the hierarchy synchronization operator to perform the
     * specified synchronization transactions on the specified patch hierarchy.
     */
    void initializeOperatorState(const SynchronizationTransactionComponent& transaction_comp,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Setup the hierarchy synchronization operator to perform the
     * specified collection of synchronization transactions on the specified
     * patch hierarchy.
     */
    void initializeOperatorState(const std::vector<SynchronizationTransactionComponent>& transaction_comps,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Reset transaction component with the synchronization operator.
     */
    void resetTransactionComponent(const SynchronizationTransactionComponent& transaction_comps);

    /*!
     * \brief Reset transaction components with the synchronization operator.
     */
    void resetTransactionComponents(const std::vector<SynchronizationTransactionComponent>& transaction_comps);

    /*!
     * \brief Clear all cached data.
     */
    void deallocateOperatorState();

    /*!
     * \brief Synchronize the data on all levels of the patch hierarchy.
     */
    void synchronizeData(double fill_time);

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FaceDataSynchronization(const FaceDataSynchronization& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FaceDataSynchronization& operator=(const FaceDataSynchronization& that) = delete;

    // Boolean indicating whether the operator is initialized.
    bool d_is_initialized = false;

    // The component synchronization operations to perform.
    std::vector<SynchronizationTransactionComponent> d_transaction_comps;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;
    int d_coarsest_ln = -1, d_finest_ln = -1;

    // Cached communications algorithms and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_coarsen_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_coarsen_scheds;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_refine_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_refine_scheds;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_FaceDataSynchronization
