// Filename: EdgeDataSynchronization.h
// Created on 02 Feb 2011 by Boyce Griffith
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

#ifndef included_EdgeDataSynchronization
#define included_EdgeDataSynchronization

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "boost/array.hpp"



namespace SAMRAI
{
namespace xfer
{

class CoarsenSchedule;

class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class EdgeDataSynchronization encapsulates the operations required to
 * "synchronize" edge-centered values defined at patch boundaries.
 */
class EdgeDataSynchronization
{
public:
    /*!
     * \brief Class EdgeDataSynchronization::SynchronizationTransactionComponent
     * encapsulates options for filling ghost cell values via class
     * EdgeDataSynchronization.
     */
    class SynchronizationTransactionComponent
    {
    public:
        friend class EdgeDataSynchronization;

        /*!
         * \brief Default constructor.
         */
        inline SynchronizationTransactionComponent(int data_idx = -1, const std::string& coarsen_op_name = "NONE")
            : d_data_idx(data_idx), d_coarsen_op_name(coarsen_op_name)
        {
            // intentionally blank
            return;
        }

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
        }

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
        }

        /*!
         * \brief Destructor.
         */
        inline ~SynchronizationTransactionComponent()
        {
            // intentionally blank
            return;
        }

    private:
        // Data.
        int d_data_idx;
        std::string d_coarsen_op_name;
    };

    /*!
     * \brief Default constructor.
     */
    EdgeDataSynchronization();

    /*!
     * \brief Destructor.
     */
    ~EdgeDataSynchronization();

    /*!
     * \brief Setup the hierarchy synchronization operator to perform the
     * specified synchronization transactions on the specified patch hierarchy.
     */
    void initializeOperatorState(const SynchronizationTransactionComponent& transaction_comp,
                                 boost::shared_ptr<SAMRAI::hier::PatchHierarchy > hierarchy);

    /*!
     * \brief Setup the hierarchy synchronization operator to perform the
     * specified collection of synchronization transactions on the specified
     * patch hierarchy.
     */
    void initializeOperatorState(const std::vector<SynchronizationTransactionComponent>& transaction_comps,
                                 boost::shared_ptr<SAMRAI::hier::PatchHierarchy > hierarchy);

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
    EdgeDataSynchronization(const EdgeDataSynchronization& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    EdgeDataSynchronization& operator=(const EdgeDataSynchronization& that);

    // Boolean indicating whether the operator is initialized.
    bool d_is_initialized;

    // The component synchronization operations to perform.
    std::vector<SynchronizationTransactionComponent> d_transaction_comps;

    // Hierarchy configuration.
    boost::shared_ptr<SAMRAI::hier::PatchHierarchy > d_hierarchy;
    boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry > d_grid_geom;
    int d_coarsest_ln, d_finest_ln;

    // Cached communications algorithms and schedules.
    boost::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm > d_coarsen_alg;
    std::vector<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule > > d_coarsen_scheds;

    boost::array<boost::shared_ptr<SAMRAI::xfer::RefineAlgorithm >, NDIM> d_refine_alg;
    boost::array<std::vector<boost::shared_ptr<SAMRAI::xfer::RefineSchedule > >, NDIM> d_refine_scheds;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_EdgeDataSynchronization
