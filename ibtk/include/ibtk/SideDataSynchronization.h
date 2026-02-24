// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
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

#ifndef included_IBTK_SideDataSynchronization
#define included_IBTK_SideDataSynchronization

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/ibtk_utilities.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICoarsenAlgorithm.h"
#include "SAMRAICoarsenSchedule.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRefineAlgorithm.h"
#include "SAMRAIRefineSchedule.h"
#include "tbox/DescribedClass.h"

#include <string>
#include <vector>

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
 * \brief Class SideDataSynchronization encapsulates the operations required to
 * "synchronize" side-centered values defined at patch boundaries.
 */
class SideDataSynchronization : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Class SideDataSynchronization::SynchronizationTransactionComponent
     * encapsulates options for filling ghost cell values via class
     * SideDataSynchronization.
     */
    class SynchronizationTransactionComponent
    {
    public:
        friend class SideDataSynchronization;

        /*!
         * \brief Default constructor.
         */
        inline SynchronizationTransactionComponent(int data_idx = invalid_index,
                                                   const std::string& coarsen_op_name = "NONE")
            : d_data_idx(data_idx), d_coarsen_op_name(coarsen_op_name)
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
    SideDataSynchronization() = default;

    /*!
     * \brief Destructor.
     */
    virtual ~SideDataSynchronization();

    /*!
     * \brief Setup the hierarchy synchronization operator to perform the
     * specified synchronization transactions on the specified patch hierarchy.
     */
    void initializeOperatorState(const SynchronizationTransactionComponent& transaction_comp,
                                 SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy);

    /*!
     * \brief Setup the hierarchy synchronization operator to perform the
     * specified collection of synchronization transactions on the specified
     * patch hierarchy.
     */
    void initializeOperatorState(const std::vector<SynchronizationTransactionComponent>& transaction_comps,
                                 SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy);

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
    SideDataSynchronization(const SideDataSynchronization& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SideDataSynchronization& operator=(const SideDataSynchronization& that) = delete;

    // Boolean indicating whether the operator is initialized.
    bool d_is_initialized = false;

    // The component synchronization operations to perform.
    std::vector<SynchronizationTransactionComponent> d_transaction_comps;

    // Hierarchy configuration.
    SAMRAIPointer<SAMRAIPatchHierarchy> d_hierarchy;
    SAMRAIPointer<SAMRAICartesianGridGeometry> d_grid_geom;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Cached communications algorithms and schedules.
    SAMRAIPointer<SAMRAICoarsenAlgorithm> d_coarsen_alg;
    std::vector<SAMRAIPointer<SAMRAICoarsenSchedule> > d_coarsen_scheds;

    SAMRAIPointer<SAMRAIRefineAlgorithm> d_refine_alg;
    std::vector<SAMRAIPointer<SAMRAIRefineSchedule> > d_refine_scheds;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SideDataSynchronization
