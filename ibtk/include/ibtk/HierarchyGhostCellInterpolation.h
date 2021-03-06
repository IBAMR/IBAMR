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

#ifndef included_IBTK_HierarchyGhostCellInterpolation
#define included_IBTK_HierarchyGhostCellInterpolation

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/ibtk_utilities.h"

#include "BoxGeometryFillPattern.h"
#include "CartesianGridGeometry.h"
#include "CoarseFineBoundary.h"
#include "CoarsenAlgorithm.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "RefineAlgorithm.h"
#include "VariableFillPattern.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace IBTK
{
class CartCellRobinPhysBdryOp;
class CartExtrapPhysBdryOp;
class CartSideRobinPhysBdryOp;
} // namespace IBTK
namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace xfer
{
template <int DIM>
class CoarsenPatchStrategy;
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefinePatchStrategy;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class HierarchyGhostCellInterpolation encapsulates the operations
 * required to set ghost cell values at physical and coarse-fine boundaries
 * across a range of levels of a locally refined patch hierarchy.
 *
 * \note In cases where physical boundary conditions are set via extrapolation
 * from interior values, setting ghost cell values may require both coarsening
 * and refining.
 */
class HierarchyGhostCellInterpolation : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Class
     * HierarchyGhostCellInterpolation::InterpolationTransactionComponent
     * encapsulates options for filling ghost cell values via class
     * HierarchyGhostCellInterpolation.
     */
    class InterpolationTransactionComponent
    {
    public:
        friend class HierarchyGhostCellInterpolation;

        /*!
         * \brief Default constructor.
         */
        inline InterpolationTransactionComponent(
            int data_idx = -1,
            const std::string& refine_op_name = "NONE",
            bool use_cf_bdry_interpolation = false,
            const std::string& coarsen_op_name = "NONE",
            const std::string& phys_bdry_extrap_type = "NONE",
            bool consistent_type_2_bdry = false,
            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* robin_bc_coef = NULL,
            SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > fill_pattern = NULL,
            const std::string& phys_bdry_type = "LINEAR")
            : d_dst_data_idx(data_idx),
              d_src_data_idx(data_idx),
              d_refine_op_name(refine_op_name),
              d_use_cf_bdry_interpolation(use_cf_bdry_interpolation),
              d_coarsen_op_name(coarsen_op_name),
              d_phys_bdry_extrap_type(phys_bdry_extrap_type),
              d_consistent_type_2_bdry(consistent_type_2_bdry),
              d_robin_bc_coefs(robin_bc_coef ? std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(1, robin_bc_coef) :
                                               std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>()),
              d_fill_pattern(fill_pattern ? fill_pattern : new SAMRAI::xfer::BoxGeometryFillPattern<NDIM>()),
              d_phys_bdry_type(phys_bdry_type)
        {
            // intentionally blank
            return;
        } // InterpolationTransactionComponent

        /*!
         * \brief Alternate constructor.
         */
        inline InterpolationTransactionComponent(
            int data_idx,
            const std::string& refine_op_name,
            bool use_cf_bdry_interpolation,
            const std::string& coarsen_op_name,
            const std::string& phys_bdry_extrap_type,
            bool consistent_type_2_bdry,
            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& robin_bc_coefs,
            SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > fill_pattern = NULL,
            const std::string& phys_bdry_type = "LINEAR")
            : d_dst_data_idx(data_idx),
              d_src_data_idx(data_idx),
              d_refine_op_name(refine_op_name),
              d_use_cf_bdry_interpolation(use_cf_bdry_interpolation),
              d_coarsen_op_name(coarsen_op_name),
              d_phys_bdry_extrap_type(phys_bdry_extrap_type),
              d_consistent_type_2_bdry(consistent_type_2_bdry),
              d_robin_bc_coefs(robin_bc_coefs),
              d_fill_pattern(fill_pattern ? fill_pattern : new SAMRAI::xfer::BoxGeometryFillPattern<NDIM>()),
              d_phys_bdry_type(phys_bdry_type)
        {
            // intentionally blank
            return;
        } // InterpolationTransactionComponent

        /*!
         * \brief Alternate constructor.
         */
        inline InterpolationTransactionComponent(
            int dst_data_idx,
            int src_data_idx,
            const std::string& refine_op_name,
            bool use_cf_bdry_interpolation,
            const std::string& coarsen_op_name,
            const std::string& phys_bdry_extrap_type,
            bool consistent_type_2_bdry,
            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* robin_bc_coef,
            SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > fill_pattern = NULL,
            const std::string& phys_bdry_type = "LINEAR")
            : d_dst_data_idx(dst_data_idx),
              d_src_data_idx(src_data_idx),
              d_refine_op_name(refine_op_name),
              d_use_cf_bdry_interpolation(use_cf_bdry_interpolation),
              d_coarsen_op_name(coarsen_op_name),
              d_phys_bdry_extrap_type(phys_bdry_extrap_type),
              d_consistent_type_2_bdry(consistent_type_2_bdry),
              d_robin_bc_coefs(robin_bc_coef ? std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(1, robin_bc_coef) :
                                               std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>()),
              d_fill_pattern(fill_pattern ? fill_pattern : new SAMRAI::xfer::BoxGeometryFillPattern<NDIM>()),
              d_phys_bdry_type(phys_bdry_type)
        {
            // intentionally blank
            return;
        } // InterpolationTransactionComponent

        /*!
         * \brief Alternate constructor.
         */
        inline InterpolationTransactionComponent(
            int dst_data_idx,
            int src_data_idx,
            const std::string& refine_op_name,
            bool use_cf_bdry_interpolation,
            const std::string& coarsen_op_name,
            const std::string& phys_bdry_extrap_type,
            bool consistent_type_2_bdry,
            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& robin_bc_coefs,
            SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > fill_pattern = NULL,
            const std::string& phys_bdry_type = "LINEAR")
            : d_dst_data_idx(dst_data_idx),
              d_src_data_idx(src_data_idx),
              d_refine_op_name(refine_op_name),
              d_use_cf_bdry_interpolation(use_cf_bdry_interpolation),
              d_coarsen_op_name(coarsen_op_name),
              d_phys_bdry_extrap_type(phys_bdry_extrap_type),
              d_consistent_type_2_bdry(consistent_type_2_bdry),
              d_robin_bc_coefs(robin_bc_coefs),
              d_fill_pattern(fill_pattern ? fill_pattern : new SAMRAI::xfer::BoxGeometryFillPattern<NDIM>()),
              d_phys_bdry_type(phys_bdry_type)
        {
            // intentionally blank
            return;
        } // InterpolationTransactionComponent

        /*!
         * \brief Copy constructor.
         *
         * \param from The value to copy to this object.
         */
        inline InterpolationTransactionComponent(const InterpolationTransactionComponent& from)
            : d_dst_data_idx(from.d_dst_data_idx),
              d_src_data_idx(from.d_src_data_idx),
              d_refine_op_name(from.d_refine_op_name),
              d_use_cf_bdry_interpolation(from.d_use_cf_bdry_interpolation),
              d_coarsen_op_name(from.d_coarsen_op_name),
              d_phys_bdry_extrap_type(from.d_phys_bdry_extrap_type),
              d_consistent_type_2_bdry(from.d_consistent_type_2_bdry),
              d_robin_bc_coefs(from.d_robin_bc_coefs),
              d_fill_pattern(from.d_fill_pattern),
              d_phys_bdry_type(from.d_phys_bdry_type)
        {
            // intentionally blank
            return;
        } // InterpolationTransactionComponent

        /*!
         * \brief Assignment operator.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        inline InterpolationTransactionComponent& operator=(const InterpolationTransactionComponent& that)
        {
            if (this != &that)
            {
                d_dst_data_idx = that.d_dst_data_idx;
                d_src_data_idx = that.d_src_data_idx;
                d_refine_op_name = that.d_refine_op_name;
                d_use_cf_bdry_interpolation = that.d_use_cf_bdry_interpolation;
                d_coarsen_op_name = that.d_coarsen_op_name;
                d_phys_bdry_extrap_type = that.d_phys_bdry_extrap_type;
                d_consistent_type_2_bdry = that.d_consistent_type_2_bdry;
                d_robin_bc_coefs = that.d_robin_bc_coefs;
                d_fill_pattern = that.d_fill_pattern;
                d_phys_bdry_type = that.d_phys_bdry_type;
            }
            return *this;
        } // operator=

        /*!
         * \brief Destructor.
         */
        inline ~InterpolationTransactionComponent()
        {
            // intentionally blank
            return;
        } // ~InterpolationTransactionComponent

        // Data.
        int d_dst_data_idx, d_src_data_idx;
        std::string d_refine_op_name;
        bool d_use_cf_bdry_interpolation;
        std::string d_coarsen_op_name;
        std::string d_phys_bdry_extrap_type;
        bool d_consistent_type_2_bdry;
        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_robin_bc_coefs;
        SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_fill_pattern;
        std::string d_phys_bdry_type;
    };

    /*!
     * \brief Default constructor.
     */
    HierarchyGhostCellInterpolation();

    /*!
     * \brief Destructor.
     */
    virtual ~HierarchyGhostCellInterpolation();

    /*!
     * \brief Specify whether the boundary conditions are homogeneous.
     */
    void setHomogeneousBc(bool homogeneous_bc);

    /*!
     * \brief Setup the hierarchy ghost cell interpolation operator to perform
     * the specified interpolation transactions on the specified patch
     * hierarchy.
     */
    void initializeOperatorState(InterpolationTransactionComponent transaction_comp,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1);

    /*!
     * \brief Setup the hierarchy ghost cell interpolation operator to perform
     * the specified collection of interpolation transactions on the specified
     * patch hierarchy.
     */
    void initializeOperatorState(const std::vector<InterpolationTransactionComponent>& transaction_comps,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1);

    /*!
     * \brief Reset transaction component with the interpolation operator.
     */
    void resetTransactionComponent(const InterpolationTransactionComponent& transaction_comps);

    /*!
     * \brief Reset transaction components with the interpolation operator.
     */
    void resetTransactionComponents(const std::vector<InterpolationTransactionComponent>& transaction_comps);

    /*!
     * \brief Reinitialize operator state following, e.g., a regridding operation.
     */
    void reinitializeOperatorState(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Clear all cached data.
     */
    void deallocateOperatorState();

    /*!
     * \brief Fill coarse-fine boundary and physical boundary ghost cells on all
     * levels of the patch hierarchy.
     */
    void fillData(double fill_time);

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    HierarchyGhostCellInterpolation(const HierarchyGhostCellInterpolation& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    HierarchyGhostCellInterpolation& operator=(const HierarchyGhostCellInterpolation& that) = delete;

    // Boolean indicating whether the operator is initialized.
    bool d_is_initialized = false;

    // Boolean indicating whether the operator should use homogeneous Robin
    // boundary conditions (when applicable).
    bool d_homogeneous_bc = false;

    // The component interpolation operations to perform.
    std::vector<InterpolationTransactionComponent> d_transaction_comps;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Cached communications algorithms and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_coarsen_alg;
    std::unique_ptr<SAMRAI::xfer::CoarsenPatchStrategy<NDIM> > d_coarsen_strategy;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_coarsen_scheds;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_refine_alg;
    std::unique_ptr<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_refine_strategy;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_refine_scheds;

    // Cached coarse-fine boundary and physical boundary condition handlers.
    std::vector<SAMRAI::tbox::Pointer<CoarseFineBoundaryRefinePatchStrategy> > d_cf_bdry_ops;
    std::vector<SAMRAI::tbox::Pointer<CartExtrapPhysBdryOp> > d_extrap_bc_ops;
    std::vector<SAMRAI::tbox::Pointer<CartCellRobinPhysBdryOp> > d_cc_robin_bc_ops;
    std::vector<SAMRAI::tbox::Pointer<CartSideRobinPhysBdryOp> > d_sc_robin_bc_ops;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_HierarchyGhostCellInterpolation
