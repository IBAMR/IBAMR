// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2026 by the IBAMR developers
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

#include <ibtk/CoarseFineBoundaryRefinePatchStrategy.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

#include <BoxGeometryFillPattern.h>
#include <CartesianGridGeometry.h>
#include <CoarseFineBoundary.h>
#include <CoarsenAlgorithm.h>
#include <IntVector.h>
#include <PatchHierarchy.h>
#include <RefineAlgorithm.h>
#include <VariableFillPattern.h>

#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace IBTK
{
class CartCellDoubleCubicCoarsen;
class CartCellRobinPhysBdryOp;
class CartExtrapPhysBdryOp;
class CartSideDoubleCubicCoarsen;
class CartSideRobinPhysBdryOp;
class RobinPhysBdryPatchStrategy;
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
class CoarsenOperator;
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
 *
 * Coarse-fine interpolation is third-order accurate within one coarse cell of a physical boundary if the coarse data
 * are extended beyond it with the quadratic Robin ghost value of a component with Robin coefficients, or by quadratic
 * extrapolation if the component has none and all destination data have at least two ghost cells. Otherwise, as for a
 * linear extension (d_extend_coarse_data_linearly), the interpolation is second-order accurate there.
 *
 * \note A component with the coarsen operator "CUBIC_COARSEN" and the refine operator "NONE" is coarsened cubically on
 * every axis if every refinement ratio is even. The ghost values of the source data may be invalid. The covered coarse
 * values are third-order accurate, except along an axis where the finer level is one coarse cell wide, and are
 * independent of the patch layout. If some ratio is 2, every patch must be a whole number of coarse cells along each
 * axis with ratio 2. If some ratio is odd and some ratio is less than 4, the component is averaged. If some ratio is
 * odd and every ratio is at least 4, the cubic coarsen operator stops with an error. A component with another refine
 * operator is coarsened by weighted averaging at ratios less than 4.
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
            int data_idx = invalid_index,
            const std::string& refine_op_name = "NONE",
            bool use_cf_bdry_interpolation = false,
            const std::string& coarsen_op_name = "NONE",
            const std::string& phys_bdry_extrap_type = "NONE",
            bool consistent_type_2_bdry = false,
            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* robin_bc_coef = nullptr,
            SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern = nullptr,
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
            SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern = nullptr,
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
            SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern = nullptr,
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
            SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern = nullptr,
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
              d_phys_bdry_type(from.d_phys_bdry_type),
              d_sc_phys_bdry_op(from.d_sc_phys_bdry_op),
              d_extend_coarse_data_linearly(from.d_extend_coarse_data_linearly)
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
                d_sc_phys_bdry_op = that.d_sc_phys_bdry_op;
                d_extend_coarse_data_linearly = that.d_extend_coarse_data_linearly;
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
        SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> d_fill_pattern;
        std::string d_phys_bdry_type;

        /*!
         * Operator, owned by another object, that sets the physical boundary values of side-centered data in the last
         * step of fillData(). If null (the default), fillData() uses the Robin boundary conditions. fillData() sets the
         * patch data index of the operator to d_dst_data_idx and its homogeneous flag to that of the
         * HierarchyGhostCellInterpolation object.
         */
        RobinPhysBdryPatchStrategy* d_sc_phys_bdry_op = nullptr;

        /*!
         * Whether to extend the coarse data of the coarse-fine interpolation linearly beyond physical boundaries
         * instead of with the Robin ghost value. Valid only for cell-centered data, such as the pressure, that do not
         * satisfy the condition of their Robin coefficients.
         */
        bool d_extend_coarse_data_linearly = false;
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
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 int coarsest_ln = invalid_level_number,
                                 int finest_ln = invalid_level_number);

    /*!
     * \brief Setup the hierarchy ghost cell interpolation operator to perform
     * the specified collection of interpolation transactions on the specified
     * patch hierarchy.
     */
    void initializeOperatorState(const std::vector<InterpolationTransactionComponent>& transaction_comps,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 int coarsest_ln = invalid_level_number,
                                 int finest_ln = invalid_level_number);

    /*!
     * \brief Reset transaction component with the interpolation operator.
     *
     * The new component must have its destination as its source if the component given to initializeOperatorState()
     * did.
     */
    void resetTransactionComponent(const InterpolationTransactionComponent& transaction_comps);

    /*!
     * \brief Reset transaction components with the interpolation operator.
     *
     * Each new component must have its destination as its source if the corresponding component given to
     * initializeOperatorState() did.
     */
    void resetTransactionComponents(const std::vector<InterpolationTransactionComponent>& transaction_comps);

    /*!
     * \brief Reinitialize operator state following, e.g., a regridding operation.
     */
    void reinitializeOperatorState(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

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

    // Return the coarsen operator of the component for the variable. If componentUsesCubicScratch() returns true, the
    // function returns a "CUBIC_COARSEN" operator of this object, which uses the ghost layer of its source.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM>>
    getCoarsenOperator(const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>& var,
                       const InterpolationTransactionComponent& comp) const;

    // Return true if the component is coarsened cubically and has the refine operator "NONE".
    bool componentUsesCubicScratch(const InterpolationTransactionComponent& comp) const;

    // Register the scratch data and the algorithms for the cubic coarsening of levels whose refinement ratios are all
    // even and include a ratio of 2.
    void setUpCubicScratch();

    // Return true if cubic coarsening of the given level uses the ghost layer of its source.
    bool levelUsesCubicScratch(int ln) const;

    // Allocate the scratch data on the level, fill them from the source data, and give the level to the cubic coarsen
    // operators. The second function releases the level and deallocates the scratch data.
    void fillCubicScratch(int ln, double fill_time);
    void deallocateCubicScratch(int ln);

    // Boolean indicating whether the operator is initialized.
    bool d_is_initialized = false;

    // Boolean indicating whether the operator should use homogeneous Robin
    // boundary conditions (when applicable).
    bool d_homogeneous_bc = false;

    // The component interpolation operations to perform.
    std::vector<InterpolationTransactionComponent> d_transaction_comps;

    // Whether the refine schedules of each component were built with its source as its destination.
    std::vector<bool> d_refine_scheds_built_in_place;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM>> d_grid_geom;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Cached communications algorithms and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM>> d_coarsen_alg;
    std::unique_ptr<SAMRAI::xfer::CoarsenPatchStrategy<NDIM>> d_coarsen_strategy;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM>>> d_coarsen_scheds;

    // The cubic coarsen operators, and the scratch data, algorithms, and schedules for coarsening a level with ratio 2:
    // a copy of the source data with a one-cell ghost layer that contains the values of the level and NaN elsewhere.
    SAMRAI::tbox::Pointer<CartCellDoubleCubicCoarsen> d_cc_cubic_coarsen_op;
    SAMRAI::tbox::Pointer<CartSideDoubleCubicCoarsen> d_sc_cubic_coarsen_op;
    std::vector<int> d_cubic_scratch_idxs;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM>> d_cubic_scratch_coarsen_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM>> d_cubic_scratch_fill_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>> d_cubic_scratch_fill_scheds;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM>> d_refine_alg;
    std::unique_ptr<SAMRAI::xfer::RefinePatchStrategy<NDIM>> d_refine_strategy;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>> d_refine_scheds;

    // Extension of the temporary coarse data beyond physical boundaries for the coarse-fine interpolation.
    class CoarseDataExtension;
    std::unique_ptr<CoarseDataExtension> d_coarse_extension;

    // Cached coarse-fine boundary and physical boundary condition handlers.
    std::vector<SAMRAI::tbox::Pointer<CoarseFineBoundaryRefinePatchStrategy>> d_cf_bdry_ops;
    std::vector<SAMRAI::tbox::Pointer<CartExtrapPhysBdryOp>> d_extrap_bc_ops;
    std::vector<SAMRAI::tbox::Pointer<CartCellRobinPhysBdryOp>> d_cc_robin_bc_ops;
    std::vector<SAMRAI::tbox::Pointer<CartSideRobinPhysBdryOp>> d_sc_robin_bc_ops;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_HierarchyGhostCellInterpolation
