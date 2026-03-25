// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_INSAveragingTurbulenceStatistics
#define included_IBAMR_INSAveragingTurbulenceStatistics

#include <ibamr/config.h>

#include <ibamr/INSTurbulenceStatistics.h>

#include <ibtk/HierarchyAveragedDataManager.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <CellVariable.h>
#include <NodeVariable.h>
#include <SideVariable.h>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class GridGeometry;
} // namespace hier
} // namespace SAMRAI

namespace IBAMR
{
/*!
 * \brief Tracks mean velocity and mean velocity products for staggered INS
 * simulations.
 *
 * This class uses HierarchyAveragedDataManager for both the velocity and its
 * symmetric second moment. The retained SnapshotCache objects can then be used
 * to reconstruct Reynolds-stress snapshots on a hierarchy.
 */
class INSAveragingTurbulenceStatistics : public INSTurbulenceStatistics
{
public:
    enum class AnalysisCentering
    {
        CELL,
        NODE
    };

    /*!
     * \brief Constructor.
     *
     * Recognized input options include the averaging controls consumed by
     * HierarchyAveragedDataManager and `statistics_start_time`, which delays
     * accumulation until the specified simulation time.
     */
    INSAveragingTurbulenceStatistics(std::string object_name,
                                     SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>>,
                                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::GridGeometry<NDIM>> grid_geom,
                                     bool register_for_restart = true);

    /*!
     * \brief Update the stored first- and second-moment statistics from the
     * supplied staggered velocity field.
     */
    bool updateStatistics(int U_idx,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> U_var,
                          const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                          double data_time,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                          SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops) override;

    /*!
     * \brief Return whether all tracked averaging snapshots are at periodic
     * steady state.
     */
    bool isAtSteadyState() const override;

    /*!
     * \brief Return the centering used for analysis data, either `CELL` or
     * `NODE`.
     */
    const std::string& getAnalysisCentering() const;

    /*!
     * \brief Access the manager storing the averaged velocity field.
     */
    IBTK::HierarchyAveragedDataManager& getAveragedVelocityManager();

    /*!
     * \brief Const access to the manager storing the averaged velocity field.
     */
    const IBTK::HierarchyAveragedDataManager& getAveragedVelocityManager() const;

    /*!
     * \brief Access the manager storing the averaged velocity products.
     */
    IBTK::HierarchyAveragedDataManager& getAveragedVelocityProductManager();

    /*!
     * \brief Const access to the manager storing the averaged velocity
     * products.
     */
    const IBTK::HierarchyAveragedDataManager& getAveragedVelocityProductManager() const;

    /*!
     * Fill a Reynolds-stress snapshot corresponding to one of the stored
     * averaging time points.
     *
     * The destination cell variable must have depth NDIM*(NDIM+1)/2 or
     * NDIM*NDIM.
     */
    void fillReynoldsStressSnapshot(int R_idx,
                                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> R_var,
                                    double time,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    double tol = 1.0e-8) const;

    /*!
     * Fill a node-centered Reynolds-stress snapshot corresponding to one of
     * the stored averaging time points.
     *
     * The destination node variable must have depth NDIM*(NDIM+1)/2 or
     * NDIM*NDIM.
     */
    void fillReynoldsStressSnapshot(int R_idx,
                                    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> R_var,
                                    double time,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    double tol = 1.0e-8) const;

private:
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_velocity_side_scratch_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_velocity_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_velocity_node_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_velocity_product_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_velocity_product_node_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_velocity_mean_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_velocity_mean_node_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_velocity_product_mean_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_velocity_product_mean_node_var;

    int d_velocity_side_scratch_idx = IBTK::invalid_index;
    int d_velocity_cell_idx = IBTK::invalid_index;
    int d_velocity_node_idx = IBTK::invalid_index;
    int d_velocity_product_cell_idx = IBTK::invalid_index;
    int d_velocity_product_node_idx = IBTK::invalid_index;
    int d_velocity_mean_cell_idx = IBTK::invalid_index;
    int d_velocity_mean_node_idx = IBTK::invalid_index;
    int d_velocity_product_mean_cell_idx = IBTK::invalid_index;
    int d_velocity_product_mean_node_idx = IBTK::invalid_index;

    /*!
     * Refinement operator used when filling recovered snapshots on a hierarchy.
     */
    std::string d_refine_type = "CONSERVATIVE_LINEAR_REFINE";

    AnalysisCentering d_analysis_centering = AnalysisCentering::CELL;

    /*!
     * Averaging manager for the first moment \f$\langle U \rangle\f$.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyAveragedDataManager> d_velocity_average_manager;
    SAMRAI::tbox::Pointer<IBTK::HierarchyAveragedDataManager> d_velocity_product_average_manager;
};
} // namespace IBAMR

#endif
