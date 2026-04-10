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
#include <ibamr/ibamr_enums.h>

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
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

namespace IBAMR
{
enum class DataCentering
{
    CELL,
    NODE
};

/*!
 * \brief Tracks mean velocity and mean velocity second moments for staggered
 * INS simulations.
 *
 * The retained statistics are
 * \f[
 * \langle u_i \rangle
 * \f]
 * and
 * \f[
 * \langle u_i u_j \rangle.
 * \f]
 *
 * This class uses HierarchyAveragedDataManager for both the velocity and its
 * symmetric second moment. The retained SnapshotCache objects can then be used
 * to reconstruct Reynolds-stress snapshots on a hierarchy.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim
 analysis_centering = "CELL"                // "CELL" or "NODE"
 refine_type = "CONSERVATIVE_LINEAR_REFINE" // refine operator used to recover snapshots
 statistics_start_time = 0.0                // delay accumulation until this simulation time
 // averaging options used by IBTK::HierarchyAveragedDataManager
 \endverbatim
 */
class INSAveragingTurbulenceStatistics : public INSTurbulenceStatistics
{
public:
    /*!
     * \brief Constructor.
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
    DataCentering getDataCentering() const;

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
     *
     * \param R_idx Patch-data index for the destination Reynolds-stress field.
     * \param R_var Destination cell-centered variable.
     * \param time Requested snapshot time.
     * \param hierarchy Patch hierarchy on which to fill data.
     * \param hier_math_ops Hierarchy math operations helper.
     * \param tol Time-point lookup tolerance used to match \p time to stored
     * snapshot times.
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
     *
     * \param R_idx Patch-data index for the destination Reynolds-stress field.
     * \param R_var Destination node-centered variable.
     * \param time Requested snapshot time.
     * \param hierarchy Patch hierarchy on which to fill data.
     * \param hier_math_ops Hierarchy math operations helper.
     * \param tol Time-point lookup tolerance used to match \p time to stored
     * snapshot times.
     */
    void fillReynoldsStressSnapshot(int R_idx,
                                    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> R_var,
                                    double time,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    double tol = 1.0e-8) const;

    /*!
     * Fill a Reynolds-stress snapshot into either cell- or node-centered data.
     *
     * This overload dispatches based on the runtime variable type.
     *
     * \param R_idx Patch-data index for the destination Reynolds-stress field.
     * \param R_var Destination variable (cell- or node-centered).
     * \param time Requested snapshot time.
     * \param hierarchy Patch hierarchy on which to fill data.
     * \param hier_math_ops Hierarchy math operations helper.
     * \param tol Time-point lookup tolerance used to match \p time to stored
     * snapshot times.
     */
    void fillReynoldsStressSnapshot(int R_idx,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> R_var,
                                    double time,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    double tol = 1.0e-8) const;

    /*!
     * Fill a TKE snapshot, \f$k=\frac{1}{2}R_{ii}\f$, into cell-centered data.
     *
     * \param k_idx Patch-data index for the destination TKE field.
     * \param k_var Destination cell-centered variable.
     * \param time Requested snapshot time.
     * \param hierarchy Patch hierarchy on which to fill data.
     * \param hier_math_ops Hierarchy math operations helper.
     * \param tol Time-point lookup tolerance used to match \p time to stored
     * snapshot times.
     */
    void fillTurbulentKineticEnergySnapshot(int k_idx,
                                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> k_var,
                                            double time,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                            SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                            double tol = 1.0e-8) const;

    /*!
     * Fill a TKE snapshot, \f$k=\frac{1}{2}R_{ii}\f$, into node-centered data.
     *
     * \param k_idx Patch-data index for the destination TKE field.
     * \param k_var Destination node-centered variable.
     * \param time Requested snapshot time.
     * \param hierarchy Patch hierarchy on which to fill data.
     * \param hier_math_ops Hierarchy math operations helper.
     * \param tol Time-point lookup tolerance used to match \p time to stored
     * snapshot times.
     */
    void fillTurbulentKineticEnergySnapshot(int k_idx,
                                            SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> k_var,
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

    DataCentering d_data_centering = DataCentering::CELL;

    /*!
     * Averaging managers for the first and second moments,
     * \f$\langle U \rangle\f$ and \f$\langle UU \rangle\f$.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyAveragedDataManager> d_velocity_average_manager;
    SAMRAI::tbox::Pointer<IBTK::HierarchyAveragedDataManager> d_velocity_product_average_manager;
};
} // namespace IBAMR

namespace IBAMR
{
template <>
inline DataCentering
string_to_enum<DataCentering>(const std::string& val)
{
    if (strcasecmp(val.c_str(), "CELL") == 0) return DataCentering::CELL;
    if (strcasecmp(val.c_str(), "NODE") == 0) return DataCentering::NODE;
    TBOX_ERROR("unsupported analysis centering value " << val << "\n");
    return DataCentering::CELL;
}

template <>
inline std::string
enum_to_string<DataCentering>(DataCentering val)
{
    switch (val)
    {
    case DataCentering::CELL:
        return "CELL";
    case DataCentering::NODE:
        return "NODE";
    default:
        TBOX_ERROR("unsupported analysis centering enum value\n");
    }
    return "UNKNOWN";
}
} // namespace IBAMR

#endif
