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
 *
 * <h2>Parameters read from the input database</h2>
 * <ol>
 *   <li><code>analysis_centering</code>: centering used for the analyzed
 *     velocity and velocity-product data. Supported values are
 *     <code>CELL</code> and <code>NODE</code>. The default is
 *     <code>CELL</code>.</li>
 *   <li><code>refine_type</code>: refinement operator used when recovering
 *     stored snapshots onto a hierarchy. The default is
 *     <code>CONSERVATIVE_LINEAR_REFINE</code>.</li>
 *   <li><code>statistics_start_time</code>: simulation time at which
 *     statistics accumulation begins. The default is <code>0.0</code>.</li>
 *   <li><code>output_reynolds_stress</code>: whether to register Reynolds
 *     stress plot output on the shared VisIt data writer. The default is
 *     <code>false</code>.</li>
 *   <li><code>output_tke</code>: whether to register turbulent kinetic energy
 *     plot output on the shared VisIt data writer. The default is
 *     <code>false</code>.</li>
 *   <li><code>reynolds_stress_output_name</code>: plot-quantity name used for
 *     Reynolds stress output. The default is
 *     <code>Reynolds stress</code>.</li>
 *   <li><code>tke_output_name</code>: plot-quantity name used for turbulent
 *     kinetic energy output. The default is <code>TKE</code>.</li>
 *   <li><code>period_start</code>, <code>period_end</code>,
 *     <code>threshold</code>, <code>num_snapshots</code>,
 *     <code>enable_logging</code>, <code>output_data</code>,
 *     <code>dir_dump_name</code>, <code>gcw</code>, and manager-local
 *     <code>refine_type</code>: options forwarded to the underlying
 *     IBTK::HierarchyAveragedDataManager objects. These entries may be given
 *     directly in this database and overridden in the optional
 *     <code>VelocityAveraging</code> and
 *     <code>VelocityProductAveraging</code> sub-databases. When
 *     <code>VelocityProductAveraging</code> does not explicitly set
 *     <code>output_data</code>, it defaults to <code>false</code>.</li>
 * </ol>
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
     * \brief Register the shared VisIt data writer and any enabled derived
     * turbulence outputs.
     */
    void registerVisItDataWriter(SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM>> visit_writer) override;

    /*!
     * \brief Prepare any enabled Reynolds-stress/TKE plot data at the
     * requested output time.
     */
    void setupPlotData(double data_time,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                       SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops) override;

    /*!
     * \brief Deallocate any plot data owned by this object.
     */
    void deallocatePlotData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy) override;

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
     * Fill a Reynolds-stress snapshot corresponding to the requested time.
     *
     * The destination cell variable must have depth NDIM*(NDIM+1)/2 or
     * NDIM*NDIM.
     *
     * \param R_idx Patch-data index for the destination Reynolds-stress field.
     * \param R_var Destination cell-centered variable.
     * \param time Requested snapshot time.
     * \param hierarchy Patch hierarchy on which to fill data.
     * \param hier_math_ops Hierarchy math operations helper.
     * \param tol Time-point lookup tolerance used when \p time matches a
     * stored snapshot exactly. Otherwise, the result is linearly interpolated
     * in time from the surrounding stored snapshots.
     */
    void fillReynoldsStressSnapshot(int R_idx,
                                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> R_var,
                                    double time,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    double tol = 1.0e-8) const;

    /*!
     * Fill a node-centered Reynolds-stress snapshot corresponding to the
     * requested time.
     *
     * The destination node variable must have depth NDIM*(NDIM+1)/2 or
     * NDIM*NDIM.
     *
     * \param R_idx Patch-data index for the destination Reynolds-stress field.
     * \param R_var Destination node-centered variable.
     * \param time Requested snapshot time.
     * \param hierarchy Patch hierarchy on which to fill data.
     * \param hier_math_ops Hierarchy math operations helper.
     * \param tol Time-point lookup tolerance used when \p time matches a
     * stored snapshot exactly. Otherwise, the result is linearly interpolated
     * in time from the surrounding stored snapshots.
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
     * \param tol Time-point lookup tolerance used when \p time matches a
     * stored snapshot exactly. Otherwise, the result is linearly interpolated
     * in time from the surrounding stored snapshots.
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
     * \param tol Time-point lookup tolerance used when \p time matches a
     * stored snapshot exactly. Otherwise, the result is linearly interpolated
     * in time from the surrounding stored snapshots.
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
     * \param tol Time-point lookup tolerance used when \p time matches a
     * stored snapshot exactly. Otherwise, the result is linearly interpolated
     * in time from the surrounding stored snapshots.
     */
    void fillTurbulentKineticEnergySnapshot(int k_idx,
                                            SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> k_var,
                                            double time,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                            SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                            double tol = 1.0e-8) const;

private:
    bool haveStoredSnapshots() const;

    double mapToStoredTime(double time) const;

    void fillAveragedSnapshots(double time,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                               double tol) const;

    void zeroPlotData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy) const;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_velocity_side_scratch_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_velocity_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_velocity_node_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_velocity_product_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_velocity_product_node_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_velocity_mean_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_velocity_mean_node_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_velocity_product_mean_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_velocity_product_mean_node_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_reynolds_stress_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_reynolds_stress_node_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_tke_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double>> d_tke_node_var;
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM>> d_visit_writer;

    int d_velocity_side_scratch_idx = IBTK::invalid_index;
    int d_velocity_cell_idx = IBTK::invalid_index;
    int d_velocity_node_idx = IBTK::invalid_index;
    int d_velocity_product_cell_idx = IBTK::invalid_index;
    int d_velocity_product_node_idx = IBTK::invalid_index;
    int d_velocity_mean_cell_idx = IBTK::invalid_index;
    int d_velocity_mean_node_idx = IBTK::invalid_index;
    int d_velocity_product_mean_cell_idx = IBTK::invalid_index;
    int d_velocity_product_mean_node_idx = IBTK::invalid_index;
    int d_reynolds_stress_plot_idx = IBTK::invalid_index;
    int d_tke_plot_idx = IBTK::invalid_index;

    /*!
     * Refinement operator used when filling recovered snapshots on a hierarchy.
     */
    std::string d_refine_type = "CONSERVATIVE_LINEAR_REFINE";
    std::string d_reynolds_stress_output_name = "Reynolds stress";
    std::string d_tke_output_name = "TKE";

    double d_period_start = 0.0;
    double d_period_end = 0.0;
    double d_period_length = 0.0;

    DataCentering d_data_centering = DataCentering::CELL;
    bool d_has_statistics_samples = false;
    bool d_output_reynolds_stress = false;
    bool d_output_tke = false;
    bool d_plot_quantities_registered = false;

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
