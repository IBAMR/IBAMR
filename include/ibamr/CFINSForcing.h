// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
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

#ifndef included_CFINSForcing
#define included_CFINSForcing

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/CFStrategy.h"
#include "ibamr/CFUpperConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/muParserCartGridFunction.h"
#include "ibtk/muParserRobinBcCoefs.h"

#include "BasePatchHierarchy.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "FaceVariable.h"
#include "HierarchyDataOpsManager.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

#include <limits>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CFINSForcing provides an interface for specifying a viscoelastic stress to be added to the
 * Navier-Stokes equations. The class uses the advection diffusion integrator to update the viscoelastic stress.
 *
 * Users must register a strategy operator <code>CFStrategy</code> that can compute (1) the relaxation of the stress and
 * (2) the conversion from conformation tensor to stress tensor. One can choose from the pre-programmed models
 * Oldroyd-B, Giesekus, or Rolie-Poly. Optionally, you can also register a your own strategy. The fluid model is
 * specified through the database parameter "fluid_parameter". By specifying "USER_DEFINED", you can register your own
 * strategy operator. This class currently solves for the conformation tensor or the square root or logarithm of the
 * conformation tensor.
 *
 * <h2>Input Database Parameters</h2>
 * <ul>
 *     <li><code>InitialConditions</code>: Database used for specifying the initial conditions via a
 * muParserCartGridFunction.</li>
 *     <li><code>interp_type</code>: String to determine order used for filling ghost
 * cells. Possible values are:
 *     <ul>
 *       <li><code>"LINEAR"</code> (default)</li>
 *       <li><code>"QUADRATIC"</code></li>
 *     </ul></li>
 *     <li><code>evolution_type</code>: Parameter specifying the method of evolving the conformation tensor. Possible
 * values are:
 *     <ul>
 *         <li><code>"STANDARD"</code>: Evolve the conformation tensor itself (default).</li>
 *         <li><code>"SQUARE_ROOT"</code>: Evolve the square root of the conformation tensor.</li>
 *         <li><code>"LOGARITHM"</code>: Evolve the logarithm of the conformation tensor.</li>
 *     </ul></li>
 *     <li><code>log_determinant</code>: Boolean on writing the determinant of the conformation tensor (default
 * false).</li>
 *     <li><code>log_divergence</code>: Boolean on writing the divergence of the conformation tensor (default
 * false).</li> <li><code>fluid_model</code>: Specifies the fluid model. Possible values are: <ul>
 *         <li><code>"OLDROYDB"</code> (default) </li>
 *         <li><code>"GIESEKUS"</code></li>
 *         <li><code>"ROLIEPOLY"</code></li>
 *         <li><code>"USER_DEFINED"</code>: In this case, users must register a relaxation function with
 * <code>registerRelaxationOperator()</code>.</li>
 *     </ul></li>
 *     <li><code>error_on_spd</code>: Boolean that if true, kills the simulation if the conformation tensor is ever not
 * positive definite (default true).</li>
 *     <li><code>project_conformation_tensor</code>: Boolean that if true, projects the conformation tensor to nearest
 * non-negative definite tensor (default true).</li>
 *     <li><code>output_conformation_tensor</code>: Boolean that if true, writes the conformation tensor to
 * visualization files (default true).</li>
 *     <li><code>output_stress_tensor</code>: Boolean that if true, writes the resulting stress tensor to visualization
 * files (default true).</li>
 *     <li><code>output_divergence</code>: Boolean that if true, writes the divergence of the conformation tensor to
 * visualization files (default false).</li>
 *     <li><code>divergence_rel_tagging</code>: Boolean that if true, tags cells for refinement if the divergence of the
 * conformation tensor exceeds some relative threshold, specified by the array <code>divergence_rel_thresh</code>
 * (default false).</li>
 *     <li><code>divergence_abs_tagging</code>: Boolean that if true, tags cells for refinement if the divergence of the
 * conformation tensor exceeds some absolute threshold, specified by the array <code>divergence_abs_thresh</code>
 * (default false).</li>
 *     <li><code>ExtraStressBoundaryConditions_d</code>: Database used for boundary conditions for the evolved tensor.
 * Here, d is a placeholder for the component of the symmetric conformation tensor, specified using Voigt notation
 * (ignored for periodic boundary conditions).</li>
 *     <li><code>convective_operator_type</code>: String for the operator used for computing the advective term.</li>
 *     <li><code>difference_form</code>: String to determine in what form the convective operator is computed. Possible
 * values are: <ul> <li><code>"ADVECTIVE"</code> (default) </li> <li><code>"CONSERVATIVE"</code></li>
 *         <li><code>"SKEW_SYMMETRIC"</code></li>
 *     </ul></li>
 * </ul>
 *
 * There are also strategy specific items searched for in the database.
 */
class CFINSForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief This constructor creates Variable and VariableContext objects for storing the viscoleastic stresses at the
     * centers of the Cartesian grid. Sets up the advection diffusion solver to use the velocity function prescribed.
     * The initial and boundary conditions should be specified on the quantity being solved for (e.g. Conformation
     * tensor or square root or logarithm of the conformation tensor).
     */
    CFINSForcing(const std::string& object_name,
                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                 SAMRAI::tbox::Pointer<IBTK::CartGridFunction> u_fcn,
                 SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry,
                 SAMRAI::tbox::Pointer<IBAMR::AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                 SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_data_writer);

    /*!
     * \brief This constructor creates Variable and VariableContext objects for storing the viscoleastic stresses at the
     * centers of the Cartesian grid. Sets up the advection diffusion solver to use the fluid velocity from the fluid
     * solver. Note that this function must be registered with the fluid solver as a forcing function for the
     * viscoelastic stress to effect the fluid velocity. The initial and boundary conditions should be specified on the
     * quantity being solved for (e.g. Conformation tensor or square root or logarithm of the conformation tensor).
     */
    CFINSForcing(const std::string& object_name,
                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> app_initializer,
                 const SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> fluid_solver,
                 SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry,
                 SAMRAI::tbox::Pointer<IBAMR::AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator,
                 SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_data_writer);

    /*!
     * \brief Deallocates draw data and deletes boundary conditions.
     */
    ~CFINSForcing();

    /*!
     * \brief This function returns a pointer to the cell variable that stores the viscoelastic stress.
     */
    inline SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getVariable()
    {
        return d_C_cc_var;
    }

    /*!
     * \brief This function returns the patch data index used to store the viscoelastic stress.
     */
    inline int getVariableIdx()
    {
        auto var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        return var_db->mapVariableAndContextToIndex(d_C_cc_var, d_adv_diff_integrator->getCurrentContext());
    }

    /*!
     * \brief This function registers a strategy with the class.
     *
     * \note This function only should be called when <code>fluid_model</code> is set to <code>USER_DEFINED</code> in
     * the input file.
     */
    void registerCFStrategy(SAMRAI::tbox::Pointer<CFStrategy> rhs);

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete INSStaggeredStochasticForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Compute the divergence of the stress tensor. Also sets up requested visualizations.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = IBTK::invalid_level_number,
                                 const int finest_ln = IBTK::invalid_level_number) override;

    /*!
     * \brief Evaluate the divergence on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(nullptr)) override;

    /*!
     * \brief Evaluate the divergence on the specified patch level.
     */
    void setDataOnPatchLevel(const int data_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
                             const double data_time,
                             const bool initial_time) override;
    //\}

    /*!
     * \brief Check whether the provided patch index stores a positive definite tensor.
     */
    void checkPositiveDefinite(const int data_idx,
                               const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                               const double data_time,
                               const bool initial_time);

    /*!
     * \brief Tag cells based on the specifications provided in the input database.
     */
    void applyGradientDetector(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double error_data_time,
                               int tag_index,
                               bool initial_time,
                               bool /*richardson_extrapolation_too*/);

    /*!
     * \brief Projects the symmetric tensor stored in data_idx to the nearest non negative matrix in the L2 norm.
     */
    void projectTensor(const int data_idx,
                       const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                       const double data_time,
                       const bool initial_time,
                       const bool extended_box);

    static void
    apply_gradient_detector_callback(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int level_number,
                                     double error_data_time,
                                     int tag_index,
                                     bool initial_time,
                                     bool richardson_extrapolation_too,
                                     void* ctx);

    static void apply_project_tensor_callback(double current_time, double new_time, int cycle_num, void* ctx);

    /*!
     * \brief Return the advection diffusion integrator used to evolve the conformation tensor.
     */
    inline SAMRAI::tbox::Pointer<AdvDiffSemiImplicitHierarchyIntegrator> getAdvDiffHierarchyIntegrator()
    {
        return d_adv_diff_integrator;
    }

private:
    void commonConstructor(const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_data_writer,
                           SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry,
                           std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> vel_bcs);

    /*!
     * \brief Compute the determinant of the symmetric tensor stored in data_idx. Fills in d_max_det and d_min_det.
     */
    void findDeterminant(const int data_idx,
                         const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                         const double data_time,
                         const bool initial_time);

    /*!
     * \brief Square the symmetric tensor stored in data_idx in place.
     */
    void squareMatrix(const int data_idx,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                      const double data_time,
                      const bool initial_time,
                      const int coarsest_ln,
                      const int finest_ln,
                      const bool extended_box);

    /*!
     * \brief Exponentiate the symmetric tensor stored in data_idx in place.
     */
    void exponentiateMatrix(const int data_idx,
                            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                            const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                            const double data_time,
                            const bool initial_time,
                            const int coarsest_ln,
                            const int finest_ln,
                            const bool extended_box);

    /*!
     * \brief Given the conformation tensor, set up any requested drawing variables.
     */
    void setupPlotConformationTensor(int C_cc_idx);

    // Scratch variables
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_C_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    int d_C_scratch_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<IBTK::muParserCartGridFunction> d_init_conds;

    // Draw Variables
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_conform_var_draw, d_stress_var_draw,
        d_div_sig_var_draw;
    int d_conform_idx_draw = IBTK::invalid_index, d_stress_idx_draw = IBTK::invalid_index,
        d_div_sig_idx_draw = IBTK::invalid_index;
    bool d_conform_draw = true, d_stress_draw = true, d_div_sig_draw = false;

    // Extra parameters
    std::string d_fluid_model = "OLDROYDB", d_interp_type = "LINEAR";
    bool d_project_conform = true;
    TensorEvolutionType d_evolve_type = STANDARD;
    SAMRAI::tbox::Pointer<AdvDiffSemiImplicitHierarchyIntegrator> d_adv_diff_integrator;
    SAMRAI::tbox::Pointer<CFUpperConvectiveOperator> d_convec_oper;

    /**
     * Boundary conditions.
     */
    std::vector<std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM> > > d_conc_bc_coefs;

    /**
     * Pointers to the previous objects (required by some APIs).
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_conc_bc_coefs_ptrs;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;

    // Logging parameters
    double d_max_det = std::numeric_limits<double>::quiet_NaN(), d_min_det = std::numeric_limits<double>::quiet_NaN();
    bool d_log_det = false, d_log_div_sig = false;
    bool d_positive_def = true, d_error_on_spd = false;
    double d_min_norm = std::numeric_limits<double>::quiet_NaN(), d_max_norm = std::numeric_limits<double>::quiet_NaN();

    // AMR tagging
    SAMRAI::tbox::Array<double> d_div_sig_rel_thresh, d_div_sig_abs_thresh;
    bool d_div_sig_rel_tag = false, d_div_sig_abs_tag = false;

    // Velocity information
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_u_fcn;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_u_var;

    // Strategy for computing stress and relaxation.
    SAMRAI::tbox::Pointer<CFStrategy> d_cf_strategy;
};
} // Namespace IBAMR
#endif
