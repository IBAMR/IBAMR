// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
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

#ifndef included_IBAMR_TemperatureSemiImplicitHierarchyIntegrator
#define included_IBAMR_TemperatureSemiImplicitHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "HierarchyFaceDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <set>
#include <string>

namespace IBAMR
{
class ConvectiveOperator;
class INSVCStaggeredHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
namespace pdat
{
template <int DIM, class TYPE>
class CellVariable;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class TemperatureSemiImplicitHierarchyIntegrator manages the spatial
 * discretization and time integration of scalar- and vector-valued quantities
 * whose dynamics are governed by the temperature equation.
 *
 * Each quantity \f$ Q \f$ managed by the integrator may have a unique diffusion
 * coefficient \f$ \kappa \f$ and damping coefficient \f$ \lambda \f$, and may
 * optionally have a forcing term \f$ F \f$.  Additionally, a different
 * advection velocity may be used with each quantity registered with the
 * integrator.
 *
 * This hierarchy integrator advances all levels of the patch hierarchy
 * synchronously in time.  In particular, subcycling in time is \em not
 * performed.
 *
 * Various options are available for the spatial and temporal discretizations.
 *
 * \see HierarchyIntegrator
 * \see SAMRAI::mesh::StandardTagAndInitStrategy
 * \see SAMRAI::algs::TimeRefinementIntegrator
 * \see SAMRAI::algs::TimeRefinementLevelStrategy
 */
class TemperatureSemiImplicitHierarchyIntegrator : public AdvDiffSemiImplicitHierarchyIntegrator
{
public:
    /*!
     * The constructor for class AdvDiffSemiImplicitHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    TemperatureSemiImplicitHierarchyIntegrator(const std::string& object_name,
                                               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                               bool register_for_restart = true);

    /*!
     * The destructor for class AdvDiffSemiImplicitHierarchyIntegrator
     * unregisters the integrator object with the restart manager when the
     * object is so registered.
     */
    ~TemperatureSemiImplicitHierarchyIntegrator() = default;

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

    /*!
     * Register the specific heat variable.
     */
    void registerSpecificHeatVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Cp_var,
                                      const bool output_Cp = false);

    /*!
     * Register the density variable.
     */
    void registerDensityVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_var,
                                 const bool output_rho = false);

    /*!
     * \brief Function to reset fluid density or viscosity if they are
     * maintained by this integrator.
     */
    using ResetFluidPropertiesFcnPtr = void (*)(int property_idx,
                                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > property_var,
                                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                int cycle_num,
                                                double time,
                                                double current_time,
                                                double new_time,
                                                void* ctx);

    /*!
     * \brief Register function to reset fluid density.
     */
    void registerResetFluidDensityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset fluid viscosity.
     */
    void registerResetSpecificHeatFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register function to reset fluid viscosity.
     */
    void registerResetDiffusionCoefficientFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * Set the cell-centered density to be used with a particular
     * cell-centered quantity.
     *
     * \note The specified source term must have been already registered with
     * the hierarchy integrator.
     */
    void setDensityVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > rho_var);

    /*!
     * Set the cell-centered density to be used with a particular
     * cell-centered quantity.
     *
     * \note The specified source term must have been already registered with
     * the hierarchy integrator.
     */
    void setSpecificHeatVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Cp_var);

    /*!
     * Register INSVCStaggeredHierarchyIntegrator class which will be
     * used to get the variables maintained by this hierarchy integrator.
     */
    void registerINSVCStaggeredHierarchyIntegrator(
        SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> ins_hier_integrator);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    TemperatureSemiImplicitHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    TemperatureSemiImplicitHierarchyIntegrator(const AdvDiffSemiImplicitHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    TemperatureSemiImplicitHierarchyIntegrator& operator=(const AdvDiffSemiImplicitHierarchyIntegrator& that) = delete;

    /*!
     * Interpolating the side centered density to cell-centered density
     * using bilinear interpolation.
     */
    void interpolateSCMassDensityToCC(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx);

    /*!
     * Additional variables required.
     */
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > >
        d_Q_rho_map, d_Q_Cp_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, bool> d_rho_output;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, bool> d_Cp_output;

    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > > d_Cp_var, d_rho_var;

    std::vector<ResetFluidPropertiesFcnPtr> d_reset_rho_fcns, d_reset_Cp_fcns, d_reset_kappa_fcns;
    std::vector<void*> d_reset_rho_fcns_ctx, d_reset_Cp_fcns_ctx, d_reset_kappa_fcns_ctx;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_C_var, d_rho_vec_cc_var;

    int d_C_scratch_idx = IBTK::invalid_index, d_C_current_idx = IBTK::invalid_index, d_C_new_idx = IBTK::invalid_index,
        d_C_rhs_scratch_idx = IBTK::invalid_index, d_rho_vec_cc_current_idx = IBTK::invalid_index,
        d_rho_vec_cc_scratch_idx = IBTK::invalid_index, d_rho_vec_cc_new_idx = IBTK::invalid_index;

    SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> d_ins_hierarchy_integrator;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_TemperatureSemiImplicitHierarchyIntegrator
