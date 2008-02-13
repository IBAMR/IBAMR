// Filename: AdvectHypPatchOps.C
// Last modified: <12.Feb.2008 21:19:07 griffith@box221.cims.nyu.edu>
// Created on 12 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "AdvectHypPatchOps.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// STOOLS INCLUDES
#include <stools/ExtendedRobinBcCoefStrategy.h>
#include <stools/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <ArrayDataBasicOps.h>
#include <BoundaryBox.h>
#include <BoxArray.h>
#include <CartesianPatchGeometry.h>
#include <CellDataFactory.h>
#include <CellIndex.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <FaceIndex.h>
#include <FaceVariable.h>
#include <Index.h>
#include <LoadBalancer.h>
#include <PatchCellDataOpsReal.h>
#include <PatchFaceDataOpsReal.h>
#include <PatchData.h>
#include <VariableDatabase.h>
#include <tbox/Array.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 1)
#define ADVECT_CONSDIFF_F77 F77_FUNC_(advect_consdiff1d, ADVECT_CONSDIFF1D)
#define ADVECT_CONSDIFFWITHDIVSOURCE_F77 F77_FUNC_(advect_consdiffwithdivsource1d, ADVECT_CONSDIFFWITHDIVSOURCE1D)
#define ADVECT_DETECTGRAD_F77 F77_FUNC_(advect_detectgrad1d, ADVECT_DETECTGRAD1D)
#endif

#if (NDIM == 2)
#define ADVECT_CONSDIFF_F77 F77_FUNC_(advect_consdiff2d, ADVECT_CONSDIFF2D)
#define ADVECT_CONSDIFFWITHDIVSOURCE_F77 F77_FUNC_(advect_consdiffwithdivsource2d, ADVECT_CONSDIFFWITHDIVSOURCE2D)
#define ADVECT_DETECTGRAD_F77 F77_FUNC_(advect_detectgrad2d, ADVECT_DETECTGRAD2D)
#endif

#if (NDIM == 3)
#define ADVECT_CONSDIFF_F77 F77_FUNC_(advect_consdiff3d, ADVECT_CONSDIFF3D)
#define ADVECT_CONSDIFFWITHDIVSOURCE_F77 F77_FUNC_(advect_consdiffwithdivsource3d, ADVECT_CONSDIFFWITHDIVSOURCE3D)
#define ADVECT_DETECTGRAD_F77 F77_FUNC_(advect_detectgrad3d, ADVECT_DETECTGRAD3D)
#endif

extern "C"
{
    void
    ADVECT_CONSDIFF_F77(
        const double* ,
#if (NDIM == 1)
        const int& , const int& ,
        const int& ,
        const int& ,
        const double* ,
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
#endif
        double*
                             );

    void
    ADVECT_CONSDIFFWITHDIVSOURCE_F77(
        const double* , const double& ,
#if (NDIM == 1)
        const int& , const int& ,
        const int& ,
        const int& ,
        const int& ,
        const int& ,
        const double* ,
        const double* ,
        const double* ,
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
#endif
        double*);

    void
    ADVECT_DETECTGRAD_F77(
#if (NDIM == 1)
        const int& , const int& ,
        const int& , const int& , const int&,
#endif
#if (NDIM == 2)
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& , const int&,
        const int& , const int& , const int&,
#endif
#if (NDIM == 3)
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& , const int&,
        const int& , const int& , const int&,
        const int& , const int& , const int&,
#endif
        const double* ,
        const double& ,
        const int&,
        const double*,
        int* , int*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_register_model_vars;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_data_on_patch;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_stable_dt_on_patch;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_fluxes_on_patch;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_conservative_difference_on_patch;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_preprocess_advance_level_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_postprocess_advance_level_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_tag_richardson_extrapolation_cells;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_tag_gradient_detector_cells;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_set_physical_boundary_conditions;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 4;
static const int FLUXG = 1;

// Values for cell tagging routines.
static const int RICHARDSON_NEWLY_TAGGED   = -10;
static const int RICHARDSON_ALREADY_TAGGED = -11;
static const int TRUE_VAL  = 1;
static const int FALSE_VAL = 0;

// Version of AdvectHypPatchOps restart file data.
static const int ADVECT_HYP_PATCH_OPS_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvectHypPatchOps::AdvectHypPatchOps(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<GodunovAdvector> godunov_advector,
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom,
    bool register_for_restart)
    : d_integrator(NULL),
      d_godunov_advector(godunov_advector),
      d_Q_vars(),
      d_F_vars(),
      d_grad_vars(),
      d_Q_conservation_form(),
      d_flux_integral_vars(),
      d_q_integral_vars(),
      d_u_integral_var(),
      d_Q_inits(),
      d_Q_bc_coefs(),
      d_F_sets(),
      d_u_var(NULL),
      d_u_set(NULL),
      d_u_is_div_free(false),
      d_u_is_registered(false),
      d_compute_init_velocity(true),
      d_compute_half_velocity(true),
      d_compute_final_velocity(true),
      d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_grid_geometry(grid_geom),
      d_visit_writer(NULL),
      d_extrap_bc_helper(),
      d_coarse_fine_bdry_op(NULL),
      d_ghosts(CELLG),
      d_flux_ghosts(FLUXG),
      d_extrap_type("CONSTANT"),
      d_refinement_criteria(),
      d_dev_tol(),
      d_dev(),
      d_dev_time_max(),
      d_dev_time_min(),
      d_grad_tol(),
      d_grad_time_max(),
      d_grad_time_min(),
      d_rich_tol(),
      d_rich_time_max(),
      d_rich_time_min()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!godunov_advector.isNull());
    TBOX_ASSERT(!grid_geom.isNull());
#endif

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from given input/restart databases.
    bool is_from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (is_from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, is_from_restart);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_register_model_vars = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::registerModelVariables()");
        t_initialize_data_on_patch = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::initializeDataOnPatch()");
        t_compute_stable_dt_on_patch = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::computeStableDtOnPatch()");
        t_compute_fluxes_on_patch = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::computeFluxesOnPatch()");
        t_conservative_difference_on_patch = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::conservativeDifferenceOnPatch()");
        t_preprocess_advance_level_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::preprocessAdvanceLevelState()");
        t_postprocess_advance_level_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::postprocessAdvanceLevelState()");
        t_tag_richardson_extrapolation_cells = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::tagRichardsonExtrapolationCells()");
        t_tag_gradient_detector_cells = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::tagGradientDetectorCells()");
        t_set_physical_boundary_conditions = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::setPhysicalBoundaryConditions()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvectHypPatchOps::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// AdvectHypPatchOps

AdvectHypPatchOps::~AdvectHypPatchOps()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
    return;
}// ~AdvectHypPatchOps

const std::string&
AdvectHypPatchOps::getName() const
{
    return d_object_name;
}// getName

void
AdvectHypPatchOps::registerCoarseFineBoundaryRefinePatchStrategy(
    SAMRAI::tbox::Pointer<STOOLS::CoarseFineBoundaryRefinePatchStrategy> coarse_fine_bdry_op)
{
    d_coarse_fine_bdry_op = coarse_fine_bdry_op;
    return;
}// registerCoarseFineBoundaryRefinePatchStrategy

///
///  The following routines:
///
///      registerAdvectedQuantity(),
///      registerAdvectedQuantityWithSourceTerm(),
///      registerAdvectionVelocity(),
///      registerVisItDataWriter()
///
///  allow the AdvectHypPatchOps to be used as a generic advection scheme.
///

void
AdvectHypPatchOps::registerAdvectedQuantity(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
    const bool conservation_form,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> Q_init,
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const Q_bc_coef,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var)
{
    registerAdvectedQuantity(
        Q_var, conservation_form, Q_init,
        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(1,Q_bc_coef),
        grad_var);
    return;
}// registerAdvectedQuantity

void
AdvectHypPatchOps::registerAdvectedQuantity(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
    const bool conservation_form,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> Q_init,
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_u_is_registered);
    TBOX_ASSERT(!Q_var.isNull());
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > Q_factory =
        Q_var->getPatchDataFactory();
    const int Q_depth = Q_factory->getDefaultDepth();

    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> Q_bc_coefs_local = Q_bc_coefs;
    if (Q_bc_coefs_local.empty())
    {
        Q_bc_coefs_local = std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(
            Q_depth,static_cast<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NULL));
    }

    if (Q_depth != static_cast<int>(Q_bc_coefs_local.size()))
    {
        TBOX_ERROR(d_object_name << "::registerAdvectedQuantity():\n"
                   << "  data depth for variable " << Q_var->getName() << " is " << Q_depth << "\n"
                   << "  but " << Q_bc_coefs_local.size() << " boundary condition coefficient objects were provided to the class constructor." << std::endl);
    }

    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > flux_integral_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > q_integral_var;

    if (conservation_form)
    {
        flux_integral_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(
            d_object_name+"::"+Q_var->getName()+" advective flux time integral",
            Q_factory->getDefaultDepth());
    }

    if (!d_u_is_div_free || !conservation_form)
    {
        q_integral_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(
            d_object_name+"::"+Q_var->getName()+" time integral",
            Q_factory->getDefaultDepth());
    }

    d_Q_vars    .push_back(Q_var);
    d_Q_inits   .push_back(Q_init);
    d_Q_bc_coefs.push_back(Q_bc_coefs_local);

    d_Q_conservation_form.push_back(conservation_form);

    d_grad_vars.push_back(grad_var);

    d_F_vars.push_back(NULL);
    d_F_sets.push_back(NULL);

    d_flux_integral_vars.push_back(flux_integral_var);
    d_q_integral_vars.push_back(q_integral_var);

    if (!conservation_form && d_u_integral_var.isNull())
    {
        d_u_integral_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(
            d_object_name+"::"+d_u_var->getName()+" time integral");
    }

    return;
}// registerAdvectedQuantity

void
AdvectHypPatchOps::registerAdvectedQuantityWithSourceTerm(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var,
    const bool conservation_form,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> Q_init,
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const Q_bc_coef,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> F_set,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var)
{
    registerAdvectedQuantityWithSourceTerm(
        Q_var, F_var, conservation_form, Q_init,
        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(1,Q_bc_coef),
        F_set, grad_var);
    return;
}// registerAdvectedQuantityWithSourceTerm

void
AdvectHypPatchOps::registerAdvectedQuantityWithSourceTerm(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var,
    const bool conservation_form,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> Q_init,
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> F_set,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_u_is_registered);
    TBOX_ASSERT(!Q_var.isNull());
    TBOX_ASSERT(!F_var.isNull());
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > Q_factory =
        Q_var->getPatchDataFactory();
    const int Q_depth = Q_factory->getDefaultDepth();

    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> Q_bc_coefs_local = Q_bc_coefs;
    if (Q_bc_coefs_local.empty())
    {
        Q_bc_coefs_local = std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(
            Q_depth,static_cast<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NULL));
    }

    if (Q_depth != static_cast<int>(Q_bc_coefs_local.size()))
    {
        TBOX_ERROR(d_object_name << "::registerAdvectedQuantityWithSourceTerm():\n"
                   << "  data depth for variable " << Q_var->getName() << " is " << Q_depth << "\n"
                   << "  but " << Q_bc_coefs_local.size() << " boundary condition coefficient objects were provided to the class constructor." << std::endl);
    }

    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > flux_integral_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > q_integral_var;

    if (conservation_form)
    {
        flux_integral_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(
            d_object_name+"::"+Q_var->getName()+" advective flux time integral",
            Q_factory->getDefaultDepth());
    }

    if (!d_u_is_div_free || !conservation_form)
    {
        q_integral_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(
            d_object_name+"::"+Q_var->getName()+" time integral",
            Q_factory->getDefaultDepth());
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > F_factory =
        F_var->getPatchDataFactory();
    TBOX_ASSERT(Q_factory->getDefaultDepth() == F_factory->getDefaultDepth());
#endif
    d_Q_vars    .push_back(Q_var);
    d_Q_inits   .push_back(Q_init);
    d_Q_bc_coefs.push_back(Q_bc_coefs_local);

    d_Q_conservation_form.push_back(conservation_form);

    d_grad_vars.push_back(grad_var);

    d_F_vars.push_back(F_var);
    d_F_sets.push_back(F_set);

    d_flux_integral_vars.push_back(flux_integral_var);
    d_q_integral_vars.push_back(q_integral_var);

    if (!conservation_form && d_u_integral_var.isNull())
    {
        d_u_integral_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(
            d_object_name+"::"+d_u_var->getName()+" time integral");
    }

    return;
}// registerAdvectedQuantityWithSourceTerm

void
AdvectHypPatchOps::registerAdvectionVelocity(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > u_var,
    const bool u_is_div_free,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> u_set)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_u_is_registered);
    TBOX_ASSERT(!u_var.isNull());
#endif
    d_u_var = u_var;
    d_u_set = u_set;
    d_u_is_div_free = u_is_div_free;

    if (!d_u_is_div_free)
    {
        d_u_integral_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(
            d_object_name+"::"+u_var->getName()+" time integral");
    }

    d_u_is_registered = true;
    return;
}// registerAdvectionVelocity

#if (NDIM>1)
void
AdvectHypPatchOps::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_visit_writer = visit_writer;
    return;
}// registerVisItDataWriter
#endif

///
///  The following routines:
///
///      registerModelVariables(),
///      initializeDataOnPatch(),
///      computeStableDtOnPatch(),
///      computeFluxesOnPatch(),
///      conservativeDifferenceOnPatch(),
///      preprocessAdvanceLevelState(),
///      postprocessAdvanceLevelState(),
///      tagGradientDetectorCells(),
///      tagRichardsonExtrapolationCells()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::algs::HyperbolicPatchStrategy abstract base class.
///

void
AdvectHypPatchOps::registerModelVariables(
    SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>* integrator)
{
    t_register_model_vars->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(integrator != static_cast<SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>*>(NULL));
#endif
    d_integrator = integrator;

    typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;

    for (CellVariableVector::iterator it = d_Q_vars.begin(); it != d_Q_vars.end(); ++it)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = *it;

        d_integrator->registerVariable(
            Q_var, d_ghosts,
            SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
            d_grid_geometry,
            "CONSERVATIVE_COARSEN",
            "CONSERVATIVE_LINEAR_REFINE");

#if (NDIM > 1)
        if (!d_visit_writer.isNull())
        {
            const int Q_idx = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase()->
                mapVariableAndContextToIndex(
                    Q_var, d_integrator->getPlotContext());

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > Q_factory =
                Q_var->getPatchDataFactory();

            const int depth = Q_factory->getDefaultDepth();

            if (depth == 1)
            {
                d_visit_writer->registerPlotQuantity(
                    Q_var->getName(), "SCALAR", Q_idx);
            }
            else
            {
                if (depth == NDIM)
                {
                    d_visit_writer->registerPlotQuantity(
                        Q_var->getName(), "VECTOR", Q_idx);
                }

                for (int d = 0; d < depth; ++d)
                {
                    std::ostringstream stream;
                    stream << d;
                    d_visit_writer->registerPlotQuantity(
                        Q_var->getName()+"_"+stream.str(), "SCALAR", Q_idx, d);
                }
            }
        }
#endif
    }

    for (CellVariableVector::iterator it = d_F_vars.begin(); it != d_F_vars.end(); ++it)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var = *it;

        // Advected quantities do not necessarily have source terms!
        if (!F_var.isNull())
        {
            d_integrator->registerVariable(
                F_var, d_ghosts,
                SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
        }
    }

    typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > > FaceVariableVector;
    for (FaceVariableVector::iterator it = d_flux_integral_vars.begin();
         it != d_flux_integral_vars.end(); ++it)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > flux_integral_var = *it;

        // Not all quantities are conservatively differenced!
        if (!flux_integral_var.isNull())
        {
            d_integrator->registerVariable(
                flux_integral_var, d_flux_ghosts,
                SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>::FLUX,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "NO_REFINE");
        }
    }

    for (FaceVariableVector::iterator it = d_q_integral_vars.begin();
         it != d_q_integral_vars.end(); ++it)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > q_integral_var = *it;

        // Not all quantities are non-conservatively differenced!
        if (!q_integral_var.isNull())
        {
            d_integrator->registerVariable(
                q_integral_var, d_flux_ghosts,
                SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>::FLUX,
                d_grid_geometry,
                "CONSERVATIVE_COARSEN",
                "NO_REFINE");
        }
    }

    d_integrator->registerVariable(
        d_u_var, d_ghosts,
        SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
        d_grid_geometry,
        "CONSERVATIVE_COARSEN",
        "CONSERVATIVE_LINEAR_REFINE");

    // Note that this quantity is only required when non-conservative
    // differencing is employed, or when the velocity field is not discretely
    // divergence free.
    if (!d_u_integral_var.isNull())
    {
        d_integrator->registerVariable(
            d_u_integral_var, d_flux_ghosts,
            SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>::FLUX,
            d_grid_geometry,
            "CONSERVATIVE_COARSEN",
            "NO_REFINE");
    }

    t_register_model_vars->stop();
    return;
}// registerModelVariables

void
AdvectHypPatchOps::initializeDataOnPatch(
    SAMRAI::hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    t_initialize_data_on_patch->start();

    if (initial_time)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;

        // We try to use the STOOLS::SetDataStrategy associated with each
        // advected quantity.  If there is no strategy associated with a given
        // quantity, initialize its value to zero.
        for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
            const int Q_idx = var_db->mapVariableAndContextToIndex(
                Q_var, getDataContext());
            SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> Q_init = d_Q_inits[l];

            if (!Q_init.isNull())
            {
                Q_init->setDataOnPatch(
                    Q_idx, Q_var, patch, data_time, initial_time);
            }
            else
            {
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Q_data =
                    patch.getPatchData(Q_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(!Q_data.isNull());
#endif
                Q_data->fillAll(0.0);
            }
        }

        // We try to use the STOOLS::SetDataStrategy associated with the source
        // term.  If there is no strategy associated with the source term,
        // initialize its value to zero.
        for (CellVariableVector::size_type l = 0; l < d_F_vars.size(); ++l)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var = d_F_vars[l];

            if (!F_var.isNull())
            {
                const int F_idx = var_db->mapVariableAndContextToIndex(
                    F_var, getDataContext());
                SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy>  F_set = d_F_sets[l];

                if (!F_set.isNull())
                {
                    F_set->setDataOnPatch(
                        F_idx, F_var, patch, data_time, initial_time);
                }
                else
                {
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_data =
                        patch.getPatchData(F_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(!F_data.isNull());
#endif
                    F_data->fillAll(0.0);
                }
            }
        }

        // We try to use the STOOLS::SetDataStrategy associated with the
        // advection velocirt.  If there is no strategy associated with the
        // advection velocity, initialize its value to zero.
        const int u_idx = var_db->mapVariableAndContextToIndex(
            d_u_var, getDataContext());

        if (!d_u_set.isNull())
        {
            d_u_set->setDataOnPatch(
                u_idx, d_u_var, patch, data_time, initial_time);
        }
        else
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_data =
                patch.getPatchData(u_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!u_data.isNull());
#endif
            u_data->fillAll(0.0);
        }
    }

    t_initialize_data_on_patch->stop();
    return;
}// initializeDataOnPatch

double
AdvectHypPatchOps::computeStableDtOnPatch(
    SAMRAI::hier::Patch<NDIM>& patch,
    const bool initial_time,
    const double dt_time)
{
    t_compute_stable_dt_on_patch->start();

    (void) initial_time;
    (void) dt_time;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_data =
        patch.getPatchData(d_u_var, getDataContext());

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_data.isNull());
#endif

    const double stable_dt = d_godunov_advector->
        computeStableDtOnPatch(*u_data,patch);

    t_compute_stable_dt_on_patch->stop();
    return stable_dt;
}// computeStableDtOnPatch

void
AdvectHypPatchOps::computeFluxesOnPatch(
    SAMRAI::hier::Patch<NDIM>& patch,
    const double time,
    const double dt)
{
    t_compute_fluxes_on_patch->start();

    (void) time;

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_vars.size() == d_F_vars.size());
    TBOX_ASSERT(d_Q_vars.size() == d_flux_integral_vars.size());
    TBOX_ASSERT(d_Q_vars.size() == d_q_integral_vars.size());
#endif

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const SAMRAI::hier::Box<NDIM>& patch_box = patch.getBox();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_data =
        patch.getPatchData(d_u_var, getDataContext());
    SAMRAI::math::PatchFaceDataOpsReal<NDIM,double> patch_fc_data_ops;

    typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;

    // Predict time- and face-centered values.
    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Q_data =
            patch.getPatchData(d_Q_vars[l], getDataContext());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > q_integral_data =
            getQIntegralData(l, patch, getDataContext());
        if (!d_F_vars[l].isNull())
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_data =
                patch.getPatchData(d_F_vars[l], getDataContext());
            d_godunov_advector->predictValueWithSourceTerm(
                *q_integral_data, *u_data, *Q_data, *F_data, patch, dt);
        }
        else
        {
            d_godunov_advector->predictValue(
                *q_integral_data, *u_data, *Q_data, patch, dt);
        }
    }

    // For incompressible flow problems, we allow for the specification of an
    // auxiliary gradient that is used to enforce the incompressibility
    // constraint in an extremely approximate manner.
    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        if (!d_grad_vars[l].isNull())
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > q_integral_data =
                getQIntegralData(l, patch, getDataContext());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > grad_data =
                patch.getPatchData(d_grad_vars[l], getDataContext());
            d_godunov_advector->enforceIncompressibility(
                *q_integral_data, *u_data, *grad_data, patch);
        }
    }

    // Set physical boundary conditions for the face-centered predicted values
    // at inflow boundaries only.
    if (pgeom->getTouchesRegularBoundary())
    {
        setInflowBoundaryConditions(patch, time+0.5*dt);
    }

    // Update the advection velocity.
    if (!d_u_set.isNull() && d_u_set->isTimeDependent() && d_compute_half_velocity)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        const int u_idx = var_db->mapVariableAndContextToIndex(
            d_u_var, getDataContext());
        d_u_set->setDataOnPatch(u_idx, d_u_var, patch, time+0.5*dt);
    }

    // Compute fluxes for those quantities that are to be conservatively
    // differenced.
    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        if (d_Q_conservation_form[l])
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > flux_integral_data =
                getFluxIntegralData(l, patch, getDataContext());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > q_integral_data =
                getQIntegralData(l, patch, getDataContext());
            d_godunov_advector->computeFlux(
                *flux_integral_data, *u_data, *q_integral_data, patch, dt);
        }
    }

    // When u is not discretely divergence free or when Q is not conservatively
    // differenced, we maintain the time integral of the predicted face values.
    // These values are used in computing the proper source terms to maintain
    // consistency between the conservative and non-conservative forms of the
    // advection equation.
    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        if (!d_u_is_div_free || !d_Q_conservation_form[l])
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > q_integral_data =
                getQIntegralData(l, patch, getDataContext());
            patch_fc_data_ops.scale(q_integral_data, // dst
                                    dt,              // alpha
                                    q_integral_data, // src
                                    patch_box);
        }
    }

    // When u is not discretely divergence free or when any quantity is not
    // conservatively differenced, we maintain the time integral of the
    // advection velocity.  This value is used in computing the proper source
    // terms to maintain consistency between the conservative and
    // non-conservative forms of the advection equation.
    if (!d_u_integral_var.isNull())
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_integral_data =
            patch.getPatchData(d_u_integral_var, getDataContext());
        patch_fc_data_ops.scale(u_integral_data, // dst
                                dt,              // alpha
                                u_data,          // src
                                patch_box);
    }


    t_compute_fluxes_on_patch->stop();
    return;
}// computeFluxesOnPatch

void
AdvectHypPatchOps::conservativeDifferenceOnPatch(
    SAMRAI::hier::Patch<NDIM>& patch,
    const double time,
    const double dt,
    bool at_synchronization)
{
    t_conservative_difference_on_patch->start();

    (void) time;
    (void) at_synchronization;

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_vars.size() == d_flux_integral_vars.size());
    TBOX_ASSERT(d_Q_vars.size() == d_q_integral_vars.size());
#endif

    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Box<NDIM>& patch_box = patch.getBox();
    const SAMRAI::hier::Index<NDIM>& ilower = patch_box.lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch_box.upper();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_integral_data =
        (!d_u_integral_var.isNull()
         ? patch.getPatchData(d_u_integral_var, getDataContext())
         : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL));
    const SAMRAI::hier::IntVector<NDIM>& u_integral_data_ghost_cells =
        (!u_integral_data.isNull()
         ? u_integral_data->getGhostCellWidth()
         : 0);
    SAMRAI::math::PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
    typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;

    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Q_data =
            patch.getPatchData(d_Q_vars[l], getDataContext());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > flux_integral_data =
            getFluxIntegralData(l, patch, getDataContext());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > q_integral_data =
            getQIntegralData(l, patch, getDataContext());

        const SAMRAI::hier::IntVector<NDIM>& Q_data_ghost_cells = Q_data->getGhostCellWidth();
        const SAMRAI::hier::IntVector<NDIM>& flux_integral_data_ghost_cells =
            (!flux_integral_data.isNull()
             ? flux_integral_data->getGhostCellWidth()
             : 0);
        const SAMRAI::hier::IntVector<NDIM>& q_integral_data_ghost_cells =
            (!q_integral_data.isNull()
             ? q_integral_data->getGhostCellWidth()
             : 0);

        if (d_Q_conservation_form[l])
        {
            for (int depth = 0; depth < Q_data->getDepth(); ++depth)
            {
                if (d_u_is_div_free)
                {
#if (NDIM == 1)
                    ADVECT_CONSDIFF_F77(
                        dx,
                        ilower(0),iupper(0),
                        flux_integral_data_ghost_cells(0),
                        Q_data_ghost_cells(0),
                        flux_integral_data->getPointer(0,depth),
                        Q_data->getPointer(depth));
#endif
#if (NDIM == 2)
                    ADVECT_CONSDIFF_F77(
                        dx,
                        ilower(0),iupper(0),ilower(1),iupper(1),
                        flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),
                        Q_data_ghost_cells(0),Q_data_ghost_cells(1),
                        flux_integral_data->getPointer(0,depth),
                        flux_integral_data->getPointer(1,depth),
                        Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                    ADVECT_CONSDIFF_F77(
                        dx,
                        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                        flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),flux_integral_data_ghost_cells(2),
                        Q_data_ghost_cells(0),Q_data_ghost_cells(1),Q_data_ghost_cells(2),
                        flux_integral_data->getPointer(0,depth),
                        flux_integral_data->getPointer(1,depth),
                        flux_integral_data->getPointer(2,depth),
                        Q_data->getPointer(depth));
#endif
                }
                else
                {
#if (NDIM == 1)
                    ADVECT_CONSDIFFWITHDIVSOURCE_F77(
                        dx,dt,
                        ilower(0),iupper(0),
                        flux_integral_data_ghost_cells(0),
                        q_integral_data_ghost_cells(0),
                        u_integral_data_ghost_cells(0),
                        Q_data_ghost_cells(0),
                        flux_integral_data->getPointer(0,depth),
                        q_integral_data->getPointer(0,depth),
                        u_integral_data->getPointer(0),
                        Q_data->getPointer(depth));
#endif
#if (NDIM == 2)
                    ADVECT_CONSDIFFWITHDIVSOURCE_F77(
                        dx,dt,
                        ilower(0),iupper(0),ilower(1),iupper(1),
                        flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),
                        q_integral_data_ghost_cells(0),q_integral_data_ghost_cells(1),
                        u_integral_data_ghost_cells(0),u_integral_data_ghost_cells(1),
                        Q_data_ghost_cells(0),Q_data_ghost_cells(1),
                        flux_integral_data->getPointer(0,depth),
                        flux_integral_data->getPointer(1,depth),
                        q_integral_data->getPointer(0,depth),
                        q_integral_data->getPointer(1,depth),
                        u_integral_data->getPointer(0),
                        u_integral_data->getPointer(1),
                        Q_data->getPointer(depth));
#endif
#if (NDIM == 3)
                    ADVECT_CONSDIFFWITHDIVSOURCE_F77(
                        dx,dt,
                        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                        flux_integral_data_ghost_cells(0),flux_integral_data_ghost_cells(1),flux_integral_data_ghost_cells(2),
                        q_integral_data_ghost_cells(0),q_integral_data_ghost_cells(1),q_integral_data_ghost_cells(2),
                        u_integral_data_ghost_cells(0),u_integral_data_ghost_cells(1),u_integral_data_ghost_cells(2),
                        Q_data_ghost_cells(0),Q_data_ghost_cells(1),Q_data_ghost_cells(2),
                        flux_integral_data->getPointer(0,depth),
                        flux_integral_data->getPointer(1,depth),
                        flux_integral_data->getPointer(2,depth),
                        q_integral_data->getPointer(0,depth),
                        q_integral_data->getPointer(1,depth),
                        q_integral_data->getPointer(2,depth),
                        u_integral_data->getPointer(0),
                        u_integral_data->getPointer(1),
                        u_integral_data->getPointer(2),
                        Q_data->getPointer(depth));
#endif
                }
            }
        }
        else
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > N_data =
                new SAMRAI::pdat::CellData<NDIM,double>(patch_box,Q_data->getDepth(),0);
            d_godunov_advector->computeAdvectiveDerivative(
                *N_data, *u_integral_data, *q_integral_data, patch);
            patch_cc_data_ops.axpy(Q_data,  // dst
                                   -1.0/dt, // alpha
                                   N_data,  // src1
                                   Q_data,  // src2
                                   patch_box);
        }
    }

    t_conservative_difference_on_patch->stop();
    return;
}// conservativeDifferenceOnPatch

void
AdvectHypPatchOps::preprocessAdvanceLevelState(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
    double current_time,
    double dt,
    bool first_step,
    bool last_step,
    bool regrid_advance)
{
    t_preprocess_advance_level_state->start();

    (void) dt;
    (void) first_step;
    (void) last_step;
    (void) regrid_advance;

    // Update the advection velocity.
    if (!d_u_set.isNull() && d_u_set->isTimeDependent() && d_compute_init_velocity)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        const int u_idx = var_db->mapVariableAndContextToIndex(
            d_u_var, d_integrator->getScratchContext());
        d_u_set->setDataOnPatchLevel(u_idx, d_u_var, level, current_time);
    }

    if (!d_coarse_fine_bdry_op.isNull() && level->inHierarchy())
    {
        // Determine the patch data indices which require further modification
        // at coarse-fine interfaces.
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        SAMRAI::hier::ComponentSelector patch_data_indices;

        typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;
        for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
            const int Q_data_idx = var_db->mapVariableAndContextToIndex(
                Q_var, d_integrator->getScratchContext());
            patch_data_indices.setFlag(Q_data_idx);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var = d_F_vars[l];
            if (!F_var.isNull())
            {
                const int F_data_idx = var_db->mapVariableAndContextToIndex(
                    F_var, d_integrator->getScratchContext());
                patch_data_indices.setFlag(F_data_idx);
            }
        }

        d_coarse_fine_bdry_op->setPatchDataIndices(patch_data_indices);

        const SAMRAI::hier::IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill = CELLG;
            d_coarse_fine_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
        }
    }

    t_preprocess_advance_level_state->stop();
    return;
}// preprocessAdvanceLevelState

void
AdvectHypPatchOps::postprocessAdvanceLevelState(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
    double current_time,
    double dt,
    bool first_step,
    bool last_step,
    bool regrid_advance)
{
    t_postprocess_advance_level_state->start();

    (void) first_step;
    (void) last_step;
    (void) regrid_advance;

    SAMRAI::tbox::Pointer<SAMRAI::math::PatchCellDataOpsReal<NDIM,double> > patch_cc_data_ops =
        new SAMRAI::math::PatchCellDataOpsReal<NDIM,double>();
    typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;

    // Update the values of any time-dependent source terms and add the values
    // of all source terms to the advected quantities.
    for (CellVariableVector::size_type l = 0; l < d_F_vars.size(); ++l)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> new_context     = d_integrator->getNewContext();
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> scratch_context = d_integrator->getScratchContext();

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var   = d_Q_vars[l];
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var = d_F_vars[l];

        SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> F_set = d_F_sets[l];

        if (!F_var.isNull() && !F_set.isNull() && F_set->isTimeDependent())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Q_data =
                    patch->getPatchData(Q_var, new_context);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_data =
                    patch->getPatchData(F_var, scratch_context);

                patch_cc_data_ops->axpy(Q_data, // dst
                                        0.5*dt, // alpha
                                        F_data, // src1
                                        Q_data, // src2
                                        patch_box);
            }

            F_set->setDataOnPatchLevel(F_var, new_context,
                                       level, current_time+dt);

            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Q_data =
                    patch->getPatchData(Q_var, new_context);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_data =
                    patch->getPatchData(F_var, scratch_context);

                patch_cc_data_ops->axpy(Q_data, // dst
                                        0.5*dt, // alpha
                                        F_data, // src1
                                        Q_data, // src2
                                        patch_box);
            }
        }
        else if (!F_var.isNull())
        {
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Q_data =
                    patch->getPatchData(Q_var, new_context);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_data =
                    patch->getPatchData(F_var, scratch_context);

                patch_cc_data_ops->axpy(Q_data, // dst
                                        dt,     // alpha
                                        F_data, // src1
                                        Q_data, // src2
                                        patch_box);
            }
        }
    }

    // Update the advection velocity.
    if (!d_u_set.isNull() && d_u_set->isTimeDependent() && d_compute_final_velocity)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        const int u_idx = var_db->mapVariableAndContextToIndex(
            d_u_var, d_integrator->getNewContext());
        d_u_set->setDataOnPatchLevel(u_idx, d_u_var, level, current_time+dt);
    }

    t_postprocess_advance_level_state->stop();
    return;
}// postprocessAdvanceLevelState

void
AdvectHypPatchOps::tagRichardsonExtrapolationCells(
    SAMRAI::hier::Patch<NDIM>& patch,
    const int error_level_number,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> coarsened_fine,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> advanced_coarse,
    const double regrid_time,
    const double deltat,
    const int error_coarsen_ratio,
    const bool initial_error,
    const int tag_index,
    const bool uses_gradient_detector_too)
{
    t_tag_richardson_extrapolation_cells->start();

    (void) initial_error;

    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Box<NDIM>& patch_box = patch.getBox();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags = patch.getPatchData(tag_index);

    // Possible tagging criteria includes
    //    QVAL_RICHARDSON
    // The criteria is specified over a time interval.
    //
    // Loop over criteria provided and check to make sure we are in the
    // specified time interval.  If so, apply appropriate tagging for the level.
    for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ++ncrit)
    {
        std::string ref = d_refinement_criteria[ncrit];
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > coarsened_fine_var;
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > advanced_coarse_var;
        int size;
        double tol;
        bool time_allowed;

        if (ref == "QVAL_RICHARDSON")
        {
            size = d_rich_tol.getSize();
            tol = (error_level_number < size
                   ? d_rich_tol[error_level_number]
                   : d_rich_tol[size-1]);
            size = d_rich_time_min.getSize();
            double time_min = (error_level_number < size
                               ? d_rich_time_min[error_level_number]
                               : d_rich_time_min[size-1]);
            size = d_rich_time_max.getSize();
            double time_max = (error_level_number < size
                               ? d_rich_time_max[error_level_number]
                               : d_rich_time_max[size-1]);
            time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

            if (time_allowed)
            {
                // We tag wherever the global error > specified tolerance
                // (i.e. d_rich_tol).  The estimated global error is the local
                // truncation error times the approximate number of steps used
                // in the simulation.  Approximate the number of steps as:
                //
                //       steps = L / (s*deltat)
                //
                // where
                //
                //       L = length of problem domain
                //       s = wave speed
                //       delta t = timestep on current level
                const double* const xdomainlo = d_grid_geometry->getXLower();
                const double* const xdomainhi = d_grid_geometry->getXUpper();

                double max_dx = 0.0;
                double max_length = 0.0;

                for (int idir = 0; idir < NDIM; ++idir)
                {
                    max_dx = std::max(max_dx, dx[idir]);
                    double length = xdomainhi[idir] - xdomainlo[idir];
                    max_length = std::max(max_length, length);
                }

                double max_wave_speed = max_dx / deltat;
                double steps = max_length / (max_wave_speed * deltat);

                // Tag cells where |w_c - w_f| * (r^n -1) * steps
                //
                // where
                //       w_c = soln on coarse level
                //       w_f = soln on fine level
                //       r   = error coarsen ratio
                //       n   = spatial order of scheme
                int order = 2;
                double r = error_coarsen_ratio;
                double rnminus1 = pow(r,order) - 1;

                double diff = 0.0;
                double error = 0.0;

                typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;

                for (CellVariableVector::size_type l = 0;
                     l < d_Q_vars.size(); ++l)
                {
                    coarsened_fine_var =
                        patch.getPatchData(d_Q_vars[l], coarsened_fine);
                    advanced_coarse_var =
                        patch.getPatchData(d_Q_vars[l], advanced_coarse);

#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(!coarsened_fine_var.isNull());
                    TBOX_ASSERT(!advanced_coarse_var.isNull());
#endif

                    for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
                    {
                        // Compute error norm.
                        diff = 0.0;

                        for (int depth = 0;
                             depth < coarsened_fine_var->getDepth();
                             ++depth)
                        {
                            diff += pow((*advanced_coarse_var)(ic(),0) -
                                        (*coarsened_fine_var)(ic(),0),
                                        2.0);
                        }

                        diff = sqrt(diff);
                        error = std::abs(diff)*rnminus1*steps;

                        // Tag cell if error exceeds the prescribed threshold.
                        // Since we are operating on the actual tag values (not
                        // temporary ones) distinguish here tags that were
                        // previously set before coming into this routine and
                        // those that are set here.
                        //
                        //     RICHARDSON_ALREADY_TAGGED: tagged before coming
                        //     into this method
                        //
                        //     RICHARDSON_NEWLY_TAGGED: newly tagged in this
                        //     method
                        if (error > tol)
                        {
                            if ((*tags)(ic(),0))
                            {
                                (*tags)(ic(),0) = RICHARDSON_ALREADY_TAGGED;
                            }
                            else
                            {
                                (*tags)(ic(),0) = RICHARDSON_NEWLY_TAGGED;
                            }
                        }
                    }
                }
            } // time_allowed
        } // if QVAL_RICHARDSON
    } // loop over refinement criteria

    // If we are NOT performing gradient detector (i.e. only doing Richardson
    // extrapolation) set tags marked in this method to TRUE_VAL and all others
    // false.  Otherwise, leave tags set to the RICHARDSON_ALREADY_TAGGED and
    // RICHARDSON_NEWLY_TAGGED.  We use this information in the gradient
    // detector.
    if (!uses_gradient_detector_too)
    {
        for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            if ((*tags)(ic(),0) == RICHARDSON_ALREADY_TAGGED ||
                (*tags)(ic(),0) == RICHARDSON_NEWLY_TAGGED)
            {
                (*tags)(ic(),0) = TRUE_VAL;
            }
            else
            {
                (*tags)(ic(),0) = FALSE_VAL;
            }
        }
    }

    t_tag_richardson_extrapolation_cells->stop();
    return;
}// tagRichardsonExtrapolationCells

void
AdvectHypPatchOps::tagGradientDetectorCells(
    SAMRAI::hier::Patch<NDIM>& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_indx,
    const bool uses_richardson_extrapolation_too)
{
    t_tag_gradient_detector_cells->start();

    const int error_level_number = patch.getPatchLevelNumber();

    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Box<NDIM>& patch_box = patch.getBox();
    const SAMRAI::hier::Index<NDIM>& ilower = patch.getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch.getBox().upper();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags = patch.getPatchData(tag_indx);

    const int not_refine_tag_val = FALSE_VAL;
    const int refine_tag_val = TRUE_VAL;

    // Create a set of temporary tags and set to untagged value.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > temp_tags = new SAMRAI::pdat::CellData<NDIM,int>(
        patch_box, 1, d_ghosts);
    temp_tags->fillAll(not_refine_tag_val);

    // Possible tagging criteria includes
    //    QVAL_DEVIATION, QVAL_GRADIENT
    // The criteria is specified over a time interval.
    //
    // Loop over criteria provided and check to make sure we are in the
    // specified time interval.  If so, apply appropriate tagging for the level.
    for (int ncrit = 0; ncrit < d_refinement_criteria.getSize(); ++ncrit)
    {
        std::string ref = d_refinement_criteria[ncrit];
        const SAMRAI::hier::IntVector<NDIM>& tagghost = tags->getGhostCellWidth();

        int size = 0;
        double tol = 0.0;
        bool time_allowed = false;

        if (ref == "QVAL_DEVIATION")
        {
            size = d_dev_tol.getSize();
            tol = (error_level_number < size
                   ? d_dev_tol[error_level_number]
                   : d_dev_tol[size-1]);
            size = d_dev.getSize();
            double dev = (error_level_number < size
                          ? d_dev[error_level_number]
                          : d_dev[size-1]);
            size = d_dev_time_min.getSize();
            double time_min = (error_level_number < size
                               ? d_dev_time_min[error_level_number]
                               : d_dev_time_min[size-1]);
            size = d_dev_time_max.getSize();
            double time_max = (error_level_number < size
                               ? d_dev_time_max[error_level_number]
                               : d_dev_time_max[size-1]);
            time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

            if (time_allowed)
            {
                // Check for tags that have already been set in a previous step.
                // Do NOT consider values tagged with value
                // RICHARDSON_NEWLY_TAGGED since these were set most recently by
                // Richardson extrapolation.
                typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;
                for (CellVariableVector::size_type l = 0;
                     l < d_Q_vars.size(); ++l)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > var =
                        patch.getPatchData(d_Q_vars[l], getDataContext());
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(!var.isNull());
#endif
                    for (int depth = 0; depth < var->getDepth(); ++depth)
                    {
                        for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
                        {
                            double locden = tol;
                            int tag_val = (*tags)(ic(),0);
                            if (tag_val)
                            {
                                if (tag_val != RICHARDSON_NEWLY_TAGGED)
                                {
                                    locden *= 0.75;
                                }
                            }
                            if (std::abs((*var)(ic(),depth)-dev) > locden)
                            {
                                (*temp_tags)(ic(),0) = refine_tag_val;
                            }
                        }
                    }
                }
            }
        }

        if (ref == "QVAL_GRADIENT")
        {
            size = d_grad_tol.getSize();
            tol = (error_level_number < size
                   ? d_grad_tol[error_level_number]
                   : d_grad_tol[size-1]);
            size = d_grad_time_min.getSize();
            double time_min = (error_level_number < size
                               ? d_grad_time_min[error_level_number]
                               : d_grad_time_min[size-1]);
            size = d_grad_time_max.getSize();
            double time_max = (error_level_number < size
                               ? d_grad_time_max[error_level_number]
                               : d_grad_time_max[size-1]);
            time_allowed = (time_min <= regrid_time) && (time_max > regrid_time);

            if (time_allowed)
            {
                typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;
                for (CellVariableVector::size_type l = 0;
                     l < d_Q_vars.size(); ++l)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > var =
                        patch.getPatchData(d_Q_vars[l], getDataContext());
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(!var.isNull());
#endif
                    const SAMRAI::hier::IntVector<NDIM>& vghost = var->getGhostCellWidth();

                    for (int depth = 0; depth < var->getDepth(); ++depth)
                    {
                        ADVECT_DETECTGRAD_F77(
#if (NDIM == 1)
                            ilower(0),iupper(0),
                            vghost(0),tagghost(0),d_ghosts(0),
#endif
#if (NDIM == 2)
                            ilower(0),iupper(0), ilower(1),iupper(1),
                            vghost(0),tagghost(0),d_ghosts(0),
                            vghost(1),tagghost(1),d_ghosts(1),
#endif
#if (NDIM == 3)
                            ilower(0),iupper(0), ilower(1),iupper(1),ilower(2),iupper(2),
                            vghost(0),tagghost(0),d_ghosts(0),
                            vghost(1),tagghost(1),d_ghosts(1),
                            vghost(2),tagghost(2),d_ghosts(2),
#endif
                            dx,
                            tol,
                            refine_tag_val,
                            var->getPointer(depth),
                            tags->getPointer(),temp_tags->getPointer());
                    }
                }
            }
        }
    } // loop over criteria

    // Adjust temp_tags from those tags set in Richardson extrapolation.  Here,
    // we just reset any tags that were set in Richardson extrapolation to be
    // the designated "refine_tag_val".
    if (uses_richardson_extrapolation_too)
    {
        for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            if ((*tags)(ic(),0) == RICHARDSON_ALREADY_TAGGED ||
                (*tags)(ic(),0) == RICHARDSON_NEWLY_TAGGED)
            {
                (*temp_tags)(ic(),0) = refine_tag_val;
            }
        }
    }

    // Update tags.
    for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
    {
        (*tags)(ic(),0) = (*temp_tags)(ic(),0);
    }

    t_tag_gradient_detector_cells->stop();
    return;
}// tagGradientDetectorCells

///
///  The following routines:
///
///      setPhysicalBoundaryConditions(),
///      getRefineOpStencilWidth(),
///      preprocessRefine(),
///      postprocessRefine()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::xfer::RefinePatchStrategy abstract base class.
///

void
AdvectHypPatchOps::setPhysicalBoundaryConditions(
    SAMRAI::hier::Patch<NDIM>& patch,
    const double fill_time,
    const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill)
{
    t_set_physical_boundary_conditions->start();

    // Extrapolate the interior data to set the ghost cell values for the state
    // variables and for any forcing terms.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    const int u_data_idx = var_db->mapVariableAndContextToIndex(
        d_u_var, d_integrator->getScratchContext());
    d_extrap_bc_helper.setExtrapolationType(d_extrap_type);
    d_extrap_bc_helper.setPatchDataIndex(u_data_idx);
    d_extrap_bc_helper.setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);

    SAMRAI::hier::ComponentSelector patch_data_indices;
    typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;
    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
        const int Q_data_idx = var_db->mapVariableAndContextToIndex(
            Q_var, d_integrator->getScratchContext());
        patch_data_indices.setFlag(Q_data_idx);

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var = d_F_vars[l];
        if (!F_var.isNull())
        {
            const int F_data_idx = var_db->mapVariableAndContextToIndex(
                F_var, d_integrator->getScratchContext());
            patch_data_indices.setFlag(F_data_idx);
        }
    }

    d_extrap_bc_helper.setExtrapolationType(d_extrap_type);
    d_extrap_bc_helper.setPatchDataIndices(patch_data_indices);
    d_extrap_bc_helper.setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);

    if (!d_coarse_fine_bdry_op.isNull())
    {
        d_coarse_fine_bdry_op->setPatchDataIndices(patch_data_indices);
        d_coarse_fine_bdry_op->setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    }

    t_set_physical_boundary_conditions->stop();
    return;
}// setPhysicalBoundaryConditions

SAMRAI::hier::IntVector<NDIM>
AdvectHypPatchOps::getRefineOpStencilWidth() const
{
    if (!d_coarse_fine_bdry_op.isNull())
    {
        return d_coarse_fine_bdry_op->getRefineOpStencilWidth();
    }
    else
    {
        return 0;
    }
}// getRefineOpStencilWidth

void
AdvectHypPatchOps::preprocessRefine(
    SAMRAI::hier::Patch<NDIM>& fine,
    const SAMRAI::hier::Patch<NDIM>& coarse,
    const SAMRAI::hier::Box<NDIM>& fine_box,
    const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    if (!d_coarse_fine_bdry_op.isNull())
    {
        // Determine the patch data indices which require further modification
        // at coarse-fine interfaces.
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        SAMRAI::hier::ComponentSelector patch_data_indices;

        typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;
        for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
            const int Q_data_idx = var_db->mapVariableAndContextToIndex(
                Q_var, getDataContext());
            patch_data_indices.setFlag(Q_data_idx);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var = d_F_vars[l];
            if (!F_var.isNull())
            {
                const int F_data_idx = var_db->mapVariableAndContextToIndex(
                    F_var, getDataContext());
                patch_data_indices.setFlag(F_data_idx);
            }
        }

        d_coarse_fine_bdry_op->setPatchDataIndices(patch_data_indices);
        d_coarse_fine_bdry_op->preprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
}// preprocessRefine

void
AdvectHypPatchOps::postprocessRefine(
    SAMRAI::hier::Patch<NDIM>& fine,
    const SAMRAI::hier::Patch<NDIM>& coarse,
    const SAMRAI::hier::Box<NDIM>& fine_box,
    const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    if (!d_coarse_fine_bdry_op.isNull())
    {
        // Determine the patch data indices which require further modification
        // at coarse-fine interfaces.
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        SAMRAI::hier::ComponentSelector patch_data_indices;

        typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;
        for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
            const int Q_data_idx = var_db->mapVariableAndContextToIndex(
                Q_var, getDataContext());
            patch_data_indices.setFlag(Q_data_idx);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var = d_F_vars[l];
            if (!F_var.isNull())
            {
                const int F_data_idx = var_db->mapVariableAndContextToIndex(
                    F_var, getDataContext());
                patch_data_indices.setFlag(F_data_idx);
            }
        }

        d_coarse_fine_bdry_op->setPatchDataIndices(patch_data_indices);
        d_coarse_fine_bdry_op->postprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
}// postprocessRefine

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::tbox::Serializable abstract base class.
///

void
AdvectHypPatchOps::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif

    db->putInteger("ADVECT_HYP_PATCH_OPS_VERSION", ADVECT_HYP_PATCH_OPS_VERSION);

    db->putIntegerArray("d_ghosts", static_cast<int*>(d_ghosts), NDIM);
    db->putIntegerArray("d_flux_ghosts", static_cast<int*>(d_flux_ghosts), NDIM);

    if (d_refinement_criteria.getSize() > 0)
    {
        db->putStringArray("d_refinement_criteria", d_refinement_criteria);
    }
    for (int i = 0; i < d_refinement_criteria.getSize(); ++i)
    {
        if (d_refinement_criteria[i] == "QVAL_DEVIATION")
        {
            db->putDoubleArray("d_dev_tol", d_dev_tol);
            db->putDoubleArray("d_dev", d_dev);
            db->putDoubleArray("d_dev_time_max", d_dev_time_max);
            db->putDoubleArray("d_dev_time_min", d_dev_time_min);
        }
        else if (d_refinement_criteria[i] == "QVAL_GRADIENT")
        {
            db->putDoubleArray("d_grad_tol", d_grad_tol);
            db->putDoubleArray("d_grad_time_max", d_grad_time_max);
            db->putDoubleArray("d_grad_time_min", d_grad_time_min);
        }
        else if (d_refinement_criteria[i] == "QVAL_RICHARDSON")
        {
            db->putDoubleArray("d_rich_tol", d_rich_tol);
            db->putDoubleArray("d_rich_time_max", d_rich_time_max);
            db->putDoubleArray("d_rich_time_min", d_rich_time_min);
        }
    }

    t_put_to_database->stop();
    return;
}// putToDatabase

///
///  The following routines:
///
///      printClassData()
///
///  are provided for your viewing pleasure.
///

void
AdvectHypPatchOps::printClassData(
    std::ostream &os) const
{
    os << "\nAdvectHypPatchOps::printClassData..." << std::endl;
    os << "this = " << const_cast<AdvectHypPatchOps*>(this) << std::endl;
    os << "d_object_name = " << d_object_name << "\n"
       << "d_registered_for_restart = " << d_registered_for_restart << std::endl;
    os << "d_grid_geometry = " << d_grid_geometry.getPointer() << std::endl;
#if (NDIM>1)
    os << "d_visit_writer = " << d_visit_writer.getPointer() << std::endl;
#endif
    os << "d_ghosts = " << d_ghosts << "\n"
       << "d_flux_ghosts = " << d_flux_ghosts << "\n"
       << "d_extrap_type = " << d_extrap_type << std::endl;
    for (int j = 0; j < d_refinement_criteria.getSize(); ++j)
    {
        os << "d_refinement_criteria[" << j << "] = "
           << d_refinement_criteria[j] << std::endl;
    }
    os << std::endl;
    for (int j = 0; j < d_dev_tol.getSize(); ++j)
    {
        os << "d_dev_tol[" << j << "] = "
           << d_dev_tol[j] << std::endl;
    }
    for (int j = 0; j < d_dev.getSize(); ++j)
    {
        os << "d_dev[" << j << "] = "
           << d_dev[j] << std::endl;
    }
    if (d_dev.getSize()) os << std::endl;
    for (int j = 0; j < d_dev_time_max.getSize(); ++j)
    {
        os << "d_dev_time_max[" << j << "] = "
           << d_dev_time_max[j] << std::endl;
    }
    if (d_dev_time_max.getSize()) os << std::endl;
    for (int j = 0; j < d_dev_time_min.getSize(); ++j)
    {
        os << "d_dev_time_min[" << j << "] = "
           << d_dev_time_min[j] << std::endl;
    }
    if (d_dev_time_min.getSize()) os << std::endl;
    for (int j = 0; j < d_grad_tol.getSize(); ++j)
    {
        os << "d_grad_tol[" << j << "] = "
           << d_grad_tol[j] << std::endl;
    }
    if (d_grad_tol.getSize()) os << std::endl;
    for (int j = 0; j < d_grad_time_max.getSize(); ++j)
    {
        os << "d_grad_time_max[" << j << "] = "
           << d_grad_time_max[j] << std::endl;
    }
    if (d_grad_time_max.getSize()) os << std::endl;
    for (int j = 0; j < d_grad_time_min.getSize(); ++j)
    {
        os << "d_grad_time_min[" << j << "] = "
           << d_grad_time_min[j] << std::endl;
    }
    if (d_grad_time_min.getSize()) os << std::endl;
    for (int j = 0; j < d_rich_tol.getSize(); ++j)
    {
        os << "d_rich_tol[" << j << "] = "
           << d_rich_tol[j] << std::endl;
    }
    if (d_rich_tol.getSize()) os << std::endl;
    for (int j = 0; j < d_rich_time_max.getSize(); ++j)
    {
        os << "d_rich_time_max[" << j << "] = "
           << d_rich_time_max[j] << std::endl;
    }
    if (d_rich_time_max.getSize()) os << std::endl;
    for (int j = 0; j < d_rich_time_min.getSize(); ++j)
    {
        os << "d_rich_time_min[" << j << "] = "
           << d_rich_time_min[j] << std::endl;
    }
    if (d_rich_time_min.getSize()) os << std::endl;

    return;
}// printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> >
AdvectHypPatchOps::getFluxIntegralData(
    const size_t l,
    SAMRAI::hier::Patch<NDIM>& patch,
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> context)
{
    if (d_Q_conservation_form[l])
    {
        return patch.getPatchData(d_flux_integral_vars[l], context);
    }
    else
    {
        return SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> >(NULL);
    }
}// getFluxIntegralData

SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> >
AdvectHypPatchOps::getQIntegralData(
    const size_t l,
    SAMRAI::hier::Patch<NDIM>& patch,
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> context)
{
    if (!d_u_is_div_free || !d_Q_conservation_form[l])
    {
        return patch.getPatchData(d_q_integral_vars[l], context);
    }
    else
    {
        return getFluxIntegralData(l, patch, context);
    }
}// getQIntegralData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
AdvectHypPatchOps::setInflowBoundaryConditions(
    SAMRAI::hier::Patch<NDIM>& patch,
    const double fill_time)
{
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

    // There is nothing to do if the patch does not touch a regular (physical)
    // boundary.
    if (!pgeom->getTouchesRegularBoundary()) return;

    // Compute the codimension one boundary boxes.
    const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> > physical_codim1_boxes =
        STOOLS::PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);

    // There is nothing to do if the patch does not have any codimension one
    // boundary boxes.
    if (physical_codim1_boxes.size() == 0) return;

    // Loop over the boundary fill boxes and set boundary conditions at inflow
    // boundaries only.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_data =
        patch.getPatchData(d_u_var, getDataContext());
    const SAMRAI::hier::Box<NDIM>& patch_box = patch.getBox();
    const double* const dx = pgeom->getDx();
    typedef std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > CellVariableVector;
    for (CellVariableVector::size_type l = 0; l < d_Q_vars.size(); ++l)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
        const int Q_data_idx = var_db->mapVariableAndContextToIndex(
            Q_var, d_integrator->getScratchContext());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > q_integral_data =
            getQIntegralData(l, patch, getDataContext());

        // Setup any extended Robin BC coef objects.
        for (int depth = 0; depth < q_integral_data->getDepth(); ++depth)
        {
            STOOLS::ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<STOOLS::ExtendedRobinBcCoefStrategy*>(d_Q_bc_coefs[l][depth]);
            if (extended_bc_coef != NULL)
            {
                extended_bc_coef->setHomogeneousBc(false);
                extended_bc_coef->setTargetPatchDataIndex(Q_data_idx);
            }
        }

        // Set the boundary conditions.
        SAMRAI::math::ArrayDataBasicOps<NDIM,double> array_ops;
        (void) array_ops;
        for (int n = 0; n < physical_codim1_boxes.size(); ++n)
        {
            const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const int location_index   = bdry_box.getLocationIndex();
            const int bdry_normal_axis = location_index/2;
            const bool is_lower        = location_index%2 == 0;

            static const SAMRAI::hier::IntVector<NDIM> gcw_to_fill = 1;
            const SAMRAI::hier::Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
            const SAMRAI::hier::BoundaryBox<NDIM> trimmed_bdry_box(
                bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
            const SAMRAI::hier::Box<NDIM> bc_coef_box = STOOLS::PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            // Loop over the boundary box indices and compute the nearest
            // interior index.
            for (int depth = 0; depth < q_integral_data->getDepth(); ++depth)
            {
                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > acoef_data =
                    new SAMRAI::pdat::ArrayData<NDIM,double>(bc_coef_box, 1);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > bcoef_data =
                    new SAMRAI::pdat::ArrayData<NDIM,double>(bc_coef_box, 1);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > gcoef_data =
                    new SAMRAI::pdat::ArrayData<NDIM,double>(bc_coef_box, 1);
#if USING_OLD_ROBIN_BC_INTERFACE
                // In the old interface, beta = (1-alpha).
                d_Q_bc_coefs[l][depth]->setBcCoefs(
                    acoef_data, gcoef_data, Q_var,
                    patch, trimmed_bdry_box, fill_time);
                array_ops.scale(*bcoef_data, -1.0, *acoef_data, bc_coef_box);
                array_ops.addScalar(*bcoef_data, *bcoef_data, 1.0, bc_coef_box);
#else
                d_Q_bc_coefs[l][depth]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, Q_var,
                    patch, trimmed_bdry_box, fill_time);
#endif

                for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = b();
                    const SAMRAI::pdat::FaceIndex<NDIM> i_f(i, bdry_normal_axis, SAMRAI::pdat::FaceIndex<NDIM>::Lower);

                    bool is_inflow_bdry = (is_lower && (*u_data)(i_f) > 0.0) || (!is_lower && (*u_data)(i_f) < 0.0);
                    if (is_inflow_bdry)
                    {
                        const double& a = (*acoef_data)(i,0);
                        const double& b = (*bcoef_data)(i,0);
                        const double& g = (*gcoef_data)(i,0);
                        const double& h = dx[bdry_normal_axis];

                        SAMRAI::hier::Index<NDIM> i_intr(i);
                        if (is_lower)
                        {
                            // intentionally blank
                        }
                        else
                        {
                            i_intr(bdry_normal_axis) -= 1;
                        }

                        const SAMRAI::pdat::FaceIndex<NDIM> i_f_bdry(i_intr, bdry_normal_axis, (is_lower ? SAMRAI::pdat::FaceIndex<NDIM>::Lower : SAMRAI::pdat::FaceIndex<NDIM>::Upper));
                        const SAMRAI::pdat::FaceIndex<NDIM> i_f_intr(i_intr, bdry_normal_axis, (is_lower ? SAMRAI::pdat::FaceIndex<NDIM>::Upper : SAMRAI::pdat::FaceIndex<NDIM>::Lower));
                        const double& q_i = (*q_integral_data)(i_f_intr,depth);
                        const double q_b = (b*q_i + g*h)/(a*h + b);
                        (*q_integral_data)(i_f_bdry,depth) = q_b;
                    }
                }
            }
        }
    }
    return;
}// setInflowBoundaryConditions

void
AdvectHypPatchOps::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    (void) is_from_restart;

    d_compute_init_velocity  = db->getBoolWithDefault("compute_init_velocity" , d_compute_init_velocity);
    d_compute_half_velocity  = db->getBoolWithDefault("compute_half_velocity" , d_compute_half_velocity);
    d_compute_final_velocity = db->getBoolWithDefault("compute_final_velocity", d_compute_final_velocity);

    d_extrap_type = db->getStringWithDefault("extrap_type", d_extrap_type);
    if (!(d_extrap_type == "CONSTANT" || d_extrap_type == "LINEAR" || d_extrap_type == "QUADRATIC"))
    {
        TBOX_ERROR("AdvectHypPatchOps::getFromInput():\n"
                   << "  unknown extrapolation type: " << d_extrap_type << "\n"
                   << "  valid selections are: CONSTANT, LINEAR, or QUADRATIC" << std::endl);
    }

    if (db->keyExists("Refinement_data"))
    {
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> refine_db = db->getDatabase("Refinement_data");
        SAMRAI::tbox::Array<std::string> refinement_keys = refine_db->getAllKeys();
        int num_keys = refinement_keys.getSize();

        if (refine_db->keyExists("refine_criteria"))
        {
            d_refinement_criteria = refine_db->getStringArray("refine_criteria");
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n"
                         << "  No key `refine_criteria' found in data for"
                         << " RefinementData. No refinement will occur.\n");
        }

        SAMRAI::tbox::Array<std::string> ref_keys_defined(num_keys);
        int def_key_cnt = 0;
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> error_db;

        for (int i = 0; i < refinement_keys.getSize(); ++i)
        {
            std::string error_key = refinement_keys[i];
            error_db.setNull();

            if (error_key != "refine_criteria")
            {
                if (!(error_key == "QVAL_DEVIATION" ||
                      error_key == "QVAL_GRADIENT"  ||
                      error_key == "QVAL_RICHARDSON"))
                {
                    TBOX_ERROR(d_object_name << ":\n"
                               << "  Unknown refinement criteria: " << error_key
                               << "\nin input.\n");
                }
                else
                {
                    error_db = refine_db->getDatabase(error_key);
                    ref_keys_defined[def_key_cnt] = error_key;
                    def_key_cnt++;
                }

                if (!error_db.isNull() && error_key == "QVAL_DEVIATION")
                {
                    if (error_db->keyExists("dev_tol"))
                    {
                        d_dev_tol = error_db->getDoubleArray("dev_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name << ":\n"
                                   << "  No key `dev_tol' found in data for "
                                   << error_key << "\n");
                    }

                    if (error_db->keyExists("qval_dev"))
                    {
                        d_dev = error_db->getDoubleArray("qval_dev");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name << ":\n"
                                   << "  No key `qval_dev' found in data for "
                                   << error_key << "\n");
                    }

                    if (error_db->keyExists("time_max"))
                    {
                        d_dev_time_max = error_db->getDoubleArray("time_max");
                    }
                    else
                    {
                        d_dev_time_max.resizeArray(1);
                        d_dev_time_max[0] = std::numeric_limits<double>::max();
                    }

                    if (error_db->keyExists("time_min"))
                    {
                        d_dev_time_min = error_db->getDoubleArray("time_min");
                    }
                    else
                    {
                        d_dev_time_min.resizeArray(1);
                        d_dev_time_min[0] = 0.0;
                    }
                }

                if (!error_db.isNull() && error_key == "QVAL_GRADIENT")
                {
                    if (error_db->keyExists("grad_tol"))
                    {
                        d_grad_tol = error_db->getDoubleArray("grad_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name << ":\n"
                                   << "  No key `grad_tol' found in data for "
                                   << error_key << "\n");
                    }

                    if (error_db->keyExists("time_max"))
                    {
                        d_grad_time_max = error_db->getDoubleArray("time_max");
                    }
                    else
                    {
                        d_grad_time_max.resizeArray(1);
                        d_grad_time_max[0] = std::numeric_limits<double>::max();
                    }

                    if (error_db->keyExists("time_min"))
                    {
                        d_grad_time_min = error_db->getDoubleArray("time_min");
                    }
                    else
                    {
                        d_grad_time_min.resizeArray(1);
                        d_grad_time_min[0] = 0.0;
                    }
                }

                if (!error_db.isNull() && error_key == "QVAL_RICHARDSON")
                {
                    if (error_db->keyExists("rich_tol"))
                    {
                        d_rich_tol = error_db->getDoubleArray("rich_tol");
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name << ":\n"
                                   << "  No key `rich_tol' found in data for "
                                   << error_key << "\n");
                    }

                    if (error_db->keyExists("time_max"))
                    {
                        d_rich_time_max = error_db->getDoubleArray("time_max");
                    }
                    else
                    {
                        d_rich_time_max.resizeArray(1);
                        d_rich_time_max[0] = std::numeric_limits<double>::max();
                    }

                    if (error_db->keyExists("time_min"))
                    {
                        d_rich_time_min = error_db->getDoubleArray("time_min");
                    }
                    else
                    {
                        d_rich_time_min.resizeArray(1);
                        d_rich_time_min[0] = 0.0;
                    }
                }
            }
        } // loop over refine criteria

        // Check that input is found for each string identifier in key list.
        for (int k0 = 0; k0 < d_refinement_criteria.getSize(); ++k0)
        {
            std::string use_key = d_refinement_criteria[k0];
            bool key_found = false;
            for (int k1 = 0; k1 < def_key_cnt; ++k1)
            {
                std::string def_key = ref_keys_defined[k1];
                if (def_key == use_key)
                {
                    key_found = true;
                }
            }

            if (!key_found)
            {
                TBOX_ERROR(d_object_name << ":\n"
                           << "  No input found for specified refine criteria: "
                           << d_refinement_criteria[k0] << "\n");
            }
        }
    } // refine db entry exists

    return;
}// getFromInput

void
AdvectHypPatchOps::getFromRestart()
{
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> root_db =
        SAMRAI::tbox::RestartManager::getManager()->getRootDatabase();

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db;

    if (root_db->isDatabase(d_object_name))
    {
        db = root_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("ADVECT_HYP_PATCH_OPS_VERSION");
    if (ver != ADVECT_HYP_PATCH_OPS_VERSION)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  Restart file version different than class version.");
    }

    db->getIntegerArray("d_ghosts", static_cast<int*>(d_ghosts), NDIM);
    if (d_ghosts != SAMRAI::hier::IntVector<NDIM>(CELLG))
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  Key data `d_ghosts' in restart file != CELLG.\n");
    }

    db->getIntegerArray("d_flux_ghosts", static_cast<int*>(d_flux_ghosts), NDIM);
    if (d_flux_ghosts != SAMRAI::hier::IntVector<NDIM>(FLUXG))
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  Key data `d_flux_ghosts' in restart file != FLUXG.\n");
    }

    if (db->keyExists("d_refinement_criteria"))
    {
        d_refinement_criteria = db->getStringArray("d_refinement_criteria");
    }
    for (int i = 0; i < d_refinement_criteria.getSize(); ++i)
    {
        if (d_refinement_criteria[i] == "QVAL_DEVIATION")
        {
            d_dev_tol = db->getDoubleArray("d_dev_tol");
            d_dev = db->getDoubleArray("d_dev");
            d_dev_time_max = db->getDoubleArray("d_dev_time_max");
            d_dev_time_min = db->getDoubleArray("d_dev_time_min");
        }
        else if (d_refinement_criteria[i] == "QVAL_GRADIENT")
        {
            d_grad_tol = db->getDoubleArray("d_grad_tol");
            d_grad_time_max = db->getDoubleArray("d_grad_time_max");
            d_grad_time_min = db->getDoubleArray("d_grad_time_min");
        }
        else if (d_refinement_criteria[i] == "QVAL_RICHARDSON")
        {
            d_rich_tol = db->getDoubleArray("d_rich_tol");
            d_rich_time_max = db->getDoubleArray("d_rich_time_max");
            d_rich_time_min = db->getDoubleArray("d_rich_time_min");
        }
    }

    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::AdvectHypPatchOps>;

//////////////////////////////////////////////////////////////////////////////
