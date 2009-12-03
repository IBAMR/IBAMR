// Filename: INSStaggeredCenteredConvectiveOperator.C
// Last modified: <14.Aug.2009 18:11:35 griffith@boyce-griffiths-mac-pro.local>
// Created on 30 Oct 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredCenteredConvectiveOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/CartExtrapPhysBdryOp.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <FaceData.h>
#include <FaceGeometry.h>
#include <SideData.h>
#include <SideGeometry.h>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_STAGGERED_DIVERGENCE_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_divergence_derivative2d,NAVIER_STOKES_STAGGERED_DIVERGENCE_DERIVATIVE2D)
#define NAVIER_STOKES_STAGGERED_ADVECTION_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_advection_derivative2d,NAVIER_STOKES_STAGGERED_ADVECTION_DERIVATIVE2D)
#define NAVIER_STOKES_STAGGERED_SKEW_SYMMETRIC_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_skew_symmetric_derivative2d,NAVIER_STOKES_STAGGERED_SKEW_SYMMETRIC_DERIVATIVE2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_STAGGERED_DIVERGENCE_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_divergence_derivative3d,NAVIER_STOKES_STAGGERED_DIVERGENCE_DERIVATIVE3D)
#define NAVIER_STOKES_STAGGERED_ADVECTION_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_advection_derivative3d,NAVIER_STOKES_STAGGERED_ADVECTION_DERIVATIVE3D)
#define NAVIER_STOKES_STAGGERED_SKEW_SYMMETRIC_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_skew_symmetric_derivative3d,NAVIER_STOKES_STAGGERED_SKEW_SYMMETRIC_DERIVATIVE3D)
#endif

extern "C"
{
    void
    NAVIER_STOKES_STAGGERED_DIVERGENCE_DERIVATIVE_FC(
        const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double* , double* , double*
#endif
                                                          );

    void
    NAVIER_STOKES_STAGGERED_ADVECTION_DERIVATIVE_FC(
        const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double* , double* , double*
#endif
                                                          );

    void
    NAVIER_STOKES_STAGGERED_SKEW_SYMMETRIC_DERIVATIVE_FC(
        const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double* , double* , double*
#endif
                                                          );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int GADVECTG = 1;

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredCenteredConvectiveOperator::INSStaggeredCenteredConvectiveOperator(
    const INSCoefs& problem_coefs,
    const std::string& differencing_form)
  : d_is_initialized(false),
    d_problem_coefs(problem_coefs),
    d_differencing_form(differencing_form),
    d_refine_alg(NULL),
    d_refine_op(NULL),
    d_refine_scheds(),
    d_hierarchy(NULL),
    d_coarsest_ln(-1),
    d_finest_ln(-1),
    d_U_var(NULL),
    d_U_scratch_idx(-1)
{
#if (NDIM != 2)
    TBOX_ERROR("INSStaggeredCenteredConvectiveOperator::INSStaggeredCenteredConvectiveOperator():\n"
               << "  incomplete implementation for NDIM != 2\n");
#endif

    if (d_differencing_form != "divergence" &&
        d_differencing_form != "advection" &&
        d_differencing_form != "skew-symmetric")
    {
        TBOX_ERROR("INSStaggeredCenteredConvectiveOperator::INSStaggeredCenteredConvectiveOperator():\n"
                   << "  invalid differencing form: " << d_differencing_form << " \n"
                   << "  valid choices are: divergence, advection, skew-symmetric\n");
    }

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> context = var_db->getContext("INSStaggeredCenteredConvectiveOperator::CONTEXT");

    const std::string U_var_name = "INSStaggeredCenteredConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var.isNull())
    {
        d_U_var = new SAMRAI::pdat::SideVariable<NDIM,double>(U_var_name);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, SAMRAI::hier::IntVector<NDIM>(GADVECTG));
    }
    else
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_U_scratch_idx >= 0);
#endif
    return;
}// INSStaggeredCenteredConvectiveOperator

INSStaggeredCenteredConvectiveOperator::~INSStaggeredCenteredConvectiveOperator()
{
    deallocateOperatorState();
    return;
}// ~INSStaggeredCenteredConvectiveOperator

void
INSStaggeredCenteredConvectiveOperator::applyConvectiveOperator(
    const int U_idx,
    const int N_idx)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSStaggeredCenteredConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to applyConvectiveOperator\n");
    }

    // Setup communications schedules.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > refine_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    refine_alg->registerRefine(d_U_scratch_idx, // destination
                               U_idx,           // source
                               d_U_scratch_idx, // temporary work space
                               d_refine_op);

    // Compute the convective derivative.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        refine_alg->resetSchedule(d_refine_scheds[ln]);
        d_refine_scheds[ln]->fillData(0.0);
        d_refine_alg->resetSchedule(d_refine_scheds[ln]);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

            const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::hier::IntVector<NDIM>& patch_lower = patch_box.lower();
            const SAMRAI::hier::IntVector<NDIM>& patch_upper = patch_box.upper();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > N_data = patch->getPatchData(N_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_data = patch->getPatchData(d_U_scratch_idx);

            const SAMRAI::hier::IntVector<NDIM>& N_ghosts = N_data->getGhostCellWidth();
            const SAMRAI::hier::IntVector<NDIM>& U_ghosts = U_data->getGhostCellWidth();

            if (d_differencing_form == "divergence")
            {
                NAVIER_STOKES_STAGGERED_DIVERGENCE_DERIVATIVE_FC(
                    dx,
#if (NDIM == 2)
                    patch_lower(0), patch_upper(0),
                    patch_lower(1), patch_upper(1),
                    U_ghosts(0), U_ghosts(1),
                    U_data->getPointer(0), U_data->getPointer(1),
                    N_ghosts(0), N_ghosts(1),
                    N_data->getPointer(0), N_data->getPointer(1)
#endif
#if (NDIM == 3)
                    patch_lower(0), patch_upper(0),
                    patch_lower(1), patch_upper(1),
                    patch_lower(2), patch_upper(2),
                    U_ghosts(0), U_ghosts(1), U_ghosts(2),
                    U_data->getPointer(0), U_data->getPointer(1), U_data->getPointer(2),
                    N_ghosts(0), N_ghosts(1), N_ghosts(2),
                    N_data->getPointer(0), N_data->getPointer(1), N_data->getPointer(2)
#endif
                                                                  );
            }
            else if (d_differencing_form == "advection")
            {
                NAVIER_STOKES_STAGGERED_ADVECTION_DERIVATIVE_FC(
                    dx,
#if (NDIM == 2)
                    patch_lower(0), patch_upper(0),
                    patch_lower(1), patch_upper(1),
                    U_ghosts(0), U_ghosts(1),
                    U_data->getPointer(0), U_data->getPointer(1),
                    N_ghosts(0), N_ghosts(1),
                    N_data->getPointer(0), N_data->getPointer(1)
#endif
#if (NDIM == 3)
                    patch_lower(0), patch_upper(0),
                    patch_lower(1), patch_upper(1),
                    patch_lower(2), patch_upper(2),
                    U_ghosts(0), U_ghosts(1), U_ghosts(2),
                    U_data->getPointer(0), U_data->getPointer(1), U_data->getPointer(2),
                    N_ghosts(0), N_ghosts(1), N_ghosts(2),
                    N_data->getPointer(0), N_data->getPointer(1), N_data->getPointer(2)
#endif
                                                                 );
            }
            else if (d_differencing_form == "skew-symmetric")
            {
                NAVIER_STOKES_STAGGERED_SKEW_SYMMETRIC_DERIVATIVE_FC(
                    dx,
#if (NDIM == 2)
                    patch_lower(0), patch_upper(0),
                    patch_lower(1), patch_upper(1),
                    U_ghosts(0), U_ghosts(1),
                    U_data->getPointer(0), U_data->getPointer(1),
                    N_ghosts(0), N_ghosts(1),
                    N_data->getPointer(0), N_data->getPointer(1)
#endif
#if (NDIM == 3)
                    patch_lower(0), patch_upper(0),
                    patch_lower(1), patch_upper(1),
                    patch_lower(2), patch_upper(2),
                    U_ghosts(0), U_ghosts(1), U_ghosts(2),
                    U_data->getPointer(0), U_data->getPointer(1), U_data->getPointer(2),
                    N_ghosts(0), N_ghosts(1), N_ghosts(2),
                    N_data->getPointer(0), N_data->getPointer(1), N_data->getPointer(2)
#endif
                                                                      );
            }
        }
    }
    return;
}// applyConvectiveOperator

void
INSStaggeredCenteredConvectiveOperator::apply(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y)
{
    // Initialize the operator (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeOperatorState(x,y);

    // Get the vector components.
    const int U_idx = x.getComponentDescriptorIndex(0);
    const int N_idx = y.getComponentDescriptorIndex(0);

    // Compute the action of the operator.
    applyConvectiveOperator(U_idx, N_idx);

    // Deallocate the operator (if necessary).
    if (deallocate_at_completion) deallocateOperatorState();
    return;
}// apply

void
INSStaggeredCenteredConvectiveOperator::initializeOperatorState(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out)
{
    if (d_is_initialized) deallocateOperatorState();

    // Get the hierarchy configuration.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#endif
    // Setup the refine algorithm, operator, patch strategy, and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    d_refine_op = grid_geom->lookupRefineOperator(d_U_var, "CONSERVATIVE_LINEAR_REFINE");
    d_refine_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    d_refine_alg->registerRefine(d_U_scratch_idx,                   // destination
                                 in.getComponentDescriptorIndex(0), // source
                                 d_U_scratch_idx,                   // temporary work space
                                 d_refine_op);
    d_refine_strategy = new IBTK::CartExtrapPhysBdryOp(d_U_scratch_idx, BDRY_EXTRAP_TYPE);
    d_refine_scheds.resize(d_finest_ln+1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_refine_scheds[ln] = d_refine_alg->createSchedule(level, ln-1, d_hierarchy, d_refine_strategy);
    }

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_U_scratch_idx))
        {
            level->allocatePatchData(d_U_scratch_idx);
        }
    }
    d_is_initialized = true;
    return;
}// initializeOperatorState

void
INSStaggeredCenteredConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_U_scratch_idx))
        {
            level->deallocatePatchData(d_U_scratch_idx);
        }
    }

    // Deallocate the refine algorithm, operator, patch strategy, and schedules.
    d_refine_op.setNull();
    d_refine_alg.setNull();
    d_refine_strategy.setNull();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_refine_scheds[ln].setNull();
    }
    d_refine_scheds.clear();

    d_is_initialized = false;
    return;
}// deallocateOperatorState

void
INSStaggeredCenteredConvectiveOperator::enableLogging(
    bool enabled)
{
    // intentionally blank
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredCenteredConvectiveOperator>;

//////////////////////////////////////////////////////////////////////////////
