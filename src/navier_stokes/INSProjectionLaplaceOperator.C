// Filename: INSProjectionLaplaceOperator.C
// Last modified: <29.Aug.2007 01:57:39 griffith@box221.cims.nyu.edu>
// Created on 24 Oct 2003 by Boyce Griffith (griffith@mstu1.cims.nyu.edu)

#include "INSProjectionLaplaceOperator.h"

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
#include <stools/RefinePatchStrategySet.h>

// SAMRAI INCLUDES
#include <CellDataFactory.h>
#include <CellVariable.h>
#include <RefineOperator.h>
#include <Variable.h>
#include <tbox/Database.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghost cells required for scratch data.
static const int SCRATCH_GHOSTS = 1;

// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_apply;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_apply_add;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_operator_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSProjectionLaplaceOperator::INSProjectionLaplaceOperator(
    const std::string& object_name,
    INSProjectionBcCoef* const bc_coef,
    const bool homogeneous_bc)
    : LinearOperator(true),
      d_object_name(object_name),
      d_is_initialized(false),
      d_apply_time(0.0),
      d_homogeneous_bc_fill_alg(NULL),
      d_homogeneous_bc_fill_scheds(),
      d_inhomogeneous_bc_fill_alg(NULL),
      d_inhomogeneous_bc_fill_scheds(),
      d_no_fill(),
      d_bc_refine_strategy(NULL),
      d_bc_op(NULL),
      d_scratch_idx(-1),
      d_x(NULL),
      d_b(NULL),
      d_correcting_rhs(false),
      d_poisson_spec(d_object_name+"::Poisson spec"),
      d_bc_coef(NULL),
      d_homogeneous_bc(false),
      d_hier_math_ops(),
      d_hier_math_ops_external(false),
      d_hierarchy(),
      d_grid_geom(),
      d_coarsest_ln(-1),
      d_finest_ln(-1)
{
    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Initialize the boundary conditions objects.
    setHomogeneousBc(homogeneous_bc);
    setPhysicalBcCoef(bc_coef);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_apply = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSProjectionLaplaceOperator::apply()");
        t_apply_add = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSProjectionLaplaceOperator::applyAdd()");
        t_initialize_operator_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSProjectionLaplaceOperator::initializeOperatorState()");
        t_deallocate_operator_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::INSProjectionLaplaceOperator::deallocateOperatorState()");
        timers_need_init = false;
    }
    return;
}// INSProjectionLaplaceOperator()

INSProjectionLaplaceOperator::~INSProjectionLaplaceOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}// ~INSProjectionLaplaceOperator()

void
INSProjectionLaplaceOperator::setPhysicalBcCoef(
    INSProjectionBcCoef* const bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(bc_coef != NULL);
#endif
    d_bc_coef = bc_coef;
    return;
}// setPhysicalBcCoef

void
INSProjectionLaplaceOperator::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSProjectionLaplaceOperator::setTime(
    const double time)
{
    d_apply_time = time;
    return;
}// setTime

void
INSProjectionLaplaceOperator::setHierarchyMathOps(
    SAMRAI::tbox::Pointer<STOOLS::HierarchyMathOps> hier_math_ops)
{
    d_hier_math_ops = hier_math_ops;
    d_hier_math_ops_external = !d_hier_math_ops.isNull();
    return;
}// setHierarchyMathOps

void
INSProjectionLaplaceOperator::modifyRhsForInhomogeneousBc(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y)
{
    // Set y := y - A*0, i.e., shift the right-hand-side vector to account for
    // inhomogeneous boundary conditions.
    if (!d_homogeneous_bc)
    {
        d_correcting_rhs = true;
        d_x->setToScalar(0.0);
        apply(*d_x,*d_b);
        y.subtract(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >(&y, false), d_b);
        d_correcting_rhs = false;
    }
    return;
}// modifyRhsForInhomogeneousBc

void
INSProjectionLaplaceOperator::apply(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y)
{
    t_apply->start();

    static const int comp = 0;
    static const int depth = 0;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_is_initialized);
#endif
    // Reset the communications algorithm.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > bdry_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& x_var = x.getComponentVariable(comp);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& y_var = y.getComponentVariable(comp);

    const int x_idx = x.getComponentDescriptorIndex(comp);
    const int y_idx = y.getComponentDescriptorIndex(comp);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > x_cc_var = x_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > y_cc_var = y_var;

    if (x_cc_var.isNull() || y_cc_var.isNull())
    {
        TBOX_ERROR(d_object_name << "::apply()\n"
                   << "  encountered non-cell centered vector components" << endl);
    }

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > x_factory =
        x_cc_var->getPatchDataFactory();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > y_factory =
        y_cc_var->getPatchDataFactory();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!x_factory.isNull());
    assert(!y_factory.isNull());
#endif

    const unsigned x_depth = x_factory->getDefaultDepth();
    const unsigned y_depth = y_factory->getDefaultDepth();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(x_depth == y_depth);
#endif
    if (x_depth != 1 || y_depth != 1)
    {
        TBOX_ERROR(d_object_name << "::apply()\n"
                   << "  each vector component must have data depth == 1" << std::endl);
    }

    const bool using_homogeneous_bc = d_correcting_rhs ? d_homogeneous_bc : true;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > x_bc_fill_op = d_grid_geom->
        lookupRefineOperator(x_var, "CONSTANT_REFINE");

    bdry_fill_alg->registerRefine(d_scratch_idx,  // destination
                                  x_idx,          // source
                                  d_scratch_idx,  // temporary work space
                                  x_bc_fill_op);

    if (!using_homogeneous_bc)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        int w_idx = d_bc_coef->getIntermediateVelocityPatchDataIndex();
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > w_var;
        var_db->mapIndexToVariable(w_idx, w_var);

        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > w_bc_fill_op = d_grid_geom->
            lookupRefineOperator(w_var, "CONSTANT_REFINE");

        bdry_fill_alg->registerRefine(w_idx,  // destination
                                      w_idx,  // source
                                      w_idx,  // temporary work space
                                      w_bc_fill_op);
    }

    d_bc_op->setPhysicalBcCoef(d_bc_coef);
    d_bc_op->setHomogeneousBc(using_homogeneous_bc);
    d_bc_coef->setHomogeneousBc(using_homogeneous_bc);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (using_homogeneous_bc)
        {
            bdry_fill_alg->resetSchedule(d_homogeneous_bc_fill_scheds[ln]);
            d_homogeneous_bc_fill_scheds[ln]->fillData(d_apply_time);
            d_homogeneous_bc_fill_alg->resetSchedule(d_homogeneous_bc_fill_scheds[ln]);
        }
        else
        {
            bdry_fill_alg->resetSchedule(d_inhomogeneous_bc_fill_scheds[ln]);
            d_inhomogeneous_bc_fill_scheds[ln]->fillData(d_apply_time);
            d_inhomogeneous_bc_fill_alg->resetSchedule(d_inhomogeneous_bc_fill_scheds[ln]);
        }
    }

    d_bc_coef->setHomogeneousBc(true);

    // Compute the action of the operator.
    d_hier_math_ops->laplace(
        y_idx, y_cc_var,
        d_poisson_spec, d_scratch_idx, x_cc_var,
        d_no_fill, d_apply_time,
        0.0, -1, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >(NULL),
        depth, depth);

    t_apply->stop();
    return;
}// apply

void
INSProjectionLaplaceOperator::applyAdd(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& z)
{
    t_apply_add->start();

    static const int comp = 0;
    static const int depth = 0;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_is_initialized);
#endif
    // Reset the communications algorithm.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > bdry_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& x_var = x.getComponentVariable(comp);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& y_var = y.getComponentVariable(comp);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& z_var = z.getComponentVariable(comp);

    const int x_idx = x.getComponentDescriptorIndex(comp);
    const int y_idx = y.getComponentDescriptorIndex(comp);
    const int z_idx = z.getComponentDescriptorIndex(comp);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > x_cc_var = x_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > y_cc_var = y_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > z_cc_var = z_var;

    if (x_cc_var.isNull() || y_cc_var.isNull() || z_cc_var.isNull())
    {
        TBOX_ERROR(d_object_name << "::applyAdd()\n"
                   << "  encountered non-cell centered vector components" << endl);
    }

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > x_factory =
        x_cc_var->getPatchDataFactory();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > y_factory =
        y_cc_var->getPatchDataFactory();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > z_factory =
        z_cc_var->getPatchDataFactory();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!x_factory.isNull());
    assert(!y_factory.isNull());
    assert(!z_factory.isNull());
#endif

    const unsigned x_depth = x_factory->getDefaultDepth();
    const unsigned y_depth = y_factory->getDefaultDepth();
    const unsigned z_depth = z_factory->getDefaultDepth();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(x_depth == y_depth);
    assert(x_depth == z_depth);
#endif
    if (x_depth != 1 || y_depth != 1 || z_depth != 1)
    {
        TBOX_ERROR(d_object_name << "::applyAdd()\n"
                   << "  each vector component must have data depth == 1" << std::endl);
    }

    const bool using_homogeneous_bc = d_correcting_rhs ? d_homogeneous_bc : true;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > x_bc_fill_op = d_grid_geom->
        lookupRefineOperator(x_var, "CONSTANT_REFINE");

    bdry_fill_alg->registerRefine(d_scratch_idx,  // destination
                                  x_idx,          // source
                                  d_scratch_idx,  // temporary work space
                                  x_bc_fill_op);

    if (!using_homogeneous_bc)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        int w_idx = d_bc_coef->getIntermediateVelocityPatchDataIndex();
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > w_var;
        var_db->mapIndexToVariable(w_idx, w_var);

        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > w_bc_fill_op = d_grid_geom->
            lookupRefineOperator(w_var, "CONSTANT_REFINE");

        bdry_fill_alg->registerRefine(w_idx,  // destination
                                      w_idx,  // source
                                      w_idx,  // temporary work space
                                      w_bc_fill_op);
    }

    d_bc_op->setPhysicalBcCoef(d_bc_coef);
    d_bc_op->setHomogeneousBc(using_homogeneous_bc);
    d_bc_coef->setHomogeneousBc(using_homogeneous_bc);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (using_homogeneous_bc)
        {
            bdry_fill_alg->resetSchedule(d_homogeneous_bc_fill_scheds[ln]);
            d_homogeneous_bc_fill_scheds[ln]->fillData(d_apply_time);
            d_homogeneous_bc_fill_alg->resetSchedule(d_homogeneous_bc_fill_scheds[ln]);
        }
        else
        {
            bdry_fill_alg->resetSchedule(d_inhomogeneous_bc_fill_scheds[ln]);
            d_inhomogeneous_bc_fill_scheds[ln]->fillData(d_apply_time);
            d_inhomogeneous_bc_fill_alg->resetSchedule(d_inhomogeneous_bc_fill_scheds[ln]);
        }
    }

    d_bc_coef->setHomogeneousBc(true);

    // Compute the action of the operator.
    d_hier_math_ops->laplace(
        z_idx, z_cc_var,
        d_poisson_spec, d_scratch_idx, x_cc_var,
        d_no_fill, d_apply_time,
        1.0, y_idx, y_cc_var,
        depth, depth, depth);

    t_apply_add->stop();
    return;
}// applyAdd

void
INSProjectionLaplaceOperator::initializeOperatorState(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out)
{
    t_initialize_operator_state->start();

    static const int comp = 0;

    // Deallocate the operator state if the solver is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_x = in .cloneVector(in .getName());
    d_b = out.cloneVector(out.getName());

    d_x->allocateVectorData();
    d_b->allocateVectorData();

    // Setup operator state.
    d_hierarchy   = in.getPatchHierarchy();
    d_grid_geom   = d_hierarchy->getGridGeometry();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln   = in.getFinestLevelNumber();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_hierarchy == out.getPatchHierarchy());
    assert(d_coarsest_ln == out.getCoarsestLevelNumber());
    assert(d_finest_ln == out.getFinestLevelNumber());
    assert(in.getNumberOfComponents() == out.getNumberOfComponents() &&
           in.getNumberOfComponents() == 1);
#endif

    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops = new STOOLS::HierarchyMathOps(
            d_object_name+"::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);
    }
    else
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!d_hier_math_ops.isNull());
#endif
    }

    // Setup the communications algorithm.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    d_homogeneous_bc_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    d_inhomogeneous_bc_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    // Get variables and patch data descriptors for temporary data.
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& in_var = in.getComponentVariable(comp);
    const int in_idx = in.getComponentDescriptorIndex(comp);
    d_scratch_idx = var_db->registerClonedPatchDataIndex(in_var, in_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(in_idx >= 0);
    assert(d_scratch_idx >= 0);
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > scratch_factory =
        var_db->getPatchDescriptor()->getPatchDataFactory(d_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!scratch_factory.isNull());
#endif
    scratch_factory->setDefaultGhostCellWidth(SAMRAI::hier::IntVector<NDIM>(SCRATCH_GHOSTS));

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > in_bc_fill_op = d_grid_geom->
        lookupRefineOperator(in_var, "CONSTANT_REFINE");
    d_homogeneous_bc_fill_alg->registerRefine(d_scratch_idx,  // destination
                                              in_idx,         // source
                                              d_scratch_idx,  // temporary work space
                                              in_bc_fill_op);

    d_inhomogeneous_bc_fill_alg->registerRefine(d_scratch_idx,  // destination
                                                in_idx,         // source
                                                d_scratch_idx,  // temporary work space
                                                in_bc_fill_op);

    int w_idx = d_bc_coef->getIntermediateVelocityPatchDataIndex();
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > w_var;
    var_db->mapIndexToVariable(w_idx, w_var);

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > w_bc_fill_op = d_grid_geom->
        lookupRefineOperator(w_var, "CONSTANT_REFINE");

    d_inhomogeneous_bc_fill_alg->registerRefine(w_idx,  // destination
                                                w_idx,  // source
                                                w_idx,  // temporary work space
                                                w_bc_fill_op);

    d_bc_op = new STOOLS::CartRobinPhysBdryOp(d_scratch_idx, d_bc_coef, false);
    d_bc_refine_strategy = d_bc_op;

    // Setup the communications schedules.
    d_no_fill.resize(d_finest_ln+1);
    d_homogeneous_bc_fill_scheds.resize(d_finest_ln+1);
    d_inhomogeneous_bc_fill_scheds.resize(d_finest_ln+1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_homogeneous_bc_fill_scheds[ln] = d_homogeneous_bc_fill_alg->
            createSchedule(level, ln-1, d_hierarchy, d_bc_refine_strategy);
        d_inhomogeneous_bc_fill_scheds[ln] = d_inhomogeneous_bc_fill_alg->
            createSchedule(level, ln-1, d_hierarchy, d_bc_refine_strategy);
    }

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_idx, d_apply_time);
    }

    d_is_initialized = true;

    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
INSProjectionLaplaceOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    t_deallocate_operator_state->start();

    // Deallocate the operator state.
    d_homogeneous_bc_fill_alg.setNull();
    d_inhomogeneous_bc_fill_alg.setNull();
    d_bc_op.setNull();
    d_bc_refine_strategy.setNull();
    d_homogeneous_bc_fill_scheds.clear();
    d_inhomogeneous_bc_fill_scheds.clear();
    d_no_fill.clear();

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_scratch_idx);
    }

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_scratch_idx);

    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    // Delete the solution and rhs vectors.
    d_x->freeVectorComponents();
    d_b->freeVectorComponents();

    d_x.setNull();
    d_b.setNull();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

void
INSProjectionLaplaceOperator::enableLogging(
    bool enabled)
{
    TBOX_WARNING("INSProjectionLaplaceOperator::enableLogging() not supported" << endl);
    return;
}// enableLogging

void
INSProjectionLaplaceOperator::printClassData(
    std::ostream& os) const
{
    TBOX_WARNING("INSProjectionLaplaceOperator::printClassData() not supported" << endl);
    return;
}// printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSProjectionLaplaceOperator>;

//////////////////////////////////////////////////////////////////////////////
