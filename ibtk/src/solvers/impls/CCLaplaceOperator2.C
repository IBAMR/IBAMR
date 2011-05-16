// Filename: CCLaplaceOperator2.C
// Created on 25 Aug 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "CCLaplaceOperator2.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CellDataFactory.h>
#include <CellVariable.h>
#include <HierarchyDataOpsManager.h>
#include <Variable.h>
#include <tbox/Database.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCLaplaceOperator2::CCLaplaceOperator2(
    const std::string& object_name,
    const PoissonSpecifications& poisson_spec,
    RobinBcCoefStrategy<NDIM>* const bc_coef,
    const bool homogeneous_bc)
    : LinearOperator(true),
      d_object_name(object_name),
      d_is_initialized(false),
      d_ncomp(0),
      d_apply_time(0.0),
      d_hier_bdry_fill(NULL),
      d_no_fill(NULL),
      d_scratch_idxs(),
      d_x(NULL),
      d_y(NULL),
      d_correcting_rhs(false),
      d_poisson_spec(d_object_name+"::Poisson spec"),
      d_default_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(
                            d_object_name+"::default_bc_coef", Pointer<Database>(NULL))),
      d_bc_coefs(),
      d_homogeneous_bc(false),
      d_hier_cc_data_ops(),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1)
{
    // Initialize the Poisson specifications.
    setPoissonSpecifications(poisson_spec);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_bc_coef->setBoundaryValue(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setHomogeneousBc(homogeneous_bc);
    setPhysicalBcCoef(bc_coef);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_apply                     = TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator2::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator2::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator2::deallocateOperatorState()");
                 );
    return;
}// CCLaplaceOperator2()

CCLaplaceOperator2::CCLaplaceOperator2(
    const std::string& object_name,
    const PoissonSpecifications& poisson_spec,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    const bool homogeneous_bc)
    : LinearOperator(true),
      d_object_name(object_name),
      d_is_initialized(false),
      d_ncomp(0),
      d_apply_time(0.0),
      d_hier_bdry_fill(NULL),
      d_no_fill(NULL),
      d_x(NULL),
      d_y(NULL),
      d_correcting_rhs(false),
      d_poisson_spec(d_object_name+"::Poisson spec"),
      d_default_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(
                            d_object_name+"::default_bc_coef", Pointer<Database>(NULL))),
      d_bc_coefs(),
      d_homogeneous_bc(false),
      d_hier_cc_data_ops(),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1)
{
    // Initialize the Poisson specifications.
    setPoissonSpecifications(poisson_spec);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_bc_coef->setBoundaryValue(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setHomogeneousBc(homogeneous_bc);
    setPhysicalBcCoefs(bc_coefs);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_apply                     = TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator2::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator2::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator2::deallocateOperatorState()");
                 );
    return;
}// CCLaplaceOperator2()

CCLaplaceOperator2::~CCLaplaceOperator2()
{
    if (d_is_initialized) deallocateOperatorState();
    delete d_default_bc_coef;
    return;
}// ~CCLaplaceOperator2()

void
CCLaplaceOperator2::setPoissonSpecifications(
    const PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
}// setPoissonSpecifications

void
CCLaplaceOperator2::setPhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1,bc_coef));
    return;
}// setPhysicalBcCoef

void
CCLaplaceOperator2::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    d_bc_coefs.resize(bc_coefs.size());
    for (unsigned int l = 0; l < bc_coefs.size(); ++l)
    {
        if (bc_coefs[l] != NULL)
        {
            d_bc_coefs[l] = bc_coefs[l];
        }
        else
        {
            d_bc_coefs[l] = d_default_bc_coef;
        }
    }
    return;
}// setPhysicalBcCoefs

void
CCLaplaceOperator2::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
CCLaplaceOperator2::setTime(
    const double time)
{
    d_apply_time = time;
    return;
}// setTime

void
CCLaplaceOperator2::modifyRhsForInhomogeneousBc(
    SAMRAIVectorReal<NDIM,double>& y)
{
    // Set y := y - A*0, i.e., shift the right-hand-side vector to account for
    // inhomogeneous boundary conditions.
    if (!d_homogeneous_bc)
    {
        d_correcting_rhs = true;
        d_x->setToScalar(0.0);
        apply(*d_x,*d_y);
        y.subtract(Pointer<SAMRAIVectorReal<NDIM,double> >(&y, false), d_y);
        d_correcting_rhs = false;
    }
    return;
}// modifyRhsForInhomogeneousBc

void
CCLaplaceOperator2::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    t_apply->start();

    int ierr;

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_is_initialized);
#endif

    for (int comp = 0; comp < d_ncomp; ++comp)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        const Pointer<Variable<NDIM> >& x_var = x.getComponentVariable(comp);
        const Pointer<Variable<NDIM> >& y_var = y.getComponentVariable(comp);

        Pointer<CellVariable<NDIM,double> > x_cc_var = x_var;
        Pointer<CellVariable<NDIM,double> > y_cc_var = y_var;

        if (x_cc_var.isNull() || y_cc_var.isNull())
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                       << "  encountered non-cell centered vector components" << std::endl);
        }

        Pointer<CellDataFactory<NDIM,double> > x_factory = x_cc_var->getPatchDataFactory();
        TBOX_ASSERT(!x_factory.isNull());
        Pointer<CellDataFactory<NDIM,double> > y_factory = y_cc_var->getPatchDataFactory();
        TBOX_ASSERT(!y_factory.isNull());

        const unsigned int x_depth = x_factory->getDefaultDepth();
        const unsigned int y_depth = y_factory->getDefaultDepth();
        TBOX_ASSERT(x_depth == y_depth);

        if (x_depth != d_bc_coefs.size() || y_depth != d_bc_coefs.size())
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                       << "  each vector component must have data depth == " << d_bc_coefs.size() << "\n"
                       << "  since d_bc_coefs.size() == " << d_bc_coefs.size() << std::endl);
        }
#endif

        const int x_idx = x.getComponentDescriptorIndex(comp);
        const int scratch_idx = d_scratch_idxs[comp];
        d_hier_cc_data_ops->copyData(scratch_idx, x_idx);
    }

    // Fill the ghost cell data.
    const bool homogeneous_bc = d_correcting_rhs ? d_homogeneous_bc : true;
    d_hier_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_bdry_fill->fillData(d_apply_time);

    // Compute the action of the operator.
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        const Pointer<Variable<NDIM> >& x_var = x.getComponentVariable(comp);
        const Pointer<Variable<NDIM> >& y_var = y.getComponentVariable(comp);

        Pointer<CellVariable<NDIM,double> > x_cc_var = x_var;
        Pointer<CellVariable<NDIM,double> > y_cc_var = y_var;

        const int scratch_idx = d_scratch_idxs[comp];
        const int y_idx = y.getComponentDescriptorIndex(comp);

        for (int ln = d_finest_ln; ln >= d_coarsest_ln; --ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            int patch_counter = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
            {
                Mat& A_mat = d_patch_mat  [comp][ln][patch_counter];
                Vec& x_vec = d_patch_x_vec[comp][ln][patch_counter];
                Vec& y_vec = d_patch_y_vec[comp][ln][patch_counter];

                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM,double> > x_data = patch->getPatchData(scratch_idx);
                Pointer<CellData<NDIM,double> > y_data = patch->getPatchData(y_idx);

                ierr = VecPlaceArray(x_vec, x_data->getPointer());  IBTK_CHKERRQ(ierr);
                ierr = VecPlaceArray(y_vec, y_data->getPointer());  IBTK_CHKERRQ(ierr);
                ierr = MatMult(A_mat, x_vec, y_vec);  IBTK_CHKERRQ(ierr);
                ierr = VecResetArray(x_vec);  IBTK_CHKERRQ(ierr);
                ierr = VecResetArray(y_vec);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    t_apply->stop();
    return;
}// apply

void
CCLaplaceOperator2::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    t_initialize_operator_state->start();

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_x = in .cloneVector(in .getName());
    d_y = out.cloneVector(out.getName());

    d_x->allocateVectorData();
    d_y->allocateVectorData();

    // Setup operator state.
    d_hierarchy   = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln   = in.getFinestLevelNumber();

    if (d_coarsest_ln != d_finest_ln)
    {
        TBOX_ERROR("CCLaplaceOperator2::initializeOperatorState()\n"
                   << "  current implementation requires uniform grid problems" << std::endl);
    }

    d_ncomp = in.getNumberOfComponents();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
    TBOX_ASSERT(d_ncomp == out.getNumberOfComponents());
#endif

    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<CellVariable<NDIM,double> > cc_var = new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, d_hierarchy, true);
    d_hier_cc_data_ops->resetLevels(d_coarsest_ln, d_finest_ln);

    // Setup scratch data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_scratch_idxs.resize(d_ncomp);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        // Get variables and patch data descriptors for temporary data.
        const Pointer<Variable<NDIM> >& in_var = in.getComponentVariable(comp);
        const int in_idx = in.getComponentDescriptorIndex(comp);
        d_scratch_idxs[comp] = var_db->registerClonedPatchDataIndex(in_var, in_idx);
        const int scratch_idx = d_scratch_idxs[comp];
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(in_idx >= 0);
        TBOX_ASSERT(scratch_idx >= 0);
#endif
        Pointer<CellDataFactory<NDIM,double> > scratch_factory =
            var_db->getPatchDescriptor()->getPatchDataFactory(scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!scratch_factory.isNull());
        TBOX_ASSERT(scratch_factory->getGhostCellWidth() == CELLG);
#endif
    }

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (int comp = 0; comp < d_ncomp; ++comp)
        {
            const int scratch_idx = d_scratch_idxs[comp];
            level->allocatePatchData(scratch_idx, d_apply_time);
        }
    }

    // Setup the interpolation transaction information.
    Pointer<VariableFillPattern<NDIM> > fill_pattern = NULL;
    if (d_poisson_spec.dIsConstant())
    {
        fill_pattern = new CellNoCornersFillPattern(CELLG, false, true);
    }
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps;
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        InterpolationTransactionComponent component(d_scratch_idxs[comp], DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_bc_coefs, fill_pattern);
        transaction_comps.push_back(component);
    }

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(transaction_comps, d_hierarchy);

    // Initialize patch operator data.
    if (!d_poisson_spec.cIsZero() && !d_poisson_spec.cIsConstant())
    {
        TBOX_ERROR("CCLaplaceOperator2::initializeOperatorState()\n"
                   << "  current implementation requires constant coefficient problems" << std::endl);
    }
    const double C = d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant();

    if (!d_poisson_spec.dIsConstant())
    {
        TBOX_ERROR("CCLaplaceOperator2::initializeOperatorState()\n"
                   << "  current implementation requires constant coefficient problems" << std::endl);
    }
    const double D = d_poisson_spec.getDConstant();

    d_patch_mat  .resize(d_ncomp);
    d_patch_x_vec.resize(d_ncomp);
    d_patch_y_vec.resize(d_ncomp);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        d_patch_mat  [comp].resize(d_finest_ln+1);
        d_patch_x_vec[comp].resize(d_finest_ln+1);
        d_patch_y_vec[comp].resize(d_finest_ln+1);
        for (int ln = d_finest_ln; ln >= d_coarsest_ln; --ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
            d_patch_mat  [comp][ln].resize(num_local_patches);
            d_patch_x_vec[comp][ln].resize(num_local_patches);
            d_patch_y_vec[comp][ln].resize(num_local_patches);

            int patch_counter = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
            {
                Mat& A_mat = d_patch_mat  [comp][ln][patch_counter];
                Vec& x_vec = d_patch_x_vec[comp][ln][patch_counter];
                Vec& y_vec = d_patch_y_vec[comp][ln][patch_counter];

                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM,double> > x_data = d_x->getComponentPatchData(comp,*patch);
                Pointer<CellData<NDIM,double> > y_data = d_y->getComponentPatchData(comp,*patch);

                PETScMatUtilities::constructPatchLaplaceOp(A_mat, C, D, *x_data, *y_data, *patch);
                PETScVecUtilities::constructPatchVecWrapper(x_vec, *x_data);
                PETScVecUtilities::constructPatchVecWrapper(y_vec, *y_data);
            }
        }
    }

    // Indicate the operator is initialized.
    d_is_initialized = true;

    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
CCLaplaceOperator2::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    t_deallocate_operator_state->start();

    // Deallocate patch operator data.
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        for (int ln = d_finest_ln; ln >= d_coarsest_ln; --ln)
        {
            for (std::vector<Mat>::iterator it = d_patch_mat[comp][ln].begin();
                 it != d_patch_mat[comp][ln].end(); ++it)
            {
                Mat& A = *it;
                int ierr = MatDestroy(A);  IBTK_CHKERRQ(ierr);
            }

            for (std::vector<Vec>::iterator it = d_patch_x_vec[comp][ln].begin();
                 it != d_patch_x_vec[comp][ln].end(); ++it)
            {
                Vec& x_vec = *it;
                int ierr = VecDestroy(x_vec);  IBTK_CHKERRQ(ierr);
            }

            for (std::vector<Vec>::iterator it = d_patch_y_vec[comp][ln].begin();
                 it != d_patch_y_vec[comp][ln].end(); ++it)
            {
                Vec& y_vec = *it;
                int ierr = VecDestroy(y_vec);  IBTK_CHKERRQ(ierr);
            }
        }
    }
    d_patch_mat  .clear();
    d_patch_x_vec.clear();
    d_patch_y_vec.clear();

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (ln <= d_hierarchy->getFinestLevelNumber())
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (int comp = 0; comp < d_ncomp; ++comp)
            {
                level->deallocatePatchData(d_scratch_idxs[comp]);
            }
        }
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        const int scratch_idx = d_scratch_idxs[comp];
        var_db->removePatchDataIndex(scratch_idx);
    }

    // Delete the solution and rhs vectors.
    d_x->resetLevels(d_x->getCoarsestLevelNumber(), std::min(d_x->getFinestLevelNumber(),d_hierarchy->getFinestLevelNumber()));
    d_x->freeVectorComponents();

    d_y->resetLevels(d_y->getCoarsestLevelNumber(), std::min(d_y->getFinestLevelNumber(),d_hierarchy->getFinestLevelNumber()));
    d_y->freeVectorComponents();

    d_x.setNull();
    d_y.setNull();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

void
CCLaplaceOperator2::enableLogging(
    bool enabled)
{
    TBOX_WARNING("CCLaplaceOperator2::enableLogging() not supported" << std::endl);
    return;
}// enableLogging

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::CCLaplaceOperator2>;

//////////////////////////////////////////////////////////////////////////////
