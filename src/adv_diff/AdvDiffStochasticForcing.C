// Filename: AdvDiffStochasticForcing.C
// Created on 29 Apr 2011 by Boyce Griffith
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

#include "AdvDiffStochasticForcing.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>
#include <ibamr/rng.h>

// IBTK INCLUDES
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/SideDataSynchronization.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>
#include <HierarchyCellDataOpsReal.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
void
genrandn(
    ArrayData<NDIM,double>& data,
    const Box<NDIM>& box)
{
    for (int depth = 0; depth < data.getDepth(); ++depth)
    {
        for (Box<NDIM>::Iterator i(box); i; i++)
        {
            ::genrandn(&data(i(),depth));
        }
    }
    return;
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffStochasticForcing::AdvDiffStochasticForcing(
    Pointer<CellVariable<NDIM,double> > C_var,
    const INSStaggeredHierarchyIntegrator* const fluid_solver,
    const AdvDiffSourceTermType& eval_type,
    const std::string& f_expression,
    const Array<double>& g)
    : d_object_name("AdvDiffStochasticForcing"),
      d_C_var(C_var),
      d_C_bc_coef(),
      d_f_parser(),
      d_g(g),
      d_fluid_solver(fluid_solver),
      d_eval_type(eval_type),
      d_context(NULL),
      d_W_sc_var(NULL),
      d_W_sc_idx(-1),
      d_C_cc_idx(-1),
      d_kappa(0.0),
      d_dt(0.0),
      d_std(0.0),
      d_regen_rand_cycle()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_g.size() == 0 || d_g.size() == NDIM);
#endif
    // Setup the expression parser.
    d_f_parser.SetExpr(f_expression);

    // Determine the number of components that need to be allocated.
    Pointer<CellDataFactory<NDIM,double> > C_factory = C_var->getPatchDataFactory();
    const int C_depth = C_factory->getDefaultDepth();

    // Setup variables and variable context objects.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const std::string& name = C_var->getName();
    d_context = var_db->getContext("AdvDiffStochasticForcing::CONTEXT::"+name);
    d_W_sc_var = new SideVariable<NDIM,double>("AdvDiffStochasticForcing::W_sc::"+name, C_depth);
    static const IntVector<NDIM> ghosts_sc = 0;
    d_W_sc_idx = var_db->registerVariableAndContext(d_W_sc_var, d_context, ghosts_sc);
    static const IntVector<NDIM> ghosts_cc = 1;
    d_C_cc_idx = var_db->registerVariableAndContext(C_var, d_context, ghosts_cc);

    // Set a default value for regen_rand_cycle.
    d_regen_rand_cycle.resizeArray(2);
    d_regen_rand_cycle[0] = 1;
    d_regen_rand_cycle[1] = 0;
    return;
}// AdvDiffStochasticForcing

AdvDiffStochasticForcing::~AdvDiffStochasticForcing()
{
    // intentionally blank
    return;
}// ~AdvDiffStochasticForcing

void
AdvDiffStochasticForcing::setPhysicalBcCoefs(
    RobinBcCoefStrategy<NDIM>* C_bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1,C_bc_coef));
    return;
}// setPhysicalBcCoefs

void
AdvDiffStochasticForcing::setPhysicalBcCoefs(
    std::vector<RobinBcCoefStrategy<NDIM>*> C_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    Pointer<CellDataFactory<NDIM,double> > C_factory = d_C_var->getPatchDataFactory();
    const int C_depth = C_factory->getDefaultDepth();
    TBOX_ASSERT(C_depth == int(C_bc_coef.size()));
#endif
    d_C_bc_coef = C_bc_coef;
    return;
}// setPhysicalBcCoefs

bool
AdvDiffStochasticForcing::isTimeDependent() const
{
    return true;
}// isTimeDependent

void
AdvDiffStochasticForcing::setDiffusionCoefficient(
    const double kappa)
{
    d_kappa = kappa;
    return;
}// setDiffusionCoefficient

void
AdvDiffStochasticForcing::setDt(
    const double dt)
{
    d_dt = dt;
    return;
}// setDt

void
AdvDiffStochasticForcing::setStd(
    const double std)
{
    d_std = std;
    return;
}// setStd

void
AdvDiffStochasticForcing::setRegenRandCycle(
    const Array<int>& regen_rand_cycle)
{
    d_regen_rand_cycle = regen_rand_cycle;
    return;
}// setRegenRandCycle

void
AdvDiffStochasticForcing::setFExpression(
    const std::string& f_expression)
{
    d_f_parser.SetExpr(f_expression);
    return;
}// setFExpression

void
AdvDiffStochasticForcing::setG(
    const Array<double>& g)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(g.size() == 0 || g.size() == NDIM);
#endif
    d_g = g;
    return;
}// setG

void
AdvDiffStochasticForcing::setAdvDiffSourceTermType(
    const AdvDiffSourceTermType& eval_type)
{
    d_eval_type = eval_type;
    return;
}// setAdvDiffSourceTermType

void
AdvDiffStochasticForcing::setDataOnPatchHierarchy(
    const int data_idx,
    Pointer<Variable<NDIM> > var,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    const double data_time,
    const bool initial_time,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    if (cycle_num >= int(d_regen_rand_cycle.size()))
    {
        TBOX_ERROR("AdvDiffStochasticForcing::setDataOnPatchHierarchy():\n"
                   << "  cycle number = " << cycle_num << "\n"
                   << "  but regen_rand_cycle.size() = " << d_regen_rand_cycle.size() << "\n");
    }

    const int coarsest_ln =
        (coarsest_ln_in == -1
         ? 0
         : coarsest_ln_in);
    const int finest_ln =
        (finest_ln_in == -1
         ? hierarchy->getFinestLevelNumber()
         : finest_ln_in);

    // Allocate data to store components of the stochastic stress components.
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
        if (!level->checkAllocated(d_W_sc_idx)) level->allocatePatchData(d_W_sc_idx);
        if (!level->checkAllocated(d_C_cc_idx)) level->allocatePatchData(d_C_cc_idx);
    }

    if (!initial_time && d_std != 0.0 && cycle_num >= 0 && d_regen_rand_cycle[cycle_num] == 1)
    {
        // Generate random components.
        Pointer<CellDataFactory<NDIM,double> > C_factory = d_C_var->getPatchDataFactory();
        const int C_depth = C_factory->getDefaultDepth();
        for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM,double> > W_sc_data = patch->getPatchData(d_W_sc_idx);
                for (int d = 0; d < NDIM; ++d)
                {
                    genrandn(W_sc_data->getArrayData(d), SideGeometry<NDIM>::toSideBox(W_sc_data->getBox(),d));
                }

                const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                if (!pgeom->getTouchesRegularBoundary()) continue;

                const Box<NDIM>& patch_box = patch->getBox();
                Box<NDIM> side_boxes[NDIM];
                for (int d = 0; d < NDIM; ++d)
                {
                    side_boxes[d] = SideGeometry<NDIM>::toSideBox(patch_box,d);
                }
                const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    const IntVector<NDIM> gcw_to_fill = 1;
                    const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                    const int location_index   = bdry_box.getLocationIndex();
                    const int bdry_normal_axis = location_index / 2;
                    const BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), location_index);
                    const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                    Pointer<ArrayData<NDIM,double> > acoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM,double> > bcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM,double> > gcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);

                    // Set the boundary condition coefficients and use them to
                    // rescale the stochastic fluxes.
                    for (int d = 0; d < C_depth; ++d)
                    {
                        RobinBcCoefStrategy<NDIM>* bc_coef = d_C_bc_coef[d];
                        bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, *patch, trimmed_bdry_box, data_time);
                        for (Box<NDIM>::Iterator it(bc_coef_box*side_boxes[bdry_normal_axis]); it; it++)
                        {
                            const Index<NDIM>& i = it();
                            const double& alpha = (*acoef_data)(i,0);
                            const double& beta  = (*bcoef_data)(i,0);
                            const bool dirichlet_bc = (alpha != 0.0 && beta == 0.0);

                            SideIndex<NDIM> s_i(i,bdry_normal_axis,0);
                            if (dirichlet_bc)
                            {
                                (*W_sc_data)(s_i,d) *= sqrt(2.0);
                            }
                            else
                            {
                                (*W_sc_data)(s_i,d) = 0.0;
                            }
                        }
                    }
                }
            }
        }

        // Synchronize side-centered values.
        typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
        SynchronizationTransactionComponent synch_component(d_W_sc_idx);
        SideDataSynchronization synch_data_op;
        synch_data_op.initializeOperatorState(synch_component, hierarchy);
        synch_data_op.synchronizeData(data_time);

        // Compute C_half.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int C_current_idx = var_db->mapVariableAndContextToIndex(d_C_var, d_fluid_solver->getCurrentContext());
        const int C_new_idx = var_db->mapVariableAndContextToIndex(d_C_var, d_fluid_solver->getNewContext());
        bool use_C_current = true;
        bool use_C_new = false;
        if (cycle_num > 0)
        {
            if (d_eval_type == TRAPEZOIDAL_RULE)
            {
                use_C_current = false;
                use_C_new = true;
            }
            else if (d_eval_type == MIDPOINT_RULE)
            {
                use_C_current = true;
                use_C_new = true;
            }
        }

        HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyCellDataOpsReal<NDIM,double> > hier_cc_data_ops = hier_ops_manager->getOperationsDouble(d_C_var, hierarchy);
        if (use_C_current && use_C_new)
        {
            hier_cc_data_ops->linearSum(d_C_cc_idx, 0.5, C_current_idx, 0.5, C_new_idx);
        }
        else if (use_C_current)
        {
            hier_cc_data_ops->copyData(d_C_cc_idx, C_current_idx);
        }
        else if (use_C_new)
        {
            hier_cc_data_ops->copyData(d_C_cc_idx, C_new_idx);
        }
        else
        {
            hier_cc_data_ops->setToScalar(d_C_cc_idx, 0.0);
        }
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_fill_components(1);
        ghost_fill_components[0] = InterpolationTransactionComponent(d_C_cc_idx, "NONE", "NONE", false, d_C_bc_coef);
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_fill_components, hierarchy);
        ghost_fill_op.fillData(data_time);
    }

    // Compute div W on each patch level.
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }
    return;
}// setDataOnPatchHierarchy

void
AdvDiffStochasticForcing::setDataOnPatchLevel(
    const int data_idx,
    Pointer<Variable<NDIM> > var,
    Pointer<PatchLevel<NDIM> > level,
    const double data_time,
    const bool initial_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!level.isNull());
#endif
    // Compute div W on each patch.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        setDataOnPatch(data_idx, var, level->getPatch(p()), data_time, initial_time, level);
    }
    return;
}// setDataOnPatchLevel

void
AdvDiffStochasticForcing::setDataOnPatch(
    const int data_idx,
    Pointer<Variable<NDIM> > var,
    Pointer<Patch<NDIM> > patch,
    const double data_time,
    const bool initial_time,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    Pointer<CellData<NDIM,double> > divW_cc_data = patch->getPatchData(data_idx);
    divW_cc_data->fillAll(0.0);
    if (d_g.size() != 0)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(d_g.size() == NDIM);
#endif
        const Box<NDIM>& patch_box = patch->getBox();

        Pointer<SideData<NDIM,double> > u_current_data = patch->getPatchData(
            d_fluid_solver->getVelocityVar(), d_fluid_solver->getCurrentContext());
        Pointer<SideData<NDIM,double> > u_new_data = patch->getPatchData(
            d_fluid_solver->getVelocityVar(), d_fluid_solver->getNewContext());

        const bool use_u_current_data = !u_current_data.isNull();
        const bool use_u_new_data = !u_new_data.isNull();

        const double u_current_factor = (use_u_current_data && use_u_new_data ? 0.5 : use_u_current_data ? 1.0 : 0.0);
        const double u_new_factor     = (use_u_current_data && use_u_new_data ? 0.5 : use_u_current_data ? 0.0 : 1.0);

        Pointer<CellDataFactory<NDIM,double> > C_factory = d_C_var->getPatchDataFactory();
        const int C_depth = C_factory->getDefaultDepth();
        for (int d = 0; d < C_depth; ++d)
        {
            for (BoxIterator<NDIM> i(patch_box); i; i++)
            {
                const Index<NDIM>& ic = i();
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    SideIndex<NDIM> is_lower(ic, axis, 0);
                    SideIndex<NDIM> is_upper(ic, axis, 1);
                    if (u_current_factor != 0.0)
                    {
                        const double& u_upper = (*u_current_data)(is_upper);
                        const double& u_lower = (*u_current_data)(is_lower);
                        (*divW_cc_data)(ic,d) -= d_g[axis]*u_current_factor*0.5*(u_lower+u_upper);
                    }
                    if (u_new_factor != 0.0)
                    {
                        const double& u_upper = (*u_new_data)(is_upper);
                        const double& u_lower = (*u_new_data)(is_lower);
                        (*divW_cc_data)(ic,d) -= d_g[axis]*u_new_factor*0.5*(u_lower+u_upper);
                    }
                }
            }
        }
    }

    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    if (cycle_num >= int(d_regen_rand_cycle.size()))
    {
        TBOX_ERROR("AdvDiffStochasticForcing::setDataOnPatch():\n"
                   << "  cycle number = " << cycle_num << "\n"
                   << "  but regen_rand_cycle.size() = " << d_regen_rand_cycle.size() << "\n");
    }

    if (initial_time || d_std == 0.0 || (cycle_num >= 0 && d_regen_rand_cycle[cycle_num] == -1)) return;

    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    double dV = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        dV *= dx[d];
    }
    const double scale = d_std*sqrt(2.0*d_kappa/(d_dt*dV));

    double C, f_C;
    d_f_parser.DefineVar("c", &C);
    d_f_parser.DefineVar("C", &C);

    Pointer<SideData<NDIM,double> > W_sc_data = patch->getPatchData(d_W_sc_idx);
    Pointer<CellData<NDIM,double> > C_cc_data = patch->getPatchData(d_C_cc_idx);
    Pointer<CellDataFactory<NDIM,double> > C_factory = d_C_var->getPatchDataFactory();
    const int C_depth = C_factory->getDefaultDepth();
    for (int d = 0; d < C_depth; ++d)
    {
        for (BoxIterator<NDIM> i(patch_box); i; i++)
        {
            const Index<NDIM>& ic = i();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                Index<NDIM> ic_lower(ic); ic_lower(axis) -= 1;
                Index<NDIM> ic_upper(ic); ic_upper(axis) += 1;
                SideIndex<NDIM> is_lower(ic, axis, 0);
                SideIndex<NDIM> is_upper(ic, axis, 1);

                C = 0.5*((*C_cc_data)(ic_upper,d) + (*C_cc_data)(ic,d));
                f_C = d_f_parser.Eval();
                if (f_C < 0.0)
                {
                    pout << "WARNING: f(C) = " << f_C << " < 0; C = " << C << "\n";
                    f_C = 0.0;
                }
                const double scale_upper = sqrt(f_C)*scale;

                C = 0.5*((*C_cc_data)(ic_lower,d) + (*C_cc_data)(ic,d));
                f_C = d_f_parser.Eval();
                if (f_C < 0.0)
                {
                    pout << "WARNING: f(C) = " << f_C << " < 0; C = " << C << "\n";
                    f_C = 0.0;
                }
                const double scale_lower = sqrt(f_C)*scale;

                (*divW_cc_data)(ic,d) += (scale_upper*(*W_sc_data)(is_upper,d) - scale_lower*(*W_sc_data)(is_lower,d))/dx[axis];
            }
        }
    }
    return;
}// setDataOnPatch

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::AdvDiffStochasticForcing>;

//////////////////////////////////////////////////////////////////////////////
