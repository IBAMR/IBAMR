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
#include <ibamr/RNG.h>
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/SideDataSynchronization.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

// C++ STDLIB INCLUDES
#include <limits>

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
            RNG::genrandn(&data(i(),depth));
        }
    }
    return;
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffStochasticForcing::AdvDiffStochasticForcing(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<CellVariable<NDIM,double> > C_var,
    const AdvDiffHierarchyIntegrator* const adv_diff_solver)
    : d_object_name(object_name),
      d_C_var(C_var),
      d_adv_diff_solver(adv_diff_solver),
      d_std(std::numeric_limits<double>::quiet_NaN()),
      d_num_rand_vals(0),
      d_weights(),
      d_dirichlet_bc_scaling(sqrt(2.0)),
      d_neumann_bc_scaling(0.0),
      d_context(NULL),
      d_F_sc_var(NULL),
      d_F_sc_idx(-1),
      d_F_sc_idxs()
{
    if (!input_db.isNull())
    {
        if (input_db->keyExists("std")) d_std = input_db->getDouble("std");
        if (input_db->keyExists("num_rand_vals")) d_num_rand_vals = input_db->getInteger("num_rand_vals");
        int k = 0;
        std::string key_name = "weights_0";
        while (input_db->keyExists(key_name))
        {
            d_weights.push_back(input_db->getDoubleArray(key_name));
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_weights.back().size() == d_num_rand_vals);
#endif
            ++k;
            std::ostringstream stream;
            stream << "weights_" << k;
            key_name = stream.str();
        }
        if (input_db->keyExists("dirichlet_bc_scaling")) d_dirichlet_bc_scaling = input_db->getDouble("dirichlet_bc_scaling");
        if (input_db->keyExists("neumann_bc_scaling")) d_neumann_bc_scaling = input_db->getDouble("neumann_bc_scaling");
    }

    // Determine the number of components that need to be allocated.
    Pointer<CellDataFactory<NDIM,double> > C_factory = d_C_var->getPatchDataFactory();
    const int C_depth = C_factory->getDefaultDepth();

    // Setup variables and variable context objects.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const std::string& name = C_var->getName();
    d_context = var_db->getContext("AdvDiffStochasticForcing::CONTEXT::"+name);
    d_F_sc_var = new SideVariable<NDIM,double>("AdvDiffStochasticForcing::F_sc::"+name, C_depth);
    static const IntVector<NDIM> ghosts_sc = 0;
    d_F_sc_idx = var_db->registerVariableAndContext(d_F_sc_var, d_context, ghosts_sc);
    for (int k = 0; k < d_num_rand_vals; ++k) d_F_sc_idxs.push_back(var_db->registerClonedPatchDataIndex(d_F_sc_var, d_F_sc_idx));
    return;
}// AdvDiffStochasticForcing

AdvDiffStochasticForcing::~AdvDiffStochasticForcing()
{
    // intentionally blank
    return;
}// ~AdvDiffStochasticForcing

bool
AdvDiffStochasticForcing::isTimeDependent() const
{
    return true;
}// isTimeDependent

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
        if (!level->checkAllocated(d_F_sc_idx)) level->allocatePatchData(d_F_sc_idx);
        for (int k = 0; k < d_num_rand_vals; ++k) if (!level->checkAllocated(d_F_sc_idxs[k])) level->allocatePatchData(d_F_sc_idxs[k]);
    }

    // Generate random components.
    const int cycle_num = d_adv_diff_solver->getCurrentCycleNumber();
    if (!initial_time && cycle_num == 0)
    {
        for (int k = 0; k < d_num_rand_vals; ++k)
        {
            for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
            {
                Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<SideData<NDIM,double> > F_sc_data = patch->getPatchData(d_F_sc_idxs[k]);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        genrandn(F_sc_data->getArrayData(d), SideGeometry<NDIM>::toSideBox(F_sc_data->getBox(),d));
                    }
                }
            }
        }
    }

    // Set random values for the present cycle as weighted combinations of the
    // generated random values.
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(cycle_num < static_cast<int>(d_weights.size()));
#endif
    const Array<double>& weights = d_weights[cycle_num];
    HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM,double> > hier_sc_data_ops = hier_data_ops_manager->getOperationsDouble(d_F_sc_var, hierarchy, /*get_unique*/ true);
    hier_sc_data_ops->setToScalar(d_F_sc_idx, 0.0);
    for (int k = 0; k < d_num_rand_vals; ++k) hier_sc_data_ops->axpy(d_F_sc_idx, weights[k], d_F_sc_idxs[k], d_F_sc_idx);

    // Modify the flux values (if necessary).
    Pointer<CellDataFactory<NDIM,double> > C_factory = d_C_var->getPatchDataFactory();
    const int C_depth = C_factory->getDefaultDepth();
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs = d_adv_diff_solver->d_Q_bc_coef.find(d_C_var)->second;
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM,double> > F_sc_data = patch->getPatchData(d_F_sc_idx);

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
                    RobinBcCoefStrategy<NDIM>* bc_coef = bc_coefs[d];
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
                            (*F_sc_data)(s_i,d) *= d_dirichlet_bc_scaling;
                        }
                        else
                        {
                            (*F_sc_data)(s_i,d) = d_neumann_bc_scaling;
                        }
                    }
                }
            }
        }
    }

    // Synchronize side-centered values.
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent synch_component(d_F_sc_idx);
    SideDataSynchronization synch_data_op;
    synch_data_op.initializeOperatorState(synch_component, hierarchy);
    synch_data_op.synchronizeData(data_time);

    // Compute div F on each patch level.
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
    // Compute div F on each patch.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        setDataOnPatch(data_idx, var, level->getPatch(p()), data_time, initial_time, level);
    }
    return;
}// setDataOnPatchLevel

void
AdvDiffStochasticForcing::setDataOnPatch(
    const int data_idx,
    Pointer<Variable<NDIM> > /*var*/,
    Pointer<Patch<NDIM> > patch,
    const double /*data_time*/,
    const bool initial_time,
    Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    Pointer<CellData<NDIM,double> > divF_cc_data = patch->getPatchData(data_idx);
    divF_cc_data->fillAll(0.0);
    if (initial_time) return;
    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    double dV = 1.0; for (unsigned int d = 0; d < NDIM; ++d) dV *= dx[d];
    const double kappa = d_adv_diff_solver->d_Q_diffusion_coef.find(d_C_var)->second;
    const double dt = d_adv_diff_solver->getCurrentTimeStepSize();
    const double scale = d_std*sqrt(2.0*kappa/(dt*dV));
    Pointer<SideData<NDIM,double> > F_sc_data = patch->getPatchData(d_F_sc_idx);
    Pointer<CellDataFactory<NDIM,double> > C_factory = d_C_var->getPatchDataFactory();
    const int C_depth = C_factory->getDefaultDepth();
    for (int d = 0; d < C_depth; ++d)
    {
        for (BoxIterator<NDIM> i(patch_box); i; i++)
        {
            const Index<NDIM>& ic = i();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                SideIndex<NDIM> is_lower(ic, axis, SideIndex<NDIM>::Lower);
                SideIndex<NDIM> is_upper(ic, axis, SideIndex<NDIM>::Upper);
                (*divF_cc_data)(ic,d) += scale*((*F_sc_data)(is_upper,d) - (*F_sc_data)(is_lower,d))/dx[axis];
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
