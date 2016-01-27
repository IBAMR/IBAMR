// Filename: AdvDiffStochasticForcing.cpp
// Created on 29 Apr 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <math.h>
#include <stddef.h>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/AdvDiffStochasticForcing.h"
#include "ibamr/RNG.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/SideDataSynchronization.h"
#include "muParser.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
void
genrandn(ArrayData<NDIM, double>& data, const Box<NDIM>& box)
{
    for (int depth = 0; depth < data.getDepth(); ++depth)
    {
        for (Box<NDIM>::Iterator i(box); i; i++)
        {
            RNG::genrandn(&data(i(), depth));
        }
    }
    return;
} // genrandn
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffStochasticForcing::AdvDiffStochasticForcing(const std::string& object_name,
                                                   Pointer<Database> input_db,
                                                   Pointer<CellVariable<NDIM, double> > C_var,
                                                   const AdvDiffSemiImplicitHierarchyIntegrator* const adv_diff_solver)
    : d_object_name(object_name),
      d_C_var(C_var),
      d_f_parser(),
      d_adv_diff_solver(adv_diff_solver),
      d_std(std::numeric_limits<double>::quiet_NaN()),
      d_num_rand_vals(0),
      d_weights(),
      d_dirichlet_bc_scaling(sqrt(2.0)),
      d_neumann_bc_scaling(0.0),
      d_context(NULL),
      d_C_cc_var(NULL),
      d_C_current_cc_idx(-1),
      d_C_half_cc_idx(-1),
      d_C_new_cc_idx(-1),
      d_F_sc_var(NULL),
      d_F_sc_idx(-1),
      d_F_sc_idxs()
{
    std::string f_expression = "1.0";
    if (input_db)
    {
        if (input_db->keyExists("std")) d_std = input_db->getDouble("std");
        if (input_db->keyExists("num_rand_vals")) d_num_rand_vals = input_db->getInteger("num_rand_vals");
        int k = 0;
        std::string key_name = "weights_0";
        while (input_db->keyExists(key_name))
        {
            d_weights.push_back(input_db->getDoubleArray(key_name));
#if !defined(NDEBUG)
            TBOX_ASSERT(d_weights.back().size() == d_num_rand_vals);
#endif
            ++k;
            std::ostringstream stream;
            stream << "weights_" << k;
            key_name = stream.str();
        }
        if (input_db->keyExists("dirichlet_bc_scaling"))
            d_dirichlet_bc_scaling = input_db->getDouble("dirichlet_bc_scaling");
        if (input_db->keyExists("neumann_bc_scaling")) d_neumann_bc_scaling = input_db->getDouble("neumann_bc_scaling");
        if (input_db->keyExists("f_expression")) f_expression = input_db->getString("f_expression");
    }
    d_f_parser.SetExpr(f_expression);

    // Determine the number of components that need to be allocated.
    Pointer<CellDataFactory<NDIM, double> > C_factory = d_C_var->getPatchDataFactory();
    const int C_depth = C_factory->getDefaultDepth();

    // Setup variables and variable context objects.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    d_C_cc_var = new CellVariable<NDIM, double>(d_object_name + "::C_cc", C_depth);
    static const IntVector<NDIM> ghosts_cc = 1;
    d_C_current_cc_idx = var_db->registerVariableAndContext(d_C_cc_var, d_context, ghosts_cc);
    d_C_half_cc_idx = var_db->registerClonedPatchDataIndex(d_C_cc_var, d_C_current_cc_idx);
    d_C_new_cc_idx = var_db->registerClonedPatchDataIndex(d_C_cc_var, d_C_current_cc_idx);
    d_F_sc_var = new SideVariable<NDIM, double>(d_object_name + "::F_sc", C_depth);
    static const IntVector<NDIM> ghosts_sc = 0;
    d_F_sc_idx = var_db->registerVariableAndContext(d_F_sc_var, d_context, ghosts_sc);
    for (int k = 0; k < d_num_rand_vals; ++k)
        d_F_sc_idxs.push_back(var_db->registerClonedPatchDataIndex(d_F_sc_var, d_F_sc_idx));
    return;
} // AdvDiffStochasticForcing

AdvDiffStochasticForcing::~AdvDiffStochasticForcing()
{
    // intentionally blank
    return;
} // ~AdvDiffStochasticForcing

bool
AdvDiffStochasticForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
AdvDiffStochasticForcing::setDataOnPatchHierarchy(const int data_idx,
                                                  Pointer<Variable<NDIM> > var,
                                                  Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                  const double data_time,
                                                  const bool initial_time,
                                                  const int coarsest_ln_in,
                                                  const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    const int cycle_num = d_adv_diff_solver->getCurrentCycleNumber();
    if (!initial_time)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(cycle_num >= 0);
#endif
        // Allocate data to store components of the stochastic stress components.
        for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
            level->allocatePatchData(d_C_current_cc_idx);
            level->allocatePatchData(d_C_half_cc_idx);
            level->allocatePatchData(d_C_new_cc_idx);
            level->allocatePatchData(d_F_sc_idx);
            for (int k = 0; k < d_num_rand_vals; ++k) level->allocatePatchData(d_F_sc_idxs[k]);
        }

        // Set concentration value used to compute concentration-dependent flux
        // scaling.
        const double dt = d_adv_diff_solver->getCurrentTimeStepSize();
        const double current_time = d_adv_diff_solver->getIntegratorTime();
        const double half_time = current_time + 0.5 * dt;
        const double new_time = current_time + dt;
        HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_cc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_C_cc_var,
                                                       hierarchy,
                                                       /*get_unique*/ true);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int C_current_idx = var_db->mapVariableAndContextToIndex(d_C_var, d_adv_diff_solver->getCurrentContext());
        const int C_new_idx = var_db->mapVariableAndContextToIndex(d_C_var, d_adv_diff_solver->getNewContext());
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_fill_components(1);
        HierarchyGhostCellInterpolation ghost_fill_op;
        const std::vector<RobinBcCoefStrategy<NDIM>*>& C_bc_coef = d_adv_diff_solver->getPhysicalBcCoefs(d_C_var);
        const TimeSteppingType convective_time_stepping_type =
            d_adv_diff_solver->getConvectiveTimeSteppingType(d_C_var);
        switch (convective_time_stepping_type)
        {
        case FORWARD_EULER:
            hier_cc_data_ops->copyData(d_C_current_cc_idx, C_current_idx);
            ghost_fill_components[0] =
                InterpolationTransactionComponent(d_C_current_cc_idx, "NONE", false, "NONE", "NONE", false, C_bc_coef);
            ghost_fill_op.initializeOperatorState(ghost_fill_components, hierarchy);
            ghost_fill_op.fillData(current_time);
            break;
        case MIDPOINT_RULE:
            hier_cc_data_ops->linearSum(d_C_half_cc_idx, 0.5, C_current_idx, 0.5, C_new_idx);
            ghost_fill_components[0] =
                InterpolationTransactionComponent(d_C_half_cc_idx, "NONE", false, "NONE", "NONE", false, C_bc_coef);
            ghost_fill_op.initializeOperatorState(ghost_fill_components, hierarchy);
            ghost_fill_op.fillData(half_time);
            break;
        case TRAPEZOIDAL_RULE:
            if (cycle_num == 0)
            {
                hier_cc_data_ops->copyData(d_C_current_cc_idx, C_current_idx);
                ghost_fill_components[0] = InterpolationTransactionComponent(
                    d_C_current_cc_idx, "NONE", false, "NONE", "NONE", false, C_bc_coef);
                ghost_fill_op.initializeOperatorState(ghost_fill_components, hierarchy);
                ghost_fill_op.fillData(current_time);
            }
            else
            {
                hier_cc_data_ops->copyData(d_C_new_cc_idx, C_new_idx);
                ghost_fill_components[0] =
                    InterpolationTransactionComponent(d_C_new_cc_idx, "NONE", false, "NONE", "NONE", false, C_bc_coef);
                ghost_fill_op.initializeOperatorState(ghost_fill_components, hierarchy);
                ghost_fill_op.fillData(new_time);
            }
            break;
        default:
            TBOX_ERROR(d_object_name << "::setDataOnPatchHierarchy():\n"
                                     << "  unsupported default convective time stepping type: "
                                     << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                                     << " \n"
                                     << "  valid choices are: FORWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }

        // Generate random components.
        if (cycle_num == 0)
        {
            for (int k = 0; k < d_num_rand_vals; ++k)
            {
                for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
                {
                    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<SideData<NDIM, double> > F_sc_data = patch->getPatchData(d_F_sc_idxs[k]);
                        for (int d = 0; d < NDIM; ++d)
                        {
                            genrandn(F_sc_data->getArrayData(d), SideGeometry<NDIM>::toSideBox(F_sc_data->getBox(), d));
                        }
                    }
                }
            }
        }

// Set random values for the present cycle as weighted combinations of
// the generated random values.
#if !defined(NDEBUG)
        TBOX_ASSERT(cycle_num >= 0 && cycle_num < static_cast<int>(d_weights.size()));
#endif
        const Array<double>& weights = d_weights[cycle_num];
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_sc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_F_sc_var,
                                                       hierarchy,
                                                       /*get_unique*/ true);
        hier_sc_data_ops->setToScalar(d_F_sc_idx, 0.0);
        for (int k = 0; k < d_num_rand_vals; ++k)
            hier_sc_data_ops->axpy(d_F_sc_idx, weights[k], d_F_sc_idxs[k], d_F_sc_idx);

        // Modify the flux values (if necessary).
        Pointer<CellDataFactory<NDIM, double> > C_factory = d_C_var->getPatchDataFactory();
        const int C_depth = C_factory->getDefaultDepth();
        const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs = d_adv_diff_solver->getPhysicalBcCoefs(d_C_var);
        for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > F_sc_data = patch->getPatchData(d_F_sc_idx);

                const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                if (!pgeom->getTouchesRegularBoundary()) continue;

                const Box<NDIM>& patch_box = patch->getBox();
                Box<NDIM> side_boxes[NDIM];
                for (int d = 0; d < NDIM; ++d)
                {
                    side_boxes[d] = SideGeometry<NDIM>::toSideBox(patch_box, d);
                }
                const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
                    PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    const IntVector<NDIM> gcw_to_fill = 1;
                    const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                    const int location_index = bdry_box.getLocationIndex();
                    const int bdry_normal_axis = location_index / 2;
                    const BoundaryBox<NDIM> trimmed_bdry_box(
                        bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);
                    const Box<NDIM> bc_coef_box =
                        PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                    Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                    // Set the boundary condition coefficients and use them to
                    // rescale the stochastic fluxes.
                    for (int d = 0; d < C_depth; ++d)
                    {
                        RobinBcCoefStrategy<NDIM>* bc_coef = bc_coefs[d];
                        bc_coef->setBcCoefs(
                            acoef_data, bcoef_data, gcoef_data, var, *patch, trimmed_bdry_box, data_time);
                        for (Box<NDIM>::Iterator it(bc_coef_box * side_boxes[bdry_normal_axis]); it; it++)
                        {
                            const Index<NDIM>& i = it();
                            const double& alpha = (*acoef_data)(i, 0);
                            const double& beta = (*bcoef_data)(i, 0);
                            const bool dirichlet_bc = (alpha != 0.0 && beta == 0.0);

                            SideIndex<NDIM> s_i(i, bdry_normal_axis, 0);
                            if (dirichlet_bc)
                            {
                                (*F_sc_data)(s_i, d) *= d_dirichlet_bc_scaling;
                            }
                            else
                            {
                                (*F_sc_data)(s_i, d) = d_neumann_bc_scaling;
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
    }

    // Compute div F on each patch level.
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }

    // Deallocate data to store components of the stochastic stress components.
    if (!initial_time)
    {
        for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
            level->deallocatePatchData(d_C_current_cc_idx);
            level->deallocatePatchData(d_C_half_cc_idx);
            level->deallocatePatchData(d_C_new_cc_idx);
            level->deallocatePatchData(d_F_sc_idx);
            for (int k = 0; k < d_num_rand_vals; ++k) level->deallocatePatchData(d_F_sc_idxs[k]);
        }
    }
    return;
} // setDataOnPatchHierarchy

void
AdvDiffStochasticForcing::setDataOnPatch(const int data_idx,
                                         Pointer<Variable<NDIM> > /*var*/,
                                         Pointer<Patch<NDIM> > patch,
                                         const double /*data_time*/,
                                         const bool initial_time,
                                         Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    Pointer<CellData<NDIM, double> > divF_cc_data = patch->getPatchData(data_idx);
    divF_cc_data->fillAll(0.0);
    if (initial_time) return;
    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    double dV = 1.0;
    for (unsigned int d = 0; d < NDIM; ++d) dV *= dx[d];
    const double kappa = d_adv_diff_solver->getDiffusionCoefficient(d_C_var);
    const double dt = d_adv_diff_solver->getCurrentTimeStepSize();
    const double scale = d_std * sqrt(2.0 * kappa / (dt * dV));
    double C;
    d_f_parser.DefineVar("c", &C);
    d_f_parser.DefineVar("C", &C);
    Pointer<CellData<NDIM, double> > C_current_cc_data = patch->getPatchData(d_C_current_cc_idx);
    Pointer<CellData<NDIM, double> > C_half_cc_data = patch->getPatchData(d_C_half_cc_idx);
    Pointer<CellData<NDIM, double> > C_new_cc_data = patch->getPatchData(d_C_new_cc_idx);
    Pointer<SideData<NDIM, double> > F_sc_data = patch->getPatchData(d_F_sc_idx);
    Pointer<CellDataFactory<NDIM, double> > C_factory = d_C_var->getPatchDataFactory();
    const int C_depth = C_factory->getDefaultDepth();
    SideData<NDIM, double> f_scale_sc_data(patch_box, C_depth, IntVector<NDIM>(0));
    const TimeSteppingType convective_time_stepping_type = d_adv_diff_solver->getConvectiveTimeSteppingType(d_C_var);
    const int cycle_num = d_adv_diff_solver->getCurrentCycleNumber();
    for (int d = 0; d < C_depth; ++d)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (BoxIterator<NDIM> i(SideGeometry<NDIM>::toSideBox(patch_box, axis)); i; i++)
            {
                const Index<NDIM>& ic = i();
                Index<NDIM> ic_lower(ic);
                ic_lower(axis) -= 1;
                SideIndex<NDIM> is(ic, axis, SideIndex<NDIM>::Lower);
                double f;
                switch (convective_time_stepping_type)
                {
                case FORWARD_EULER:
                {
                    C = 0.5 * ((*C_current_cc_data)(ic, d) + (*C_current_cc_data)(ic_lower, d));
                    f = d_f_parser.Eval();
                    break;
                }
                case MIDPOINT_RULE:
                {
                    C = 0.5 * ((*C_half_cc_data)(ic, d) + (*C_half_cc_data)(ic_lower, d));
                    f = d_f_parser.Eval();
                    break;
                }
                case TRAPEZOIDAL_RULE:
                {
                    if (cycle_num == 0)
                    {
                        C = 0.5 * ((*C_current_cc_data)(ic, d) + (*C_current_cc_data)(ic_lower, d));
                        f = d_f_parser.Eval();
                    }
                    else
                    {
                        C = 0.5 * ((*C_current_cc_data)(ic, d) + (*C_current_cc_data)(ic_lower, d));
                        f = 0.5 * d_f_parser.Eval();
                        C = 0.5 * ((*C_new_cc_data)(ic, d) + (*C_new_cc_data)(ic_lower, d));
                        f += 0.5 * d_f_parser.Eval();
                    }
                    break;
                }
                default:
                    f = std::numeric_limits<double>::quiet_NaN();
                    TBOX_ERROR(d_object_name << "::setDataOnPatch():\n"
                                             << "  unsupported default convective time stepping type: "
                                             << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                                             << " \n"
                                             << "  valid choices are: FORWARD_EULER, MIDPOINT_RULE, "
                                                "TRAPEZOIDAL_RULE\n");
                }
                f_scale_sc_data(is, d) = sqrt(f) * scale;
            }
            for (BoxIterator<NDIM> i(patch_box); i; i++)
            {
                const Index<NDIM>& ic = i();
                SideIndex<NDIM> is_lower(ic, axis, SideIndex<NDIM>::Lower);
                SideIndex<NDIM> is_upper(ic, axis, SideIndex<NDIM>::Upper);
                const double scale_lower = f_scale_sc_data(is_lower, d);
                const double scale_upper = f_scale_sc_data(is_upper, d);
                (*divF_cc_data)(ic, d) +=
                    (scale_upper * (*F_sc_data)(is_upper, d) - scale_lower * (*F_sc_data)(is_lower, d)) / dx[axis];
            }
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
