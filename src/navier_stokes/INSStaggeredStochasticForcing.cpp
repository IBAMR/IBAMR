// Filename: INSStaggeredStochasticForcing.cpp
// Created on 02 Feb 2011 by Boyce Griffith
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
#include "CellIndex.h"
#include "CellVariable.h"
#include "EdgeData.h"     // IWYU pragma: keep
#include "EdgeGeometry.h" // IWYU pragma: keep
#include "EdgeIndex.h"    // IWYU pragma: keep
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "NodeData.h"     // IWYU pragma: keep
#include "NodeGeometry.h" // IWYU pragma: keep
#include "NodeIndex.h"    // IWYU pragma: keep
#include "NodeVariable.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/INSStaggeredStochasticForcing.h"
#include "ibamr/RNG.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#if (NDIM == 2)
#include "ibtk/NodeDataSynchronization.h"
#endif
#if (NDIM == 3)
#include "ibtk/EdgeDataSynchronization.h"
#endif

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_STOCHASTIC_STRESS_DIV_FC                                                                         \
    IBAMR_FC_FUNC_(navier_stokes_stochastic_stress_div2d, NAVIER_STOKES_STOCHASTIC_STRESS_DIV2D)
#endif
#if (NDIM == 3)
#define NAVIER_STOKES_STOCHASTIC_STRESS_DIV_FC                                                                         \
    IBAMR_FC_FUNC_(navier_stokes_stochastic_stress_div3d, NAVIER_STOKES_STOCHASTIC_STRESS_DIV3D)
#endif

extern "C" {
void NAVIER_STOKES_STOCHASTIC_STRESS_DIV_FC(
#if (NDIM == 2)
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const double&,
    const int&,
    const int&,
    const double*,
    const int&,
    const int&,
    const double*,
    const int&,
    const int&,
    double*,
    double*
#endif
#if (NDIM == 3)
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double&,
    const int&,
    const int&,
    const int&,
    const double*,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*,
    const int&,
    const int&,
    const int&,
    double*,
    double*,
    double*
#endif
    );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline Box<NDIM>
compute_tangential_extension(const Box<NDIM>& box, const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
} // compute_tangential_extension

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

INSStaggeredStochasticForcing::INSStaggeredStochasticForcing(const std::string& object_name,
                                                             Pointer<Database> input_db,
                                                             const INSStaggeredHierarchyIntegrator* const fluid_solver)
    : d_object_name(object_name),
      d_fluid_solver(fluid_solver),
      d_stress_tensor_type(UNCORRELATED),
      d_std(std::numeric_limits<double>::quiet_NaN()),
      d_num_rand_vals(0),
      d_weights(),
      d_velocity_bc_scaling(NDIM == 2 ? 2.0 : 5.0 / 3.0),
      d_traction_bc_scaling(0.0),
      d_context(NULL),
      d_W_cc_var(NULL),
      d_W_cc_idx(-1),
      d_W_cc_idxs(),
#if (NDIM == 2)
      d_W_nc_var(NULL),
      d_W_nc_idx(-1),
      d_W_nc_idxs()
#endif
#if (NDIM == 3)
          d_W_ec_var(NULL),
      d_W_ec_idx(-1),
      d_W_ec_idxs()
#endif
{
    if (input_db)
    {
        if (input_db->keyExists("stress_tensor_type"))
            d_stress_tensor_type =
                string_to_enum<StochasticStressTensorType>(input_db->getString("stress_tensor_type"));
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
        if (input_db->keyExists("velocity_bc_scaling"))
            d_velocity_bc_scaling = input_db->getDouble("velocity_bc_scaling");
        if (input_db->keyExists("traction_bc_scaling"))
            d_traction_bc_scaling = input_db->getDouble("traction_bc_scaling");
    }

    // Setup variables and variable context objects.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    d_W_cc_var = new CellVariable<NDIM, double>(d_object_name + "::W_cc", NDIM);
    static const IntVector<NDIM> ghosts_cc = 1;
    d_W_cc_idx = var_db->registerVariableAndContext(d_W_cc_var, d_context, ghosts_cc);
    for (int k = 0; k < d_num_rand_vals; ++k)
        d_W_cc_idxs.push_back(var_db->registerClonedPatchDataIndex(d_W_cc_var, d_W_cc_idx));
#if (NDIM == 2)
    d_W_nc_var = new NodeVariable<NDIM, double>(d_object_name + "::W_nc", 2);
    static const IntVector<NDIM> ghosts_nc = 0;
    d_W_nc_idx = var_db->registerVariableAndContext(d_W_nc_var, d_context, ghosts_nc);
    for (int k = 0; k < d_num_rand_vals; ++k)
        d_W_nc_idxs.push_back(var_db->registerClonedPatchDataIndex(d_W_nc_var, d_W_nc_idx));
#endif
#if (NDIM == 3)
    d_W_ec_var = new EdgeVariable<NDIM, double>(d_object_name + "::W_ec", 2);
    static const IntVector<NDIM> ghosts_ec = 0;
    d_W_ec_idx = var_db->registerVariableAndContext(d_W_ec_var, d_context, ghosts_ec);
    for (int k = 0; k < d_num_rand_vals; ++k)
        d_W_ec_idxs.push_back(var_db->registerClonedPatchDataIndex(d_W_ec_var, d_W_ec_idx));
#endif
    return;
} // INSStaggeredStochasticForcing

INSStaggeredStochasticForcing::~INSStaggeredStochasticForcing()
{
    // intentionally blank
    return;
} // ~INSStaggeredStochasticForcing

bool
INSStaggeredStochasticForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
INSStaggeredStochasticForcing::setDataOnPatchHierarchy(const int data_idx,
                                                       Pointer<Variable<NDIM> > var,
                                                       Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                       const double data_time,
                                                       const bool initial_time,
                                                       const int coarsest_ln_in,
                                                       const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    if (!initial_time)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(cycle_num >= 0);
#endif
        // Allocate data to store components of the stochastic stress components.
        for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
            level->allocatePatchData(d_W_cc_idx);
            for (int k = 0; k < d_num_rand_vals; ++k) level->allocatePatchData(d_W_cc_idxs[k]);
#if (NDIM == 2)
            level->allocatePatchData(d_W_nc_idx);
            for (int k = 0; k < d_num_rand_vals; ++k) level->allocatePatchData(d_W_nc_idxs[k]);
#endif
#if (NDIM == 3)
            level->allocatePatchData(d_W_ec_idx);
            for (int k = 0; k < d_num_rand_vals; ++k) level->allocatePatchData(d_W_ec_idxs[k]);
#endif
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
                        Pointer<CellData<NDIM, double> > W_cc_data = patch->getPatchData(d_W_cc_idxs[k]);
                        genrandn(W_cc_data->getArrayData(), W_cc_data->getBox());
#if (NDIM == 2)
                        Pointer<NodeData<NDIM, double> > W_nc_data = patch->getPatchData(d_W_nc_idxs[k]);
                        genrandn(W_nc_data->getArrayData(), NodeGeometry<NDIM>::toNodeBox(W_nc_data->getBox()));
#endif
#if (NDIM == 3)
                        Pointer<EdgeData<NDIM, double> > W_ec_data = patch->getPatchData(d_W_ec_idxs[k]);
                        for (int d = 0; d < NDIM; ++d)
                        {
                            genrandn(W_ec_data->getArrayData(d), EdgeGeometry<NDIM>::toEdgeBox(W_ec_data->getBox(), d));
                        }
#endif
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
        HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_cc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_W_cc_var,
                                                       hierarchy,
                                                       /*get_unique*/ true);
        hier_cc_data_ops->setToScalar(d_W_cc_idx, 0.0);
        for (int k = 0; k < d_num_rand_vals; ++k)
            hier_cc_data_ops->axpy(d_W_cc_idx, weights[k], d_W_cc_idxs[k], d_W_cc_idx);
#if (NDIM == 2)
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_nc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_W_nc_var,
                                                       hierarchy,
                                                       /*get_unique*/ true);
        hier_nc_data_ops->setToScalar(d_W_nc_idx, 0.0);
        for (int k = 0; k < d_num_rand_vals; ++k)
            hier_nc_data_ops->axpy(d_W_nc_idx, weights[k], d_W_nc_idxs[k], d_W_nc_idx);
#endif
#if (NDIM == 3)
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_ec_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_W_ec_var,
                                                       hierarchy,
                                                       /*get_unique*/ true);
        hier_ec_data_ops->setToScalar(d_W_ec_idx, 0.0);
        for (int k = 0; k < d_num_rand_vals; ++k)
            hier_ec_data_ops->axpy(d_W_ec_idx, weights[k], d_W_ec_idxs[k], d_W_ec_idx);
#endif

        // Modify the stress tensor values (if necessary).
        const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs =
            d_fluid_solver->getIntermediateVelocityBoundaryConditions();
        for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double> > W_cc_data = patch->getPatchData(d_W_cc_idx);
#if (NDIM == 2)
                Pointer<NodeData<NDIM, double> > W_nc_data = patch->getPatchData(d_W_nc_idx);
#endif
#if (NDIM == 3)
                Pointer<EdgeData<NDIM, double> > W_ec_data = patch->getPatchData(d_W_ec_idx);
#endif
                if (d_stress_tensor_type == SYMMETRIC || d_stress_tensor_type == SYMMETRIC_TRACELESS)
                {
// Symmetrize the stress tensor.
//
// NOTE: By averaging random variates instead of just using
// one of the two, we do more work than necessary.
#if (NDIM == 2)
                    for (Box<NDIM>::Iterator it(NodeGeometry<NDIM>::toNodeBox(patch_box)); it; it++)
                    {
                        const NodeIndex<NDIM> i_n(it(), 0);
                        double avg = 0.5 * ((*W_nc_data)(i_n, 0) + (*W_nc_data)(i_n, 1));
                        (*W_nc_data)(i_n, 0) = sqrt(2.0) * avg;
                        (*W_nc_data)(i_n, 1) = sqrt(2.0) * avg;
                    }
#endif
#if (NDIM == 3)
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator it(EdgeGeometry<NDIM>::toEdgeBox(patch_box, axis)); it; it++)
                        {
                            const EdgeIndex<NDIM> i_e(it(), axis, 0);
                            double avg = 0.5 * ((*W_ec_data)(i_e, 0) + (*W_ec_data)(i_e, 1));
                            (*W_ec_data)(i_e, 0) = sqrt(2.0) * avg;
                            (*W_ec_data)(i_e, 1) = sqrt(2.0) * avg;
                        }
                    }
#endif
                    if (d_stress_tensor_type == SYMMETRIC)
                    {
                        // Multiply the diagonal by sqrt(2) to make the variance
                        // 2.
                        for (Box<NDIM>::Iterator it(patch_box); it; it++)
                        {
                            const Index<NDIM>& i = it();
                            for (int d = 0; d < NDIM; ++d)
                            {
                                (*W_cc_data)(i, d) *= sqrt(2.0);
                            }
                        }
                    }
                    else if (d_stress_tensor_type == SYMMETRIC_TRACELESS)
                    {
                        // Subtract the trace from the diagonal and multiply the
                        // diagonal by sqrt(2) to make the variance 2.
                        for (Box<NDIM>::Iterator it(patch_box); it; it++)
                        {
                            const Index<NDIM>& i = it();
                            double trace = 0.0;
                            for (int d = 0; d < NDIM; ++d)
                            {
                                trace += (*W_cc_data)(i, d);
                            }
                            for (int d = 0; d < NDIM; ++d)
                            {
                                (*W_cc_data)(i, d) = sqrt(2.0) * ((*W_cc_data)(i, d) - trace / double(NDIM));
                            }
                        }
                    }
                }
                else if (d_stress_tensor_type != UNCORRELATED)
                {
                    TBOX_ERROR(d_object_name << "::setDataOnPatchHierarchy():\n"
                                             << "  unrecognized stress tensor type: "
                                             << enum_to_string<StochasticStressTensorType>(d_stress_tensor_type)
                                             << "."
                                             << std::endl);
                }

                const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                if (!pgeom->getTouchesRegularBoundary()) continue;
                const double* const dx = pgeom->getDx();
                const double* const patch_x_lower = pgeom->getXLower();
                const double* const patch_x_upper = pgeom->getXUpper();
                const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
                Array<Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    touches_regular_bdry[axis].resizeArray(2);
                    touches_periodic_bdry[axis].resizeArray(2);
                    for (int upperlower = 0; upperlower < 2; ++upperlower)
                    {
                        touches_regular_bdry[axis][upperlower] = pgeom->getTouchesRegularBoundary(axis, upperlower);
                        touches_periodic_bdry[axis][upperlower] = pgeom->getTouchesPeriodicBoundary(axis, upperlower);
                    }
                }

                const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
                    PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();

#if (NDIM == 2)
                const Box<NDIM> node_box = NodeGeometry<NDIM>::toNodeBox(patch_box);
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    const IntVector<NDIM> gcw_to_fill = 1;
                    const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                    const int location_index = bdry_box.getLocationIndex();
                    const int bdry_normal_axis = location_index / 2;
                    const int bdry_tangent_axis = (bdry_normal_axis + 1) % 2; // NOTE: NDIM == 2
                    const BoundaryBox<NDIM> trimmed_bdry_box(
                        bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);
                    const Box<NDIM> bc_coef_box = compute_tangential_extension(
                        PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), bdry_tangent_axis);
                    Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                    // Temporarily reset the patch geometry object associated
                    // with the patch so that boundary conditions are set at the
                    // correct spatial locations.
                    double shifted_patch_x_lower[NDIM], shifted_patch_x_upper[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        shifted_patch_x_lower[d] = patch_x_lower[d];
                        shifted_patch_x_upper[d] = patch_x_upper[d];
                    }
                    shifted_patch_x_lower[bdry_tangent_axis] -= 0.5 * dx[bdry_tangent_axis];
                    shifted_patch_x_upper[bdry_tangent_axis] -= 0.5 * dx[bdry_tangent_axis];
                    patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                             touches_regular_bdry,
                                                                             touches_periodic_bdry,
                                                                             dx,
                                                                             shifted_patch_x_lower,
                                                                             shifted_patch_x_upper));

                    // Set the boundary condition coefficients and use them to
                    // rescale the stochastic fluxes.
                    for (int d = 0; d < NDIM; ++d)
                    {
                        RobinBcCoefStrategy<NDIM>* bc_coef = u_bc_coefs[d];
                        bc_coef->setBcCoefs(
                            acoef_data, bcoef_data, gcoef_data, var, *patch, trimmed_bdry_box, data_time);
                        for (Box<NDIM>::Iterator it(bc_coef_box * node_box); it; it++)
                        {
                            const Index<NDIM>& i = it();
                            const NodeIndex<NDIM> n_i(i, 0);
                            const double& alpha = (*acoef_data)(i, 0);
                            const double& beta = (*bcoef_data)(i, 0);
                            const bool velocity_bc = (alpha != 0.0 && beta == 0.0);
                            if (velocity_bc)
                            {
                                (*W_nc_data)(n_i, d) *= d_velocity_bc_scaling;
                            }
                            else
                            {
                                (*W_nc_data)(n_i, d) *= d_traction_bc_scaling;
                            }
                        }
                    }

                    // Restore the original patch geometry object.
                    patch->setPatchGeometry(pgeom);
                }
#endif
#if (NDIM == 3)
                Box<NDIM> edge_boxes[NDIM];
                for (int d = 0; d < NDIM; ++d)
                {
                    edge_boxes[d] = EdgeGeometry<NDIM>::toEdgeBox(patch_box, d);
                }
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    const IntVector<NDIM> gcw_to_fill = 1;
                    const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                    const int location_index = bdry_box.getLocationIndex();
                    const int bdry_normal_axis = location_index / 2;
                    const BoundaryBox<NDIM> trimmed_bdry_box(
                        bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);
                    for (int edge_axis = 0; edge_axis < NDIM; ++edge_axis)
                    {
                        if (edge_axis == bdry_normal_axis)
                            continue; // we only care about edges that are on the boundary

                        const Box<NDIM> bc_coef_box = compute_tangential_extension(
                            PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), edge_axis);
                        Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                        Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                        Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                        // Temporarily reset the patch geometry object
                        // associated with the patch so that boundary conditions
                        // are set at the correct spatial locations.
                        double shifted_patch_x_lower[NDIM], shifted_patch_x_upper[NDIM];
                        for (int d = 0; d < NDIM; ++d)
                        {
                            shifted_patch_x_lower[d] = patch_x_lower[d];
                            shifted_patch_x_upper[d] = patch_x_upper[d];
                        }
                        shifted_patch_x_lower[edge_axis] -= 0.5 * dx[edge_axis];
                        shifted_patch_x_upper[edge_axis] -= 0.5 * dx[edge_axis];
                        patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                                 touches_regular_bdry,
                                                                                 touches_periodic_bdry,
                                                                                 dx,
                                                                                 shifted_patch_x_lower,
                                                                                 shifted_patch_x_upper));

                        // Set the boundary condition coefficients and use them
                        // to rescale the stochastic fluxes.
                        for (int d = 0; d < NDIM; ++d)
                        {
                            if (d == edge_axis) continue;
                            const int data_depth = ((d == 1 && edge_axis == 2) || (d == 2)) ? 1 : 0;

                            RobinBcCoefStrategy<NDIM>* bc_coef = u_bc_coefs[d];
                            bc_coef->setBcCoefs(
                                acoef_data, bcoef_data, gcoef_data, var, *patch, trimmed_bdry_box, data_time);
                            for (Box<NDIM>::Iterator it(bc_coef_box * edge_boxes[edge_axis]); it; it++)
                            {
                                const Index<NDIM>& i = it();
                                const EdgeIndex<NDIM> e_i(i, edge_axis, 0);
                                const double& alpha = (*acoef_data)(i, 0);
                                const double& beta = (*bcoef_data)(i, 0);
                                const bool velocity_bc = (alpha != 0.0 && beta == 0.0);
                                if (velocity_bc)
                                {
                                    (*W_ec_data)(e_i, data_depth) *= d_velocity_bc_scaling;
                                }
                                else
                                {
                                    (*W_ec_data)(e_i, data_depth) *= d_traction_bc_scaling;
                                }
                            }
                        }

                        // Restore the original patch geometry object.
                        patch->setPatchGeometry(pgeom);
                    }
                }
#endif
            }
        }

#if (NDIM == 2)
        // Synchronize node-centered values.
        typedef NodeDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
        SynchronizationTransactionComponent synch_component(d_W_nc_idx);
        NodeDataSynchronization synch_data_op;
        synch_data_op.initializeOperatorState(synch_component, hierarchy);
        synch_data_op.synchronizeData(data_time);
#endif
#if (NDIM == 3)
        // Synchronize edge-centered values.
        typedef EdgeDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
        SynchronizationTransactionComponent synch_component(d_W_ec_idx);
        EdgeDataSynchronization synch_data_op;
        synch_data_op.initializeOperatorState(synch_component, hierarchy);
        synch_data_op.synchronizeData(data_time);
#endif

        // Communicate ghost-cell data.
        LocationIndexRobinBcCoefs<NDIM> bc_coef(d_object_name + "::bc_coef", Pointer<Database>(NULL));
        for (int d = 0; d < NDIM; ++d)
        {
            bc_coef.setBoundarySlope(2 * d, 0.0);
            bc_coef.setBoundarySlope(2 * d + 1, 0.0);
        }
        std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM, &bc_coef);
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_fill_components(1);
        ghost_fill_components[0] =
            InterpolationTransactionComponent(d_W_cc_idx, "NONE", false, "NONE", "NONE", false, bc_coefs);
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_fill_components, hierarchy);
        ghost_fill_op.fillData(data_time);
    }

    // Compute div W on each patch level.
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
            level->allocatePatchData(d_W_cc_idx);
            for (int k = 0; k < d_num_rand_vals; ++k) level->allocatePatchData(d_W_cc_idxs[k]);
#if (NDIM == 2)
            level->allocatePatchData(d_W_nc_idx);
            for (int k = 0; k < d_num_rand_vals; ++k) level->allocatePatchData(d_W_nc_idxs[k]);
#endif
#if (NDIM == 3)
            level->allocatePatchData(d_W_ec_idx);
            for (int k = 0; k < d_num_rand_vals; ++k) level->allocatePatchData(d_W_ec_idxs[k]);
#endif
        }
    }
    return;
} // computeStochasticForcingOnPatchHierarchy

void
INSStaggeredStochasticForcing::setDataOnPatch(const int data_idx,
                                              Pointer<Variable<NDIM> > /*var*/,
                                              Pointer<Patch<NDIM> > patch,
                                              const double /*data_time*/,
                                              const bool initial_time,
                                              Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    Pointer<SideData<NDIM, double> > divW_sc_data = patch->getPatchData(data_idx);
    const IntVector<NDIM> divW_sc_ghosts = divW_sc_data->getGhostCellWidth();
    divW_sc_data->fillAll(0.0);
    if (initial_time) return;
    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    double dV = 1.0;
    for (unsigned int d = 0; d < NDIM; ++d) dV *= dx[d];
    const double mu = d_fluid_solver->getStokesSpecifications()->getMu();
    const double dt = d_fluid_solver->getCurrentTimeStepSize();
    // NOTE: We are solving the momentum equation, not the velocity equation.
    const double scale = d_std * sqrt(2.0 * mu / (dt * dV));
    Pointer<CellData<NDIM, double> > W_cc_data = patch->getPatchData(d_W_cc_idx);
    const IntVector<NDIM> W_cc_ghosts = W_cc_data->getGhostCellWidth();
#if (NDIM == 2)
    Pointer<NodeData<NDIM, double> > W_nc_data = patch->getPatchData(d_W_nc_idx);
    const IntVector<NDIM> W_nc_ghosts = W_nc_data->getGhostCellWidth();
    double* const divW_sc0 = divW_sc_data->getPointer(0);
    double* const divW_sc1 = divW_sc_data->getPointer(1);
    const double* const W_cc = W_cc_data->getPointer();
    const double* const W_nc = W_nc_data->getPointer();
    NAVIER_STOKES_STOCHASTIC_STRESS_DIV_FC(dx,
                                           patch_box.lower(0),
                                           patch_box.upper(0),
                                           patch_box.lower(1),
                                           patch_box.upper(1),
                                           scale,
                                           W_cc_ghosts(0),
                                           W_cc_ghosts(1),
                                           W_cc,
                                           W_nc_ghosts(0),
                                           W_nc_ghosts(1),
                                           W_nc,
                                           divW_sc_ghosts(0),
                                           divW_sc_ghosts(1),
                                           divW_sc0,
                                           divW_sc1);
#endif
#if (NDIM == 3)
    Pointer<EdgeData<NDIM, double> > W_ec_data = patch->getPatchData(d_W_ec_idx);
    const IntVector<NDIM> W_ec_ghosts = W_ec_data->getGhostCellWidth();
    double* const divW_sc0 = divW_sc_data->getPointer(0);
    double* const divW_sc1 = divW_sc_data->getPointer(1);
    double* const divW_sc2 = divW_sc_data->getPointer(2);
    const double* const W_cc = W_cc_data->getPointer();
    const double* const W_ec0 = W_ec_data->getPointer(0);
    const double* const W_ec1 = W_ec_data->getPointer(1);
    const double* const W_ec2 = W_ec_data->getPointer(2);
    NAVIER_STOKES_STOCHASTIC_STRESS_DIV_FC(dx,
                                           patch_box.lower(0),
                                           patch_box.upper(0),
                                           patch_box.lower(1),
                                           patch_box.upper(1),
                                           patch_box.lower(2),
                                           patch_box.upper(2),
                                           scale,
                                           W_cc_ghosts(0),
                                           W_cc_ghosts(1),
                                           W_cc_ghosts(2),
                                           W_cc,
                                           W_ec_ghosts(0),
                                           W_ec_ghosts(1),
                                           W_ec_ghosts(2),
                                           W_ec0,
                                           W_ec1,
                                           W_ec2,
                                           divW_sc_ghosts(0),
                                           divW_sc_ghosts(1),
                                           divW_sc_ghosts(2),
                                           divW_sc0,
                                           divW_sc1,
                                           divW_sc2);
#endif
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
