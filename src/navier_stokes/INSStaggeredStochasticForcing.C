// Filename: INSStaggeredStochasticForcing.C
// Created on 02 Feb 2011 by Boyce Griffith
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

#include "INSStaggeredStochasticForcing.h"

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
#if (NDIM == 3)
#include <ibtk/EdgeDataSynchronization.h>
#endif
#include <ibtk/HierarchyGhostCellInterpolation.h>
#if (NDIM == 2)
#include <ibtk/NodeDataSynchronization.h>
#endif
#include <ibtk/PhysicalBoundaryUtilities.h>

// C++ STDLIB INCLUDES
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_STOCHASTIC_STRESS_DIV_FC FC_FUNC_(navier_stokes_stochastic_stress_div2d,NAVIER_STOKES_STOCHASTIC_STRESS_DIV2D)
#endif
#if (NDIM == 3)
#define NAVIER_STOKES_STOCHASTIC_STRESS_DIV_FC FC_FUNC_(navier_stokes_stochastic_stress_div3d,NAVIER_STOKES_STOCHASTIC_STRESS_DIV3D)
#endif

extern "C"
{
    void
    NAVIER_STOKES_STOCHASTIC_STRESS_DIV_FC(
#if (NDIM == 2)
        const double* ,
        const int& , const int& ,
        const int& , const int& ,
        const double& ,
        const int& , const int& ,
        const double* ,
        const int& , const int& ,
        const double* ,
        const int& , const int& ,
        double* , double*
#endif
#if (NDIM == 3)
        const double* ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double& ,
        const int& , const int& , const int& ,
        const double* ,
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
inline Box<NDIM>
compute_tangential_extension(
    const Box<NDIM>& box,
    const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
}// compute_tangential_extension

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
}// genrandn
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredStochasticForcing::INSStaggeredStochasticForcing(
    const INSStaggeredHierarchyIntegrator* const fluid_solver)
    : d_object_name("INSStaggeredStochasticForcing"),
      d_u_bc_coef(),
      d_dirichlet_scaling(NDIM == 2 ? 2.0 : 5.0/3.0),
      d_neumann_scaling(0.0),
      d_fluid_solver(fluid_solver),
      d_context(NULL),
      d_W_cc_var(NULL),
      d_W_cc_idx(-1),
#if (NDIM == 2)
      d_W_nc_var(NULL),
      d_W_nc_idx(-1),
#endif
#if (NDIM == 3)
      d_W_ec_var(NULL),
      d_W_ec_idx(-1),
#endif
      d_stress_tensor_type(UNCORRELATED),
      d_rho(std::numeric_limits<double>::quiet_NaN()),
      d_mu(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_std(std::numeric_limits<double>::quiet_NaN()),
      d_regen_rand_cycle()
{
    // Setup variables and variable context objects.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext("INSStaggeredStochasticForcing::CONTEXT");

    d_W_cc_var = new CellVariable<NDIM,double>("INSStaggeredStochasticForcing::W_cc", NDIM);
    static const IntVector<NDIM> ghosts_cc = 1;
    d_W_cc_idx = var_db->registerVariableAndContext(d_W_cc_var, d_context, ghosts_cc);
#if (NDIM == 2)
    d_W_nc_var = new NodeVariable<NDIM,double>("INSStaggeredStochasticForcing::W_nc", 2);
    static const IntVector<NDIM> ghosts_nc = 0;
    d_W_nc_idx = var_db->registerVariableAndContext(d_W_nc_var, d_context, ghosts_nc);
#endif
#if (NDIM == 3)
    d_W_ec_var = new EdgeVariable<NDIM,double>("INSStaggeredStochasticForcing::W_ec", 2);
    static const IntVector<NDIM> ghosts_ec = 0;
    d_W_ec_idx = var_db->registerVariableAndContext(d_W_ec_var, d_context, ghosts_ec);
#endif

    // Set a default value for regen_rand_cycle.
    d_regen_rand_cycle.resizeArray(2);
    d_regen_rand_cycle[0] = 1;
    d_regen_rand_cycle[1] = 0;
    return;
}// INSStaggeredStochasticForcing

INSStaggeredStochasticForcing::~INSStaggeredStochasticForcing()
{
    // intentionally blank
    return;
}// ~INSStaggeredStochasticForcing

void
INSStaggeredStochasticForcing::setPhysicalBcCoefs(
    blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM> u_bc_coef)
{
    d_u_bc_coef = u_bc_coef;
    return;
}// setPhysicalBcCoefs

bool
INSStaggeredStochasticForcing::isTimeDependent() const
{
    return true;
}// isTimeDependent

void
INSStaggeredStochasticForcing::setStochasticStressTensorType(
    const StochasticStressTensorType stress_tensor_type)
{
    d_stress_tensor_type = stress_tensor_type;
    return;
}// setStochasticStressTensorType

void
INSStaggeredStochasticForcing::setFluidDensity(
    const double rho)
{
    d_rho = rho;
    return;
}// setFluidDensity

void
INSStaggeredStochasticForcing::setFluidViscosity(
    const double mu)
{
    d_mu = mu;
    return;
}// setFluidViscosity

void
INSStaggeredStochasticForcing::setDt(
    const double dt)
{
    d_dt = dt;
    return;
}// setDt

void
INSStaggeredStochasticForcing::setStd(
    const double std)
{
    d_std = std;
    return;
}// setStd

void
INSStaggeredStochasticForcing::setRegenRandCycle(
    const Array<int>& regen_rand_cycle)
{
    d_regen_rand_cycle = regen_rand_cycle;
    return;
}// setRegenRandCycle

void
INSStaggeredStochasticForcing::setDirichletBcScaling(
    const double dirichlet_scaling)
{
    d_dirichlet_scaling = dirichlet_scaling;
    return;
}// setDirichletBcScaling

void
INSStaggeredStochasticForcing::setNeumannBcScaling(
    const double neumann_scaling)
{
    d_neumann_scaling = neumann_scaling;
    return;
}// setNeumannBcScaling

void
INSStaggeredStochasticForcing::setDataOnPatchHierarchy(
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
        TBOX_ERROR("INSStaggeredStochasticForcing::setDataOnPatchHierarchy():\n"
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
        if (!level->checkAllocated(d_W_cc_idx)) level->allocatePatchData(d_W_cc_idx);
#if (NDIM == 2)
        if (!level->checkAllocated(d_W_nc_idx)) level->allocatePatchData(d_W_nc_idx);
#endif
#if (NDIM == 3)
        if (!level->checkAllocated(d_W_ec_idx)) level->allocatePatchData(d_W_ec_idx);
#endif
    }

    if (!initial_time && d_std != 0.0 && cycle_num >= 0 && d_regen_rand_cycle[cycle_num] == 1)
    {
        // Generate random components.
        for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_num);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM,double> > W_cc_data = patch->getPatchData(d_W_cc_idx);
                genrandn(W_cc_data->getArrayData(), W_cc_data->getBox());
#if (NDIM == 2)
                Pointer<NodeData<NDIM,double> > W_nc_data = patch->getPatchData(d_W_nc_idx);
                genrandn(W_nc_data->getArrayData(), NodeGeometry<NDIM>::toNodeBox(W_nc_data->getBox()));
#endif
#if (NDIM == 3)
                Pointer<EdgeData<NDIM,double> > W_ec_data = patch->getPatchData(d_W_ec_idx);
                for (int d = 0; d < NDIM; ++d)
                {
                    genrandn(W_ec_data->getArrayData(d), EdgeGeometry<NDIM>::toEdgeBox(W_ec_data->getBox(),d));
                }
#endif

                // Modify the stress tensor values (if necessary).
                const Box<NDIM>& patch_box = patch->getBox();
                if (d_stress_tensor_type == SYMMETRIC || d_stress_tensor_type == SYMMETRIC_TRACELESS)
                {
                    // Symmetrize the stress tensor.
                    //
                    // NOTE: By averaging random variates instead of just using
                    // one of the two, we do more work than necessary.
#if (NDIM == 2)
                    for (Box<NDIM>::Iterator it(NodeGeometry<NDIM>::toNodeBox(patch_box)); it; it++)
                    {
                        const NodeIndex<NDIM> i_n(it(),0);
                        double avg = 0.5*((*W_nc_data)(i_n,0) + (*W_nc_data)(i_n,1));
                        (*W_nc_data)(i_n,0) = sqrt(2.0)*avg;
                        (*W_nc_data)(i_n,1) = sqrt(2.0)*avg;
                    }
#endif
#if (NDIM == 3)
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator it(EdgeGeometry<NDIM>::toEdgeBox(patch_box,axis)); it; it++)
                        {
                            const EdgeIndex<NDIM> i_e(it(),axis,0);
                            double avg = 0.5*((*W_ec_data)(i_e,0) + (*W_ec_data)(i_e,1));
                            (*W_ec_data)(i_e,0) = sqrt(2.0)*avg;
                            (*W_ec_data)(i_e,1) = sqrt(2.0)*avg;
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
                                (*W_cc_data)(i,d) *= sqrt(2.0);
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
                                trace += (*W_cc_data)(i,d);
                            }
                            for (int d = 0; d < NDIM; ++d)
                            {
                                (*W_cc_data)(i,d) = sqrt(2.0)*((*W_cc_data)(i,d) - trace/double(NDIM));
                            }
                        }
                    }
                }
                else if (d_stress_tensor_type != UNCORRELATED)
                {
                    TBOX_ERROR(d_object_name << "::setDataOnPatchHierarchy():\n"
                               << "  unrecognized stress tensor type: " << enum_to_string<StochasticStressTensorType>(d_stress_tensor_type) << "." << std::endl);
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
                    touches_regular_bdry [axis].resizeArray(2);
                    touches_periodic_bdry[axis].resizeArray(2);
                    for (int upperlower = 0; upperlower < 2; ++upperlower)
                    {
                        touches_regular_bdry [axis][upperlower] = pgeom->getTouchesRegularBoundary( axis,upperlower);
                        touches_periodic_bdry[axis][upperlower] = pgeom->getTouchesPeriodicBoundary(axis,upperlower);
                    }
                }

                const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();

#if (NDIM == 2)
                const Box<NDIM> node_box = NodeGeometry<NDIM>::toNodeBox(patch_box);
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    const IntVector<NDIM> gcw_to_fill = 1;
                    const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                    const int location_index   = bdry_box.getLocationIndex();
                    const int bdry_normal_axis = location_index / 2;
                    const int bdry_tangent_axis = (bdry_normal_axis+1)%2;  // NOTE: NDIM == 2
                    const BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), location_index);
                    const Box<NDIM> bc_coef_box = compute_tangential_extension(PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), bdry_tangent_axis);
                    Pointer<ArrayData<NDIM,double> > acoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM,double> > bcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM,double> > gcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);

                    // Temporarily reset the patch geometry object associated
                    // with the patch so that boundary conditions are set at the
                    // correct spatial locations.
                    double shifted_patch_x_lower[NDIM], shifted_patch_x_upper[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        shifted_patch_x_lower[d] = patch_x_lower[d];
                        shifted_patch_x_upper[d] = patch_x_upper[d];
                    }
                    shifted_patch_x_lower[bdry_tangent_axis] -= 0.5*dx[bdry_tangent_axis];
                    shifted_patch_x_upper[bdry_tangent_axis] -= 0.5*dx[bdry_tangent_axis];
                    patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero, touches_regular_bdry, touches_periodic_bdry, dx, shifted_patch_x_lower, shifted_patch_x_upper));

                    // Set the boundary condition coefficients and use them to
                    // rescale the stochastic fluxes.
                    for (int d = 0; d < NDIM; ++d)
                    {
                        RobinBcCoefStrategy<NDIM>* bc_coef = d_u_bc_coef[d];
                        bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, *patch, trimmed_bdry_box, data_time);
                        for (Box<NDIM>::Iterator it(bc_coef_box*node_box); it; it++)
                        {
                            const Index<NDIM>& i = it();
                            const NodeIndex<NDIM> n_i(i,0);
                            const double& alpha = (*acoef_data)(i,0);
                            const double& beta  = (*bcoef_data)(i,0);
                            const bool dirichlet_bc = (alpha != 0.0 && beta == 0.0);
                            if (dirichlet_bc)
                            {
                                (*W_nc_data)(n_i,d) *= d_dirichlet_scaling;
                            }
                            else
                            {
                                (*W_nc_data)(n_i,d) *= d_neumann_scaling;
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
                    const int location_index   = bdry_box.getLocationIndex();
                    const int bdry_normal_axis = location_index / 2;
                    const BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), location_index);
                    for (int edge_axis = 0; edge_axis < NDIM; ++edge_axis)
                    {
                        if (edge_axis == bdry_normal_axis) continue;  // we only care about edges that are on the boundary

                        const Box<NDIM> bc_coef_box = compute_tangential_extension(PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), edge_axis);
                        Pointer<ArrayData<NDIM,double> > acoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                        Pointer<ArrayData<NDIM,double> > bcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                        Pointer<ArrayData<NDIM,double> > gcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);

                        // Temporarily reset the patch geometry object
                        // associated with the patch so that boundary conditions
                        // are set at the correct spatial locations.
                        double shifted_patch_x_lower[NDIM], shifted_patch_x_upper[NDIM];
                        for (int d = 0; d < NDIM; ++d)
                        {
                            shifted_patch_x_lower[d] = patch_x_lower[d];
                            shifted_patch_x_upper[d] = patch_x_upper[d];
                        }
                        shifted_patch_x_lower[edge_axis] -= 0.5*dx[edge_axis];
                        shifted_patch_x_upper[edge_axis] -= 0.5*dx[edge_axis];
                        patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero, touches_regular_bdry, touches_periodic_bdry, dx, shifted_patch_x_lower, shifted_patch_x_upper));

                        // Set the boundary condition coefficients and use them
                        // to rescale the stochastic fluxes.
                        for (int d = 0; d < NDIM; ++d)
                        {
                            if (d == edge_axis) continue;
                            const int data_depth = ((d == 1 && edge_axis == 2) || (d == 2)) ? 1 : 0;

                            RobinBcCoefStrategy<NDIM>* bc_coef = d_u_bc_coef[d];
                            bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, var, *patch, trimmed_bdry_box, data_time);
                            for (Box<NDIM>::Iterator it(bc_coef_box*edge_boxes[edge_axis]); it; it++)
                            {
                                const Index<NDIM>& i = it();
                                const EdgeIndex<NDIM> e_i(i, edge_axis, 0);
                                const double& alpha = (*acoef_data)(i,0);
                                const double& beta  = (*bcoef_data)(i,0);
                                const bool dirichlet_bc = (alpha != 0.0 && beta == 0.0);
                                if (dirichlet_bc)
                                {
                                    (*W_ec_data)(e_i,data_depth) *= d_dirichlet_scaling;
                                }
                                else
                                {
                                    (*W_ec_data)(e_i,data_depth) *= d_neumann_scaling;
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
        LocationIndexRobinBcCoefs<NDIM> bc_coef(d_object_name+"::bc_coef", Pointer<Database>(NULL));
        for (int d = 0; d < NDIM; ++d)
        {
            bc_coef.setBoundarySlope(2*d  ,0.0);
            bc_coef.setBoundarySlope(2*d+1,0.0);
        }
        std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM,&bc_coef);
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_fill_components(1);
        ghost_fill_components[0] = InterpolationTransactionComponent(d_W_cc_idx, "NONE", "NONE", false, bc_coefs);
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
}// computeStochasticForcingOnPatchHierarchy

void
INSStaggeredStochasticForcing::setDataOnPatchLevel(
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
INSStaggeredStochasticForcing::setDataOnPatch(
    const int data_idx,
    Pointer<Variable<NDIM> > /*var*/,
    Pointer<Patch<NDIM> > patch,
    const double /*data_time*/,
    const bool initial_time,
    Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    Pointer<SideData<NDIM,double> > divW_sc_data = patch->getPatchData(data_idx);
    const IntVector<NDIM> divW_sc_ghosts = divW_sc_data->getGhostCellWidth();
    divW_sc_data->fillAll(0.0);

    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    if (cycle_num >= int(d_regen_rand_cycle.size()))
    {
        TBOX_ERROR("INSStaggeredStochasticForcing::setDataOnPatch():\n"
                   << "  cycle number = " << cycle_num << "\n"
                   << "  but regen_rand_cycle.size() = " << d_regen_rand_cycle.size() << "\n");
    }

    if (initial_time || d_std == 0.0 || (cycle_num >= 0 && d_regen_rand_cycle[cycle_num] == -1)) return;

    Pointer<CellData<NDIM,double> > W_cc_data = patch->getPatchData(d_W_cc_idx);
    const IntVector<NDIM> W_cc_ghosts = W_cc_data->getGhostCellWidth();

#if (NDIM == 2)
    Pointer<NodeData<NDIM,double> > W_nc_data = patch->getPatchData(d_W_nc_idx);
    const IntVector<NDIM> W_nc_ghosts = W_nc_data->getGhostCellWidth();
#endif

#if (NDIM == 3)
    Pointer<EdgeData<NDIM,double> > W_ec_data = patch->getPatchData(d_W_ec_idx);
    const IntVector<NDIM> W_ec_ghosts = W_ec_data->getGhostCellWidth();
#endif

    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    double dV = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        dV *= dx[d];
    }
    // Note that this is really solving the *momentum*, rather than the velocity equation:
    const double scale = d_std*sqrt(2.0*d_mu/(d_dt*dV));
#if (NDIM == 2)
    double* const divW_sc0 = divW_sc_data->getPointer(0);
    double* const divW_sc1 = divW_sc_data->getPointer(1);
    const double* const W_cc = W_cc_data->getPointer();
    const double* const W_nc = W_nc_data->getPointer();
    NAVIER_STOKES_STOCHASTIC_STRESS_DIV_FC(
        dx,
        patch_box.lower(0), patch_box.upper(0),
        patch_box.lower(1), patch_box.upper(1),
        scale,
        W_cc_ghosts(0), W_cc_ghosts(1),
        W_cc,
        W_nc_ghosts(0), W_nc_ghosts(1),
        W_nc,
        divW_sc_ghosts(0), divW_sc_ghosts(1),
        divW_sc0, divW_sc1
                                           );
#endif
#if (NDIM == 3)
    double* const divW_sc0 = divW_sc_data->getPointer(0);
    double* const divW_sc1 = divW_sc_data->getPointer(1);
    double* const divW_sc2 = divW_sc_data->getPointer(2);
    const double* const W_cc = W_cc_data->getPointer();
    const double* const W_ec0 = W_ec_data->getPointer(0);
    const double* const W_ec1 = W_ec_data->getPointer(1);
    const double* const W_ec2 = W_ec_data->getPointer(2);
    NAVIER_STOKES_STOCHASTIC_STRESS_DIV_FC(
        dx,
        patch_box.lower(0), patch_box.upper(0),
        patch_box.lower(1), patch_box.upper(1),
        patch_box.lower(2), patch_box.upper(2),
        scale,
        W_cc_ghosts(0), W_cc_ghosts(1), W_cc_ghosts(2),
        W_cc,
        W_ec_ghosts(0), W_ec_ghosts(1), W_ec_ghosts(2),
        W_ec0, W_ec1, W_ec2,
        divW_sc_ghosts(0), divW_sc_ghosts(1), divW_sc_ghosts(2),
        divW_sc0, divW_sc1, divW_sc2
                                           );
#endif
    return;
}// setDataOnPatch

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::INSStaggeredStochasticForcing>;

//////////////////////////////////////////////////////////////////////////////
