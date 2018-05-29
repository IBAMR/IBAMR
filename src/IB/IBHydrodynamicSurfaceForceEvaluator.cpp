// Filename: IBHydrodynamicSurfaceForceEvaluator.cpp
// Created on 11 May 2018 by Nishant Nangia
//
// Copyright (c) 2002-2018, Nishant Nangia and Amneet Bhalla
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

#include "ibamr/IBHydrodynamicSurfaceForceEvaluator.h"
#include "ArrayDataBasicOps.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CoarseFineBoundary.h"
#include "HierarchyDataOpsManager.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "SideData.h"
#include "SideIndex.h"
#include "boost/array.hpp"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/INSStaggeredPressureBcCoef.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredPressureBcCoef.h"
#include "ibamr/namespaces.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/IndexUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
static const int GLEVELSETG = 1;
static const int GVELOCITYG = 1;
static const int GPRESSUREG = 1;
static const int GVISCOSITYG = 1;

inline int
sign(const double X)
{
    return ((X > 0) ? 1 : ((X < 0) ? -1 : 0));
} // sign
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBHydrodynamicSurfaceForceEvaluator::IBHydrodynamicSurfaceForceEvaluator(
    const std::string& object_name,
    Pointer<CellVariable<NDIM, double> > ls_solid_var,
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
    Pointer<INSHierarchyIntegrator> fluid_solver,
    Pointer<Database> db)
    : d_object_name(object_name),
      d_ls_solid_var(ls_solid_var),
      d_adv_diff_solver(adv_diff_solver),
      d_fluid_solver(fluid_solver),
      d_ls_solid_idx(-1),
      d_u_idx(-1),
      d_p_idx(-1),
      d_mu_idx(-1),
      d_mu(std::numeric_limits<double>::quiet_NaN())

{
    // Set some default values
    d_surface_contour_value = 0.0;

    // Get from database
    if (!db.isNull()) getFromInput(db);

#if !defined(NDEBUG)
    TBOX_ASSERT(d_adv_diff_solver);
    TBOX_ASSERT(d_fluid_solver);
#endif

    // Registered required patch data
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

#if !defined(NDEBUG)
    TBOX_ASSERT(d_ls_solid_var);
#endif
    Pointer<VariableContext> ls_ctx = var_db->getContext(d_object_name + "::ls_ctx");
    d_ls_solid_idx = var_db->registerVariableAndContext(d_ls_solid_var, ls_ctx, IntVector<NDIM>(GLEVELSETG));

    Pointer<SideVariable<NDIM, double> > u_var = d_fluid_solver->getVelocityVariable();
#if !defined(NDEBUG)
    TBOX_ASSERT(u_var);
#endif
    Pointer<VariableContext> u_ctx = var_db->getContext(d_object_name + "::u_ctx");
    d_u_idx = var_db->registerVariableAndContext(u_var, u_ctx, IntVector<NDIM>(GVELOCITYG));

    Pointer<CellVariable<NDIM, double> > p_var = d_fluid_solver->getPressureVariable();
#if !defined(NDEBUG)
    TBOX_ASSERT(p_var);
#endif
    Pointer<VariableContext> p_ctx = var_db->getContext(d_object_name + "::p_ctx");
    d_p_idx = var_db->registerVariableAndContext(p_var, p_ctx, IntVector<NDIM>(GPRESSUREG));

    INSStaggeredHierarchyIntegrator* p_ins_hier_integrator =
        dynamic_cast<INSStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());

    INSVCStaggeredHierarchyIntegrator* p_vc_ins_hier_integrator =
        dynamic_cast<INSVCStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());

    if (p_ins_hier_integrator)
    {
        d_mu_is_const = true;
    }
    else if (p_vc_ins_hier_integrator)
    {
        d_mu_is_const = p_vc_ins_hier_integrator->muIsConstant();
        if (!d_mu_is_const)
        {
            Pointer<CellVariable<NDIM, double> > mu_adv_diff_var = p_vc_ins_hier_integrator->getTransportedViscosityVariable();
            Pointer<CellVariable<NDIM, double> > mu_ins_var = p_vc_ins_hier_integrator->getViscosityVariable();
            Pointer<VariableContext> mu_ctx = var_db->getContext(d_object_name + "::mu_ctx");
            if (mu_adv_diff_var)
            {
                d_mu_idx = var_db->registerVariableAndContext(mu_adv_diff_var, mu_ctx, IntVector<NDIM>(GVISCOSITYG));
            }
            else if (mu_ins_var)
            {
                d_mu_idx = var_db->registerVariableAndContext(mu_ins_var, mu_ctx, IntVector<NDIM>(GVISCOSITYG));
            }
            else
            {
                TBOX_ERROR(
                    d_object_name << "::IBHydrodynamicSurfaceForceEvaluator():\n"
                                  << " no valid viscosity variable registered with the INS or AdvDiff integrators");
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << "::IBHydrodynamicSurfaceForceEvaluator():\n"
                                 << " unsupported INSHierarchyIntegrator type");
    }

    d_mu =
        d_mu_is_const ? d_fluid_solver->getStokesSpecifications()->getMu() : std::numeric_limits<double>::quiet_NaN();

    // Output data stream
    // Set up the streams for printing drag and torque
    if (SAMRAI_MPI::getRank() == 0)
    {
        std::ostringstream drag;
        drag << "Hydro_Force_" << d_ls_solid_var->getName();
        bool from_restart = RestartManager::getManager()->isFromRestart();
        if (from_restart)
        {
            d_hydro_force_stream = new std::ofstream(drag.str().c_str(), std::fstream::app);
            d_hydro_force_stream->precision(10);
        }
        else
        {
            d_hydro_force_stream = new std::ofstream(drag.str().c_str(), std::fstream::out);
            d_hydro_force_stream->precision(10);
        }
    }

    return;
} // IBHydrodynamicSurfaceForceEvaluator

IBHydrodynamicSurfaceForceEvaluator::~IBHydrodynamicSurfaceForceEvaluator()
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_ls_solid_idx);
    var_db->removePatchDataIndex(d_u_idx);
    var_db->removePatchDataIndex(d_p_idx);
    var_db->removePatchDataIndex(d_mu_idx);
    delete (d_hydro_force_stream);

    return;
} // ~IBHydrodynamicSurfaceForceEvaluator

void
IBHydrodynamicSurfaceForceEvaluator::computeHydrodynamicForce()
{
    // Get the current integrator time
    const double integrator_time = d_fluid_solver->getIntegratorTime();

    // Allocate required patch data
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_fluid_solver->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_ls_solid_idx, integrator_time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_u_idx, integrator_time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_p_idx, integrator_time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_mu_idx, integrator_time);
    }

    // Fill patch data and ghost cells
    fillPatchData(patch_hierarchy, integrator_time);

    // Object to hold net hydrodynamic force
    IBTK::Vector3d pressure_force = IBTK::Vector3d::Zero();
    IBTK::Vector3d viscous_force = IBTK::Vector3d::Zero();

    // Loop over side-centered DoFs of the computational domain to compute sum of n.(-pI + mu*(grad U + grad U)^T)
    // Note: n points outward from the solid into the fluid domain, which makes the above expression the force of the
    // fluid on the solid.
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        // Assumes that the structure is placed on the finest level
        if (ln < finest_ln) continue;

        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            double cell_vol = 1.0;
            for (unsigned int d = 0; d < NDIM; ++d) cell_vol *= patch_dx[d];

            // Get the required patch data
            Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(d_ls_solid_idx);
            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(d_u_idx);
            Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(d_p_idx);
            Pointer<CellData<NDIM, double> > mu_data;
            if (!d_mu_is_const) mu_data = patch->getPatchData(d_mu_idx);

            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                // Compute the required area element
                const double dS = cell_vol / patch_dx[axis];

                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SideIndex<NDIM> s_i(it(), axis, 0);
                    CellIndex<NDIM> c_l = s_i.toCell(0);
                    CellIndex<NDIM> c_u = s_i.toCell(1);
                    const double phi_lower = (*ls_solid_data)(c_l);
                    const double phi_upper = (*ls_solid_data)(c_u);

                    // If not within a band near the body, do not use this cell in the force calculation
                    if (sign((phi_lower - d_surface_contour_value) * (phi_upper - d_surface_contour_value)) >= 0.0)
                        continue;

                    // Compute the required unit normal
                    IBTK::Vector3d n = IBTK::Vector3d::Zero();
                    n(axis) = sign(phi_upper - phi_lower);

                    // Compute pressure on the face using simple averaging (n. -p I) * dA
                    IBTK::Vector3d pn = 0.5 * n * ((*p_data)(c_l) + (*p_data)(c_u));

                    // Compute the viscosity on the face using harmonic averaging, or setting to constant
                    const double mu_side =
                        d_mu_is_const ? d_mu :
                                        2.0 * (*mu_data)(c_l) * (*mu_data)(c_u) / ((*mu_data)(c_l) + (*mu_data)(c_u));

                    // Viscous traction force := n . mu(grad u + grad u ^ T) * dA
                    IBTK::Vector3d viscous_trac = IBTK::Vector3d::Zero();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            viscous_trac(axis) = (2.0 * mu_side) / (2.0 * patch_dx[axis]) *
                                                 ((*u_data)(SideIndex<NDIM>(c_u, axis, SideIndex<NDIM>::Upper)) -
                                                  (*u_data)(SideIndex<NDIM>(c_l, axis, SideIndex<NDIM>::Lower)));
                        }
                        else
                        {
                            CellIndex<NDIM> offset(0);
                            offset(d) = 1;

                            viscous_trac(d) =
                                mu_side / (2.0 * patch_dx[d]) *
                                    ((*u_data)(SideIndex<NDIM>(c_u + offset, axis, SideIndex<NDIM>::Lower)) -
                                     (*u_data)(SideIndex<NDIM>(c_u - offset, axis, SideIndex<NDIM>::Lower)))

                                +

                                mu_side / (2.0 * patch_dx[axis]) *
                                    ((*u_data)(SideIndex<NDIM>(c_u, d, SideIndex<NDIM>::Upper)) +
                                     (*u_data)(SideIndex<NDIM>(c_u, d, SideIndex<NDIM>::Lower)) -
                                     (*u_data)(SideIndex<NDIM>(c_l, d, SideIndex<NDIM>::Upper)) -
                                     (*u_data)(SideIndex<NDIM>(c_l, d, SideIndex<NDIM>::Lower))

                                         );
                        }
                    }

                    // Add up the pressure forces n.(-pI)dS
                    pressure_force += (-pn * dS);

                    // Add up the viscous forces n.(mu*(grad U + grad U)^T)dS
                    viscous_force += (n(axis) * viscous_trac * dS);
                }
            }
        }
    }
    // Print the hydrodynamic force to file.
    SAMRAI_MPI::sumReduction(pressure_force.data(), 3);
    SAMRAI_MPI::sumReduction(viscous_force.data(), 3);

    if (SAMRAI_MPI::getRank() == 0)
    {
        *d_hydro_force_stream << integrator_time << '\t' << pressure_force[0] << '\t' << pressure_force[1] << '\t'
                              << pressure_force[2] << '\t' << viscous_force[0] << '\t' << viscous_force[1] << '\t'
                              << viscous_force[2] << std::endl;
    }

    // Deallocate patch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_ls_solid_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_u_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_p_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_mu_idx);
    }

    return;

} // computeHydrodynamicForce

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBHydrodynamicSurfaceForceEvaluator::fillPatchData(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                   const double fill_time)
{
    // Fill ghost cells for level set
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;

    const int ls_solid_current_idx =
        var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getCurrentContext());
    InterpolationTransactionComponent ls_solid_transaction(d_ls_solid_idx,
                                                           ls_solid_current_idx,
                                                           /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                           /*USE_CF_INTERPOLATION*/ true,
                                                           /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                           /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                           /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                           d_adv_diff_solver->getPhysicalBcCoefs(d_ls_solid_var),
                                                           Pointer<VariableFillPattern<NDIM> >(NULL));
    Pointer<HierarchyGhostCellInterpolation> hier_ls_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_ls_bdry_fill->initializeOperatorState(ls_solid_transaction, patch_hierarchy);
    hier_ls_bdry_fill->setHomogeneousBc(false);
    hier_ls_bdry_fill->fillData(fill_time);

    // Fill ghost cells for velocity
    Pointer<SideVariable<NDIM, double> > u_var = d_fluid_solver->getVelocityVariable();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, d_fluid_solver->getCurrentContext());
    InterpolationTransactionComponent u_transaction(d_u_idx,
                                                    u_current_idx,
                                                    /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                    /*USE_CF_INTERPOLATION*/ true,
                                                    /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                    /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                    /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                    d_fluid_solver->getVelocityBoundaryConditions(),
                                                    Pointer<VariableFillPattern<NDIM> >(NULL));

    Pointer<HierarchyGhostCellInterpolation> hier_vel_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_vel_bdry_fill->initializeOperatorState(u_transaction, patch_hierarchy);
    hier_vel_bdry_fill->setHomogeneousBc(false);
    hier_vel_bdry_fill->fillData(fill_time);

    // Fill in ghost cells for viscosity, when necessary
    if (!d_mu_is_const)
    {
        INSVCStaggeredHierarchyIntegrator* p_vc_ins_hier_integrator =
            dynamic_cast<INSVCStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());
#if !defined(NDEBUG)
        TBOX_ASSERT(p_vc_ins_hier_integrator);
#endif
        Pointer<CellVariable<NDIM, double> > mu_adv_diff_var = p_vc_ins_hier_integrator->getTransportedViscosityVariable();
        Pointer<CellVariable<NDIM, double> > mu_ins_var = p_vc_ins_hier_integrator->getViscosityVariable();
        RobinBcCoefStrategy<NDIM>* mu_bc_coef = NULL;
        int mu_current_idx = -1;
        if (mu_adv_diff_var)
        {
            mu_current_idx =
                var_db->mapVariableAndContextToIndex(mu_adv_diff_var, d_adv_diff_solver->getCurrentContext());
            mu_bc_coef = (d_adv_diff_solver->getPhysicalBcCoefs(mu_adv_diff_var)).front();
        }
        else if (mu_ins_var)
        {
            mu_current_idx = var_db->mapVariableAndContextToIndex(mu_ins_var, d_fluid_solver->getCurrentContext());
            INSVCStaggeredHierarchyIntegrator* p_vc_ins_hier_integrator =
                dynamic_cast<INSVCStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());
            mu_bc_coef = p_vc_ins_hier_integrator->getViscosityBoundaryConditions();
        }
        else
        {
            TBOX_ERROR("This statement should not be reached");
        }

        InterpolationTransactionComponent mu_transaction_comp(d_mu_idx,
                                                              mu_current_idx,
                                                              /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                              /*USE_CF_INTERPOLATION*/ true,
                                                              /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                              /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                              /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                              mu_bc_coef,
                                                              Pointer<VariableFillPattern<NDIM> >(NULL));
        Pointer<HierarchyGhostCellInterpolation> hier_mu_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_mu_bdry_fill->initializeOperatorState(mu_transaction_comp, patch_hierarchy);
        hier_mu_bdry_fill->setHomogeneousBc(false);
        hier_mu_bdry_fill->fillData(fill_time);

        // Fill ghost cells for pressure
        Pointer<CellVariable<NDIM, double> > p_var = d_fluid_solver->getPressureVariable();
        const int p_current_idx = var_db->mapVariableAndContextToIndex(p_var, d_fluid_solver->getCurrentContext());
        INSStaggeredPressureBcCoef* p_ins_bc_coef =
            dynamic_cast<INSStaggeredPressureBcCoef*>(d_fluid_solver->getPressureBoundaryConditions());
        INSVCStaggeredPressureBcCoef* p_vc_ins_bc_coef =
            dynamic_cast<INSVCStaggeredPressureBcCoef*>(d_fluid_solver->getPressureBoundaryConditions());
        InterpolationTransactionComponent p_transaction_comp;
        if (p_ins_bc_coef)
        {
            p_ins_bc_coef->setTargetVelocityPatchDataIndex(d_u_idx);
            InterpolationTransactionComponent p_transaction_comp(d_p_idx,
                                                                 p_current_idx,
                                                                 /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                                 /*USE_CF_INTERPOLATION*/ true,
                                                                 /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                                 /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                 /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                 p_ins_bc_coef,
                                                                 Pointer<VariableFillPattern<NDIM> >(NULL));
            Pointer<HierarchyGhostCellInterpolation> hier_p_bdry_fill = new HierarchyGhostCellInterpolation();
            hier_p_bdry_fill->initializeOperatorState(p_transaction_comp, patch_hierarchy);
            hier_p_bdry_fill->setHomogeneousBc(false);
            hier_p_bdry_fill->fillData(fill_time);
        }
        else if (p_vc_ins_bc_coef)
        {
            p_vc_ins_bc_coef->setTargetVelocityPatchDataIndex(d_u_idx);
            InterpolationTransactionComponent p_transaction_comp(d_p_idx,
                                                                 p_current_idx,
                                                                 /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                                 /*USE_CF_INTERPOLATION*/ true,
                                                                 /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                                 /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                 /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                 p_vc_ins_bc_coef,
                                                                 Pointer<VariableFillPattern<NDIM> >(NULL));

            Pointer<HierarchyGhostCellInterpolation> hier_p_bdry_fill = new HierarchyGhostCellInterpolation();
            hier_p_bdry_fill->initializeOperatorState(p_transaction_comp, patch_hierarchy);
            hier_p_bdry_fill->setHomogeneousBc(false);
            hier_p_bdry_fill->fillData(fill_time);
        }
        else
        {
            TBOX_ERROR(d_object_name << "::IBHydrodynamicSurfaceForceEvaluator():\n"
                                     << " no valid pressure boundary condition object registered with INS integrator.\n"
                                     << " This statement should not have been reached");
        }
    }

    return;
} // fillPatchData

void
IBHydrodynamicSurfaceForceEvaluator::getFromInput(Pointer<Database> input_db)
{
    if (input_db->keyExists("surface_contour_value"))
    {
        d_surface_contour_value = input_db->getDouble("surface_contour_value");
    }

    return;
} // getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
