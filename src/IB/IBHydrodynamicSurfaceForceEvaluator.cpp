// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBHydrodynamicSurfaceForceEvaluator.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/INSStaggeredPressureBcCoef.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredPressureBcCoef.h"
#include "ibamr/StokesSpecifications.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_utilities.h"

#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "IntVector.h"
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
#include "VariableFillPattern.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "Eigen/Core"

#include <fstream>
#include <utility>
#include <vector>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

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
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBHydrodynamicSurfaceForceEvaluator::IBHydrodynamicSurfaceForceEvaluator(
    std::string object_name,
    Pointer<CellVariable<NDIM, double> > ls_solid_var,
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
    Pointer<INSHierarchyIntegrator> fluid_solver,
    Pointer<Database> db)
    : d_object_name(std::move(object_name)),
      d_ls_solid_var(ls_solid_var),
      d_adv_diff_solver(adv_diff_solver),
      d_fluid_solver(fluid_solver)

{
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

    auto p_ins_hier_integrator = dynamic_cast<INSStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());

    auto p_vc_ins_hier_integrator = dynamic_cast<INSVCStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());

    if (p_ins_hier_integrator)
    {
        d_mu_is_const = true;
    }
    else if (p_vc_ins_hier_integrator)
    {
        d_mu_is_const = p_vc_ins_hier_integrator->muIsConstant();
        if (!d_mu_is_const)
        {
            Pointer<CellVariable<NDIM, double> > mu_adv_diff_var =
                p_vc_ins_hier_integrator->getTransportedViscosityVariable();
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

    return;
} // IBHydrodynamicSurfaceForceEvaluator

IBHydrodynamicSurfaceForceEvaluator::~IBHydrodynamicSurfaceForceEvaluator()
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_ls_solid_idx);
    var_db->removePatchDataIndex(d_u_idx);
    var_db->removePatchDataIndex(d_p_idx);
    if (!d_mu_is_const)
    {
        var_db->removePatchDataIndex(d_mu_idx);
    }
    return;
} // ~IBHydrodynamicSurfaceForceEvaluator

void
IBHydrodynamicSurfaceForceEvaluator::computeHydrodynamicForceTorque(IBTK::Vector3d& pressure_force,
                                                                    IBTK::Vector3d& viscous_force,
                                                                    IBTK::Vector3d& pressure_torque,
                                                                    IBTK::Vector3d& viscous_torque,
                                                                    const IBTK::Vector3d& X0)
{
    const double time = d_fluid_solver->getIntegratorTime();
    computeHydrodynamicForceTorque(
        pressure_force, viscous_force, pressure_torque, viscous_torque, X0, time, time, time);

    if (d_write_to_file && IBTK_MPI::getRank() == 0)
    {
        *d_hydro_force_stream << time << '\t' << pressure_force[0] << '\t' << pressure_force[1] << '\t'
                              << pressure_force[2] << '\t' << viscous_force[0] << '\t' << viscous_force[1] << '\t'
                              << viscous_force[2] << std::endl;

        *d_hydro_torque_stream << time << '\t' << pressure_torque[0] << '\t' << pressure_torque[1] << '\t'
                               << pressure_torque[2] << '\t' << viscous_torque[0] << '\t' << viscous_torque[1] << '\t'
                               << viscous_torque[2] << std::endl;
    }
    return;

} // computeHydrodynamicForce

void
IBHydrodynamicSurfaceForceEvaluator::computeHydrodynamicForceTorque(IBTK::Vector3d& pressure_force,
                                                                    IBTK::Vector3d& viscous_force,
                                                                    IBTK::Vector3d& pressure_torque,
                                                                    IBTK::Vector3d& viscous_torque,
                                                                    const IBTK::Vector3d& X0,
                                                                    double time,
                                                                    double current_time,
                                                                    double new_time)
{
    bool use_current_ctx = false, use_new_ctx = false;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        use_current_ctx = true;
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        use_new_ctx = true;
    }
    else
    {
        TBOX_ERROR("IBHydrodynamicSurfaceForceEvaluator::computeHydrodynamicForceTorque()"
                   << " Forces are evalauted at only current or new time. \n");
    }

    // Allocate required patch data
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_fluid_solver->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_ls_solid_idx, time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_u_idx, time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_p_idx, time);
        if (!d_mu_is_const)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_mu_idx, time);
        }
    }

    // Fill patch data and ghost cells
    fillPatchData(patch_hierarchy, time, use_current_ctx, use_new_ctx);

    // Zero out the vectors.
    pressure_force.setZero();
    viscous_force.setZero();
    pressure_torque.setZero();
    viscous_torque.setZero();

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

            auto signof = [](const double x) { return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0); };

            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                // Compute the required area element
                const double dS = cell_vol / patch_dx[axis];

                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);
                    CellIndex<NDIM> c_l = s_i.toCell(SideIndex<NDIM>::Lower);
                    CellIndex<NDIM> c_u = s_i.toCell(SideIndex<NDIM>::Upper);
                    const double phi_lower = (*ls_solid_data)(c_l);
                    const double phi_upper = (*ls_solid_data)(c_u);

                    // If not within a band near the body, do not use this cell in the force calculation
                    if ((phi_lower - d_surface_contour_value) * (phi_upper - d_surface_contour_value) >= 0.0) continue;

                    // Compute the required unit normal
                    IBTK::Vector3d n = IBTK::Vector3d::Zero();
                    n(axis) = signof(phi_upper - phi_lower);

                    // Get the relative coordinate from X0
                    const IBTK::Vector3d r_vec = IBTK::IndexUtilities::getSideCenter<IBTK::Vector3d>(*patch, s_i) - X0;

                    // Compute pressure on the face using simple averaging (n. -p I) * dA
                    const IBTK::Vector3d pn = 0.5 * n * ((*p_data)(c_l) + (*p_data)(c_u));

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
                                     (*u_data)(SideIndex<NDIM>(c_l, d, SideIndex<NDIM>::Lower)));
                        }
                    }

                    // Add up the pressure forces n.(-pI)dS
                    // and pressure torques r X n.(-pI)dS
                    pressure_force += (-pn * dS);
                    pressure_torque += (r_vec.cross(-pn)) * dS;

                    // Add up the viscous forces n.(mu*(grad U + grad U)^T)dS
                    // and viscous torque r X n.(mu*(grad U + grad U)^T)dS
                    viscous_force += (n(axis) * viscous_trac * dS);
                    viscous_torque += r_vec.cross(n(axis) * viscous_trac) * dS;
                }
            }
        }
    }
    // Sum the net force and torque across processors.
    IBTK_MPI::sumReduction(pressure_force.data(), pressure_force.size());
    IBTK_MPI::sumReduction(viscous_force.data(), viscous_force.size());
    IBTK_MPI::sumReduction(pressure_torque.data(), pressure_torque.size());
    IBTK_MPI::sumReduction(viscous_torque.data(), viscous_force.size());

    // Deallocate patch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_ls_solid_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_u_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_p_idx);
        if (!d_mu_is_const)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_mu_idx);
        }
    }

    return;
} //  computeHydrodynamicForce

void
IBHydrodynamicSurfaceForceEvaluator::setSurfaceContourLevel(double s)
{
    d_surface_contour_value = s;
    return;
} // setSurfaceContourLevel

void
IBHydrodynamicSurfaceForceEvaluator::writeToFile(bool write_to_file)
{
    d_write_to_file = write_to_file;

    // Set up the streams for printing force and torque
    if (d_write_to_file && IBTK_MPI::getRank() == 0)
    {
        std::string force;
        force = "Hydro_Force_" + d_ls_solid_var->getName();
        bool from_restart = RestartManager::getManager()->isFromRestart();
        if (from_restart)
        {
            d_hydro_force_stream.reset(new std::ofstream(force.c_str(), std::fstream::app));
            d_hydro_force_stream->precision(10);
        }
        else
        {
            d_hydro_force_stream.reset(new std::ofstream(force.c_str(), std::fstream::out));
            d_hydro_force_stream->precision(10);
        }

        std::string torque;
        torque = "Hydro_Torque_" + d_ls_solid_var->getName();
        if (from_restart)
        {
            d_hydro_torque_stream.reset(new std::ofstream(torque.c_str(), std::fstream::app));
            d_hydro_torque_stream->precision(10);
        }
        else
        {
            d_hydro_torque_stream.reset(new std::ofstream(torque.c_str(), std::fstream::out));
            d_hydro_torque_stream->precision(10);
        }
    }
    return;
} // writeToFile

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBHydrodynamicSurfaceForceEvaluator::fillPatchData(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                   const double fill_time,
                                                   bool use_current_ctx,
                                                   bool use_new_ctx)
{
    // Fill ghost cells for level set
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;

    const int ls_solid_idx =
        use_current_ctx ? var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getCurrentContext()) :
        use_new_ctx     ? var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getNewContext()) :
                          IBTK::invalid_index;
    InterpolationTransactionComponent ls_solid_transaction(d_ls_solid_idx,
                                                           ls_solid_idx,
                                                           /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                           /*USE_CF_INTERPOLATION*/ true,
                                                           /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                           /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                           /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                           d_adv_diff_solver->getPhysicalBcCoefs(d_ls_solid_var),
                                                           Pointer<VariableFillPattern<NDIM> >(nullptr));
    Pointer<HierarchyGhostCellInterpolation> hier_ls_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_ls_bdry_fill->initializeOperatorState(ls_solid_transaction, patch_hierarchy);
    hier_ls_bdry_fill->setHomogeneousBc(false);
    hier_ls_bdry_fill->fillData(fill_time);

    // Fill ghost cells for velocity
    Pointer<SideVariable<NDIM, double> > u_var = d_fluid_solver->getVelocityVariable();
    const int u_idx = use_current_ctx ?
                                    var_db->mapVariableAndContextToIndex(u_var, d_fluid_solver->getCurrentContext()) :
                      use_new_ctx ? var_db->mapVariableAndContextToIndex(u_var, d_fluid_solver->getNewContext()) :
                                    IBTK::invalid_index;
    InterpolationTransactionComponent u_transaction(d_u_idx,
                                                    u_idx,
                                                    /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                    /*USE_CF_INTERPOLATION*/ true,
                                                    /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                    /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                    /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                    d_fluid_solver->getVelocityBoundaryConditions(),
                                                    Pointer<VariableFillPattern<NDIM> >(nullptr));

    Pointer<HierarchyGhostCellInterpolation> hier_vel_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_vel_bdry_fill->initializeOperatorState(u_transaction, patch_hierarchy);
    hier_vel_bdry_fill->setHomogeneousBc(false);
    hier_vel_bdry_fill->fillData(fill_time);

    // Fill in ghost cells for viscosity, when necessary
    if (!d_mu_is_const)
    {
        auto p_vc_ins_hier_integrator = dynamic_cast<INSVCStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());
#if !defined(NDEBUG)
        TBOX_ASSERT(p_vc_ins_hier_integrator);
#endif
        Pointer<CellVariable<NDIM, double> > mu_adv_diff_var =
            p_vc_ins_hier_integrator->getTransportedViscosityVariable();
        Pointer<CellVariable<NDIM, double> > mu_ins_var = p_vc_ins_hier_integrator->getViscosityVariable();
        RobinBcCoefStrategy<NDIM>* mu_bc_coef = nullptr;
        int mu_idx = IBTK::invalid_index;
        if (mu_adv_diff_var)
        {
            mu_idx = use_current_ctx ?
                         var_db->mapVariableAndContextToIndex(mu_adv_diff_var, d_adv_diff_solver->getCurrentContext()) :
                     use_new_ctx ?
                         var_db->mapVariableAndContextToIndex(mu_adv_diff_var, d_adv_diff_solver->getNewContext()) :
                         IBTK::invalid_index;
            mu_bc_coef = (d_adv_diff_solver->getPhysicalBcCoefs(mu_adv_diff_var)).front();
        }
        else if (mu_ins_var)
        {
            mu_idx = use_current_ctx ?
                                   var_db->mapVariableAndContextToIndex(mu_ins_var, d_fluid_solver->getCurrentContext()) :
                     use_new_ctx ? var_db->mapVariableAndContextToIndex(mu_ins_var, d_fluid_solver->getNewContext()) :
                                   IBTK::invalid_index;
            mu_bc_coef = p_vc_ins_hier_integrator->getViscosityBoundaryConditions();
        }
        else
        {
            TBOX_ERROR("This statement should not be reached");
        }

        InterpolationTransactionComponent mu_transaction_comp(d_mu_idx,
                                                              mu_idx,
                                                              /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                              /*USE_CF_INTERPOLATION*/ true,
                                                              /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                              /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                              /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                              mu_bc_coef,
                                                              Pointer<VariableFillPattern<NDIM> >(nullptr));
        Pointer<HierarchyGhostCellInterpolation> hier_mu_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_mu_bdry_fill->initializeOperatorState(mu_transaction_comp, patch_hierarchy);
        hier_mu_bdry_fill->setHomogeneousBc(false);
        hier_mu_bdry_fill->fillData(fill_time);
    }

    // Fill ghost cells for pressure
    Pointer<CellVariable<NDIM, double> > p_var = d_fluid_solver->getPressureVariable();
    const int p_idx = use_current_ctx ?
                                    var_db->mapVariableAndContextToIndex(p_var, d_fluid_solver->getCurrentContext()) :
                      use_new_ctx ? var_db->mapVariableAndContextToIndex(p_var, d_fluid_solver->getNewContext()) :
                                    IBTK::invalid_index;
    auto p_ins_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(d_fluid_solver->getPressureBoundaryConditions());
    auto p_vc_ins_bc_coef =
        dynamic_cast<INSVCStaggeredPressureBcCoef*>(d_fluid_solver->getPressureBoundaryConditions());
    InterpolationTransactionComponent p_transaction_comp;
    if (p_ins_bc_coef)
    {
        p_ins_bc_coef->setTargetVelocityPatchDataIndex(d_u_idx);
        InterpolationTransactionComponent p_transaction_comp(d_p_idx,
                                                             p_idx,
                                                             /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                             /*USE_CF_INTERPOLATION*/ true,
                                                             /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                             /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                             /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                             p_ins_bc_coef,
                                                             Pointer<VariableFillPattern<NDIM> >(nullptr));
        Pointer<HierarchyGhostCellInterpolation> hier_p_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_p_bdry_fill->initializeOperatorState(p_transaction_comp, patch_hierarchy);
        hier_p_bdry_fill->setHomogeneousBc(false);
        hier_p_bdry_fill->fillData(fill_time);
    }
    else if (p_vc_ins_bc_coef)
    {
        p_vc_ins_bc_coef->setTargetVelocityPatchDataIndex(d_u_idx);
        InterpolationTransactionComponent p_transaction_comp(d_p_idx,
                                                             p_idx,
                                                             /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                             /*USE_CF_INTERPOLATION*/ true,
                                                             /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                             /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                             /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                             p_vc_ins_bc_coef,
                                                             Pointer<VariableFillPattern<NDIM> >(nullptr));

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

    return;
} // fillPatchData

void
IBHydrodynamicSurfaceForceEvaluator::getFromInput(Pointer<Database> input_db)
{
    if (input_db->keyExists("surface_contour_value"))
    {
        d_surface_contour_value = input_db->getDouble("surface_contour_value");
    }

    if (input_db->keyExists("write_to_file"))
    {
        d_write_to_file = input_db->getBool("write_to_file");
    }

    return;
} // getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
