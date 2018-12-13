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
#include <utility>

#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "HierarchyDataOpsManager.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "SideData.h"
#include "SideIndex.h"
#include "ibamr/IBHydrodynamicSurfaceForceEvaluator.h"
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
static const int GLEVELSETG = 2;
static const int GVELOCITYG = 2;
static const int GPRESSUREG = 2;
static const int GVISCOSITYG = 2;
static const int GDENSITYG = 2;

inline int
sign(const double X)
{
    return ((X > 0) ? 1 : ((X < 0) ? -1 : 0));
} // sign

inline void
get_physical_coordinate(IBTK::Vector3d& side_coord, Pointer<Patch<NDIM> > patch, const SideIndex<NDIM>& side_idx)
{
    const int axis = side_idx.getAxis();
    Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* patch_X_lower = patch_geom->getXLower();
    const Box<NDIM>& patch_box = patch->getBox();
    const Index<NDIM>& patch_lower_idx = patch_box.lower();
    const double* const patch_dx = patch_geom->getDx();

    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis)
        {
            side_coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(side_idx(d) - patch_lower_idx(d)));
        }
        else
        {
            side_coord[d] =
                patch_X_lower[d] + patch_dx[d] * (static_cast<double>(side_idx(d) - patch_lower_idx(d)) + 0.5);
        }
    }
    return;
} // get_physical_coordinate

inline void
get_physical_coordinate(IBTK::Vector3d& cell_coord, Pointer<Patch<NDIM> > patch, const CellIndex<NDIM>& cell_idx)
{
    Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* patch_X_lower = patch_geom->getXLower();
    const Box<NDIM>& patch_box = patch->getBox();
    const Index<NDIM>& patch_lower_idx = patch_box.lower();
    const double* const patch_dx = patch_geom->getDx();

    for (int d = 0; d < NDIM; ++d)
    {
        cell_coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(cell_idx(d) - patch_lower_idx(d)) + 0.5);
    }
    return;
} // get_physical_coordinate

inline double
get_cell_average(Pointer<SideData<NDIM, double> > side_data, const CellIndex<NDIM>& cell_idx)
{
    double cell_average = 0.0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        SideIndex<NDIM> s_i_lower(cell_idx, axis, SideIndex<NDIM>::Lower);
        SideIndex<NDIM> s_i_upper(cell_idx, axis, SideIndex<NDIM>::Upper);
        cell_average += (*side_data)(s_i_lower) + (*side_data)(s_i_upper);
    }
    cell_average /= (NDIM * 2);
    return cell_average;
} // get_cell_average

inline IBTK::Vector3d
get_cell_acceleration(Pointer<SideData<NDIM, double> > u_data,
                      Pointer<SideData<NDIM, double> > u_old_data,
                      const CellIndex<NDIM>& cell_idx,
                      const double dt)
{
    IBTK::Vector3d accn = IBTK::Vector3d::Zero();
    for (int axis = 0; axis < NDIM; ++axis)
    {
        SideIndex<NDIM> s_i_lower(cell_idx, axis, SideIndex<NDIM>::Lower);
        SideIndex<NDIM> s_i_upper(cell_idx, axis, SideIndex<NDIM>::Upper);
        accn(axis) =
            0.5 * ((*u_data)(s_i_lower) + (*u_data)(s_i_upper) - (*u_old_data)(s_i_lower) - (*u_old_data)(s_i_upper));
    }

    accn /= dt;
    accn(1) += 9.81;

    return accn;
} // get_cell_acceleration
}

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

    auto p_ins_hier_integrator = dynamic_cast<INSStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());
    auto p_vc_ins_hier_integrator = dynamic_cast<INSVCStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());

#if !defined(NDEBUG)
    TBOX_ASSERT(d_ls_solid_var);
#endif
    Pointer<VariableContext> ls_ctx = var_db->getContext(d_object_name + "::ls_ctx");
    d_ls_solid_idx = var_db->registerVariableAndContext(d_ls_solid_var, ls_ctx, IntVector<NDIM>(GLEVELSETG));

    Pointer<SideVariable<NDIM, double> > u_var = d_fluid_solver->getVelocityVariable();
    Pointer<SideVariable<NDIM, double> > u_old_var = p_vc_ins_hier_integrator->getOldVelocityVariable();
#if !defined(NDEBUG)
    TBOX_ASSERT(u_var);
    TBOX_ASSERT(u_old_var);
#endif
    Pointer<VariableContext> u_ctx = var_db->getContext(d_object_name + "::u_ctx");
    d_u_idx = var_db->registerVariableAndContext(u_var, u_ctx, IntVector<NDIM>(GVELOCITYG));
    d_u_old_idx = var_db->registerVariableAndContext(u_old_var, u_ctx, IntVector<NDIM>(GVELOCITYG));

    Pointer<CellVariable<NDIM, double> > p_var = d_fluid_solver->getPressureVariable();
#if !defined(NDEBUG)
    TBOX_ASSERT(p_var);
#endif
    Pointer<VariableContext> p_ctx = var_db->getContext(d_object_name + "::p_ctx");
    d_p_idx = var_db->registerVariableAndContext(p_var, p_ctx, IntVector<NDIM>(GPRESSUREG));


    if (p_ins_hier_integrator)
    {
        d_mu_is_const = true;
        d_rho_is_const = true;
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

        Pointer<Variable<NDIM> > rho_ins_var = p_vc_ins_hier_integrator->getMassDensityVariable();
        Pointer<VariableContext> rho_ctx = var_db->getContext(d_object_name + "::rho_ctx");
        d_rho_idx = var_db->registerVariableAndContext(rho_ins_var, rho_ctx, IntVector<NDIM>(GDENSITYG));
    }
    else
    {
        TBOX_ERROR(d_object_name << "::IBHydrodynamicSurfaceForceEvaluator():\n"
                                 << " unsupported INSHierarchyIntegrator type");
    }

    d_mu =
        d_mu_is_const ? d_fluid_solver->getStokesSpecifications()->getMu() : std::numeric_limits<double>::quiet_NaN();
    d_rho =
        d_rho_is_const ? d_fluid_solver->getStokesSpecifications()->getRho() : std::numeric_limits<double>::quiet_NaN();

    return;
} // IBHydrodynamicSurfaceForceEvaluator

IBHydrodynamicSurfaceForceEvaluator::~IBHydrodynamicSurfaceForceEvaluator()
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_ls_solid_idx);
    var_db->removePatchDataIndex(d_u_idx);
    var_db->removePatchDataIndex(d_p_idx);
    var_db->removePatchDataIndex(d_mu_idx);
    delete d_hydro_force_stream;
    delete d_hydro_torque_stream;

    return;
} // ~IBHydrodynamicSurfaceForceEvaluator

void
IBHydrodynamicSurfaceForceEvaluator::computeHydrodynamicForceTorque(IBTK::Vector3d& pressure_force,
                                                                    IBTK::Vector3d& viscous_force,
                                                                    IBTK::Vector3d& pressure_torque,
                                                                    IBTK::Vector3d& viscous_torque,
                                                                    const IBTK::Vector3d& X0)
{
    double time = d_fluid_solver->getIntegratorTime();
    computeHydrodynamicForceTorque(
        pressure_force, viscous_force, pressure_torque, viscous_torque, X0, time, time, time);

    if (d_write_to_file && SAMRAI_MPI::getRank() == 0)
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
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        use_current_ctx = true;
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        use_new_ctx = true;
    }
    else
    {
        TBOX_ERROR("IBHydrodynamicSurfaceForceEvaluator::computeHydrodynamicForceTorque()"
                   << " Forces are evalauted at only current or new time. \n");
    }

    const double dt = new_time - current_time;

    // Allocate required patch data
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_fluid_solver->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_ls_solid_idx, time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_u_idx, time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_u_old_idx, time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_p_idx, time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_mu_idx, time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(d_rho_idx, time);
    }

    // Fill patch data and ghost cells
    fillPatchData(patch_hierarchy, time, use_current_ctx, use_new_ctx);

    // Zero out the vectors.
    pressure_force.setZero();
    viscous_force.setZero();
    pressure_torque.setZero();
    viscous_torque.setZero();
    IBTK::Vector3d r_vec = IBTK::Vector3d::Zero();
    IBTK::Vector3d dr = IBTK::Vector3d::Zero();

    // Loop over side-centered DoFs of the computational domain to compute sum of n.(-pI + mu*(grad U + grad U)^T)
    // Note: n points outward from the solid into the fluid domain, which makes the above expression the force of the
    // fluid on the solid.
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        // Assumes that the structure is placed on the finest level
        if (ln < finest_ln) continue;

        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
        const IntVector<NDIM>& ratio = level->getRatio();
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
            Pointer<SideData<NDIM, double> > u_old_data = patch->getPatchData(d_u_old_idx);
            Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(d_p_idx);
            Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(d_rho_idx);
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
                    // IBTK::Vector3d pn = 0.5 * n * ((*p_data)(c_l) + (*p_data)(c_u));

                    // Extract the pressure on the face using momentum equation in the normal direction.
                    const CellIndex<NDIM> c_f = phi_upper > 0 ? c_u : c_l;
                    IBTK::Vector3d n_ls = IBTK::Vector3d::Zero();
                    double h_interp = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        CellIndex<NDIM> offset;
                        offset(d) = 1;
                        n_ls(d) = ((*ls_solid_data)(c_f + offset) - (*ls_solid_data)(c_f - offset)) / (2 * patch_dx[d]);
                        h_interp += patch_dx[d] * patch_dx[d];
                    }
                    n_ls /= n_ls.norm();
                    h_interp = std::sqrt(h_interp);
                    get_physical_coordinate(r_vec, patch, c_f);
                    r_vec += h_interp * n_ls;
                    const CellIndex<NDIM> c_interp = IndexUtilities::getCellIndex(&r_vec[0], grid_geom, ratio);
                    const double p_interp = (*p_data)(c_interp);
                    const double rho_fluid_cell = get_cell_average(rho_data, c_f);
                    IBTK::Vector3d acc_fluid_cell = get_cell_acceleration(u_data, u_old_data, c_f, dt);
                    const double p_extracted_cell = p_interp + h_interp * rho_fluid_cell * (n_ls.dot(acc_fluid_cell));
                    const double p_extracted_face =
                        p_extracted_cell + (patch_dx[axis] / 2.0) * rho_fluid_cell * (n.dot(acc_fluid_cell));
                    IBTK::Vector3d pn = n * p_extracted_face;

                    // Get the relative coordinate from X0
                    get_physical_coordinate(r_vec, patch, s_i);
                    r_vec -= X0;

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
    SAMRAI_MPI::sumReduction(pressure_force.data(), 3);
    SAMRAI_MPI::sumReduction(viscous_force.data(), 3);
    SAMRAI_MPI::sumReduction(pressure_torque.data(), 3);
    SAMRAI_MPI::sumReduction(viscous_torque.data(), 3);

    // Deallocate patch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_ls_solid_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_u_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_u_old_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_p_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_mu_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_rho_idx);
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
    if (d_write_to_file && SAMRAI_MPI::getRank() == 0)
    {
        std::ostringstream force;
        force << "Hydro_Force_" << d_ls_solid_var->getName();
        bool from_restart = RestartManager::getManager()->isFromRestart();
        if (from_restart)
        {
            d_hydro_force_stream = new std::ofstream(force.str().c_str(), std::fstream::app);
            d_hydro_force_stream->precision(10);
        }
        else
        {
            d_hydro_force_stream = new std::ofstream(force.str().c_str(), std::fstream::out);
            d_hydro_force_stream->precision(10);
        }

        std::ostringstream torque;
        torque << "Hydro_Torque_" << d_ls_solid_var->getName();
        if (from_restart)
        {
            d_hydro_torque_stream = new std::ofstream(torque.str().c_str(), std::fstream::app);
            d_hydro_torque_stream->precision(10);
        }
        else
        {
            d_hydro_torque_stream = new std::ofstream(torque.str().c_str(), std::fstream::out);
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
    auto p_vc_ins_hier_integrator = dynamic_cast<INSVCStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());

    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;

    const int ls_solid_idx =
        use_current_ctx ?
            var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getCurrentContext()) :
            use_new_ctx ? var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getNewContext()) : -1;
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
    const int u_idx =
        use_current_ctx ?
            var_db->mapVariableAndContextToIndex(u_var, d_fluid_solver->getCurrentContext()) :
            use_new_ctx ? var_db->mapVariableAndContextToIndex(u_var, d_fluid_solver->getNewContext()) : -1;
    Pointer<SideVariable<NDIM, double> > u_old_var = p_vc_ins_hier_integrator->getOldVelocityVariable();
    const int u_old_idx =
        use_current_ctx ?
            var_db->mapVariableAndContextToIndex(u_old_var, p_vc_ins_hier_integrator->getCurrentContext()) :
            use_new_ctx ? var_db->mapVariableAndContextToIndex(u_old_var, p_vc_ins_hier_integrator->getNewContext()) :
                          -1;
    std::vector<InterpolationTransactionComponent> u_transaction(2);
    u_transaction[0] = InterpolationTransactionComponent(d_u_idx,
                                                         u_idx,
                                                         /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                         /*USE_CF_INTERPOLATION*/ true,
                                                         /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                         /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                         /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                         d_fluid_solver->getVelocityBoundaryConditions(),
                                                         Pointer<VariableFillPattern<NDIM> >(nullptr));
    u_transaction[1] = InterpolationTransactionComponent(d_u_old_idx,
                                                         u_old_idx,
                                                         /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                         /*USE_CF_INTERPOLATION*/ true,
                                                         /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                         /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                         /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                         p_vc_ins_hier_integrator->getVelocityBoundaryConditions(),
                                                         Pointer<VariableFillPattern<NDIM> >(nullptr));

    Pointer<HierarchyGhostCellInterpolation> hier_vel_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_vel_bdry_fill->initializeOperatorState(u_transaction, patch_hierarchy);
    hier_vel_bdry_fill->setHomogeneousBc(false);
    hier_vel_bdry_fill->fillData(fill_time);

    // Fill in ghost cells for viscosity, when necessary
    if (!d_mu_is_const)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(p_vc_ins_hier_integrator);
#endif
        Pointer<CellVariable<NDIM, double> > mu_adv_diff_var = p_vc_ins_hier_integrator->getTransportedViscosityVariable();
        Pointer<CellVariable<NDIM, double> > mu_ins_var = p_vc_ins_hier_integrator->getViscosityVariable();
        RobinBcCoefStrategy<NDIM>* mu_bc_coef = nullptr;
        int mu_idx = -1;
        if (mu_adv_diff_var)
        {
            mu_idx = use_current_ctx ?
                         var_db->mapVariableAndContextToIndex(mu_adv_diff_var, d_adv_diff_solver->getCurrentContext()) :
                         use_new_ctx ?
                         var_db->mapVariableAndContextToIndex(mu_adv_diff_var, d_adv_diff_solver->getNewContext()) :
                         -1;
            mu_bc_coef = (d_adv_diff_solver->getPhysicalBcCoefs(mu_adv_diff_var)).front();
        }
        else if (mu_ins_var)
        {
            mu_idx = use_current_ctx ?
                         var_db->mapVariableAndContextToIndex(mu_ins_var, d_fluid_solver->getCurrentContext()) :
                         use_new_ctx ?
                         var_db->mapVariableAndContextToIndex(mu_ins_var, d_fluid_solver->getNewContext()) :
                         -1;
            auto p_vc_ins_hier_integrator =
                dynamic_cast<INSVCStaggeredHierarchyIntegrator*>(d_fluid_solver.getPointer());
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
    const int p_idx =
        use_current_ctx ?
            var_db->mapVariableAndContextToIndex(p_var, d_fluid_solver->getCurrentContext()) :
            use_new_ctx ? var_db->mapVariableAndContextToIndex(p_var, d_fluid_solver->getNewContext()) : -1;
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

    // Fill ghost cells of density
    Pointer<Variable<NDIM> > rho_var = p_vc_ins_hier_integrator->getMassDensityVariable();
    const int rho_idx =
        use_current_ctx ?
            var_db->mapVariableAndContextToIndex(rho_var, p_vc_ins_hier_integrator->getCurrentContext()) :
            use_new_ctx ? var_db->mapVariableAndContextToIndex(rho_var, p_vc_ins_hier_integrator->getNewContext()) : -1;
    InterpolationTransactionComponent rho_transaction(
        d_rho_idx,
        rho_idx,
        /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
        /*USE_CF_INTERPOLATION*/ true,
        /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
        /*BDRY_EXTRAP_TYPE*/ "LINEAR",
        /*CONSISTENT_TYPE_2_BDRY*/ false,
        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr),
        Pointer<VariableFillPattern<NDIM> >(nullptr));

    Pointer<HierarchyGhostCellInterpolation> hier_rho_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_rho_bdry_fill->initializeOperatorState(rho_transaction, patch_hierarchy);
    hier_rho_bdry_fill->setHomogeneousBc(true);
    hier_rho_bdry_fill->fillData(fill_time);

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
