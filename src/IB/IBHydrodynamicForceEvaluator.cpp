// Filename: IBHydrodynamicForceEvaluator.cpp
// Created on 22 Oct 2016 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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

#include "ibamr/IBHydrodynamicForceEvaluator.h"
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
#include "ibamr/INSStaggeredPressureBcCoef.h"
#include "ibamr/namespaces.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/IndexUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBHydrodynamicForceEvaluator::IBHydrodynamicForceEvaluator(const std::string& object_name,
                                                           double rho,
                                                           double mu,
                                                           double current_time,
                                                           bool register_for_restart)
{
    d_object_name = object_name;
    d_rho = rho;
    d_mu = mu;
    d_current_time = current_time;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_u_var = new SideVariable<NDIM, double>(d_object_name + "::u_var", 1);
    d_p_var = new CellVariable<NDIM, double>(d_object_name + "::p_var", 1);
    Pointer<VariableContext> u_ctx = var_db->getContext(d_object_name + "::u_ctx");
    Pointer<VariableContext> p_ctx = var_db->getContext(d_object_name + "::p_ctx");
    d_u_idx = var_db->registerVariableAndContext(d_u_var, u_ctx, /*ghost_width*/ 1);
    d_p_idx = var_db->registerVariableAndContext(d_p_var, p_ctx, /*ghost_width*/ 1);

    Pointer<SideVariable<NDIM, double> > wgt_var = new SideVariable<NDIM, double>(d_object_name + "::wgt_var", 1);
    Pointer<VariableContext> face_wgt_ctx = var_db->getContext(d_object_name + "::face_wgt_ctx");
    Pointer<VariableContext> vol_wgt_ctx = var_db->getContext(d_object_name + "::vol_wgt_ctx");
    d_face_wgt_sc_idx = var_db->registerVariableAndContext(wgt_var, face_wgt_ctx, /*ghost_width*/ 0);
    d_vol_wgt_sc_idx = var_db->registerVariableAndContext(wgt_var, vol_wgt_ctx, /*ghost_width*/ 0);

    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }
    return;
} // IBHydrodynamicForceEvaluator

IBHydrodynamicForceEvaluator::~IBHydrodynamicForceEvaluator()
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_u_idx);
    var_db->removePatchDataIndex(d_p_idx);
    var_db->removePatchDataIndex(d_face_wgt_sc_idx);
    var_db->removePatchDataIndex(d_vol_wgt_sc_idx);

    return;
} // ~IBHydrodynamicForceEvaluator

void
IBHydrodynamicForceEvaluator::registerStructure(IBTK::Vector3d& box_X_lower,
                                                IBTK::Vector3d& box_X_upper,
                                                Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                const IBTK::Vector3d& box_vel,
                                                int strct_id)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hydro_objs.find(strct_id) == d_hydro_objs.end());
#endif

    IBHydrodynamicForceObject force_obj;
    force_obj.strct_id = strct_id;

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (!from_restart)
    {
        // Ensure that box is aligned to grid cell sides at coarsest level
	const int coarsest_ln = 0;
	Pointer<PatchLevel<NDIM> > coarsest_level = patch_hierarchy->getPatchLevel(coarsest_ln);
	const Pointer<CartesianGridGeometry<NDIM> > coarsest_grid_geom = coarsest_level->getGridGeometry();
	const double* const dx_coarsest = coarsest_grid_geom->getDx();
	const double* const grid_X_lower = coarsest_grid_geom->getXLower();
	bool modified_box = false;
        const IBTK::Vector3d box_X_lower_old = box_X_lower;
        const IBTK::Vector3d box_X_upper_old = box_X_upper;

        for (int d = 0; d < NDIM; ++d)
	{
	    double num_cells_lower = (box_X_lower[d]-grid_X_lower[d])/dx_coarsest[d];
	    double num_cells_upper = (box_X_upper[d]-grid_X_lower[d])/dx_coarsest[d];
	    
	    if (!MathUtilities<double>::equalEps(num_cells_lower, floor(num_cells_lower)))
	    {
                TBOX_WARNING("Lower side of integration box is not aligned with sides on coarsest level in dimension "
                             << d
                             << ". Modifying coordinate to nearest box side\n");
                const int N = floor(num_cells_lower);
		box_X_lower[d] = grid_X_lower[d] + dx_coarsest[d] * N;
		modified_box = true;
	    }
	    
	    if (!MathUtilities<double>::equalEps(num_cells_upper, floor(num_cells_upper)))
	    {
                TBOX_WARNING("Upper side of integration box is not aligned with sides on coarsest level in dimension "
                             << d
                             << ". Modifying coordinate to nearest box side\n");
                const int N = ceil(num_cells_upper);
		box_X_upper[d] = grid_X_lower[d] + dx_coarsest[d] * N;
		modified_box = true;
	    }
	}
	
	if (modified_box)
	{
            pout << "IBHydrodynamicForceEvaluator::registerStructure: integration box modified from\n"
                 << "[" << box_X_lower_old[0] << ", " << box_X_upper_old[0] << "]"
                 << " x [" << box_X_lower_old[1] << ", " << box_X_upper_old[1] << "]"
                 << " x [" << box_X_lower_old[2] << ", " << box_X_upper_old[2] << "]\n"
                 << "to\n"
                 << "[" << box_X_lower[0] << ", " << box_X_upper[0] << "]"
                 << " x [" << box_X_lower[1] << ", " << box_X_upper[1] << "]"
                 << " x [" << box_X_lower[2] << ", " << box_X_upper[2] << "]\n";
        }

        force_obj.box_u_current = box_vel;
        force_obj.box_X_lower_current = box_X_lower;
        force_obj.box_X_upper_current = box_X_upper;
	force_obj.box_vol_current = (box_X_upper[0] - box_X_lower[0]) * (box_X_upper[1] - box_X_lower[1])
#if (NDIM == 3)
				    * (box_X_upper[2] - box_X_lower[2])
#endif
				    ;

        force_obj.F_current.setZero();
        force_obj.T_current.setZero();
        force_obj.P_current.setZero();
        force_obj.L_current.setZero();
        force_obj.P_box_current.setZero();
        force_obj.L_box_current.setZero();
        force_obj.r0.setZero();
    }
    else
    {
        Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
        Pointer<Database> db;
        if (restart_db->isDatabase(d_object_name))
        {
            db = restart_db->getDatabase(d_object_name);
        }
        else
        {
            TBOX_ERROR("IBHydrodynamicForceEvaluator::registerStructure(). Restart database corresponding to "
                       << d_object_name
                       << " not found in restart file.\n");
        }

        std::ostringstream F, T, P, L, P_box, L_box, X_lo, X_hi, r_or, vol_curr;
        F << "F_" << strct_id;
        T << "T_" << strct_id;
        P << "P_" << strct_id;
        L << "L_" << strct_id;
        P_box << "P_box_" << strct_id;
        L_box << "L_box_" << strct_id;
        X_lo << "X_lo_" << strct_id;
        X_hi << "X_hi_" << strct_id;
        r_or << "r_or_" << strct_id;
        vol_curr << "vol_curr_" << strct_id;

        db->getDoubleArray(F.str(), force_obj.F_current.data(), 3);
        db->getDoubleArray(T.str(), force_obj.T_current.data(), 3);
        db->getDoubleArray(P.str(), force_obj.P_current.data(), 3);
        db->getDoubleArray(L.str(), force_obj.L_current.data(), 3);
        db->getDoubleArray(P_box.str(), force_obj.P_box_current.data(), 3);
        db->getDoubleArray(L_box.str(), force_obj.L_box_current.data(), 3);
        db->getDoubleArray(X_lo.str(), force_obj.box_X_lower_current.data(), 3);
        db->getDoubleArray(X_hi.str(), force_obj.box_X_upper_current.data(), 3);
        db->getDoubleArray(r_or.str(), force_obj.r0.data(), 3);
        force_obj.box_vol_current = db->getDouble(vol_curr.str());
    }

    // Set up the streams for printing drag and torque
    if (SAMRAI_MPI::getRank() == 0)
    {
        std::ostringstream drag, torque;
        drag << "Drag_CV_strct_id_" << strct_id;
        torque << "Torque_CV_strct_id_" << strct_id;
        if (from_restart)
        {
            force_obj.drag_CV_stream = new std::ofstream(drag.str().c_str(), std::fstream::app);
            force_obj.torque_CV_stream = new std::ofstream(torque.str().c_str(), std::fstream::app);
            (force_obj.drag_CV_stream)->precision(10);
            (force_obj.torque_CV_stream)->precision(10);
        }
        else
        {
            force_obj.drag_CV_stream = new std::ofstream(drag.str().c_str(), std::fstream::out);
            force_obj.torque_CV_stream = new std::ofstream(torque.str().c_str(), std::fstream::out);
            (force_obj.drag_CV_stream)->precision(10);
            (force_obj.torque_CV_stream)->precision(10);
        }
    }

    d_hydro_objs[strct_id] = force_obj;
    return;

} // registerStructure

void
IBHydrodynamicForceEvaluator::updateStructureDomain(const IBTK::Vector3d& box_vel_new,
                                                    double dt,
                                                    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                    int strct_id)

{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hydro_objs.find(strct_id) != d_hydro_objs.end());
#endif

    IBHydrodynamicForceEvaluator::IBHydrodynamicForceObject& fobj = d_hydro_objs[strct_id];
    fobj.box_u_new = box_vel_new;
    fobj.box_X_lower_new = fobj.box_X_lower_current + box_vel_new * dt;
    fobj.box_X_upper_new = fobj.box_X_upper_current + box_vel_new * dt;

    // Assert that the volume of the box has not changed
    fobj.box_vol_new = (fobj.box_X_upper_new[0] - fobj.box_X_lower_new[0]) *
                       (fobj.box_X_upper_new[1] - fobj.box_X_lower_new[1])
#if (NDIM == 3)
                       * (fobj.box_X_upper_new[2] - fobj.box_X_lower_new[2])
#endif
        ;

    TBOX_ASSERT(MathUtilities<double>::equalEps(fobj.box_vol_current, fobj.box_vol_new));

    return;

} // updateStructureDomain

void
IBHydrodynamicForceEvaluator::setTorqueOrigin(const IBTK::Vector3d& X0, int strct_id)
{
// Set the origin of the position vector about which torques are computed
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hydro_objs.find(strct_id) != d_hydro_objs.end());
#endif

    IBHydrodynamicForceEvaluator::IBHydrodynamicForceObject& fobj = d_hydro_objs[strct_id];
    fobj.r0 = X0;
}

void
IBHydrodynamicForceEvaluator::preprocessIntegrateData(double /*current_time*/, double /*new_time*/)
{
    pout << "WARNING:: IBHydrodynamicForceEvaluator::preprocessIntegrateData() not implemented.\n";

    return;
} // preprocessIntegrateData

const IBHydrodynamicForceEvaluator::IBHydrodynamicForceObject&
IBHydrodynamicForceEvaluator::getHydrodynamicForceObject(int strct_id)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hydro_objs.find(strct_id) != d_hydro_objs.end());
#endif

    const IBHydrodynamicForceEvaluator::IBHydrodynamicForceObject& force_obj = d_hydro_objs[strct_id];
    return force_obj;

} // getHydrodynamicForceObject

void
IBHydrodynamicForceEvaluator::computeLaggedMomentumIntegral(
    int u_old_idx,
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_src_bc_coef)
{
    resetFaceAreaWeight(patch_hierarchy);
    resetFaceVolWeight(patch_hierarchy);
    fillPatchData(u_old_idx, -1, patch_hierarchy, u_src_bc_coef, NULL, d_current_time);

    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Whether or not the simulation has adaptive mesh refinement
    const bool amr_case = (coarsest_ln != finest_ln);

    for (std::map<int, IBHydrodynamicForceObject>::iterator it = d_hydro_objs.begin(); it != d_hydro_objs.end(); ++it)
    {
        IBHydrodynamicForceObject& fobj = it->second;

        // Compute the momentum integral:= (rho * u * dv) for the previous time step (integral is over new control
        // volume)
        fobj.P_box_current.setZero();

        // Compute the rotational momentum integral:= (rho * r x u * dv) for the previous time step (integral is over
        // new control volume)
        fobj.L_box_current.setZero();

        // Coordinate of the side index and r vector needed for cross product
        IBTK::Vector3d side_coord, r_vec;

        for (int ln = finest_ln; ln >= coarsest_ln; --ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            Box<NDIM> integration_box(
                IndexUtilities::getCellIndex(fobj.box_X_lower_new.data(), level->getGridGeometry(), level->getRatio()),
                IndexUtilities::getCellIndex(fobj.box_X_upper_new.data(), level->getGridGeometry(), level->getRatio()));

            // Shorten the integration box so it only includes the control volume
            integration_box.upper() -= 1;

            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const bool boxes_intersect = patch_box.intersects(integration_box);
                if (!boxes_intersect) continue;

                // Part of the box on this patch.
                Box<NDIM> trim_box = patch_box * integration_box;

                // Loop over the box and compute momentum.
                Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(d_u_idx);
                Pointer<SideData<NDIM, double> > vol_sc_data = patch->getPatchData(d_vol_wgt_sc_idx);

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(trim_box, axis)); b; b++)
                    {
                        const CellIndex<NDIM>& cell_idx = *b;
                        const SideIndex<NDIM> side_idx(cell_idx, axis, SideIndex<NDIM>::Lower);
                        const double& u_axis = (*u_data)(side_idx);
                        const double& vol = (*vol_sc_data)(side_idx);
                        double dV;

                        // Check if cell is a CV boundary
                        const bool lower_bdry_vel = (cell_idx(axis) == (integration_box.lower())(axis));
                        const bool upper_bdry_vel = (cell_idx(axis) == (integration_box.upper())(axis) + 1);

                        // Check if CV boundary intersects a patch boundary
                        const bool lower_patch_bdry_eq_box_bdry =
                            ((patch_box.lower())(axis) == (integration_box.lower())(axis));
                        const bool upper_patch_bdry_eq_box_bdry =
                            ((patch_box.upper())(axis) + 1 == (integration_box.upper())(axis) + 1);

                        if (!amr_case)
                        {
                            /* Uniform mesh scaling correction
                             * If the velocity is on the CV boundary, scale the volume element by 1/2
                             * If the patch boundary equals the CV boundary, then volume element is correct (dx * dy)/2
                             */

                            const bool scale_dV = (lower_bdry_vel && !lower_patch_bdry_eq_box_bdry) ||
                                                  (upper_bdry_vel && !upper_patch_bdry_eq_box_bdry);

                            dV = scale_dV ? 0.5 * vol : vol;
                        }
                        else
                        {
                            /* Adaptive mesh scaling correction
                             * If on a CV boundary, set dV to (dx * dy)/2, using the patch grid spacing
                             * If vol == 0, don't change anything
                             */

                            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                            const double* const patch_dx = patch_geom->getDx();
                            const double box_edge_dV = 0.5 * patch_dx[0] * patch_dx[1]
#if (NDIM == 3)
                                                       * patch_dx[2]
#endif
                                ;

                            const bool modify_dV = (lower_bdry_vel || upper_bdry_vel) && vol > 0;
                            dV = modify_dV ? box_edge_dV : vol;
                        }

                        fobj.P_box_current(axis) += d_rho * u_axis * dV;

                        // Compute angular momentum by looping over all the sides in one axis direction

                        if (axis == 0)
                        {
                            // Get the coordinate of the side index and r vector
                            side_coord.setZero();
                            getPhysicalCoordinateFromSideIndex(side_coord, level, patch, side_idx, axis);
                            r_vec = side_coord - fobj.r0;
                            IBTK::Vector3d u_vec = IBTK::Vector3d::Zero();
                            u_vec(axis) = u_axis;

                            for (int d = 0; d < NDIM; ++d)
                            {
                                if (d == axis) continue;

                                CellIndex<NDIM> cell_left_idx = cell_idx;
                                cell_left_idx(axis) -= 1;
                                u_vec(d) =
                                    0.25 * ((*u_data)(SideIndex<NDIM>(cell_left_idx, d, SideIndex<NDIM>::Lower)) +
                                            (*u_data)(SideIndex<NDIM>(cell_left_idx, d, SideIndex<NDIM>::Upper)) +
                                            (*u_data)(SideIndex<NDIM>(cell_idx, d, SideIndex<NDIM>::Lower)) +
                                            (*u_data)(SideIndex<NDIM>(cell_idx, d, SideIndex<NDIM>::Upper)));
                            }

                            fobj.L_box_current += d_rho * r_vec.cross(u_vec) * dV;
                        }
                    }
                }
            }
        }

        SAMRAI_MPI::sumReduction(fobj.P_box_current.data(), 3);
        SAMRAI_MPI::sumReduction(fobj.L_box_current.data(), 3);
    }

    return;

} // computeLaggedMomentumIntegral

void
IBHydrodynamicForceEvaluator::updateStructureMomentum(const IBTK::Vector3d& P_strct_new,
                                                      const IBTK::Vector3d& L_strct_new,
                                                      int strct_id)

{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hydro_objs.find(strct_id) != d_hydro_objs.end());
#endif

    IBHydrodynamicForceEvaluator::IBHydrodynamicForceObject& force_obj = d_hydro_objs[strct_id];

    force_obj.P_new = P_strct_new;
    force_obj.L_new = L_strct_new;

    return;

} // updateStructureMomentum

void
IBHydrodynamicForceEvaluator::computeHydrodynamicForce(int u_idx,
                                                       int p_idx,
                                                       int /*f_idx*/,
                                                       Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                       double dt,
                                                       const std::vector<RobinBcCoefStrategy<NDIM>*>& u_src_bc_coef,
                                                       RobinBcCoefStrategy<NDIM>* p_src_bc_coef)
{
    resetFaceAreaWeight(patch_hierarchy);
    resetFaceVolWeight(patch_hierarchy);
    fillPatchData(u_idx, p_idx, patch_hierarchy, u_src_bc_coef, p_src_bc_coef, d_current_time + dt);

    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Whether or not the simulation has adaptive mesh refinement
    const bool amr_case = (coarsest_ln != finest_ln);

    for (std::map<int, IBHydrodynamicForceObject>::iterator it = d_hydro_objs.begin(); it != d_hydro_objs.end(); ++it)
    {
        IBHydrodynamicForceObject& fobj = it->second;

        // Compute the momentum integral:= (rho * u * dv)
        fobj.P_box_new.setZero();

        // Compute the rotational momentum integral:= (rho * r x u * dv) for the new time step (integral is over new
        // control volume)
        fobj.L_box_new.setZero();

        // Coordinate of the side index and r vector needed for cross product
        IBTK::Vector3d side_coord, r_vec;

        for (int ln = finest_ln; ln >= coarsest_ln; --ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            Box<NDIM> integration_box(
                IndexUtilities::getCellIndex(fobj.box_X_lower_new.data(), level->getGridGeometry(), level->getRatio()),
                IndexUtilities::getCellIndex(fobj.box_X_upper_new.data(), level->getGridGeometry(), level->getRatio()));

            // Shorten the integration box so it only includes the control volume
            integration_box.upper() -= 1;
	    
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const bool boxes_intersect = patch_box.intersects(integration_box);
                if (!boxes_intersect) continue;

                // Part of the box on this patch.
                Box<NDIM> trim_box = patch_box * integration_box;
		
                // Loop over the box and compute momentum.
                Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(d_u_idx);
                Pointer<SideData<NDIM, double> > vol_sc_data = patch->getPatchData(d_vol_wgt_sc_idx);
		
		for (int axis = 0; axis < NDIM; ++axis)
		{
		    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(trim_box, axis)); b; b++)
		    {
                        const CellIndex<NDIM>& cell_idx = *b;
                        const SideIndex<NDIM> side_idx(cell_idx, axis, SideIndex<NDIM>::Lower);
                        const double& u_axis = (*u_data)(side_idx);
                        const double& vol = (*vol_sc_data)(side_idx);
                        double dV;

                        // Check if cell is a CV boundary
                        const bool lower_bdry_vel = (cell_idx(axis) == (integration_box.lower())(axis));
                        const bool upper_bdry_vel = (cell_idx(axis) == (integration_box.upper())(axis) + 1);

                        // Check if CV boundary intersects a patch boundary
                        const bool lower_patch_bdry_eq_box_bdry =
                            ((patch_box.lower())(axis) == (integration_box.lower())(axis));
                        const bool upper_patch_bdry_eq_box_bdry =
                            ((patch_box.upper())(axis) + 1 == (integration_box.upper())(axis) + 1);

                        if (!amr_case)
                        {
                            /* Uniform mesh scaling correction
                             * If the velocity is on the CV boundary, scale the volume element by 1/2
                             * If the patch boundary equals the CV boundary, then volume element is correct (dx * dy)/2
                             */
                            const bool scale_dV = (lower_bdry_vel && !lower_patch_bdry_eq_box_bdry) ||
                                                  (upper_bdry_vel && !upper_patch_bdry_eq_box_bdry);

                            dV = scale_dV ? 0.5 * vol : vol;
                        }
                        else
                        {
                            /* Adaptive mesh scaling correction
                             * If on a CV boundary, set dV to (dx * dy)/2, using the patch grid spacing
                             * If vol == 0, don't change anything
                             */

                            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                            const double* const patch_dx = patch_geom->getDx();
                            const double box_edge_dV = 0.5 * patch_dx[0] * patch_dx[1]
#if (NDIM == 3)
                                                       * patch_dx[2]
#endif
                                ;

                            const bool modify_dV = (lower_bdry_vel || upper_bdry_vel) && vol > 0;
                            dV = modify_dV ? box_edge_dV : vol;
                        }

                        fobj.P_box_new(axis) += d_rho * u_axis * dV;

                        // Compute angular momentum by looping over all the sides in one axis direction

                        if (axis == 0)
                        {
                            // Get the coordinate of the side index and r vector
                            side_coord.setZero();
                            getPhysicalCoordinateFromSideIndex(side_coord, level, patch, side_idx, axis);
                            r_vec = side_coord - fobj.r0;
                            IBTK::Vector3d u_vec = IBTK::Vector3d::Zero();
                            u_vec(axis) = u_axis;

                            for (int d = 0; d < NDIM; ++d)
                            {
                                if (d == axis) continue;

                                CellIndex<NDIM> cell_left_idx = cell_idx;
                                cell_left_idx(axis) -= 1;
                                u_vec(d) =
                                    0.25 * ((*u_data)(SideIndex<NDIM>(cell_left_idx, d, SideIndex<NDIM>::Lower)) +
                                            (*u_data)(SideIndex<NDIM>(cell_left_idx, d, SideIndex<NDIM>::Upper)) +
                                            (*u_data)(SideIndex<NDIM>(cell_idx, d, SideIndex<NDIM>::Lower)) +
                                            (*u_data)(SideIndex<NDIM>(cell_idx, d, SideIndex<NDIM>::Upper)));
                            }

                            fobj.L_box_new += d_rho * r_vec.cross(u_vec) * dV;
                        }
                    }
                }
            }
        }

        SAMRAI_MPI::sumReduction(fobj.P_box_new.data(), 3);
        SAMRAI_MPI::sumReduction(fobj.L_box_new.data(), 3);

        // Compute surface integral term.
        IBTK::Vector3d trac, torque_trac;
        trac.setZero();
        torque_trac.setZero();

        for (int ln = finest_ln; ln >= coarsest_ln; --ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            Box<NDIM> integration_box(
                IndexUtilities::getCellIndex(fobj.box_X_lower_new.data(), level->getGridGeometry(), level->getRatio()),
                IndexUtilities::getCellIndex(fobj.box_X_upper_new.data(), level->getGridGeometry(), level->getRatio()));
	    
	    // Shorten the integration box so it only includes the control volume
	    integration_box.upper() -= 1;
	    
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const patch_dx = patch_geom->getDx();
                const bool boxes_intersect = patch_box.intersects(integration_box);
                if (!boxes_intersect) continue;

                // Store boxes corresponding to integration domain boundaries.
                boost::array<boost::array<Box<NDIM>, 2>, NDIM> bdry_boxes;
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    Box<NDIM> bdry_box;

                    static const int lower_side = 0;
                    bdry_box = integration_box;
                    bdry_box.upper()(axis) = bdry_box.lower()(axis);
                    bdry_boxes[axis][lower_side] = bdry_box;

                    static const int upper_side = 1;
                    bdry_box = integration_box;
                    bdry_box.lower()(axis) = bdry_box.upper()(axis);
                    bdry_boxes[axis][upper_side] = bdry_box;
                }

                // Integrate over boundary boxes.
                Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(d_p_idx);
                Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(d_u_idx);
                Pointer<SideData<NDIM, double> > face_sc_data = patch->getPatchData(d_face_wgt_sc_idx);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (int upperlower = 0; upperlower <= 1; ++upperlower)
                    {
                        const Box<NDIM>& side_box = bdry_boxes[axis][upperlower];
                        if (!patch_box.intersects(side_box)) continue;

                        Box<NDIM> trim_box = patch_box * side_box;
                        IBTK::Vector3d n = IBTK::Vector3d::Zero();
                        n(axis) = upperlower ? 1 : -1;
                        for (Box<NDIM>::Iterator b(trim_box); b; b++)
                        {
                            const CellIndex<NDIM>& cell_idx = *b;
                            CellIndex<NDIM> cell_nbr_idx = cell_idx;
                            cell_nbr_idx(axis) += n(axis);

                            SideIndex<NDIM> bdry_idx(
                                cell_idx, axis, upperlower ? SideIndex<NDIM>::Upper : SideIndex<NDIM>::Lower);
                            const double& dA = (*face_sc_data)(bdry_idx);

                            // Get the coordinate of the side index and r vector
                            side_coord.setZero();
                            getPhysicalCoordinateFromSideIndex(side_coord, level, patch, bdry_idx, axis);
                            r_vec = side_coord - fobj.r0;

                            IBTK::Vector3d pn = 0.5 * n * ((*p_data)(cell_idx) + (*p_data)(cell_nbr_idx));

                            // Pressure force := (n. -p I) * dA
                            trac += -pn * dA;

                            // Pressure torque := r x (-p n I) * dA
                            torque_trac += r_vec.cross(-pn) * dA;

                            // Momentum force := (n. -rho*(u)u) * dA
                            IBTK::Vector3d u = IBTK::Vector3d::Zero();
                            for (int d = 0; d < NDIM; ++d)
                            {
                                if (d == axis)
                                {
                                    u(d) = (*u_data)(bdry_idx);
                                }
                                else
                                {
                                    u(d) =
                                        0.25 *
                                        ((*u_data)(SideIndex<NDIM>(cell_idx, d, SideIndex<NDIM>::Lower)) +
                                         (*u_data)(SideIndex<NDIM>(cell_idx, d, SideIndex<NDIM>::Upper)) +
                                         (*u_data)(SideIndex<NDIM>(cell_nbr_idx, d, SideIndex<NDIM>::Lower)) +
                                         (*u_data)(SideIndex<NDIM>(cell_nbr_idx, d, SideIndex<NDIM>::Upper)));
                                }
                            }
                            trac += -d_rho * n.dot(u) * u * dA;

                            // Momentum torque := -(n. u) * rho * (r x u) * dA
                            torque_trac += -n.dot(u) * d_rho * r_vec.cross(u) * dA;

                            // Viscous traction force := n . mu(grad u + grad u ^ T) * dA
                            IBTK::Vector3d viscous_force = IBTK::Vector3d::Zero();
                            for (int d = 0; d < NDIM; ++d)
                            {
                                if (d == axis)
                                {
                                    viscous_force(axis) =
                                        n(axis) * (2.0 * d_mu) / (2.0 * patch_dx[axis]) *
                                        ((*u_data)(SideIndex<NDIM>(cell_nbr_idx, axis, upperlower ? SideIndex<NDIM>::Upper : SideIndex<NDIM>::Lower)) -
                                         (*u_data)(SideIndex<NDIM>(cell_idx, axis, upperlower ? SideIndex<NDIM>::Lower : SideIndex<NDIM>::Upper)));
                                }
                                else
                                {
                                    CellIndex<NDIM> offset(0);
                                    offset(d) = 1;

                                    viscous_force(d) =
                                        d_mu / (2.0 * patch_dx[d]) *
                                            ((*u_data)(
                                                 SideIndex<NDIM>(cell_idx + offset, axis, upperlower ? SideIndex<NDIM>::Upper : SideIndex<NDIM>::Lower)) -
                                             (*u_data)(
                                                 SideIndex<NDIM>(cell_idx - offset, axis, upperlower ? SideIndex<NDIM>::Upper : SideIndex<NDIM>::Lower)))

                                        +

                                        d_mu * n(axis) / (2.0 * patch_dx[axis]) *
                                            ((*u_data)(SideIndex<NDIM>(cell_nbr_idx, d, SideIndex<NDIM>::Lower)) +
                                             (*u_data)(
                                                 SideIndex<NDIM>(cell_nbr_idx + offset, d, SideIndex<NDIM>::Lower)) -
                                             (*u_data)(SideIndex<NDIM>(cell_idx, d, SideIndex<NDIM>::Lower)) -
                                             (*u_data)(SideIndex<NDIM>(cell_idx + offset, d, SideIndex<NDIM>::Lower))

                                                 );
                                }
                            }
                            IBTK::Vector3d n_dot_T = n(axis) * viscous_force;

                            trac += n_dot_T * dA;

                            // Viscous traction torque r x ( n . mu(grad u + grad u ^ T) * dA
                            torque_trac += r_vec.cross(n_dot_T) * dA;
                        }
                    }
                }
            }
        }
        SAMRAI_MPI::sumReduction(trac.data(), 3);
        SAMRAI_MPI::sumReduction(torque_trac.data(), 3);

        // Compute hydrodynamic force on the body : -integral_{box_new} (rho du/dt) + d/dt(rho u)_body + trac
        fobj.F_new = -(fobj.P_box_new - fobj.P_box_current) / dt + (fobj.P_new - fobj.P_current) / dt + trac;

        // Compute hydrodynamic torque on the body : -integral_{box_new} (rho d (r x u)/dt) + d/dt(rho r x u)_body +
        // torque_trac
        fobj.T_new = -(fobj.L_box_new - fobj.L_box_current) / dt + (fobj.L_new - fobj.L_current) / dt + torque_trac;
    }

    return;

} // computeHydrodynamicForce

void
IBHydrodynamicForceEvaluator::postprocessIntegrateData(double /*current_time*/, double new_time)
{
    for (std::map<int, IBHydrodynamicForceObject>::iterator it = d_hydro_objs.begin(); it != d_hydro_objs.end(); ++it)
    {
        IBHydrodynamicForceObject& force_obj = it->second;

        // Output drag and torque to stream
        if (SAMRAI_MPI::getRank() == 0)
        {
            *force_obj.drag_CV_stream << new_time << '\t' << force_obj.F_new(0) << '\t' << force_obj.F_new(1) << '\t'
                                      << force_obj.F_new(2) << std::endl;
            *force_obj.torque_CV_stream << new_time << '\t' << force_obj.T_new(0) << '\t' << force_obj.T_new(1) << '\t'
                                        << force_obj.T_new(2) << std::endl;
        }
        d_current_time = new_time;
        force_obj.box_u_current = force_obj.box_u_new;
        force_obj.box_X_lower_current = force_obj.box_X_lower_new;
        force_obj.box_X_upper_current = force_obj.box_X_upper_new;
	force_obj.box_vol_current = force_obj.box_vol_new;
        force_obj.F_current = force_obj.F_new;
        force_obj.T_current = force_obj.T_new;
        force_obj.P_current = force_obj.P_new;
        force_obj.L_current = force_obj.L_new;
        force_obj.P_box_current = force_obj.P_box_new;
        force_obj.L_box_current = force_obj.L_box_new;
    }

    return;

} // postprocessIntegrateData

void
IBHydrodynamicForceEvaluator::putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    for (std::map<int, IBHydrodynamicForceObject>::const_iterator it = d_hydro_objs.begin(); it != d_hydro_objs.end();
         ++it)
    {
        int strct_id = it->first;
        const IBHydrodynamicForceObject& force_obj = it->second;

        std::ostringstream F, T, P, L, P_box, L_box, X_lo, X_hi, r_or, vol_curr;
        F << "F_" << strct_id;
        T << "T_" << strct_id;
        P << "P_" << strct_id;
        L << "L_" << strct_id;
        P_box << "P_box_" << strct_id;
        L_box << "L_box_" << strct_id;
        X_lo << "X_lo_" << strct_id;
        X_hi << "X_hi_" << strct_id;
        r_or << "r_or_" << strct_id;
        vol_curr << "vol_curr_" << strct_id;

        db->putDoubleArray(F.str(), force_obj.F_current.data(), 3);
        db->putDoubleArray(T.str(), force_obj.T_current.data(), 3);
        db->putDoubleArray(P.str(), force_obj.P_current.data(), 3);
        db->putDoubleArray(L.str(), force_obj.L_current.data(), 3);
        db->putDoubleArray(P_box.str(), force_obj.P_box_current.data(), 3);
        db->putDoubleArray(L_box.str(), force_obj.L_box_current.data(), 3);
        db->putDoubleArray(X_lo.str(), force_obj.box_X_lower_current.data(), 3);
        db->putDoubleArray(X_hi.str(), force_obj.box_X_upper_current.data(), 3);
        db->putDoubleArray(r_or.str(), force_obj.r0.data(), 3);
        db->putDouble(vol_curr.str(), force_obj.box_vol_current);
    }

    return;

} // putToDatabase

void
IBAMR::IBHydrodynamicForceEvaluator::registerStructurePlotData(Pointer<VisItDataWriter<NDIM> > visit_data_writer,
                                                               Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                               int strct_id)

{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hydro_objs.find(strct_id) != d_hydro_objs.end());
#endif

    IBHydrodynamicForceEvaluator::IBHydrodynamicForceObject& fobj = d_hydro_objs[strct_id];

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Create a variable that is strct_id + 1 within control volume and 0 outside
    std::stringstream strct_id_stream;
    strct_id_stream << strct_id;
    std::string struct_no = strct_id_stream.str();
    Pointer<CellVariable<NDIM, double> > inside_strct_var = new CellVariable<NDIM, double>("box" + struct_no, 1);
    Pointer<VariableContext> ctx = var_db->getContext("box" + struct_no);
    fobj.inside_strct_idx = var_db->registerVariableAndContext(inside_strct_var, ctx, (IntVector<NDIM>)0);

    int coarsest_ln = 0;
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(fobj.inside_strct_idx)) level->allocatePatchData(fobj.inside_strct_idx);
    }

    // Register the indicator variable with the VisitDataWriter
    visit_data_writer->registerPlotQuantity("box" + struct_no, "SCALAR", fobj.inside_strct_idx);

    // Set the plot data for the initial box
    HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double> > hier_data_ops =
        hier_data_ops_manager->getOperationsDouble(inside_strct_var, patch_hierarchy, true);
    hier_data_ops->setToScalar(fobj.inside_strct_idx, 0.0, /*interior_only*/ true);

    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        Box<NDIM> integration_box(
            IndexUtilities::getCellIndex(fobj.box_X_lower_current.data(), level->getGridGeometry(), level->getRatio()),
            IndexUtilities::getCellIndex(fobj.box_X_upper_current.data(), level->getGridGeometry(), level->getRatio()));

        // Shorten the integration box so it only includes the control volume
        integration_box.upper() -= 1;

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const bool boxes_intersect = patch_box.intersects(integration_box);
            if (!boxes_intersect) continue;

            // Part of the box on this patch.
            Box<NDIM> trim_box = patch_box * integration_box;

            // Set plot variable to strct_id + 1
            Pointer<CellData<NDIM, double> > inside_strct_data = patch->getPatchData(fobj.inside_strct_idx);
            inside_strct_data->fillAll((double)strct_id + 1, trim_box);
        }
    }

    return;

} // registerStructurePlotData

void
IBAMR::IBHydrodynamicForceEvaluator::updateStructurePlotData(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                             int strct_id)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hydro_objs.find(strct_id) != d_hydro_objs.end());
#endif

    IBHydrodynamicForceEvaluator::IBHydrodynamicForceObject& fobj = d_hydro_objs[strct_id];

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Get the plot variable
    std::stringstream strct_id_stream;
    strct_id_stream << strct_id;
    std::string struct_no = strct_id_stream.str();
    Pointer<CellVariable<NDIM, double> > inside_strct_var = var_db->getVariable("box" + struct_no);

    int coarsest_ln = 0;
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(fobj.inside_strct_idx)) level->allocatePatchData(fobj.inside_strct_idx);
    }

    // Set the plot data for the new box to 0
    HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double> > hier_data_ops =
        hier_data_ops_manager->getOperationsDouble(inside_strct_var, patch_hierarchy, true);
    hier_data_ops->setToScalar(fobj.inside_strct_idx, 0.0, /*interior_only*/ true);

    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        Box<NDIM> integration_box(
            IndexUtilities::getCellIndex(fobj.box_X_lower_new.data(), level->getGridGeometry(), level->getRatio()),
            IndexUtilities::getCellIndex(fobj.box_X_upper_new.data(), level->getGridGeometry(), level->getRatio()));

        // Shorten the integration box so it only includes the control volume
        integration_box.upper() -= 1;

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const bool boxes_intersect = patch_box.intersects(integration_box);
            if (!boxes_intersect) continue;

            // Part of the box on this patch.
            Box<NDIM> trim_box = patch_box * integration_box;

            // Set plot variable to strct_id + 1
            Pointer<CellData<NDIM, double> > inside_strct_data = patch->getPatchData(fobj.inside_strct_idx);
            inside_strct_data->fillAll((double)strct_id + 1, trim_box);
        }
    }

    return;

} // updateStructurePlotData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBHydrodynamicForceEvaluator::resetFaceAreaWeight(Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
{
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_face_wgt_sc_idx)) level->allocatePatchData(d_face_wgt_sc_idx);
    }

    // Each cell's face weight is set to its face area, unless the cell is refined
    // on a finer level, in which case the weight is set to zero.  This ensures
    // that no part of the physical domain is counted twice when discrete norms
    // and surface integrals are calculated on the entire hierarchy.
    //
    // Away from coarse-fine interfaces and boundaries of the computational
    // domain, each cell face's weight is set to the face area associated with
    // the level of the patch hierarchy.
    ArrayDataBasicOps<NDIM, double> array_ops;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        BoxArray<NDIM> refined_region_boxes;
        if (ln < finest_ln)
        {
            Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
        }

        const IntVector<NDIM> max_gcw(1);
        const CoarseFineBoundary<NDIM> cf_bdry(*patch_hierarchy, ln, max_gcw);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

            const double* const dx = pgeom->getDx();
            const double cell_vol = dx[0] * dx[1]
#if (NDIM > 2)
                                    * dx[2]
#endif
                ;
            Pointer<SideData<NDIM, double> > face_wgt_sc_data = patch->getPatchData(d_face_wgt_sc_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                ArrayData<NDIM, double>& axis_data = face_wgt_sc_data->getArrayData(axis);
                axis_data.fill(cell_vol / dx[axis]);
            }

            // Zero-out weights within the refined region.
            if (ln < finest_ln)
            {
                const IntVector<NDIM>& periodic_shift = level->getGridGeometry()->getPeriodicShift(level->getRatio());
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        if (periodic_shift(axis) != 0)
                        {
                            for (int sgn = -1; sgn <= 1; sgn += 2)
                            {
                                IntVector<NDIM> periodic_offset = 0;
                                periodic_offset(axis) = sgn * periodic_shift(axis);
                                const Box<NDIM> refined_box =
                                    Box<NDIM>::shift(refined_region_boxes[i], periodic_offset);
                                const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                                if (!intersection.empty())
                                {
                                    face_wgt_sc_data->fillAll(0.0, intersection);
                                }
                            }
                        }
                    }
                    const Box<NDIM>& refined_box = refined_region_boxes[i];
                    const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        face_wgt_sc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }
    }
    return;

} // resetFaceAreaWeight

void
IBHydrodynamicForceEvaluator::resetFaceVolWeight(Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
{
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_vol_wgt_sc_idx))
        {
            level->allocatePatchData(d_vol_wgt_sc_idx);
        }
    }
    // Each cell's weight is set to its cell volume, unless the cell is refined
    // on a finer level, in which case the weight is set to zero.  This ensures
    // that no part of the physical domain is counted twice when discrete norms
    // and integrals are calculated on the entire hierarchy.
    //
    // Away from coarse-fine interfaces and boundaries of the computational
    // domain, each cell face's weight is set to the cell volume associated with
    // the level of the patch hierarchy.  Along coarse-fine interfaces or
    // physical boundaries, the weights associated with the cell faces are
    // modified so that the sum of the weights equals to volume of the
    // computational domain.
    ArrayDataBasicOps<NDIM, double> array_ops;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        BoxArray<NDIM> refined_region_boxes;
        if (ln < finest_ln)
        {
            Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
        }
        const IntVector<NDIM> max_gcw(1);
        const CoarseFineBoundary<NDIM> cf_bdry(*patch_hierarchy, ln, max_gcw);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double cell_vol = dx[0] * dx[1]
#if (NDIM > 2)
            * dx[2]
#endif
            ;
            Pointer<SideData<NDIM, double> > wgt_sc_data = patch->getPatchData(d_vol_wgt_sc_idx);
            wgt_sc_data->fillAll(cell_vol);
            // Rescale values along the edges of the patches.
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                Box<NDIM> side_lower_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                side_lower_box.upper()(axis) = side_lower_box.lower()(axis);
                array_ops.scale(wgt_sc_data->getArrayData(axis), 0.5, wgt_sc_data->getArrayData(axis), side_lower_box);
            }
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                Box<NDIM> side_upper_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                side_upper_box.lower()(axis) = side_upper_box.upper()(axis);
                array_ops.scale(wgt_sc_data->getArrayData(axis), 0.5, wgt_sc_data->getArrayData(axis), side_upper_box);
            }
            // Correct the values along coarse-fine interfaces.
            if (ln > coarsest_ln)
            {
                const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
                const int bdry_type = 1;
                const Array<BoundaryBox<NDIM> >& cf_bdry_boxes = cf_bdry.getBoundaries(p(), bdry_type);
                for (int k = 0; k < cf_bdry_boxes.getSize(); ++k)
                {
                    const Box<NDIM>& bdry_box = cf_bdry_boxes[k].getBox();
                    const unsigned int axis = cf_bdry_boxes[k].getLocationIndex() / 2;
                    const int lower_upper = cf_bdry_boxes[k].getLocationIndex() % 2;
                    if (!pgeom->getTouchesRegularBoundary(axis, lower_upper))
                    {
                        const double extra_vol = 0.5 * static_cast<double>(ratio(axis)) * cell_vol;
                        Box<NDIM> side_bdry_box = SideGeometry<NDIM>::toSideBox(bdry_box, axis);
                        array_ops.addScalar(
                                            wgt_sc_data->getArrayData(axis), wgt_sc_data->getArrayData(axis), extra_vol, side_bdry_box);
                    }
                }
            }
            // Zero-out weights within the refined region.
            if (ln < finest_ln)
            {
                const IntVector<NDIM>& periodic_shift = level->getGridGeometry()->getPeriodicShift(level->getRatio());
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        if (periodic_shift(axis) != 0)
                        {
                            for (int sgn = -1; sgn <= 1; sgn += 2)
                            {
                                IntVector<NDIM> periodic_offset = 0;
                                periodic_offset(axis) = sgn * periodic_shift(axis);
                                const Box<NDIM> refined_box =
                                Box<NDIM>::shift(refined_region_boxes[i], periodic_offset);
                                const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                                if (!intersection.empty())
                                {
                                    wgt_sc_data->fillAll(0.0, intersection);
                                }
                            }
                        }
                    }
                    const Box<NDIM>& refined_box = refined_region_boxes[i];
                    const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        wgt_sc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }
    }
    return;
} // resetFaceVolWeight

void
IBHydrodynamicForceEvaluator::fillPatchData(const int u_src_idx,
                                            const int p_src_idx,
                                            Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                            const std::vector<RobinBcCoefStrategy<NDIM>*>& u_src_bc_coef,
                                            RobinBcCoefStrategy<NDIM>* p_src_bc_coef,
                                            const double fill_time)
{
    // Whether or not to fill u and p
    const bool fill_velocity = (u_src_idx > 0);
    const bool fill_pressure = (p_src_idx > 0);

    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_u_idx) && fill_velocity) level->allocatePatchData(d_u_idx);
        if (!level->checkAllocated(d_p_idx) && fill_pressure) level->allocatePatchData(d_p_idx);
    }

    if (fill_velocity)
    {
        // Fill velocity data from integrator index.
        HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_sc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_u_var, patch_hierarchy, true);
        hier_sc_data_ops->copyData(d_u_idx, u_src_idx, true);

        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> transaction_comp(1);
        transaction_comp[0] = InterpolationTransactionComponent(d_u_idx,
                                                                u_src_idx,
                                                                /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ true,
                                                                /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                u_src_bc_coef,
                                                                Pointer<VariableFillPattern<NDIM> >(NULL));

        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(transaction_comp, patch_hierarchy);
        hier_bdry_fill->setHomogeneousBc(false);
        hier_bdry_fill->fillData(fill_time);
    }

    if (fill_pressure)
    {
        // Fill pressure data from integrator index
        HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();

        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_cc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_p_var, patch_hierarchy, true);
        hier_cc_data_ops->copyData(d_p_idx, p_src_idx, true);

        INSStaggeredPressureBcCoef* p_ins_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(p_src_bc_coef);
#if !defined(NDEBUG)
        TBOX_ASSERT(p_ins_bc_coef);
#endif
        p_ins_bc_coef->setTargetVelocityPatchDataIndex(d_u_idx);

        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> transaction_comp(1);
        transaction_comp[0] = InterpolationTransactionComponent(d_p_idx,
                                                                p_src_idx,
                                                                /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ true,
                                                                /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                p_ins_bc_coef,
                                                                Pointer<VariableFillPattern<NDIM> >(NULL));
        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(transaction_comp, patch_hierarchy);
        hier_bdry_fill->setHomogeneousBc(false);
        hier_bdry_fill->fillData(fill_time);
    }

    return;
} // fillPatchData

void
IBHydrodynamicForceEvaluator::getPhysicalCoordinateFromSideIndex(IBTK::Vector3d& side_coord,
                                                                 Pointer<PatchLevel<NDIM> > patch_level,
                                                                 Pointer<Patch<NDIM> > patch,
                                                                 const SideIndex<NDIM> side_idx,
                                                                 const int axis)
{
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
} // getPhysicalCoordinateFromSideIndex

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
