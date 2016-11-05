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

#include "ArrayDataBasicOps.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CoarseFineBoundary.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "SideData.h"
#include "SideIndex.h"
#include "boost/array.hpp"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "ibamr/IBHydrodynamicForceEvaluator.h"
#include "ibamr/namespaces.h"
#include "ibtk/IndexUtilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBHydrodynamicForceEvaluator::IBHydrodynamicForceEvaluator(const std::string& object_name,
                                                           double rho,
                                                           double mu,
                                                           bool register_for_restart)
{
    d_object_name = object_name;
    d_rho = rho;
    d_mu = mu;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
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
    var_db->removePatchDataIndex(d_face_wgt_sc_idx);

    return;
} // ~IBHydrodynamicForceEvaluator

void
IBHydrodynamicForceEvaluator::registerStructure(int strct_id,
                                                int strct_ln,
                                                const Eigen::Vector3d& box_vel,
                                                Eigen::Vector3d& box_X_lower,
                                                Eigen::Vector3d& box_X_upper,
						Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hydro_objs.find(strct_id) == d_hydro_objs.end());
#endif

    IBHydrodynamicForceObject force_obj;
    force_obj.strct_id = strct_id;
    force_obj.strct_ln = strct_ln;

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
	const Eigen::Vector3d box_X_lower_old = box_X_lower;
	const Eigen::Vector3d box_X_upper_old = box_X_upper;
	
	for (int d = 0; d < NDIM; ++d)
	{
	    double num_cells_lower = (box_X_lower[d]-grid_X_lower[d])/dx_coarsest[d];
	    double num_cells_upper = (box_X_upper[d]-grid_X_lower[d])/dx_coarsest[d];
	    
	    if (num_cells_lower != (int) num_cells_lower)
	    {
		TBOX_WARNING("Lower side of integration box is not aligned with sides on coarsest level in dimension "
			      << d << ". Modifying coordinate to nearest box side" << std::endl);
		const int N = floor(num_cells_lower);
		box_X_lower[d] = grid_X_lower[d] + dx_coarsest[d] * N;
		modified_box = true;
	    }
	    
	    if (num_cells_upper != (int) num_cells_upper)
	    {
		TBOX_WARNING("Upper side of integration box is not aligned with sides on coarsest level in dimension "
			      << d << ". Modifying coordinate to nearest box side" << std::endl);
		const int N = ceil(num_cells_upper);
		box_X_upper[d] = grid_X_lower[d] + dx_coarsest[d] * N;
		modified_box = true;
	    }
	}
	
	if (modified_box)
	{
	    pout << "IBHydrodynamicForceEvaluator::registerStructure: integration box modified from\n"
	         <<    "[" << box_X_lower_old[0] << ", " << box_X_lower_old[1] << ", " << box_X_lower_old[2] << "]"
		 << " x [" << box_X_upper_old[0] << ", " << box_X_upper_old[1] << ", " << box_X_upper_old[2] << "]\n"
		 << "to\n"
		 <<    "[" << box_X_lower[0] << ", " << box_X_lower[1] << ", " << box_X_lower[2] << "]"
		 << " x [" << box_X_upper[0] << ", " << box_X_upper[1] << ", " << box_X_upper[2] << "]"
		 << std::endl;
	}
      
        force_obj.box_u_current = box_vel;
        force_obj.box_X_lower_current = box_X_lower;
        force_obj.box_X_upper_current = box_X_upper;
        force_obj.F_current.setZero();
        force_obj.T_current.setZero();
        force_obj.P_current.setZero();
        force_obj.L_current.setZero();
        force_obj.P_box_current.setZero();
        force_obj.L_box_current.setZero();
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

        std::ostringstream F, T, P, L, P_box, L_box, X_lo, X_hi;
        F << "F_" << strct_id;
        T << "T_" << strct_id;
        P << "P_" << strct_id;
        L << "L_" << strct_id;
        P_box << "P_box_" << strct_id;
        L_box << "L_box_" << strct_id;
        X_lo << "X_lo_" << strct_id;
        X_hi << "X_hi_" << strct_id;
        db->getDoubleArray(F.str(), force_obj.F_current.data(), 3);
        db->getDoubleArray(T.str(), force_obj.T_current.data(), 3);
        db->getDoubleArray(P.str(), force_obj.P_current.data(), 3);
        db->getDoubleArray(L.str(), force_obj.L_current.data(), 3);
        db->getDoubleArray(P_box.str(), force_obj.P_box_current.data(), 3);
        db->getDoubleArray(L_box.str(), force_obj.L_box_current.data(), 3);
        db->getDoubleArray(X_lo.str(), force_obj.box_X_lower_current.data(), 3);
        db->getDoubleArray(X_hi.str(), force_obj.box_X_upper_current.data(), 3);
    }

    d_hydro_objs[strct_id] = force_obj;
    return;

} // registerStructure

void
IBHydrodynamicForceEvaluator::updateStructureDomain(int strct_id,
                                                    int /*strct_ln*/,
                                                    double current_time,
                                                    double new_time,
                                                    const Eigen::Vector3d& box_vel_new,
                                                    const Eigen::Vector3d& P_strct_new,
                                                    const Eigen::Vector3d& L_strct_new)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hydro_objs.find(strct_id) != d_hydro_objs.end());
#endif

    IBHydrodynamicForceEvaluator::IBHydrodynamicForceObject& force_obj = d_hydro_objs[strct_id];
    const double dt = new_time - current_time;
    force_obj.box_u_new = box_vel_new;
    force_obj.box_X_lower_new = force_obj.box_X_lower_current + box_vel_new * dt;
    force_obj.box_X_upper_new = force_obj.box_X_upper_current + box_vel_new * dt;
    force_obj.P_new = P_strct_new;
    force_obj.L_new = L_strct_new;

    return;

} // updateStructureDomain

void
IBHydrodynamicForceEvaluator::preprocessIntegrateData(double /*current_time*/, double /*new_time*/)
{
    pout << "WARNING:: IBHydrodynamicForceEvaluator::preprocessIntegrateData() not implemented.\n";

    return;
} // preprocessIntegrateData

const IBHydrodynamicForceEvaluator::IBHydrodynamicForceObject&
IBHydrodynamicForceEvaluator::getHydrodynamicForceObject(int strct_id, int /*strct_ln*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hydro_objs.find(strct_id) != d_hydro_objs.end());
#endif

    const IBHydrodynamicForceEvaluator::IBHydrodynamicForceObject& force_obj = d_hydro_objs[strct_id];
    return force_obj;

} // getHydrodynamicForceObject

void
IBHydrodynamicForceEvaluator::computeHydrodynamicForce(int u_idx,
                                                       int p_idx,
                                                       int /*f_idx*/,
                                                       int /*vol_sc_idx*/,
                                                       int /*vol_cc_idx*/,
                                                       Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                       int coarsest_ln,
                                                       int finest_ln,
                                                       double current_time,
                                                       double new_time)
{
    resetFaceAreaWeight(patch_hierarchy);
    resetFaceVolWeight(patch_hierarchy);
    const double dt = new_time - current_time;

    // Get the grid spacing at the coarsest level 
    // and lower left corner of the computational domain (required for computing AMR offset)
    Pointer<PatchLevel<NDIM> > coarsest_level = patch_hierarchy->getPatchLevel(coarsest_ln);
    const Pointer<CartesianGridGeometry<NDIM> > coarsest_grid_geom = coarsest_level->getGridGeometry();
    const double* const dx_coarsest           = coarsest_grid_geom->getDx();
    const double* const X_lower_left_coarsest = coarsest_grid_geom->getXLower();
    
    for (std::map<int, IBHydrodynamicForceObject>::iterator it = d_hydro_objs.begin(); it != d_hydro_objs.end(); ++it)
    {
        IBHydrodynamicForceObject& fobj = it->second;

        // Compute the momentum integral:= (rho * u * dv)
        fobj.P_box_new.setZero();
        for (int ln = finest_ln; ln >= coarsest_ln; --ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            Box<NDIM> integration_box(
                IndexUtilities::getCellIndex(fobj.box_X_lower_new.data(), level->getGridGeometry(), level->getRatio()),
                IndexUtilities::getCellIndex(fobj.box_X_upper_new.data(), level->getGridGeometry(), level->getRatio()));
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
		const Index<NDIM>& patch_lower_idx = patch_box.lower();
		const double* patch_X_lower = patch_geom->getXLower();
		const double* patch_dx = patch_geom->getDx();
                const bool boxes_intersect = patch_box.intersects(integration_box);
                if (!boxes_intersect) continue;

                // Part of the box on this patch.
                Box<NDIM> trim_box = patch_box * integration_box;
		
                // Loop over the box and compute momentum.
                Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
                Pointer<SideData<NDIM, double> > vol_sc_data = patch->getPatchData(d_vol_wgt_sc_idx);
		
		for (int axis = 0; axis < NDIM; ++axis)
		{
		    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(trim_box, axis)); b; b++)
		    {
			const CellIndex<NDIM>& cell_idx = *b;
		    
// 			const bool inside_domain = insideBoxDomain(cell_idx, axis, 
// 								   patch_lower_idx, patch_X_lower, patch_dx,
// 								   fobj.box_X_lower_new, fobj.box_X_upper_new);
			
			// Check whether the lower sides of the grid cell is within the integration box
// 			if (!withinIntegrationBox(cell_idx, fobj.box_X_lower_new, fobj.box_X_upper_new, patch_geom, dx_coarsest, X_lower_left_coarsest))
// 			    continue;
			
			// If a cell is outside of the integration box domain, then the grid cells are not aligned
			// Only need to check this once for each cell.
// 			if (!inside_domain)
// 			    TBOX_ERROR("Non aligned box with level " << ln << "\n");
                        
                        const SideIndex<NDIM> side_idx(cell_idx, axis, SideIndex<NDIM>::Lower);
                        const double& u_axis = (*u_data)(side_idx);
                        const double& vol = (*vol_sc_data)(side_idx);
                        fobj.P_box_new(axis) += d_rho * vol * u_axis;
                    }
                }
            }
        }
        SAMRAI_MPI::sumReduction(fobj.P_box_new.data(), 3);

        // Compute surface integral term.
        Eigen::Vector3d trac;
        trac.setZero();
        for (int ln = finest_ln; ln >= coarsest_ln; --ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            Box<NDIM> integration_box(
                IndexUtilities::getCellIndex(fobj.box_X_lower_new.data(), level->getGridGeometry(), level->getRatio()),
                IndexUtilities::getCellIndex(fobj.box_X_upper_new.data(), level->getGridGeometry(), level->getRatio()));
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
                Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(p_idx);
                Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
                Pointer<SideData<NDIM, double> > face_sc_data = patch->getPatchData(d_face_wgt_sc_idx);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (int upperlower = 0; upperlower <= 1; ++upperlower)
                    {
                        const Box<NDIM>& side_box = bdry_boxes[axis][upperlower];
                        if (!patch_box.intersects(side_box)) continue;

                        Box<NDIM> trim_box = patch_box * side_box;
                        Eigen::Vector3d n = Eigen::Vector3d::Zero();
                        n(axis) = upperlower ? 1 : -1;
                        for (Box<NDIM>::Iterator b(trim_box); b; b++)
                        {
                            const CellIndex<NDIM>& cell_idx = *b;
                            CellIndex<NDIM> cell_nbr_idx = cell_idx;
                            cell_nbr_idx(axis) += n(axis);

                            SideIndex<NDIM> bdry_idx(
                                cell_idx, axis, upperlower ? SideIndex<NDIM>::Upper : SideIndex<NDIM>::Lower);
                            const double& dA = (*face_sc_data)(bdry_idx);

                            // Pressure force := (n. -p I) * dA
                            trac += -0.5 * n * ((*p_data)(cell_idx) + (*p_data)(cell_nbr_idx)) * dA;

                            // Momentum force := (n. -rho*(u - u_b)u) * dA
                            Eigen::Vector3d u = Eigen::Vector3d::Zero();
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
                            trac += -d_rho * n.dot(u - fobj.box_u_new) * u * dA;

                            // Viscous traction force := n . mu(grad u + grad u ^ T) * ds
                            Eigen::Vector3d viscous_force = Eigen::Vector3d::Zero();
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
                            trac += n(axis) * viscous_force * dA;
                        }
                    }
                }
            }
        }
        SAMRAI_MPI::sumReduction(trac.data(), 3);

        // Compute hydrodynamic force on the body : -d/dt(rho u)_box + d/dt(rho u)_body + trac
        fobj.F_new = (fobj.P_box_current - fobj.P_box_new + fobj.P_new - fobj.P_current) / dt + trac;
    }

    return;

} // computeHydrodynamicForce

void
IBHydrodynamicForceEvaluator::postprocessIntegrateData(double /*current_time*/, double /*new_time*/)
{
    for (std::map<int, IBHydrodynamicForceObject>::iterator it = d_hydro_objs.begin(); it != d_hydro_objs.end(); ++it)
    {
        IBHydrodynamicForceObject& force_obj = it->second;

        force_obj.box_u_current = force_obj.box_u_new;
        force_obj.box_X_lower_current = force_obj.box_X_lower_new;
        force_obj.box_X_upper_current = force_obj.box_X_upper_new;
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

        std::ostringstream F, T, P, L, P_box, L_box, X_lo, X_hi;
        F << "F_" << strct_id;
        T << "T_" << strct_id;
        P << "P_" << strct_id;
        L << "L_" << strct_id;
        P_box << "P_box_" << strct_id;
        L_box << "L_box_" << strct_id;
        X_lo << "X_lo_" << strct_id;
        X_hi << "X_hi_" << strct_id;

        db->putDoubleArray(F.str(), force_obj.F_current.data(), 3);
        db->putDoubleArray(T.str(), force_obj.T_current.data(), 3);
        db->putDoubleArray(P.str(), force_obj.P_current.data(), 3);
        db->putDoubleArray(L.str(), force_obj.L_current.data(), 3);
        db->putDoubleArray(P_box.str(), force_obj.P_box_current.data(), 3);
        db->putDoubleArray(L_box.str(), force_obj.L_box_current.data(), 3);
        db->putDoubleArray(X_lo.str(), force_obj.box_X_lower_current.data(), 3);
        db->putDoubleArray(X_hi.str(), force_obj.box_X_upper_current.data(), 3);
    }

    return;

} // putToDatabase

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
    // on a finer level, in which case the weight is set to zero.  This insures
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
    // on a finer level, in which case the weight is set to zero.  This insures
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

bool
IBHydrodynamicForceEvaluator::withinIntegrationBox(const CellIndex<NDIM>& cc_idx,
						   Eigen::Vector3d& X_lower,
						   Eigen::Vector3d& X_upper,
						   const Pointer<CartesianPatchGeometry<NDIM> > patch_geom,
						   const double* const dx_coarsest,
						   const double* const X_lower_left_coarsest)
{
    // Returns true if the grid cell in question is within the integration box.
    // This ensures that if the indices returned by getCellIndex correspond
    // to physical locations outside of the integration box (defined by X_lower
    // and X_upper), then they are not counted towards the drag computation
  
    // You only need to check whether the lower side of each box is within the 
    // integration box boundary. Computation of the physical location is done
    // by counting the number of grid cells away a side is from the left corner
    // of the given patch
    
    // Get grid spacing and coordinate of bottom left corner of the patch
    const double* const patch_dx           = patch_geom->getDx();
    const double* const patch_X_lower_left = patch_geom->getXLower();
    
    // Get the number of grid cells to between lower left patch corner and 
    // lower left corner of the computational domain
    int N_coarsest[NDIM];
    for (int axis = 0; axis < NDIM; ++axis)
    {
	N_coarsest[axis] = (patch_X_lower_left[axis] - X_lower_left_coarsest[axis])/dx_coarsest[axis];
    }
    
    const int* const ratio_to_coarsest = patch_geom->getRatio();
    
    for (int axis = 0; axis < NDIM; ++axis)
    {
	const SideIndex<NDIM> lower_sc_idx(cc_idx, axis, SideIndex<NDIM>::Lower);
	
	// Get number of cells between lower left patch corner and sides of the grid cell
	const int AMR_offset = N_coarsest[axis] > 1 ? N_coarsest[axis] * ratio_to_coarsest[axis] : 0;
	const int N_lower = lower_sc_idx[axis] - AMR_offset;
	
	// Calculate physical location of the side indexes
	double lower_phys_loc = patch_X_lower_left[axis] + patch_dx[axis] * N_lower;

	// Only need to check the lower position of the box
	if (lower_phys_loc < X_lower(axis) || lower_phys_loc > X_upper(axis))
	{
	    return false;
	}
    }
    return true;
} // withinIntegrationBox

bool
IBHydrodynamicForceEvaluator::insideBoxDomain(const SAMRAI::pdat::CellIndex<NDIM>& cell_idx,
					      const int axis,
					      const Index<NDIM>& patch_lower_idx,
					      const double* patch_X_lower,
					      const double* patch_dx,
					      const Eigen::Vector3d& X_lower, 
					      const Eigen::Vector3d& X_upper)
{
    Eigen::Vector3d X_loc = Eigen::Vector3d::Zero();
    bool inside;
    for (int d = 0; d < NDIM; ++d)
    {
	
	if (d != axis)
	{
	    X_loc(d) = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(cell_idx(d) - patch_lower_idx(d)) + 0.5);
	}
	else
	{
	    X_loc(d) = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(cell_idx(d) - patch_lower_idx(d)));
	}
	
	
	if (X_loc[0] < X_lower[0] || X_loc[1] < X_lower[1] || X_loc[2] < X_lower[2])
	{
	    inside =  false;
	}
	else if (X_loc[0] > X_upper[0] || X_loc[1] > X_upper[1] || X_loc[2] > X_upper[2])
	{
	    inside =  false;
	}
	else
	{
	    inside =  true;
#if !defined (NDEBUG)  // Maybe this too much...
	    TBOX_ASSERT(X_loc[0] >= X_lower[0] && X_loc[0] <= X_upper[0] 
		    &&  X_loc[1] >= X_lower[1] && X_loc[1] <= X_upper[1] 
		    &&  X_loc[2] >= X_lower[2] && X_loc[2] <= X_upper[2]);
	
#endif
	}
  }

	  return inside;
} // insideBoxDomain


/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
