// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2020 by the IBAMR developers
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

#include "ibamr/NonbondedForceEvaluator.h"

#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LIndexSetData.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeSet.h"
#include "ibtk/LNodeSetData.h"
#include "ibtk/LSetData.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "Index.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Database.h"
#include "tbox/Utilities.h"

#include "petscvec.h"
#include <petscsys.h>

#include <assert.h>

#include <algorithm>
#include <cmath>
#include <string>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

NonbondedForceEvaluator::NonbondedForceEvaluator(Pointer<Database> input_db,
                                                 Pointer<CartesianGridGeometry<NDIM> > grid_geometry)
    : d_grid_geometry(grid_geometry)
{
    // get interaction radius
    if (input_db->keyExists("interaction_radius"))
    {
        d_interaction_radius = input_db->getDouble("interaction_radius");
    }
    else
    {
        TBOX_ERROR("Must specify interaction_radius for NonbondedForceEvaluator.");
    }

    // get regrid_alpha
    if (input_db->keyExists("regrid_alpha"))
    {
        d_regrid_alpha = input_db->getDouble("regrid_alpha");
    }
    else
    {
        TBOX_ERROR("Must specify regrid_alpha for NonbondedForceEvaluator.");
    }

    // this will only work if the domain is a single box.
    assert(d_grid_geometry->getDomainIsSingleBox());

    // get parameters for force function
    d_parameters = input_db->getDoubleArray("parameters");
}

void
NonbondedForceEvaluator::evaluateForces(int mstr_petsc_idx,
                                        int search_petsc_idx,
                                        Pointer<LData> X_data,
                                        std::vector<int> cell_offset,
                                        Pointer<LData> F_data)
{
    //   Function to add nonbonded forces from the interaction between the nodes at
    //   mstr_petsc_idx and search_petsc_idx.
    //
    //   inputs: mstr_petsc_idx:    integer petsc index of the master node
    //           search_petsc_idx:  integer petsc index of the search node
    //           X_data:            pointer to LData of location values for the particles
    //           cell_offset:       vector of periodic offsets of the search particle,
    //                    integer value in units of the entire domain.  E.G. a particle
    //                    image at position (120, -20, 10)  in a length 100^3 domain
    //                    would have cell_offset = (1, -1, 0).
    //
    //   outputs:
    //          F_data - pointer to LData object containing forces on particles.  Will
    //                   be added to by this function.
    //
    //////////////////////////////////////////////////////////////////////////////////

    // get vectors of data
    // How costly is this?  TODO: Check, maybe just do this once in the loop.
    PetscScalar* position;
    VecGetArray(X_data->getVec(), &position);
    PetscScalar* force;
    VecGetArray(F_data->getVec(), &force);

    // get domain bounds
    const double* x_lower = d_grid_geometry->getXLower();
    const double* x_upper = d_grid_geometry->getXUpper();

    double R = 0.0; // distance between particles
    double D[NDIM]; // vector connecting particles.

    for (int k = 0; k < NDIM; ++k)
    {
        D[k] = (position[mstr_petsc_idx * NDIM + k] - position[search_petsc_idx * NDIM + k] -
                cell_offset[k] * (x_upper[k] - x_lower[k]));
        R += std::pow(D[k], 2);
    }
    R = std::sqrt(R);

    double nonbdd_force[NDIM];
    (d_force_fcn_ptr)(D, d_parameters, nonbdd_force);
    for (int k = 0; k < NDIM; ++k)
    {
        force[mstr_petsc_idx * NDIM + k] += nonbdd_force[k];
        force[search_petsc_idx * NDIM + k] += -1.0 * nonbdd_force[k];
    }
    VecRestoreArray(F_data->getVec(), &force);
    return;
} // evaluateForces

void
NonbondedForceEvaluator::computeLagrangianForce(Pointer<LData> F_data,
                                                Pointer<LData> X_data,
                                                Pointer<LData> /*U_data*/,
                                                const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                const int level_number,
                                                const double /*data_time*/,
                                                LDataManager* const l_data_manager)
{
    // Get grid geometry and relevant lower and upper limits.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    if (!grid_geom->getDomainIsSingleBox()) TBOX_ERROR("physical domain must be a single box...\n");

    // These will only work if the domain is a single box.
    assert(grid_geom->getDomainIsSingleBox());
    const double* const x_lower = grid_geom->getXLower();
    const double* const x_upper = grid_geom->getXUpper();

    // we will grow the search box by interaction_radius + 2.0*regrid_alpha
    IntVector<NDIM> grow_amount(static_cast<int>(ceil(d_interaction_radius + 2.0 * d_regrid_alpha)));
    const int lag_node_idx_current_idx = l_data_manager->getLNodePatchDescriptorIndex();

    // iterate through levels.
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<LNodeSetData> current_idx_data = patch->getPatchData(lag_node_idx_current_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();

        std::vector<int> cell_offset(NDIM);
        // Loop through cells in this processors patch. For each iteration, this
        // is the "master" cell. Iterate through particles in the box, and add
        // springs.
        for (LNodeSetData::CellIterator cit(patch_box); cit; cit++)
        {
            // get list of particles in this cell
            const hier::Index<NDIM>& first_cell_idx = *cit;
            LNodeSet* const mstr_node_set = current_idx_data->getItem(first_cell_idx);
            if (mstr_node_set)
            {
                Box<NDIM> search_box(first_cell_idx, first_cell_idx);
                // loop over neighboring cells, up to interaction_radius +
                // 2*regrid_alpha away.
                for (LNodeSetData::CellIterator scit(Box<NDIM>::grow(search_box, grow_amount)); scit; scit++)
                {
                    // loop over particles in the neighbor cell, adding up
                    // forces onto the particle in the "master" cell At this
                    // point we know both cells, need to figure out periodic
                    // additions.
                    const hier::Index<NDIM>& search_cell_idx = *scit;

                    // search across periodic boundaries.
                    for (int k = 0; k < NDIM; ++k)
                    {
                        // Difference between lower boundary and this search cell.
                        double absolute_diff = search_cell_idx[k] * patch_dx[k];
                        // Periodic offset of this cell.
                        cell_offset[k] = floor(absolute_diff / (x_upper[k] - x_lower[k]));
                    }
                    LNodeSet* search_node_set = current_idx_data->getItem(search_cell_idx);
                    if (search_node_set)
                    {
                        // we have a set of nodes in the first cell and the search cell,
                        // add up forces
                        // and accumulate for the first cell.
                        for (const auto& mstr_node_idx : *mstr_node_set)
                        {
                            // master nodes
                            const int mstr_lag_idx = mstr_node_idx->getLagrangianIndex();
                            const int mstr_petsc_idx = mstr_node_idx->getLocalPETScIndex();

                            for (const auto& search_node_idx : *search_node_set)
                            {
                                const int search_lag_idx = search_node_idx->getLagrangianIndex();
                                const int search_petsc_idx = search_node_idx->getLocalPETScIndex();
                                if (mstr_lag_idx < search_lag_idx)
                                {
                                    // apply forces with the force evaluator
                                    evaluateForces(mstr_petsc_idx, search_petsc_idx, X_data, cell_offset, F_data);
                                }
                            } // search node index
                        }     // mstr node index
                    }         // if(search node set)
                }             // search cell loop
            }                 // if mastr_node_idx
        }                     // first cell
    }                         // patches
    return;
} // computeLagrangianForce

void
NonbondedForceEvaluator::registerForceFcnPtr(NonBddForceFcnPtr force_fcn_ptr)
{
    // set the nonbonded force function pointer to the given force function pointer
    d_force_fcn_ptr = force_fcn_ptr;
    return;
} // registerForceFcnPtr

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
