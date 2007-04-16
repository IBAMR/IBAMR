// Filename: IBTargetPointForceGen.C
// Last modified: <16.Apr.2007 05:39:35 boyce@bigboy.nyconnect.com>
// Created on 21 Mar 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBTargetPointForceGen.h"

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
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/LNodeIndexData.h>

// STOOLS INCLUDES
#include <stools/PETSC_SAMRAI_ERROR.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_lagrangian_force;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBTargetPointForceGen::IBTargetPointForceGen(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_lagrangian_force = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBTargetPointForceGen::computeLagrangianForce()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBTargetPointForceGen::initializeLevelData()");
        timers_need_init = false;
    }
    return;
}// IBTargetPointForceGen

IBTargetPointForceGen::~IBTargetPointForceGen()
{
    // intentionally blank
    return;
}// ~IBTargetPointForceGen

void
IBTargetPointForceGen::computeLagrangianForce(
    SAMRAI::tbox::Pointer<LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    t_compute_lagrangian_force->start();

    int ierr;

    // Extract the local arrays.
    Vec F_vec = F_data->getGlobalVec();
    double* F_arr;
    ierr = VecGetArray(F_vec, &F_arr);  PETSC_SAMRAI_ERROR(ierr);

    Vec X_vec = X_data->getGlobalVec();
    double* X_arr;
    ierr = VecGetArray(X_vec, &X_arr);  PETSC_SAMRAI_ERROR(ierr);

    // Get the grid geometry object and determine the extents of the physical
    // domain.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    if (!grid_geom->getDomainIsSingleBox())
    {
        TBOX_ERROR("IBTargetPointForceGen::computeLagrangianForce():\n"
                   << "  physical domain must be a single box.\n");
    }
    const double* const XLower = grid_geom->getXLower();
    const double* const XUpper = grid_geom->getXUpper();

    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->
        getLNodeIndexPatchDescriptorIndex();

    // Compute the penalty force associated with the Lagrangian target points.
    static double max_displacement = 0.0;
    double max_config_displacement = 0.0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data =
            patch->getPatchData(lag_node_index_idx);

        for (LNodeIndexData::Iterator it(*idx_data); it; it++)
        {
            if (patch_box.contains(it.getIndex()))
            {
                const LNodeIndexSet& node_set = *it;
                for (LNodeIndexSet::const_iterator n = node_set.begin();
                     n != node_set.end(); ++n)
                {
                    const LNodeIndexSet::value_type& node_idx = *n;
                    const std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data =
                        node_idx->getStashData();
                    for (unsigned l = 0; l < stash_data.size(); ++l)
                    {
                        SAMRAI::tbox::Pointer<IBTargetPointForceSpec> force_spec = stash_data[l];
                        if (!force_spec.isNull())
                        {
#ifdef DEBUG_CHECK_ASSERTIONS
                            const int& mastr_idx = node_idx->getLagrangianIndex();
                            assert(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                            const double& kappa_target = force_spec->getStiffness();
                            if (!SAMRAI::tbox::Utilities::deq(kappa_target,0.0))
                            {
                                const int& petsc_idx = node_idx->getLocalPETScIndex();
                                const double* const X = &X_arr[NDIM*petsc_idx];
                                const std::vector<double>& X_target = force_spec->getTargetPointPosition();

                                double* const F = &F_arr[NDIM*petsc_idx];
                                double displacement = 0.0;
                                for (int d = 0; d < NDIM; ++d)
                                {
                                    const double shift =
                                        (X[d] < XLower[d]
                                         ? XUpper[d] - XLower[d]
                                         : (X[d] > XUpper[d]
                                            ? XLower[d] - XUpper[d]
                                            : 0.0));
#ifdef DEBUG_CHECK_ASSERTIONS
                                    assert(X[d]+shift >= XLower[d]);
                                    assert(X[d]+shift <= XUpper[d]);
#endif
                                    F[d] += kappa_target*(X_target[d] - (X[d]+shift));
                                    displacement += pow(X_target[d] - (X[d]+shift),2.0);
                                }
                                displacement = sqrt(displacement);
                                if (displacement > max_config_displacement)
                                {
                                    max_config_displacement = displacement;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    ierr = VecRestoreArray(F_vec, &F_arr);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecRestoreArray(X_vec, &X_arr);  PETSC_SAMRAI_ERROR(ierr);

    max_config_displacement = SAMRAI::tbox::MPI::maxReduction(max_config_displacement);
    if (max_config_displacement > max_displacement)
    {
        max_displacement = max_config_displacement;
    }

    if (!SAMRAI::tbox::Utilities::deq(max_config_displacement,0.0))
    {
        SAMRAI::tbox::plog << "IBTargetPointForceGen::computeLagrangianForce():" << endl;
        SAMRAI::tbox::plog << "  maximum target point displacement [present configuration] = " << max_config_displacement << endl;
        SAMRAI::tbox::plog << "  maximum target point displacement [entire simulation] = " << max_displacement << endl;
    }

    t_compute_lagrangian_force->stop();
    return;
}// computeLagrangianForce

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBTargetPointForceGen::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    if (!db.isNull())
    {
        // intentionally blank
    }
    return;
}// getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBTargetPointForceGen>;

//////////////////////////////////////////////////////////////////////////////
