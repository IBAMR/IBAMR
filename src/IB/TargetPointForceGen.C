// Filename: TargetPointForceGen.C
// Created on 23 Oct 2006 by Boyce Griffith (boyce@bigboy.nyconnect.com)
// Last modified: <03.Nov.2006 12:40:06 boyce@boyce-griffiths-powerbook-g4-15.local>

#include "TargetPointForceGen.h"

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
#include <ibamr/LNodeIndexData.h>
#include <ibamr/TargetPointForceSpec.h>

// STOOLS INCLUDES
#include <stools/PETSC_SAMRAI_ERROR.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/MPI.h>
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

TargetPointForceGen::TargetPointForceGen()
{
    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_lagrangian_force = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::TargetPointForceGen::computeLagrangianForce()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::TargetPointForceGen::initializeLevelData()");
        timers_need_init = false;
    }
    return;
}// TargetPointForceGen

TargetPointForceGen::~TargetPointForceGen()
{
    // intentionally blank
    return;
}// ~TargetPointForceGen

void
TargetPointForceGen::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    const LDataManager* const lag_manager)
{
    t_initialize_level_data->start();

    // intentionally blank

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
TargetPointForceGen::computeLagrangianForce(
    SAMRAI::tbox::Pointer<LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    const LDataManager* const lag_manager)
{
    t_compute_lagrangian_force->start();

    // Extract the local arrays.
    int ierr;

    Vec F_vec = F_data->getGlobalVec();
    double* F_arr;
    ierr = VecGetArray(F_vec, &F_arr);  PETSC_SAMRAI_ERROR(ierr);

    Vec X_vec = X_data->getGlobalVec();
    double* X_arr;
    ierr = VecGetArray(X_vec, &X_arr);  PETSC_SAMRAI_ERROR(ierr);

    // Get the grid geometry object and determine the extents of the
    // physical domain.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    if (!grid_geom->getDomainIsSingleBox())
    {
        TBOX_ERROR("physical domain must be a single box...\n");
    }
    const double* const XLower = grid_geom->getXLower();
    const double* const XUpper = grid_geom->getXUpper();

    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->
        getLNodeIndexPatchDescriptorIndex();

    // Compute the penalty force associated with the Lagrangian target
    // points.
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
                    SAMRAI::tbox::Pointer<TargetPointForceSpec> force_spec =
                        node_idx->getStashData()[0];
                    if (!force_spec.isNull())
                    {
                        const int& petsc_idx = node_idx->getLocalPETScIndex();
                        const double* const X = &X_arr[NDIM*petsc_idx];
                        const vector<double>& X_target = force_spec->getPosition();
                        const double& kappa = force_spec->getStiffness();

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
                            F[d] = kappa*(X_target[d] - (X[d]+shift));
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

    ierr = VecRestoreArray(F_vec, &F_arr);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecRestoreArray(X_vec, &X_arr);  PETSC_SAMRAI_ERROR(ierr);

    max_config_displacement = SAMRAI::tbox::MPI::maxReduction(max_config_displacement);
    if (max_config_displacement > max_displacement)
    {
        max_displacement = max_config_displacement;
    }
    SAMRAI::tbox::plog << "TargetPointForceGen::computeLagrangianForce():" << endl;
    SAMRAI::tbox::plog << "  maximum target point displacement [present configuration] = " << max_config_displacement << endl;
    SAMRAI::tbox::plog << "  maximum target point displacement [entire simulation] = " << max_displacement << endl;

    t_compute_lagrangian_force->stop();
    return;
}// computeLagrangianForce

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::TargetPointForceGen>;

//////////////////////////////////////////////////////////////////////////////
