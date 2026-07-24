// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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
#include <ibamr/VelocityCartGridFunction.h>
#include <ibamr/ibamr_utilities.h>

#include <CartesianGridGeometry.h>
#include <FaceData.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>

#include <cmath>

#include <ibamr/namespaces.h>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_SIDE_TO_FACE_FC IBAMR_FC_FUNC_(navier_stokes_side_to_face2d, NAVIER_STOKES_SIDE_TO_FACE2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_SIDE_TO_FACE_FC IBAMR_FC_FUNC_(navier_stokes_side_to_face3d, NAVIER_STOKES_SIDE_TO_FACE3D)
#endif

extern "C"
{
    void NAVIER_STOKES_SIDE_TO_FACE_FC(
#if (NDIM == 2)
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        const double*,
        const int&,
        double*,
        double*,
        const int&
#endif
#if (NDIM == 3)
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        const double*,
        const double*,
        const int&,
        double*,
        double*,
        double*,
        const int&
#endif
    );
}

// Copy data from a side-centered variable to a face-centered variable.
void
copy_side_to_face(const int U_fc_idx, const int U_sc_idx, Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const hier::Index<NDIM>& ilower = patch->getBox().lower();
            const hier::Index<NDIM>& iupper = patch->getBox().upper();
            Pointer<SideData<NDIM, double>> U_sc_data = patch->getPatchData(U_sc_idx);
            Pointer<FaceData<NDIM, double>> U_fc_data = patch->getPatchData(U_fc_idx);

#if !defined(NDEBUG)
            TBOX_ASSERT(U_sc_data->getGhostCellWidth().min() == U_sc_data->getGhostCellWidth().max());
            TBOX_ASSERT(U_fc_data->getGhostCellWidth().min() == U_fc_data->getGhostCellWidth().max());
#endif
            const int U_sc_gcw = U_sc_data->getGhostCellWidth().max();
            const int U_fc_gcw = U_fc_data->getGhostCellWidth().max();
            NAVIER_STOKES_SIDE_TO_FACE_FC(ilower(0),
                                          iupper(0),
                                          ilower(1),
                                          iupper(1),
#if (NDIM == 3)
                                          ilower(2),
                                          iupper(2),
#endif
                                          U_sc_data->getPointer(0),
                                          U_sc_data->getPointer(1),
#if (NDIM == 3)
                                          U_sc_data->getPointer(2),
#endif
                                          U_sc_gcw,
                                          U_fc_data->getPointer(0),
                                          U_fc_data->getPointer(1),
#if (NDIM == 3)
                                          U_fc_data->getPointer(2),
#endif
                                          U_fc_gcw);
        }
    }
    return;
} // copy_side_to_face

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////
/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityCartGridFunction::VelocityCartGridFunction(const std::string& object_name,
                                                   Pointer<CartesianGridGeometry<NDIM>> grid_geometry,
                                                   int velocity_index,
                                                   const std::string& db_base_name)
    : CartGridFunction(object_name),
      d_grid_geometry(grid_geometry),
      d_velocity_index(velocity_index),
      d_db_base_name(db_base_name)
{
    ensureDirectoryExists(d_db_base_name);
}

bool
VelocityCartGridFunction::isTimeDependent() const
{
    return true;
}

void
VelocityCartGridFunction::setDataOnPatchHierarchy(const int data_idx,
                                                  Pointer<Variable<NDIM>> /*var*/,
                                                  Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                  const double data_time,
                                                  const bool /*initial_time*/,
                                                  int coarsest_ln,
                                                  int finest_ln)
{
    copy_side_to_face(data_idx, d_velocity_index, hierarchy);
}

void
VelocityCartGridFunction::setDataOnPatch(const int /*data_idx*/,
                                         Pointer<Variable<NDIM>> /*var*/,
                                         Pointer<Patch<NDIM>> /*patch*/,
                                         const double data_time,
                                         const bool /*initial_time*/,
                                         Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
}

void
VelocityCartGridFunction::saveVelocity(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                                       int iteration_num,
                                       double loop_time)
{
    // create <base_name>/<iteration_num>/
    ensureDirectoryExists(getIterationDir(iteration_num));

    const std::string filename = getIterationFilename(iteration_num);

    Pointer<HDFDatabase> hdf_db = new HDFDatabase(d_object_name + "::hdf_db");
    hdf_db->create(filename);

    hdf_db->putDouble("time", loop_time);

    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(*level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());

            const std::string key = "patch_data_" + std::to_string(ln) + "_" + std::to_string(p());
            patch->getPatchData(d_velocity_index)->putToDatabase(hdf_db->putDatabase(key));
        }
    }

    hdf_db->close();
}

double
VelocityCartGridFunction::readVelocity(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                                       int iteration_num)
{
    const std::string filename = getIterationFilename(iteration_num);
    TBOX_ASSERT(fileExists(filename) &&
                "VelocityCartGridFunction: cannot read, file does not exist for this iteration");

    Pointer<HDFDatabase> hdf_db = new HDFDatabase(d_object_name + "::hdf_db");
    hdf_db->open(filename);

    const double loop_time = hdf_db->getDouble("time");

    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(*level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const std::string key = "patch_data_" + std::to_string(ln) + "_" + std::to_string(p());
            TBOX_ASSERT(hdf_db->keyExists(key));
            patch->getPatchData(d_velocity_index)->getFromDatabase(hdf_db->getDatabase(key));
        }
    }

    hdf_db->close();
    return loop_time;
}

void
VelocityCartGridFunction::setVelocityIndex(const int velocity_index)
{
    d_velocity_index = velocity_index;
}

void
VelocityCartGridFunction::ensureDirectoryExists(const std::string& dir_path)
{
    if (!fileExists(dir_path))
    {
        // mode 0755: rwxr-xr-x
        int status = mkdir(dir_path.c_str(), 0755);
        TBOX_ASSERT(status == 0 && "VelocityCartGridFunction: failed to create directory");
    }
}

bool
VelocityCartGridFunction::fileExists(const std::string& filename)
{
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

std::string
VelocityCartGridFunction::getIterationDir(int iteration_num) const
{
    std::ostringstream oss;
    oss << d_db_base_name << "/velocity." << std::setfill('0') << std::setw(6) << iteration_num;
    return oss.str();
}

std::string
VelocityCartGridFunction::getIterationFilename(int iteration_num) const
{
    return getIterationDir(iteration_num) + "/velocity.hdf5";
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
