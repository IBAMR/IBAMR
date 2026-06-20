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
                                                   int velocity_index)
    : CartGridFunction(object_name), d_grid_geometry(grid_geometry), d_velocity_index(velocity_index)
{
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
VelocityCartGridFunction::save_velocity_data(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                                             const std::string& base_folder,
                                             double data_time,
                                             int iteration)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // create base folder and iteration subfolder on rank 0
    const std::string iter_folder = base_folder + "/" + std::to_string(iteration);
    if (rank == 0)
    {
        mkdir(base_folder.c_str(), 0755);
        mkdir(iter_folder.c_str(), 0755);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const std::string filename = iter_folder + "/proc." + std::to_string(rank);
    std::ofstream ofs(filename, std::ios::binary);

    const int ndim = NDIM;
    ofs.write(reinterpret_cast<const char*>(&data_time), sizeof(double));
    ofs.write(reinterpret_cast<const char*>(&ndim), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&rank), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&iteration), sizeof(int));

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(*level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> sd = patch->getPatchData(d_velocity_index);

            auto box = patch->getBox();
            auto lo = box.lower();
            auto hi = box.upper();

            for (int d = 0; d < NDIM; ++d)
            {
                const int lo_d = lo(d);
                const int hi_d = hi(d);
                ofs.write(reinterpret_cast<const char*>(&lo_d), sizeof(int));
                ofs.write(reinterpret_cast<const char*>(&hi_d), sizeof(int));
            }

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> s(box, axis); s; s++)
                {
                    double val = (*sd)(*s, 0);
                    ofs.write(reinterpret_cast<const char*>(&val), sizeof(double));
                }
            }
        }
    }
}

double
VelocityCartGridFunction::read_velocity_data(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                                             const std::string& base_folder,
                                             int iteration)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const std::string filename = base_folder + "/" + std::to_string(iteration) + "/proc." + std::to_string(rank);
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs.is_open())
    {
        TBOX_ERROR("read_side_data: cannot open file: " + filename);
    }

    double data_time;
    int ndim, file_rank, file_iter;
    ifs.read(reinterpret_cast<char*>(&data_time), sizeof(double));
    ifs.read(reinterpret_cast<char*>(&ndim), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&file_rank), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&file_iter), sizeof(int));

    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(*level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> sd = patch->getPatchData(d_velocity_index);

            auto box = patch->getBox();

            // Read (and discard) the box bounds written by save_side_data.
            int file_lo, file_hi;
            for (int d = 0; d < NDIM; ++d)
            {
                ifs.read(reinterpret_cast<char*>(&file_lo), sizeof(int));
                ifs.read(reinterpret_cast<char*>(&file_hi), sizeof(int));
            }

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> s(box, axis); s; s++)
                {
                    double val;
                    ifs.read(reinterpret_cast<char*>(&val), sizeof(double));
                    (*sd)(*s, 0) = val;
                }
            }
        }
    }

    return data_time;
}

void
VelocityCartGridFunction::setVelocityIndex(const int velocity_index)
{
    d_velocity_index = velocity_index;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
