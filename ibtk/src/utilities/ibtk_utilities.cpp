// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>

#include <petscsys.h>

#include <CartesianPatchGeometry.h>
#include <PatchLevel.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

#include <ibtk/app_namespaces.h>

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

void
load_petsc_options_file(Pointer<Database> input_db, int argc, char* argv[], const std::string& source_dir)
{
    std::string options_key;
    if (input_db->keyExists("PETSC_OPTIONS_FILE"))
    {
        options_key = "PETSC_OPTIONS_FILE";
    }
    else if (input_db->keyExists("petsc_options_file"))
    {
        options_key = "petsc_options_file";
    }
    else
    {
        return;
    }

    const std::string petsc_options_file = input_db->getString(options_key);
    std::string resolved_options_file = petsc_options_file;
    auto path_exists = [](const std::string& path) -> bool
    {
        std::ifstream stream(path.c_str());
        return stream.good();
    };
    if (!path_exists(resolved_options_file) && argc > 1)
    {
        std::filesystem::path input_filename = argv[1];
        std::error_code error_code;
        input_filename = std::filesystem::weakly_canonical(input_filename, error_code);
        if (error_code) input_filename = argv[1];
        const std::filesystem::path input_dir = input_filename.parent_path();
        if (!input_dir.empty())
        {
            const std::string candidate = (input_dir / petsc_options_file).string();
            if (path_exists(candidate)) resolved_options_file = candidate;
        }
    }
    if (!path_exists(resolved_options_file) && !source_dir.empty())
    {
        const std::string candidate = source_dir + "/" + petsc_options_file;
        if (path_exists(candidate)) resolved_options_file = candidate;
    }
    if (!path_exists(resolved_options_file))
    {
        TBOX_ERROR("could not open PETSc options file: " << petsc_options_file << "\n");
    }

    int ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, nullptr, resolved_options_file.c_str(), PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
}

double
get_min_patch_dx(const PatchLevel<NDIM>& patch_level)
{
    double result = std::numeric_limits<double>::max();

    // Some processors might not have any patches so its easier to just quit
    // after one loop operation than to check
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level.getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);
        result = std::min(result, patch_dx_min);
        break; // all patches on the same level have the same dx values
    }

    result = IBTK_MPI::minReduction(result);

    return result;
} // get_min_patch_dx
} // namespace IBTK
