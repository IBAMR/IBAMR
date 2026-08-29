// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/InputDatabase.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <petscsys.h>

#include <fstream>
#include <string>

#include <ibamr/app_namespaces.h>

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // The expected-error case must still leave an output file for attest.
    std::ofstream initial_output("output");
    initial_output << '\n';
    initial_output.close();

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "petsc_options_file.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    const std::string test_case = input_db->getString("test_case");

    // Normal process success makes attest reject this expected-error case.
    if (test_case == "missing") return 0;

    PetscInt expected_value = 0;
    if (test_case == "lowercase")
    {
        expected_value = 1101;
    }
    else if (test_case == "uppercase_precedence")
    {
        expected_value = 2202;
    }
    else if (test_case == "source_fallback")
    {
        Pointer<InputDatabase> options_db = new InputDatabase("options_db");
        options_db->putString("petsc_options_file", "petsc_options_source_fallback.opts");
        load_petsc_options_file(options_db, 1, argv, SOURCE_DIR);
        expected_value = 3303;
    }
    else
    {
        TBOX_ERROR("unknown PETSc options-file test case: " << test_case << "\n");
    }

    PetscInt loaded_value = 0;
    PetscBool value_is_set = PETSC_FALSE;
    int ierr = PetscOptionsGetInt(nullptr, nullptr, "-r00_value", &loaded_value, &value_is_set);
    IBTK_CHKERRQ(ierr);
    if (!value_is_set || loaded_value != expected_value)
    {
        TBOX_ERROR("PETSc options-file test expected " << expected_value << " but loaded " << loaded_value << "\n");
    }

    std::ofstream output("output");
    output << "loaded value = " << loaded_value << '\n';

    return 0;
}
