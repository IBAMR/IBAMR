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
#include <ibtk/IBTK_MPI.h>

#include <tbox/Database.h>
#include <tbox/Logger.h>
#include <tbox/Utilities.h>

#include <petscoptions.h>
#include <petscsys.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
constexpr int ROOT_RANK = 0;

void
verify_values(Pointer<Database> expected)
{
    const Array<std::string> keys = expected->getAllKeys();
    int failures = 0;
    for (int k = 0; k < keys.size(); ++k)
    {
        const std::string name = "-" + keys[k];
        PetscBool present = PETSC_FALSE;
        int ierr = 0;
        switch (expected->getArrayType(keys[k]))
        {
        case Database::SAMRAI_BOOL:
        {
            PetscBool value = PETSC_FALSE;
            ierr = PetscOptionsGetBool(nullptr, nullptr, name.c_str(), &value, &present);
            failures += value != expected->getBool(keys[k]);
            break;
        }
        case Database::SAMRAI_INT:
        {
            PetscInt value = 0;
            ierr = PetscOptionsGetInt(nullptr, nullptr, name.c_str(), &value, &present);
            failures += value != expected->getInteger(keys[k]);
            break;
        }
        case Database::SAMRAI_DOUBLE:
        {
            PetscReal value = 0.0;
            ierr = PetscOptionsGetReal(nullptr, nullptr, name.c_str(), &value, &present);
            failures += value != expected->getDouble(keys[k]);
            break;
        }
        case Database::SAMRAI_STRING:
        {
            const std::string expected_value = expected->getString(keys[k]);
            std::vector<char> value(expected_value.size() + 1);
            ierr = PetscOptionsGetString(nullptr, nullptr, name.c_str(), value.data(), value.size(), &present);
            failures += expected_value != value.data();
            break;
        }
        default:
            TBOX_ERROR("unsupported expected type\n");
        }
        IBTK_CHKERRQ(ierr);
        failures += !present;
    }
    if (IBTK_MPI::sumReduction(failures)) TBOX_ERROR("PETSc settings value mismatch on one or more ranks\n");
}

// Record the real abort message, without volatile source-file/line metadata.
class OptionsErrorAppender : public Logger::Appender
{
public:
    /*! \brief Write the SAMRAI abort message to the test output. */
    void logMessage(const std::string& message, const std::string&, const int) override
    {
        if (IBTK_MPI::getRank() != ROOT_RANK) return;
        std::ofstream output("output");
        // SAMRAI appends a NUL: use c_str() as in TestAppender.
        output << message.c_str() << std::flush;
    }
};
} // namespace

int
main(int argc, char* argv[])
{
    if (argc < 2) return EXIT_FAILURE;
    // Leave input symlinks for AppInitializer to resolve.
    const std::filesystem::path input_path = std::filesystem::absolute(argv[1]);
    unsetenv("PETSC_OPTIONS");
    unsetenv("PETSC_OPTIONS_YAML");
    unsetenv("PETSC_OPTIONS_YAML_DIRECTORY");

    std::vector<std::string> storage = { argv[0], input_path.string(), "-skip_petscrc" };
    if (input_path.string().find(".lowercase.") != std::string::npos)
        storage.insert(storage.end(), { "-r00_value", "4404" });
    if (input_path.string().find(".inline.") != std::string::npos)
        storage.insert(storage.end(), { "-SP_pc_type", "jacobi" });
    std::vector<char*> args;
    for (std::string& value : storage) args.push_back(value.data());
    args.push_back(nullptr);
    std::vector<char*> init_args = args;
    IBTKInit ibtk_init(static_cast<int>(storage.size()), init_args.data(), MPI_COMM_WORLD);
    Pointer<Logger::Appender> abort_appender = new OptionsErrorAppender();
    Logger::getInstance()->setAbortAppender(abort_appender);

    Pointer<AppInitializer> initializer =
        new AppInitializer(static_cast<int>(storage.size()), args.data(), "petsc_options_file.log");
    // Normal continuation must return success, so attest rejects an expected-error case.
    if (input_path.string().find(".expect_error=true.") != std::string::npos) return EXIT_SUCCESS;
    Pointer<Database> input = initializer->getInputDatabase();
    if (input->isDatabase("Expected")) verify_values(input->getDatabase("Expected"));
    if (IBTK_MPI::getRank() == ROOT_RANK)
    {
        std::ofstream output("output");
        output << "PETSc options verified\n";
    }
    return EXIT_SUCCESS;
}
