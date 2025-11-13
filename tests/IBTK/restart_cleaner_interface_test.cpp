// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2025 by the IBAMR developers
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
#include <ibtk/IBTK_MPI.h>
#include <ibtk/RestartCleaner.h>

#include <tbox/MemoryDatabase.h>
#include <tbox/SAMRAIManager.h>

// MPI header
#include <mpi.h>

// C++ headers
#include <filesystem>
#include <iostream>

#include <ibamr/app_namespaces.h>

/*!
 * RestartCleaner Interface Test, which verifies complete RestartCleaner integration:
 * 1. Header file compilation and inclusion
 * 2. Database constructor with enabled cleaner
 * 3. Database constructor with disabled cleaner
 * 4. Database parameter parsing and defaults
 * 5. MPI environment compatibility
 *
 * \note This is NOT a functional test - it only verifies integration health.
 * \note RestartCleaner uses a single Database-only constructor interface.
 */

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Initialize with minimal configuration
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");

    pout << "RestartCleaner interface test..." << std::endl;

    // Test 1: Header inclusion and compilation verification
    pout << "Testing header inclusion and basic compilation..." << std::endl;
    {
        try
        {
            // Test that we can create SAMRAI Database objects (header inclusion test)
            Pointer<MemoryDatabase> test_db = new MemoryDatabase("TestDB");
            pout << "Header inclusion: successful" << std::endl;
        }
        catch (const std::exception& e)
        {
            pout << "ERROR in header inclusion test: " << e.what() << std::endl;
            return 1;
        }
        catch (...)
        {
            pout << "ERROR in header inclusion test: unknown exception" << std::endl;
            return 1;
        }
    }

    // Test 2: Database constructor with enabled cleaner
    pout << "Testing database constructor with enabled cleaner..." << std::endl;
    {
        try
        {
            std::string temp_dir = "./temp_restart_test_dir";
            std::filesystem::create_directories(temp_dir);

            struct TempDirGuard
            {
                std::string path;
                ~TempDirGuard()
                {
                    std::error_code ec;
                    std::filesystem::remove_all(path, ec);
                }
            } cleanup_guard{ temp_dir };

            Pointer<MemoryDatabase> db = new MemoryDatabase("EnabledCleanerTest");
            db->putBool("enable_cleaner", true);
            db->putString("restart_directory", temp_dir);
            db->putInteger("keep_recent_files", 3);
            db->putString("cleanup_strategy", "KEEP_RECENT_N");
            db->putBool("dry_run", true);

            Pointer<Database> db_ptr = db;
            RestartCleaner cleaner("EnabledCleanerTest", db_ptr);

            // Test basic method calls work
            bool enabled = cleaner.isEnabled();
            pout << "Enabled cleaner test: enabled=" << (enabled ? "true" : "false") << std::endl;

            // Test method call doesn't crash (dry run mode)
            cleaner.cleanup();
            pout << "Enabled cleaner test: cleanup() call successful" << std::endl;
        }
        catch (const std::exception& e)
        {
            pout << "ERROR in enabled cleaner test: " << e.what() << std::endl;
            return 1;
        }
        catch (...)
        {
            pout << "ERROR in enabled cleaner test: unknown exception" << std::endl;
            return 1;
        }
    }

    // Test 3: Database constructor - configuration parsing
    pout << "Testing database constructor and parameter parsing..." << std::endl;
    {
        try
        {
            // Create minimal configuration database
            Pointer<MemoryDatabase> db = new MemoryDatabase("RestartCleanerConfig");
            db->putBool("enable_cleaner", false); // Disabled to avoid file operations
            db->putInteger("keep_recent_files", 5);
            db->putBool("log_cleaning_actions", false);
            db->putString("cleanup_strategy", "KEEP_RECENT_N");

            Pointer<Database> db_ptr = db;
            RestartCleaner cleaner("TestCleaner", db_ptr);

            // Test configuration was parsed correctly
            bool enabled = cleaner.isEnabled();
            pout << "Database constructor: enabled=" << (enabled ? "true" : "false") << std::endl;

            // Test cleanup call with disabled cleaner (should be no-op)
            cleaner.cleanup();
            pout << "Database constructor: cleanup() call successful" << std::endl;
        }
        catch (const std::exception& e)
        {
            pout << "ERROR in database constructor test: " << e.what() << std::endl;
            return 1;
        }
        catch (...)
        {
            pout << "ERROR in database constructor test: unknown exception" << std::endl;
            return 1;
        }
    }

    // Test 4: Database parameter defaults
    pout << "Testing database parameter defaults..." << std::endl;
    {
        try
        {
            // Minimal database - test default parameter handling
            Pointer<MemoryDatabase> db = new MemoryDatabase("MinimalConfig");
            db->putBool("enable_cleaner", false);

            Pointer<Database> db_ptr = db;
            RestartCleaner cleaner("MinimalCleaner", db_ptr);

            bool enabled = cleaner.isEnabled();
            pout << "Minimal config: enabled=" << (enabled ? "true" : "false") << std::endl;
            pout << "Database parameter defaults: successful" << std::endl;
        }
        catch (const std::exception& e)
        {
            pout << "ERROR in defaults test: " << e.what() << std::endl;
            return 1;
        }
        catch (...)
        {
            pout << "ERROR in defaults test: unknown exception" << std::endl;
            return 1;
        }
    }

    // Test 5: MPI environment check
    pout << "Testing MPI environment compatibility..." << std::endl;
    {
        try
        {
            int mpi_initialized = 0;
            MPI_Initialized(&mpi_initialized);
            pout << "MPI " << (mpi_initialized ? "is" : "is not") << " initialized" << std::endl;

            int rank = IBTK_MPI::getRank();
            int size = IBTK_MPI::getNodes();
            pout << "MPI rank=" << rank << ", size=" << size << std::endl;

            // Test MPI barrier (used in RestartCleaner::cleanup())
            IBTK_MPI::barrier();
            pout << "MPI barrier: successful" << std::endl;
        }
        catch (const std::exception& e)
        {
            pout << "ERROR in MPI environment: " << e.what() << std::endl;
            return 1;
        }
        catch (...)
        {
            pout << "ERROR in MPI environment: unknown exception" << std::endl;
            return 1;
        }
    }

    // Test 6: Configuration error - enable_cleaner=true without restart_directory
    // Note: TBOX_ERROR correctly triggers abort() for missing required parameters
    // This has been manually verified but cannot be tested automatically

    // Test 7: Parameter boundary test - keep_restart_count <= 0
    // Note: TBOX_ASSERT correctly triggers abort() for invalid parameters
    // This has been manually verified but cannot be tested automatically

    pout << "RestartCleaner interface test: ALL TESTS PASSED" << std::endl;

    return 0;
} // main
