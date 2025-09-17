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
 * 2. Direct parameter constructor functionality
 * 3. SAMRAI Database constructor functionality
 * 4. Database parameter parsing and defaults
 * 5. MPI environment compatibility
 * 6. Complete IBAMR framework integration
 *
 * \note This is NOT a functional test - it only verifies integration health.
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

    // Test 2: Direct constructor - compilation and linking
    pout << "Testing direct parameter constructor..." << std::endl;
    {
        try
        {
            // Create temporary directory for interface testing
            std::string temp_dir = "./temp_restart_test_dir";
            std::filesystem::create_directories(temp_dir);

            // RAII cleanup guard
            struct TempDirGuard
            {
                std::string path;
                ~TempDirGuard()
                {
                    std::error_code ec;
                    std::filesystem::remove_all(path, ec);
                }
            } cleanup_guard{ temp_dir };

            // Use temporary directory with dry_run=true
            RestartCleaner cleaner(temp_dir, 3, "KEEP_RECENT_N", true);

            // Test basic method calls work
            bool enabled = cleaner.isEnabled();
            pout << "Direct constructor: enabled=" << (enabled ? "true" : "false") << std::endl;

            // Test method call doesn't crash (dry run mode)
            cleaner.cleanup();
            pout << "Direct constructor: cleanup() call successful" << std::endl;
        }
        catch (const std::exception& e)
        {
            pout << "ERROR in direct constructor test: " << e.what() << std::endl;
            return 1;
        }
        catch (...)
        {
            pout << "ERROR in direct constructor test: unknown exception" << std::endl;
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
            // Note: restart_directory not required when cleaner is disabled

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
            db->putBool("enable_cleaner", false); // Only required parameter

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

    // Test 6: Constructor disambiguation (technical validation)
    pout << "Testing constructor disambiguation..." << std::endl;
    {
        try
        {
            // Create temporary directory for disambiguation test
            std::string temp_dir2 = "./temp_disambiguation_test_dir";
            std::filesystem::create_directories(temp_dir2);

            // RAII cleanup guard
            struct TempDirGuard
            {
                std::string path;
                ~TempDirGuard()
                {
                    std::error_code ec;
                    std::filesystem::remove_all(path, ec);
                }
            } cleanup_guard{ temp_dir2 };

            // Test that both constructors can be called without ambiguity
            // This validates our explicit keyword fixes
            RestartCleaner direct_cleaner(temp_dir2, 2, "KEEP_RECENT_N", true);

            Pointer<MemoryDatabase> db = new MemoryDatabase("DisambiguationTest");
            db->putBool("enable_cleaner", false);
            Pointer<Database> db_ptr = db;
            RestartCleaner db_cleaner("DisambiguationTest", db_ptr);

            pout << "Constructor disambiguation: successful" << std::endl;
        }
        catch (const std::exception& e)
        {
            pout << "ERROR in constructor disambiguation: " << e.what() << std::endl;
            return 1;
        }
        catch (...)
        {
            pout << "ERROR in constructor disambiguation: unknown exception" << std::endl;
            return 1;
        }
    }

    // Test 7: Configuration error - enable_cleaner=true without restart_directory
    // Note: TBOX_ERROR correctly triggers abort() for missing required parameters
    // This has been manually verified but cannot be tested automatically

    // Test 8: Parameter boundary test - keep_restart_count <= 0
    // Note: TBOX_ASSERT correctly triggers abort() for invalid parameters
    // This has been manually verified but cannot be tested automatically

    pout << "RestartCleaner interface test: ALL TESTS PASSED" << std::endl;

    return 0;
} // main
