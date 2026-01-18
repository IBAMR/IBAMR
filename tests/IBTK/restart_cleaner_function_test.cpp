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

#include <tbox/Database.h>
#include <tbox/MemoryDatabase.h>
#include <tbox/SAMRAIManager.h>

#include <mpi.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

/*!
 * RestartCleaner Functional Test, which verifies complete RestartCleaner functionality with real cases:
 *
 * This comprehensive test suite validates RestartCleaner's core capabilities using MemoryDatabase
 * configuration for all scenarios:
 * 1. Directory pattern parsing and filtering
 * 2. File system cleanup operations
 * 3. Dry run functionality verification
 * 4. Database constructor parameter handling
 * 5. Error handling and edge cases
 * 6. MPI environment compatibility
 *
 * \note This is an END-TO-END functional test using real file system operations.
 */

namespace
{
// RAII guard for automatic test directory cleanup
class TestDirGuard
{
public:
    explicit TestDirGuard(const std::string& path) : d_path(path)
    {
    }

    ~TestDirGuard()
    {
        if (!d_path.empty() && std::filesystem::exists(d_path))
        {
            std::error_code ec;
            std::filesystem::remove_all(d_path, ec);
            // Silently ignore errors in destructor
        }
    }

    // Disable copy to prevent double deletion
    TestDirGuard(const TestDirGuard&) = delete;
    TestDirGuard& operator=(const TestDirGuard&) = delete;

private:
    std::string d_path;
};

void
create_test_restart_dirs(const std::string& base_path, const std::vector<int>& iterations)
{
    std::filesystem::create_directories(base_path);

    for (int iter : iterations)
    {
        std::ostringstream dirname;
        dirname << "restore." << std::setfill('0') << std::setw(6) << iter;
        std::filesystem::path dir_path = std::filesystem::path(base_path) / dirname.str();
        std::filesystem::create_directories(dir_path);

        // Create realistic IBAMR restart files
        for (int i = 0; i < 3; ++i)
        {
            std::filesystem::path data_file = dir_path / ("samrai_data_" + std::to_string(i) + ".dat");
            std::ofstream file(data_file);
            file << "SAMRAI restart data file " << i << " for iteration " << iter << std::endl;
            file.close();
        }

        // Create subdirectory with hierarchy data
        std::filesystem::path sub_dir = dir_path / "hier_data";
        std::filesystem::create_directories(sub_dir);
        std::filesystem::path hier_file = sub_dir / "hierarchy.samrai.00000";
        std::ofstream hier(hier_file);
        hier << "Hierarchy data for iteration " << iter << std::endl;
        hier.close();
    }
}

void
create_invalid_dirs(const std::string& base_path, const std::vector<std::string>& dir_names)
{
    std::filesystem::create_directories(base_path);

    for (const auto& name : dir_names)
    {
        std::filesystem::path dir_path = std::filesystem::path(base_path) / name;
        std::filesystem::create_directories(dir_path);

        // Create dummy files to make directories realistic
        std::filesystem::path test_file = dir_path / "invalid_data.txt";
        std::ofstream file(test_file);
        file << "Invalid directory content for " << name << std::endl;
        file.close();
    }
}

int
count_dirs_matching_pattern(const std::string& base_path)
{
    int count = 0;
    if (!std::filesystem::exists(base_path)) return 0;

    std::regex pattern("restore\\.([0-9]{6,})");
    for (const auto& entry : std::filesystem::directory_iterator(base_path))
    {
        if (entry.is_directory())
        {
            std::string dirname = entry.path().filename().string();
            if (std::regex_match(dirname, pattern))
            {
                count++;
            }
        }
    }

    return count;
}
} // namespace

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Initialize with proper configuration
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");

    pout << "RestartCleaner functional test..." << std::endl;
    pout << "Testing complete end-to-end functionality with real file operations" << std::endl;

    int test_failures = 0;

    // Test 1: Directory pattern parsing and filtering
    pout << "\n=== Test 1: Directory Pattern Parsing and Filtering ===" << std::endl;
    {
        const std::string test_dir = "test_restart_parsing";
        TestDirGuard guard(test_dir);

        // Create mix of valid and invalid directories
        // Test coverage: 6-digit (standard), 6-digit boundary, and 7+ digit (long simulations)
        std::vector<int> valid_iterations = { 1, 100, 200, 300, 1000, 2500, 3000, 5000, 999999, 1000000, 1234567, 10000000 };
        std::vector<std::string> invalid_names = {
            "restore.12345",        // string less than 6 digits
            "restore.abc123",       // contains letters
            "restore_invalid",      // wrong format (underscore instead of dot)
            "other_directory",      // completely different
            "restore.",             // no digits
            "restore.000100_backup" // extra suffix
        };

        create_test_restart_dirs(test_dir, valid_iterations);
        create_invalid_dirs(test_dir, invalid_names);

        try
        {
            // Create database configuration for parsing test
            Pointer<MemoryDatabase> db = new MemoryDatabase("ParsingTestConfig");
                        db->putString("restart_directory", test_dir);
            db->putInteger("keep_recent_files", 20);
            db->putString("cleanup_strategy", "KEEP_RECENT_N");
            db->putBool("dry_run", true);

            RestartCleaner cleaner("ParsingTest", db);
            auto parsed_iterations = cleaner.getAvailableIterations();

            pout << "Valid iterations detected: ";
            for (int iter : parsed_iterations) pout << iter << " ";
            pout << std::endl;

            // Check count - should ignore invalid directories
            if (parsed_iterations.size() != valid_iterations.size())
            {
                pout << "FAILED: Expected " << valid_iterations.size() << " valid iterations, found "
                     << parsed_iterations.size() << std::endl;
                test_failures++;
            }
            else
            {
                // Check sorting and content
                std::sort(valid_iterations.begin(), valid_iterations.end());
                if (parsed_iterations == valid_iterations)
                {
                    pout << "Test 1 PASSED: Correctly parsed and sorted all valid directories" << std::endl;
                }
                else
                {
                    pout << "FAILED: Parsing mismatch in content or sorting" << std::endl;
                    test_failures++;
                }
            }
        }
        catch (const std::exception& e)
        {
            pout << "FAILED: Exception in parsing test: " << e.what() << std::endl;
            test_failures++;
        }
    }

    // Test 2: Cleanup functionality with real file deletion
    pout << "\n=== Test 2: Real File System Cleanup Operations ===" << std::endl;
    {
        const std::string test_dir = "test_restart_cleanup";
        TestDirGuard guard(test_dir);
        std::vector<int> iterations = { 100, 200, 300, 400, 500, 600, 700, 800 };

        create_test_restart_dirs(test_dir, iterations);

        try
        {
            int dirs_before = count_dirs_matching_pattern(test_dir);
            pout << "Directories before cleanup: " << dirs_before << std::endl;

            // Create database configuration for cleanup test
            Pointer<MemoryDatabase> db = new MemoryDatabase("CleanupTestConfig");
                        db->putString("restart_directory", test_dir);
            db->putInteger("keep_recent_files", 3);
            db->putString("cleanup_strategy", "KEEP_RECENT_N");
            db->putBool("dry_run", false);

            RestartCleaner cleaner("CleanupTest", db);
            cleaner.cleanup();

            int dirs_after = count_dirs_matching_pattern(test_dir);
            pout << "Directories after cleanup: " << dirs_after << std::endl;

            auto remaining_iterations = cleaner.getAvailableIterations();
            pout << "Remaining iterations: ";
            for (int iter : remaining_iterations) pout << iter << " ";
            pout << std::endl;

            // Should keep exactly 3 directories
            if (dirs_after == 3 && remaining_iterations.size() == 3)
            {
                // Check that we kept the highest 3
                std::vector<int> expected = { 600, 700, 800 };
                if (remaining_iterations == expected)
                {
                    pout << "Test 2 PASSED: Correctly cleaned up and kept 3 most recent" << std::endl;
                }
                else
                {
                    pout << "FAILED: Wrong directories were kept" << std::endl;
                    test_failures++;
                }
            }
            else
            {
                pout << "FAILED: Wrong number of directories after cleanup" << std::endl;
                test_failures++;
            }
        }
        catch (const std::exception& e)
        {
            pout << "FAILED: Exception in cleanup test: " << e.what() << std::endl;
            test_failures++;
        }
    }

    // Test 3: Dry run functionality verification
    pout << "\n=== Test 3: Dry Run Functionality Verification ===" << std::endl;
    {
        const std::string test_dir = "test_restart_dryrun";
        TestDirGuard guard(test_dir);
        std::vector<int> iterations = { 10, 20, 30, 40, 50, 60 };

        create_test_restart_dirs(test_dir, iterations);

        try
        {
            int dirs_before = count_dirs_matching_pattern(test_dir);
            pout << "Directories before dry run: " << dirs_before << std::endl;

            // Create database configuration for dry run test
            Pointer<MemoryDatabase> db = new MemoryDatabase("DryRunTestConfig");
                        db->putString("restart_directory", test_dir);
            db->putInteger("keep_recent_files", 2);
            db->putString("cleanup_strategy", "KEEP_RECENT_N");
            db->putBool("dry_run", true);

            RestartCleaner cleaner("DryRunTest", db);
            cleaner.cleanup();

            int dirs_after = count_dirs_matching_pattern(test_dir);
            pout << "Directories after dry run: " << dirs_after << std::endl;

            // All directories should still exist after dry run
            if (dirs_before == dirs_after && dirs_after == 6)
            {
                pout << "Test 3 PASSED: Dry run correctly preserved all directories" << std::endl;
            }
            else
            {
                pout << "FAILED: Dry run should not modify file system" << std::endl;
                test_failures++;
            }
        }
        catch (const std::exception& e)
        {
            pout << "FAILED: Exception in dry run test: " << e.what() << std::endl;
            test_failures++;
        }
    }

    // Test 4: Database constructor functionality
    pout << "\n=== Test 4: Database Constructor Functionality ===" << std::endl;
    {
        const std::string test_dir = "test_restart_database";
        TestDirGuard guard(test_dir);
        std::vector<int> iterations = { 10, 20, 30, 40, 50, 60, 70, 80 };

        create_test_restart_dirs(test_dir, iterations);

        try
        {
            // Create database configuration
            Pointer<MemoryDatabase> db = new MemoryDatabase("DatabaseTest");
                        db->putString("restart_directory", test_dir);
            db->putInteger("keep_recent_files", 4);
            db->putBool("enable_logging", true);
            db->putString("cleanup_strategy", "KEEP_RECENT_N");

            RestartCleaner cleaner("DatabaseCleaner", db);
            cleaner.cleanup();

            // Should keep 4 recent (database specified value)
            auto iterations_after = cleaner.getAvailableIterations();
            pout << "Iterations after database cleanup: ";
            for (int iter : iterations_after) pout << iter << " ";
            pout << std::endl;

            // Should keep: 50, 60, 70, 80
            std::vector<int> expected = { 50, 60, 70, 80 };
            if (iterations_after == expected)
            {
                pout << "Test 4 PASSED: Database constructor worked correctly" << std::endl;
            }
            else
            {
                pout << "FAILED: Database constructor cleanup error" << std::endl;
                test_failures++;
            }
        }
        catch (const std::exception& e)
        {
            pout << "FAILED: Exception in database test: " << e.what() << std::endl;
            test_failures++;
        }
    }

    // Test 5: Error handling and edge cases
    pout << "\n=== Test 5: Error Handling and Edge Cases ===" << std::endl;
    {
        try
        {
            // Test empty directory
            const std::string empty_dir = "test_empty_dir";
            TestDirGuard empty_guard(empty_dir);
            std::filesystem::create_directories(empty_dir);

            // Create database configuration for empty directory test
            Pointer<MemoryDatabase> db1 = new MemoryDatabase("EmptyDirTestConfig");
                        db1->putString("restart_directory", empty_dir);
            db1->putInteger("keep_recent_files", 5);
            db1->putString("cleanup_strategy", "KEEP_RECENT_N");
            db1->putBool("dry_run", true);

            RestartCleaner cleaner("EmptyDirTest", db1);
            auto result = cleaner.getAvailableIterations();
            if (result.empty())
            {
                pout << "Empty directory test: PASSED" << std::endl;
            }
            else
            {
                pout << "FAILED: Empty directory should return empty vector" << std::endl;
                test_failures++;
            }

            // Test directory with no valid restore directories
            const std::string invalid_dir = "test_invalid_content";
            TestDirGuard invalid_guard(invalid_dir);
            std::filesystem::create_directories(invalid_dir);
            std::filesystem::create_directories(invalid_dir + "/not_a_restore_dir");
            std::filesystem::create_directories(invalid_dir + "/restore_wrong_format");

            // Create database configuration for invalid content test
            Pointer<MemoryDatabase> db2 = new MemoryDatabase("InvalidContentTestConfig");
                        db2->putString("restart_directory", invalid_dir);
            db2->putInteger("keep_recent_files", 5);
            db2->putString("cleanup_strategy", "KEEP_RECENT_N");
            db2->putBool("dry_run", true);

            RestartCleaner cleaner2("InvalidContentTest", db2);
            auto result2 = cleaner2.getAvailableIterations();
            if (result2.empty())
            {
                pout << "Invalid content test: PASSED" << std::endl;
            }
            else
            {
                pout << "FAILED: Directory with no valid restores should return empty vector" << std::endl;
                test_failures++;
            }

            if (test_failures == 0)
            {
                pout << "Test 5 PASSED: Error handling works correctly" << std::endl;
            }
        }
        catch (const std::exception& e)
        {
            pout << "FAILED: Exception in error handling test: " << e.what() << std::endl;
            test_failures++;
        }
    }

    // Test 6: MPI environment compatibility
    pout << "\n=== Test 6: MPI Environment Compatibility ===" << std::endl;
    {
        try
        {
            int rank = IBTK_MPI::getRank();
            int size = IBTK_MPI::getNodes();
            pout << "MPI rank=" << rank << ", size=" << size << std::endl;

            // Test MPI barrier
            IBTK_MPI::barrier();

            pout << "Test 6 PASSED: MPI environment compatibility verified" << std::endl;
        }
        catch (const std::exception& e)
        {
            pout << "FAILED: Exception in MPI test: " << e.what() << std::endl;
            test_failures++;
        }
    }

    // Final summary
    pout << "\n=== RestartCleaner Functional Test Summary ===" << std::endl;
    if (test_failures == 0)
    {
        pout << "ALL TESTS PASSED!" << std::endl;
        return 0;
    }
    else
    {
        pout << "TEST FAILURES DETECTED: " << test_failures << " test(s) failed" << std::endl;
        pout << "Please check the implementation and test results above." << std::endl;
        return 1;
    }
} // main
