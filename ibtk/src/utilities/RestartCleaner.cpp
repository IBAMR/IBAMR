// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2025 by the IBAMR developers
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

#include "ibtk/IBTK_MPI.h"
#include "ibtk/RestartCleaner.h"

#include "tbox/PIO.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <filesystem>
#include <optional>
#include <regex>

#include "ibtk/app_namespaces.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::regex RestartCleaner::s_restart_dir_pattern("restore\\.([0-9]{6,})");

/////////////////////////////// PUBLIC ///////////////////////////////////////

RestartCleaner::RestartCleaner(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_object_name(object_name),
      d_restart_base_path(),
      d_strategy(parseStrategy(input_db->getStringWithDefault("cleanup_strategy", "KEEP_RECENT_N"))),
      d_keep_restart_count(input_db->getIntegerWithDefault("keep_recent_files", 5)),
      d_enable_logging(input_db->getBoolWithDefault("enable_logging", true)),
      d_dry_run(input_db->getBoolWithDefault("dry_run", false))
{
    if (d_keep_restart_count < 0)
    {
        TBOX_ERROR(d_object_name << "::RestartCleaner(): "
                                 << "Invalid configuration: keep_recent_files = " << d_keep_restart_count << "\n"
                                 << "keep_recent_files must be non-negative.\n"
                                 << "  - Use 0 to delete ALL restart directories\n"
                                 << "  - Use 1 or higher to keep that many recent directories\n"
                                 << "Please correct your configuration file." << std::endl);
    }
    if (d_keep_restart_count == 0)
    {
        plog << "WARNING: RestartCleaner::RestartCleaner(): "
             << "keep_recent_files is set to 0.\n"
             << "ALL restart directories will be deleted during cleanup!" << std::endl;
    }

    if (input_db->keyExists("restart_directory"))
    {
        d_restart_base_path = input_db->getString("restart_directory");
    }
    else
    {
        TBOX_ERROR(d_object_name << "::RestartCleaner(): "
                                 << "'restart_directory' must be specified" << std::endl);
    }
} // RestartCleaner

void
RestartCleaner::cleanup()
{
    // Only perform file system operations on the master processor
    if (IBTK_MPI::getRank() == 0)
    {
        executeStrategy();
    }
    IBTK_MPI::barrier();
} // cleanup

std::vector<int>
RestartCleaner::getAvailableRestartRestoreNumbers() const
{
    std::vector<int> restart_restore_numbers;

    if (IBTK_MPI::getRank() == 0)
    {
        auto dirs = getAllRestartDirs(d_restart_base_path);
        for (const auto& dir : dirs)
        {
            auto restart_restore_number = parseRestartRestoreNumber(dir.filename().string());
            if (restart_restore_number)
            {
                restart_restore_numbers.push_back(*restart_restore_number);
            }
        }
        std::sort(restart_restore_numbers.begin(), restart_restore_numbers.end());
    }

    int vec_size = static_cast<int>(restart_restore_numbers.size());
    vec_size = IBTK_MPI::bcast(vec_size, 0);

    // Resize vector on non-master processors
    if (IBTK_MPI::getRank() != 0)
    {
        restart_restore_numbers.resize(vec_size);
    }

    if (vec_size > 0)
    {
        IBTK_MPI::bcast(restart_restore_numbers.data(), vec_size, 0);
    }

    return restart_restore_numbers;
} // getAvailableRestartRestoreNumbers

/////////////////////////////// PRIVATE //////////////////////////////////////

RestartCleaner::CleanupStrategy
RestartCleaner::parseStrategy(const std::string& strategy_str) const
{
    if (strategy_str == "KEEP_RECENT_N")
    {
        return CleanupStrategy::KEEP_RECENT_N;
    }

    TBOX_ERROR(d_object_name << "::parseStrategy(): Unknown strategy: " << strategy_str << std::endl);
    return CleanupStrategy::KEEP_RECENT_N; // Never reached
} // parseStrategy

void
RestartCleaner::executeStrategy() const
{
    switch (d_strategy)
    {
    case CleanupStrategy::KEEP_RECENT_N:
        keepRecentN();
        break;
    default:
        TBOX_ERROR(d_object_name << "::executeStrategy(): Unknown strategy" << std::endl);
    }
} // executeStrategy

std::optional<int>
RestartCleaner::parseRestartRestoreNumber(const std::string& dirname) const
{
    std::smatch match;

    if (std::regex_match(dirname, match, s_restart_dir_pattern))
    {
        return std::stoi(match[1].str());
    }
    return std::nullopt;
} // parseRestartRestoreNumber

std::vector<std::filesystem::path>
RestartCleaner::getAllRestartDirs(const std::string& restart_dir) const
{
    std::vector<std::filesystem::path> restart_dirs;

    if (!std::filesystem::exists(restart_dir))
    {
        TBOX_ERROR(d_object_name << "::getAllRestartDirs(): Directory does not exist: " << restart_dir << std::endl);
    }

    if (!std::filesystem::is_directory(restart_dir))
    {
        TBOX_ERROR(d_object_name << "::getAllRestartDirs(): Path is not a directory: " << restart_dir << std::endl);
    }

    try
    {
        for (const auto& entry : std::filesystem::directory_iterator(restart_dir))
        {
            if (entry.is_directory())
            {
                std::string dirname = entry.path().filename().string();
                if (std::regex_match(dirname, s_restart_dir_pattern))
                {
                    restart_dirs.push_back(entry.path());
                }
            }
        }
    }
    catch (const std::filesystem::filesystem_error& ex)
    {
        TBOX_ERROR(d_object_name << "::getAllRestartDirs(): Filesystem error: " << ex.what() << std::endl);
    }

    return restart_dirs;
} // getAllRestartDirs

void
RestartCleaner::keepRecentN() const
{
    // Structure to cache parsed restart restore numbers with their paths
    struct DirInfo
    {
        int restart_restore_number;
        std::filesystem::path path;
    };

    // Build cache: parse once, filter invalid directories
    std::vector<DirInfo> valid_dirs;
    for (const auto& dir : getAllRestartDirs(d_restart_base_path))
    {
        auto restart_restore_number = parseRestartRestoreNumber(dir.filename().string());
        if (restart_restore_number)
        {
            valid_dirs.push_back({ *restart_restore_number, dir });
        }
        else
        {
            // Invalid directories are skipped but logged for user awareness
            TBOX_WARNING(d_object_name << "::keepRecentN(): Skipping directory with unparseable name: " << dir
                                       << "\n  This directory will not be managed by RestartCleaner." << std::endl);
        }
    }

    if (valid_dirs.size() <= static_cast<size_t>(d_keep_restart_count))
    {
        if (d_enable_logging)
        {
            plog << d_object_name << "::keepRecentN(): Found " << valid_dirs.size()
                 << " valid directories, keeping all (threshold: " << d_keep_restart_count << ")" << std::endl;
        }
        return;
    }

    // Sort by cached restart restore number (pure integer comparison, no regex)
    std::sort(valid_dirs.begin(),
              valid_dirs.end(),
              [](const DirInfo& a, const DirInfo& b) { return a.restart_restore_number < b.restart_restore_number; });

    // Delete older directories
    size_t dirs_to_delete_count = valid_dirs.size() - d_keep_restart_count;

    if (d_enable_logging)
    {
        plog << d_object_name << "::keepRecentN(): Found " << valid_dirs.size() << " valid directories, "
             << "keeping " << d_keep_restart_count << " most recent, "
             << "deleting " << dirs_to_delete_count << " oldest" << std::endl;
    }

    for (size_t i = 0; i < dirs_to_delete_count; ++i)
    {
        const auto& dir_path = valid_dirs[i].path;

        if (d_dry_run)
        {
            if (d_enable_logging)
            {
                plog << d_object_name << "::keepRecentN(): [DRY RUN] Would delete " << dir_path << std::endl;
            }
        }
        else
        {
            if (d_enable_logging)
            {
                plog << d_object_name << "::keepRecentN(): Deleting " << dir_path << std::endl;
            }

            std::error_code ec;
            std::uintmax_t removed_count = std::filesystem::remove_all(dir_path, ec);

            if (ec)
            {
                TBOX_WARNING(d_object_name << "::keepRecentN(): "
                                           << "Failed to delete " << dir_path << ": " << ec.message() << "\n"
                                           << "Stopping cleanup to maintain consistency. "
                                           << "Remaining old directories will not be deleted." << std::endl);
                return;
            }
            else if (d_enable_logging)
            {
                plog << d_object_name << "::keepRecentN(): Successfully deleted " << removed_count
                     << " files/directories from " << dir_path << std::endl;
            }
        }
    }

    if (d_enable_logging)
    {
        plog << d_object_name << "::keepRecentN(): Cleanup completed" << std::endl;
    }
} // keepRecentN

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
