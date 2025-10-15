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
#include <regex>

#include "ibtk/app_namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

RestartCleaner::RestartCleaner(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_object_name(object_name),
      d_restart_base_path(),
      d_strategy(parseStrategy(input_db->getStringWithDefault("cleanup_strategy", "KEEP_RECENT_N"))),
      d_keep_restart_count(input_db->getIntegerWithDefault("keep_recent_files", 5)),
      d_enabled(input_db->getBoolWithDefault("enable_cleaner", false)),
      d_log_actions(input_db->getBoolWithDefault("log_cleaning_actions", true)),
      d_dry_run(input_db->getBoolWithDefault("dry_run", false))
{
    if (d_enabled)
    {
        if (input_db->keyExists("restart_directory"))
        {
            d_restart_base_path = input_db->getString("restart_directory");
        }
        else
        {
            TBOX_ERROR(d_object_name << "::RestartCleaner(): "
                                     << "'restart_directory' must be specified when enable_cleaner is true"
                                     << std::endl);
        }
    }
} // RestartCleaner

void
RestartCleaner::cleanup()
{
    if (!d_enabled) return;

    // Only perform file system operations on the master processor
    if (IBTK_MPI::getRank() == 0)
    {
        executeStrategy();
    }
    IBTK_MPI::barrier();
} // cleanup

std::vector<int>
RestartCleaner::getAvailableIterations() const
{
    std::vector<int> iterations;

    if (IBTK_MPI::getRank() == 0)
    {
        auto dirs = getAllRestartDirs(d_restart_base_path);
        for (const auto& dir : dirs)
        {
            int iter = parseIterationNum(dir.filename().string());
            if (iter >= 0)
            {
                iterations.push_back(iter);
            }
        }
        std::sort(iterations.begin(), iterations.end());
    }

    int vec_size = static_cast<int>(iterations.size());
    vec_size = IBTK_MPI::bcast(vec_size, 0);

    // Resize vector on non-master processors
    if (IBTK_MPI::getRank() != 0)
    {
        iterations.resize(vec_size);
    }

    if (vec_size > 0)
    {
        IBTK_MPI::bcast(iterations.data(), vec_size, 0);
    }

    return iterations;
} // getAvailableIterations

bool
RestartCleaner::isEnabled() const
{
    return d_enabled;
} // isEnabled

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

int
RestartCleaner::parseIterationNum(const std::string& dirname) const
{
    std::regex pattern("restore\\.([0-9]{6})");
    std::smatch match;

    if (std::regex_match(dirname, match, pattern))
    {
        return std::stoi(match[1].str());
    }
    return -1;
} // parseIterationNum

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

    std::regex pattern("restore\\.([0-9]{6})");

    try
    {
        for (const auto& entry : std::filesystem::directory_iterator(restart_dir))
        {
            if (entry.is_directory())
            {
                std::string dirname = entry.path().filename().string();
                if (std::regex_match(dirname, pattern))
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
    auto dirs = getAllRestartDirs(d_restart_base_path);

    if (dirs.size() <= static_cast<size_t>(d_keep_restart_count))
    {
        if (d_log_actions)
        {
            plog << d_object_name << "::keepRecentN(): Found " << dirs.size()
                 << " directories, keeping all (threshold: " << d_keep_restart_count << ")" << std::endl;
        }
        return;
    }

    // Sort directories by iteration number (ascending order)
    std::sort(dirs.begin(),
              dirs.end(),
              [this](const std::filesystem::path& a, const std::filesystem::path& b)
              {
                  int iter_a = parseIterationNum(a.filename().string());
                  int iter_b = parseIterationNum(b.filename().string());
                  return iter_a < iter_b;
              });

    // Delete older directories
    size_t dirs_to_delete_count = dirs.size() - d_keep_restart_count;

    if (d_log_actions)
    {
        plog << d_object_name << "::keepRecentN(): Found " << dirs.size() << " directories, "
             << "keeping " << d_keep_restart_count << " most recent, "
             << "deleting " << dirs_to_delete_count << " oldest" << std::endl;
    }

    for (size_t i = 0; i < dirs_to_delete_count; ++i)
    {
        if (d_dry_run)
        {
            if (d_log_actions)
            {
                plog << d_object_name << "::keepRecentN(): [DRY RUN] Would delete " << dirs[i] << std::endl;
            }
        }
        else
        {
            if (d_log_actions)
            {
                plog << d_object_name << "::keepRecentN(): Deleting " << dirs[i] << std::endl;
            }

            std::error_code ec;
            std::uintmax_t removed_count = std::filesystem::remove_all(dirs[i], ec);

            if (ec)
            {
                TBOX_ERROR(d_object_name << "::keepRecentN(): Failed to delete " << dirs[i] << ": " << ec.message()
                                         << std::endl);
            }
            else if (d_log_actions)
            {
                plog << d_object_name << "::keepRecentN(): Successfully deleted " << removed_count
                     << " files/directories from " << dirs[i] << std::endl;
            }
        }
    }

    if (d_log_actions)
    {
        plog << d_object_name << "::keepRecentN(): Cleanup completed" << std::endl;
    }
} // keepRecentN

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
