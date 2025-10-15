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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_RestartCleaner
#define included_IBTK_RestartCleaner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/IBTK_MPI.h"

#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <filesystem>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class RestartCleaner provides functionality to manage restart directories.
 *
 * This class provides methods to automatically clean up old restart directories while
 * keeping the most recent ones. It is designed to work with IBAMR's standard restart
 * directory naming convention and integrates with the SAMRAI Database configuration system.
 *
 * The main functionalities include:
 * -# Scan a specified directory for subdirectories matching the pattern "restore.XXXXXX"
 * -# Parse iteration numbers from these directory names
 * -# Sort directories based on iteration numbers
 * -# Keep the N most recent directories and delete the rest
 *
 * \note This class assumes restart directories follow the naming pattern "restore.XXXXXX"
 * where XXXXXX is a zero-padded iteration number.
 *
 * \note This is a stateless utility class. Each cleanup operation performs a complete scan
 * of the target directory and makes decisions based on the current filesystem state. No
 * internal state is maintained between operations.
 *
 * \note In parallel environments, file operations are performed only on the master processor
 * (rank 0) to ensure consistency and avoid race conditions.
 *
 * Supported cleanup strategies:
 * - "KEEP_RECENT_N": Keep the N most recent restart directories (default)
 *
 * Sample usage:
 * \code
 * // Configuration-driven usage with SAMRAI Database
 * Pointer<Database> restart_db = input_db->getDatabase("RestartCleaner");
 * RestartCleaner cleaner("RestartCleaner", restart_db);
 * if (cleaner.isEnabled()) {
 *     cleaner.cleanup();
 * }
 *
 * // Check available iterations
 * auto iterations = cleaner.getAvailableIterations();
 * \endcode
 */
class RestartCleaner
{
public:
    /*!
     * \brief Constructor using SAMRAI Database configuration.
     *
     * This is the constructor for RestartCleaner and integrates with IBAMR's
     * standard configuration system. The cleaner can be enabled or disabled through
     * the configuration.
     *
     * The input database is searched for the following keys:
     * - 'enable_cleaner': bool (default: false)
     * - 'keep_recent_files': int (default: 5)
     * - 'restart_directory': string (required if enable_cleaner is true)
     * - 'log_cleaning_actions': bool (default: true)
     * - 'cleanup_strategy': string (default: "KEEP_RECENT_N")
     * - 'dry_run': bool (default: false)
     *
     * \param object_name Name for this object (used in error messages and logging)
     * \param input_db    Database containing configuration parameters
     */
    explicit RestartCleaner(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~RestartCleaner() = default;

    /*!
     * \brief Scan and cleanup old restart directories.
     *
     * This method performs the complete cleanup process:
     * 1. Scans the base directory for restart folders
     * 2. Parses iteration numbers from directory names
     * 3. Sorts directories by iteration number
     * 4. Keeps the N most recent directories and deletes the rest
     *
     * In parallel environments, file operations are performed only on the master
     * processor (rank 0) to ensure consistency.
     *
     * \note Call cleanup() after each restart save to maintain only the most recent directories.
     *
     * \note This method only performs cleanup if the cleaner is enabled (for Database
     * constructor) or always executes (for explicit parameter constructor).
     */
    void cleanup();

    /*!
     * \brief Get available iteration numbers.
     *
     * Scans the restart directory and returns iteration numbers for all currently
     * existing restore directories, sorted in ascending order.
     *
     * \return Vector of iteration numbers for available restore directories,
     *         sorted in ascending order
     */
    std::vector<int> getAvailableIterations() const;

    /*!
     * \brief Check if the cleaner is enabled.
     *
     * Returns the value of 'enable_cleaner' parameter from the Database.
     *
     * \return true if cleaning is enabled, false otherwise
     */
    bool isEnabled() const;

private:
    RestartCleaner() = delete;
    RestartCleaner(const RestartCleaner& from) = delete;
    RestartCleaner& operator=(const RestartCleaner& that) = delete;

    /*!
     * \brief Internal strategy enumeration.
     */
    enum class CleanupStrategy
    {
        KEEP_RECENT_N
    };

    /*!
     * \brief Parse strategy string to enum.
     */
    CleanupStrategy parseStrategy(const std::string& strategy_str) const;

    /*!
     * \brief Execute cleanup based on current strategy.
     */
    void executeStrategy() const;

    /*!
     * \brief Parse iteration number from directory name.
     *
     * Extracts the iteration number from directory names following the pattern
     * "restore.XXXXXX" where XXXXXX is a zero-padded number.
     *
     * \param dirname Directory name to parse
     * \return Iteration number, or -1 if parsing fails
     */
    int parseIterationNum(const std::string& dirname) const;

    /*!
     * \brief Get all restart directories from the base path.
     *
     * Scans the restart base directory and returns all subdirectories that
     * match the restart naming pattern.
     *
     * \param restart_dir Base directory to scan
     * \return Vector of filesystem paths to valid restart directories
     */
    std::vector<std::filesystem::path> getAllRestartDirs(const std::string& restart_dir) const;

    /*!
     * \brief KEEP_RECENT_N strategy implementation.
     */
    void keepRecentN() const;

    // Member variables
    std::string d_object_name;
    std::string d_restart_base_path;
    CleanupStrategy d_strategy;
    int d_keep_restart_count;
    bool d_enabled;
    bool d_log_actions;
    bool d_dry_run;
};

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_RestartCleaner
