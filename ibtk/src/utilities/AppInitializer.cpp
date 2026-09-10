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

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LSiloDataWriter.h>

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/InputDatabase.h>
#include <tbox/InputManager.h>
#include <tbox/NullDatabase.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/RestartManager.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

#include <petsclog.h>
#include <petscoptions.h>
#include <petscsys.h>

#include <VisItDataWriter.h>

#include <charconv>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <ostream>
#include <set>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
constexpr std::string_view PETSC_SETTINGS_KEY = "petsc_settings";
constexpr std::string_view TAGGED_PETSC_SETTINGS_PREFIX = "petsc_settings_";

std::filesystem::path
resolve_petsc_options_file(const std::filesystem::path& advertised_path, const std::filesystem::path& input_filename)
{
    constexpr int ROOT_RANK = 0;
    std::string resolved_path_string;
    if (IBTK_MPI::getRank() == ROOT_RANK)
    {
        std::error_code error_code;
        bool path_exists = std::filesystem::exists(advertised_path, error_code);
        if (error_code)
            TBOX_ERROR("IBTK_PETSC_OPTIONS_FILESYSTEM_ERROR: could not inspect PETSc options file '"
                       << advertised_path.string() << "': " << error_code.message() << '\n');

        std::filesystem::path resolved_path = advertised_path;
        if (!path_exists)
        {
            resolved_path = input_filename.parent_path() / advertised_path;
            path_exists = std::filesystem::exists(resolved_path, error_code);
            if (error_code)
                TBOX_ERROR("IBTK_PETSC_OPTIONS_FILESYSTEM_ERROR: could not inspect PETSc options file '"
                           << resolved_path.string() << "': " << error_code.message() << '\n');
        }

        if (!path_exists)
            TBOX_ERROR("IBTK_PETSC_OPTIONS_FILE_MISSING: could not open PETSc options file '"
                       << advertised_path.string() << "'\n");
        resolved_path_string = resolved_path.string();
    }

    int resolved_path_size = IBTK_MPI::getRank() == ROOT_RANK ? static_cast<int>(resolved_path_string.size()) : 0;
    resolved_path_size = IBTK_MPI::bcast(resolved_path_size, ROOT_RANK);
    resolved_path_string.resize(resolved_path_size);
    if (resolved_path_size > 0) IBTK_MPI::bcast(resolved_path_string.data(), resolved_path_size, ROOT_RANK);
    return std::filesystem::path(resolved_path_string);
}

bool
is_settings_key(const std::string& key)
{
    return key == PETSC_SETTINGS_KEY.data() ||
           key.compare(0, TAGGED_PETSC_SETTINGS_PREFIX.size(), TAGGED_PETSC_SETTINGS_PREFIX.data()) == 0;
}

template <class T>
std::string
floating_point_to_string(const T value)
{
    constexpr std::size_t FLOATING_POINT_FORMAT_OVERHEAD = 16;
    // max_digits10 preserves round trips instead of rounding small tolerances to zero as std::to_string() can.
    char buffer[std::numeric_limits<T>::max_digits10 + FLOATING_POINT_FORMAT_OVERHEAD];
    const auto result = std::to_chars(
        buffer, buffer + sizeof(buffer), value, std::chars_format::general, std::numeric_limits<T>::max_digits10);
    if (result.ec != std::errc{}) TBOX_ERROR("Could not format a floating-point PETSc setting\n");
    return std::string(buffer, result.ptr);
}

void
insert_petsc_options(const std::string& options_file, int argc, char* argv[])
{
    // PetscOptionsInsert() may change argc and argv, so work with local copies.
    std::vector<std::string> argument_storage;
    argument_storage.reserve(argc);
    for (int k = 0; k < argc; ++k) argument_storage.emplace_back(argv[k]);
    std::vector<char*> arguments;
    arguments.reserve(argc + 1);
    for (std::string& argument : argument_storage) arguments.push_back(argument.data());
    arguments.push_back(nullptr);
    int local_argc = argc;
    char** local_argv = arguments.data();

    int ierr =
        PetscOptionsInsert(nullptr, &local_argc, &local_argv, options_file.empty() ? nullptr : options_file.c_str());
    IBTK_CHKERRQ(ierr);
}

bool
insert_petsc_settings(Pointer<Database> input_db, std::set<std::string>& option_names)
{
    bool found_settings = false;
    const Array<std::string> database_keys = input_db->getAllKeys();
    for (int database_key_n = 0; database_key_n < database_keys.size(); ++database_key_n)
    {
        const std::string& database_key = database_keys[database_key_n];
        if (!is_settings_key(database_key))
        {
            if (input_db->isDatabase(database_key) &&
                insert_petsc_settings(input_db->getDatabase(database_key), option_names))
                found_settings = true;
            continue;
        }

        found_settings = true;
        if (!input_db->isDatabase(database_key))
            TBOX_ERROR("PETSc settings entry '" << database_key << "' must be a database\n");
        Pointer<Database> settings_db = input_db->getDatabase(database_key);
        const std::string prefix = database_key == PETSC_SETTINGS_KEY.data() ?
                                       "" :
                                       database_key.substr(TAGGED_PETSC_SETTINGS_PREFIX.size()) + "_";
        const Array<std::string> keys = settings_db->getAllKeys();
        for (int key_n = 0; key_n < keys.size(); ++key_n)
        {
            const std::string& key = keys[key_n];
            if (settings_db->isDatabase(key) || settings_db->getArraySize(key) != 1)
                TBOX_ERROR("PETSc setting '" << key << "' must be a scalar value\n");

            std::string value;
            switch (settings_db->getArrayType(key))
            {
            case Database::SAMRAI_BOOL:
                value = settings_db->getBool(key) ? "true" : "false";
                break;
            case Database::SAMRAI_INT:
                value = std::to_string(settings_db->getInteger(key));
                break;
            case Database::SAMRAI_FLOAT:
                value = floating_point_to_string(settings_db->getFloat(key));
                break;
            case Database::SAMRAI_DOUBLE:
                value = floating_point_to_string(settings_db->getDouble(key));
                break;
            case Database::SAMRAI_STRING:
                value = settings_db->getString(key);
                break;
            default:
                TBOX_ERROR("PETSc setting '" << key
                                             << "' must be a Boolean, integer, floating-point value, or string\n");
            }

            const std::string name = "-" + prefix + key;
            if (!option_names.insert(name).second)
                TBOX_ERROR("PETSc option '" << name << "' is defined more than once in inline settings\n");
            int ierr = PetscOptionsSetValue(nullptr, name.c_str(), value.c_str());
            IBTK_CHKERRQ(ierr);
        }
    }
    return found_settings;
}

} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

AppInitializer::AppInitializer(int argc, char* argv[], const std::string& default_log_file_name)
{
    if (argc == 1)
    {
        TBOX_ERROR("USAGE: " << argv[0] << " <input filename> <restart dir> <restore number> [options]\n"
                             << "OPTIONS: PETSc command line options; use -help for more information.\n");
    }

    // Process command line options.
    const std::string input_filename = argv[1];
    if (argc >= 4)
    {
        // Check whether this appears to be a restarted run.
        FILE* fstream = (IBTK_MPI::getRank() == 0 ? fopen(argv[2], "r") : nullptr);
        if (IBTK_MPI::bcast(fstream ? 1 : 0, 0) == 1)
        {
            d_restart_read_dirname = argv[2];
            d_restart_restore_num = atoi(argv[3]);
            d_is_from_restart = true;
        }
        if (fstream)
        {
            fclose(fstream);
        }
    }

    // Process restart data if this is a restarted run.
    if (d_is_from_restart)
    {
        RestartManager::getManager()->openRestartFile(
            d_restart_read_dirname, d_restart_restore_num, IBTK_MPI::getNodes());
    }

    // Create input database and parse all data in input file.
    d_input_db = new InputDatabase("input_db");
    InputManager::getManager()->parseInputFile(input_filename, d_input_db);

    // Process "Main" section of the input database.
    Pointer<Database> main_db = new NullDatabase();
    if (d_input_db->isDatabase("Main"))
    {
        main_db = d_input_db->getDatabase("Main");
    }

    // Configure PETSc options and then reapply normal PETSc sources so command-line options win.
    std::set<std::string> option_names;
    const bool found_settings = insert_petsc_settings(d_input_db, option_names);
    const bool has_options_file =
        d_input_db->keyExists("PETSC_OPTIONS_FILE") || d_input_db->keyExists("petsc_options_file");
    if (found_settings && has_options_file)
        TBOX_ERROR("Inline PETSc settings and a PETSc options file cannot both be specified\n");
    if (has_options_file)
    {
        const std::string key =
            d_input_db->keyExists("PETSC_OPTIONS_FILE") ? "PETSC_OPTIONS_FILE" : "petsc_options_file";
        const std::filesystem::path options_file =
            resolve_petsc_options_file(d_input_db->getString(key), input_filename);
        insert_petsc_options(options_file.string(), argc, argv);
    }
    else if (found_settings)
    {
        insert_petsc_options("", argc, argv);
    }

    // Configure logging options.
    std::string log_file_name = default_log_file_name;
    bool log_all_nodes = false;
    if (main_db->keyExists("log_file_name")) log_file_name = main_db->getString("log_file_name");
    if (main_db->keyExists("log_all_nodes")) log_all_nodes = main_db->getBool("log_all_nodes");
    if (!log_file_name.empty())
    {
        if (log_all_nodes)
        {
            PIO::logAllNodes(log_file_name);
        }
        else
        {
            PIO::logOnlyNodeZero(log_file_name);
        }
    }

    // Configure timer options.
    std::string timer_dump_interval_key_name;
    if (main_db->keyExists("timer_interval"))
    {
        timer_dump_interval_key_name = "timer_interval";
    }
    else if (main_db->keyExists("timer_dump_interval"))
    {
        timer_dump_interval_key_name = "timer_dump_interval";
    }
    else if (main_db->keyExists("timer_write_interval"))
    {
        timer_dump_interval_key_name = "timer_write_interval";
    }

    if (!timer_dump_interval_key_name.empty())
    {
        d_timer_dump_interval = main_db->getInteger(timer_dump_interval_key_name);
    }

    // Avoid some warnings by unconditionally creating the timer database, even if
    // we never use it:
    {
        Pointer<Database> timer_manager_db;
        if (d_input_db->isDatabase("TimerManager"))
        {
            timer_manager_db = d_input_db->getDatabase("TimerManager");
        }
        TimerManager::createManager(timer_manager_db);
    }

    // Configure visualization options.
    std::string viz_dump_interval_key_name;
    if (main_db->keyExists("viz_interval"))
    {
        viz_dump_interval_key_name = "viz_interval";
    }
    else if (main_db->keyExists("viz_dump_interval"))
    {
        viz_dump_interval_key_name = "viz_dump_interval";
    }
    else if (main_db->keyExists("viz_write_interval"))
    {
        viz_dump_interval_key_name = "viz_write_interval";
    }

    std::string viz_dump_dirname_key_name;
    if (main_db->keyExists("viz_dirname"))
    {
        viz_dump_dirname_key_name = "viz_dirname";
    }
    else if (main_db->keyExists("viz_dump_dirname"))
    {
        viz_dump_dirname_key_name = "viz_dump_dirname";
    }
    else if (main_db->keyExists("viz_write_dirname"))
    {
        viz_dump_dirname_key_name = "viz_write_dirname";
    }

    if (!viz_dump_interval_key_name.empty())
    {
        d_viz_dump_interval = main_db->getInteger(viz_dump_interval_key_name);
        if (!viz_dump_dirname_key_name.empty())
        {
            d_viz_dump_dirname = main_db->getString(viz_dump_dirname_key_name);
            if (d_viz_dump_dirname.empty())
            {
                pout << "WARNING: AppInitializer::AppInitializer(): " << viz_dump_interval_key_name << " > 0, but `"
                     << d_viz_dump_dirname << "' is empty\n";
            }
        }
        else
        {
            pout << "WARNING: AppInitializer::AppInitializer(): " << viz_dump_interval_key_name
                 << " > 0, but `viz_dump_dirname' is not specified in input file\n";
        }
    }

    std::string viz_writers_key_name;
    Array<std::string> viz_writers_arr;
    if (main_db->keyExists("viz_writer"))
    {
        viz_writers_key_name = "viz_writer";
    }
    else if (main_db->keyExists("viz_writers"))
    {
        viz_writers_key_name = "viz_writers";
    }

    if (!viz_writers_key_name.empty())
    {
        viz_writers_arr = main_db->getStringArray(viz_writers_key_name);
    }
    if (viz_writers_arr.size() > 0)
    {
        d_viz_writers = std::vector<std::string>(viz_writers_arr.getPointer(),
                                                 viz_writers_arr.getPointer() + viz_writers_arr.size());
    }

    if (d_viz_dump_interval == 0 && d_viz_writers.size() > 0)
    {
        if (main_db->keyExists(viz_dump_dirname_key_name))
        {
            d_viz_dump_dirname = main_db->getString(viz_dump_dirname_key_name);
            if (d_viz_dump_dirname.empty())
            {
                pout << "WARNING: AppInitializer::AppInitializer(): `" << viz_writers_key_name << "' is set, but `"
                     << viz_dump_dirname_key_name << "' is empty\n";
            }
        }
        else
        {
            pout << "WARNING: AppInitializer::AppInitializer(): `" << viz_writers_key_name
                 << "' is set, but key `viz_dump_dirname' not specifed in input file\n";
        }
    }

    for (const auto& viz_writer : d_viz_writers)
    {
        if (viz_writer == "VisIt")
        {
            int visit_number_procs_per_file = 1;
            if (main_db->keyExists("visit_number_procs_per_file"))
                visit_number_procs_per_file = main_db->getInteger("visit_number_procs_per_file");
            d_visit_data_writer =
                new VisItDataWriter<NDIM>("VisItDataWriter", d_viz_dump_dirname, visit_number_procs_per_file);
        }

        if (viz_writer == "Silo")
        {
            d_silo_data_writer = new LSiloDataWriter("LSiloDataWriter", d_viz_dump_dirname);
        }

        if (viz_writer == "ExodusII")
        {
            if (main_db->keyExists("exodus_filename")) d_exodus_filename = main_db->getString("exodus_filename");
        }

        if (viz_writer == "GMV")
        {
            if (main_db->keyExists("gmv_filename")) d_gmv_filename = main_db->getString("gmv_filename");
        }
    }

    // Configure restart options.
    std::string restart_dump_interval_key_name;
    if (main_db->keyExists("restart_interval"))
    {
        restart_dump_interval_key_name = "restart_interval";
    }
    else if (main_db->keyExists("restart_dump_interval"))
    {
        restart_dump_interval_key_name = "restart_dump_interval";
    }
    else if (main_db->keyExists("restart_write_interval"))
    {
        restart_dump_interval_key_name = "restart_write_interval";
    }

    std::string restart_dump_dirname_key_name;
    if (main_db->keyExists("restart_dirname"))
    {
        restart_dump_dirname_key_name = "restart_dirname";
    }
    else if (main_db->keyExists("restart_dump_dirname"))
    {
        restart_dump_dirname_key_name = "restart_dump_dirname";
    }
    else if (main_db->keyExists("restart_write_dirname"))
    {
        restart_dump_dirname_key_name = "restart_write_dirname";
    }

    if (!restart_dump_interval_key_name.empty())
    {
        d_restart_dump_interval = main_db->getInteger(restart_dump_interval_key_name);
        if (!restart_dump_dirname_key_name.empty())
        {
            d_restart_dump_dirname = main_db->getString(restart_dump_dirname_key_name);
            if (d_restart_dump_dirname.empty())
            {
                pout << "WARNING: AppInitializer::AppInitializer(): " << restart_dump_interval_key_name << " > 0, but `"
                     << d_restart_dump_dirname << "' is empty\n";
            }
        }
        else
        {
            pout << "WARNING: AppInitializer::AppInitializer(): " << restart_dump_interval_key_name
                 << " > 0, but `restart_dump_dirname' is not specified in input file\n";
        }
    }

    // Configure post-processing data output options.
    std::string data_dump_interval_key_name;
    if (main_db->keyExists("data_interval"))
    {
        data_dump_interval_key_name = "data_interval";
    }
    else if (main_db->keyExists("data_dump_interval"))
    {
        data_dump_interval_key_name = "data_dump_interval";
    }
    else if (main_db->keyExists("data_write_interval"))
    {
        data_dump_interval_key_name = "data_write_interval";
    }

    std::string data_dump_dirname_key_name;
    if (main_db->keyExists("data_dirname"))
    {
        data_dump_dirname_key_name = "data_dirname";
    }
    else if (main_db->keyExists("data_dump_dirname"))
    {
        data_dump_dirname_key_name = "data_dump_dirname";
    }
    else if (main_db->keyExists("data_write_dirname"))
    {
        data_dump_dirname_key_name = "data_write_dirname";
    }

    if (!data_dump_interval_key_name.empty())
    {
        d_data_dump_interval = main_db->getInteger(data_dump_interval_key_name);
        if (!data_dump_dirname_key_name.empty())
        {
            d_data_dump_dirname = main_db->getString(data_dump_dirname_key_name);
            if (d_data_dump_dirname.empty())
            {
                pout << "WARNING: AppInitializer::AppInitializer(): " << data_dump_interval_key_name << " > 0, but `"
                     << d_data_dump_dirname << "' is empty\n";
            }
        }
        else
        {
            pout << "WARNING: AppInitializer::AppInitializer(): " << data_dump_interval_key_name
                 << " > 0, but `data_dump_dirname' is not specified in input file\n";
        }
    }
    return;
} // AppInitializer

AppInitializer::~AppInitializer()
{
    InputManager::freeManager();
    return;
} // ~AppInitializer

Pointer<Database>
AppInitializer::getInputDatabase()
{
    return d_input_db;
} // getInputDatabase

bool
AppInitializer::isFromRestart() const
{
    return d_is_from_restart;
} // isFromRestart

const std::string&
AppInitializer::getRestartReadDirectory() const
{
    return d_restart_read_dirname;
}

int
AppInitializer::getRestartRestoreNumber() const
{
    return d_restart_restore_num;
}

Pointer<Database>
AppInitializer::getRestartDatabase(const bool suppress_warning)
{
    if (!d_is_from_restart && !suppress_warning)
    {
        pout << "WARNING: AppInitializer::getRestartDatabase(): Not a restarted run, restart "
                "database is empty\n";
    }
    return RestartManager::getManager()->getRootDatabase();
} // getRestartDatabase

Pointer<Database>
AppInitializer::getComponentDatabase(const std::string& component_name, const bool suppress_warning)
{
    const bool db_exists = d_input_db->isDatabase(component_name);
    if (!db_exists && !suppress_warning)
    {
        pout << "WARNING: AppInitializer::getComponentDatabase(): Database corresponding to "
                "component `"
             << component_name << "' not found in input\n";
        return new NullDatabase();
    }
    else
    {
        return d_input_db->getDatabase(component_name);
    }
} // getComponentDatabase

bool
AppInitializer::dumpVizData() const
{
    return d_viz_dump_interval > 0;
} // dumpVizData

int
AppInitializer::getVizDumpInterval() const
{
    return d_viz_dump_interval;
} // getVizDumpInterval

std::string
AppInitializer::getVizDumpDirectory() const
{
    return d_viz_dump_dirname;
} // getVizDumpDirectory

std::vector<std::string>
AppInitializer::getVizWriters() const
{
    return d_viz_writers;
} // getVizDumpDirectory

Pointer<VisItDataWriter<NDIM>>
AppInitializer::getVisItDataWriter() const
{
    return d_visit_data_writer;
} // getVisItDataWriter

Pointer<LSiloDataWriter>
AppInitializer::getLSiloDataWriter() const
{
    return d_silo_data_writer;
} // getLSiloDataWriter

std::string
AppInitializer::getExodusIIFilename(const std::string& prefix) const
{
    std::string exodus_filename;
    if (!d_exodus_filename.empty())
    {
        exodus_filename = d_viz_dump_dirname + "/" + prefix + d_exodus_filename;
    }
    return exodus_filename;
} // getExodusIIFilename

std::string
AppInitializer::getGMVFilename(const std::string& prefix) const
{
    std::string gmv_filename;
    if (!d_gmv_filename.empty())
    {
        gmv_filename = d_viz_dump_dirname + "/" + prefix + d_gmv_filename;
    }
    return gmv_filename;
} // getGMVFilename

bool
AppInitializer::dumpRestartData() const
{
    return d_restart_dump_interval > 0;
} // dumpRestartData

int
AppInitializer::getRestartDumpInterval() const
{
    return d_restart_dump_interval;
} // getRestartDumpInterval

std::string
AppInitializer::getRestartDumpDirectory() const
{
    return d_restart_dump_dirname;
} // getRestartDumpDirectory

bool
AppInitializer::dumpPostProcessingData() const
{
    return d_data_dump_interval > 0;
} // dumpPostProcessingData

int
AppInitializer::getPostProcessingDataDumpInterval() const
{
    return d_data_dump_interval;
} // getPostProcessingDataDumpInterval

std::string
AppInitializer::getPostProcessingDataDumpDirectory() const
{
    return d_data_dump_dirname;
} // getPostProcessingDataDumpDirectory

bool
AppInitializer::dumpTimerData() const
{
    return d_timer_dump_interval > 0;
} // dumpTimerData

int
AppInitializer::getTimerDumpInterval() const
{
    return d_timer_dump_interval;
} // getTimerDumpInterval

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
