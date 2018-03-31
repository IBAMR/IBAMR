// Filename: AppInitializer.cpp
// Created on 19 Aug 2011 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <ostream>
#include <string>
#include <vector>

#include "VisItDataWriter.h"
#include "ibtk/AppInitializer.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/NullDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

AppInitializer::AppInitializer(int argc, char* argv[], const std::string& default_log_file_name)
    : d_input_db(NULL),
      d_is_from_restart(false),
      d_viz_dump_interval(0),
      d_viz_dump_dirname(""),
      d_viz_writers(),
      d_visit_data_writer(NULL),
      d_silo_data_writer(NULL),
      d_exodus_filename("output.ex2"),
      d_gmv_filename("output.gmv"),
      d_restart_dump_interval(0),
      d_restart_dump_dirname(""),
      d_data_dump_interval(0),
      d_data_dump_dirname(""),
      d_timer_dump_interval(0)
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
        FILE* fstream = (SAMRAI_MPI::getRank() == 0 ? fopen(argv[2], "r") : NULL);
        if (SAMRAI_MPI::bcast(fstream ? 1 : 0, 0) == 1)
        {
            d_restart_read_dirname = argv[2];
            d_restart_restore_num = atoi(argv[3]);
            d_is_from_restart = true;
        }
        else
        {
            d_restart_read_dirname = "";
            d_restart_restore_num = 0;
            d_is_from_restart = false;
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
            d_restart_read_dirname, d_restart_restore_num, SAMRAI_MPI::getNodes());
    }

    // Create input database and parse all data in input file.
    d_input_db = new InputDatabase("input_db");
    InputManager::getManager()->parseInputFile(input_filename, d_input_db);

    // Set custom PETSc options file when one is specified.
    if (d_input_db->keyExists("petsc_options_file"))
    {
        std::string petsc_options_file = d_input_db->getString("petsc_options_file");
        PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, petsc_options_file.c_str(), PETSC_TRUE);
    }

    // Process "Main" section of the input database.
    Pointer<Database> main_db = new NullDatabase();
    if (d_input_db->isDatabase("Main"))
    {
        main_db = d_input_db->getDatabase("Main");
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
    viz_writers_arr = main_db->getStringArray(viz_writers_key_name);
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

    for (unsigned int i = 0; i < d_viz_writers.size(); i++)
    {
        if (d_viz_writers[i] == "VisIt")
        {
            int visit_number_procs_per_file = 1;
            if (main_db->keyExists("visit_number_procs_per_file"))
                visit_number_procs_per_file = main_db->getInteger("visit_number_procs_per_file");
            d_visit_data_writer =
                new VisItDataWriter<NDIM>("VisItDataWriter", d_viz_dump_dirname, visit_number_procs_per_file);
        }

        if (d_viz_writers[i] == "Silo")
        {
            d_silo_data_writer = new LSiloDataWriter("LSiloDataWriter", d_viz_dump_dirname);
        }

        if (d_viz_writers[i] == "ExodusII")
        {
            if (main_db->keyExists("exodus_filename")) d_exodus_filename = main_db->getString("exodus_filename");
        }

        if (d_viz_writers[i] == "GMV")
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

    if (d_timer_dump_interval > 0)
    {
        Pointer<Database> timer_manager_db = new NullDatabase();
        if (d_input_db->isDatabase("TimerManager"))
        {
            timer_manager_db = d_input_db->getDatabase("TimerManager");
        }
        else
        {
            pout << "WARNING: AppInitializer::AppInitializer(): " << timer_dump_interval_key_name
                 << " > 0, but `TimerManager' input entries not specifed in input file\n";
        }
        TimerManager::createManager(timer_manager_db);
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

Pointer<VisItDataWriter<NDIM> >
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
    std::string exodus_filename = "";
    if (!d_exodus_filename.empty())
    {
        std::ostringstream exodus_filename_stream;
        exodus_filename_stream << d_viz_dump_dirname << "/" << prefix << d_exodus_filename;
        exodus_filename = exodus_filename_stream.str();
    }
    return exodus_filename;
} // getExodusIIFilename

std::string
AppInitializer::getGMVFilename(const std::string& prefix) const
{
    std::string gmv_filename = "";
    if (!d_gmv_filename.empty())
    {
        std::ostringstream gmv_filename_stream;
        gmv_filename_stream << d_viz_dump_dirname << "/" << prefix << d_gmv_filename;
        gmv_filename = gmv_filename_stream.str();
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
