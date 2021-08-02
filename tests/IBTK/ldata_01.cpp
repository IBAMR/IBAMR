// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

#include <SAMRAI_config.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>

// Set up application namespace declarations
#include <boost/multi_array.hpp>

#include <ibamr/app_namespaces.h>

// Do some basic verification tests for LData to guarantee that we no longer
// leak memory.

void
write_ldata_info(LData& l_data)
{
    if (IBTK_MPI::getRank() == 0)
    {
        std::ofstream output("output", std::ios_base::app);
        output << "object name: " << l_data.getName() << '\n';
    }

    for (int rank = 0; rank < IBTK_MPI::getNodes(); ++rank)
    {
        if (IBTK_MPI::getRank() == rank)
        {
            std::ofstream output("output", std::ios_base::app);
            output << "\nrank: " << rank << '\n'
                   << "global nodes: " << l_data.getGlobalNodeCount() << " local nodes: " << l_data.getLocalNodeCount()
                   << " ghost nodes: " << l_data.getGhostNodeCount() << '\n';

            output << "test locally indexed:\n";

            if (l_data.getDepth() == 1)
            {
                boost::multi_array_ref<double, 1> local_entries = *l_data.getLocalFormArray();
                output << "local entries: ";
                auto it = local_entries.begin();
                for (; it != local_entries.end() - 1; ++it) output << *it << ", ";
                output << *it << '\n';

                boost::multi_array_ref<double, 1> ghost_entries = *l_data.getGhostedLocalFormArray();
                output << "local ghost entries: ";
                it = ghost_entries.begin();
                for (; it != ghost_entries.end() - 1; ++it) output << *it << ", ";
                output << *it << '\n';

                output << '\n';
            }

            {
                output << "with 2D vector:\n";
                boost::multi_array_ref<double, 2> local_entries = *l_data.getLocalFormVecArray();
                output << "local entries: ";
                auto it = local_entries.begin();
                for (; it != local_entries.end(); ++it)
                {
                    output << '(';
                    for (unsigned int i = 0; i < local_entries.shape()[1] - 1; ++i) output << (*it)[i] << ", ";
                    output << (*it)[local_entries.shape()[1] - 1] << ')';
                    if (it != local_entries.end() - 1)
                        output << ", ";
                    else
                        output << '\n';
                }

                boost::multi_array_ref<double, 2> ghost_entries = *l_data.getGhostedLocalFormVecArray();
                output << "local ghost entries: ";
                it = ghost_entries.begin();
                for (; it != ghost_entries.end(); ++it)
                {
                    output << '(';
                    for (unsigned int i = 0; i < ghost_entries.shape()[1] - 1; ++i) output << (*it)[i] << ", ";
                    output << (*it)[ghost_entries.shape()[1] - 1] << ')';
                    if (it != ghost_entries.end() - 1)
                        output << ", ";
                    else
                        output << '\n';
                }

                output << '\n';
            }

            output << "test globally indexed:\n";

            if (l_data.getDepth() == 1)
            {
                boost::multi_array_ref<double, 1> local_entries = *l_data.getArray();
                const int begin = IBTK_MPI::getRank() * 10;
                const int end = begin + (IBTK_MPI::getRank() == 3 ? 6 : 10);
                output << "local range: " << '[' << begin << ", " << end << ")\n";
                output << "local entries: ";
                for (int i = begin; i < end - 1; ++i) output << local_entries[i] << ", ";
                output << local_entries[end - 1] << '\n';
            }

            if (l_data.getDepth() == 2)
            {
                output << "with 2D vector:\n" << std::flush;
                boost::multi_array_ref<double, 2> local_entries = *l_data.getVecArray();
                const int begin = IBTK_MPI::getRank() * 10;
                const int end = begin + (IBTK_MPI::getRank() == 3 ? 3 : 5);
                output << "local range: " << '[' << begin << ", " << end << ")\n";
                output << "index bases: " << local_entries.index_bases()[0] << ", " << local_entries.index_bases()[1]
                       << '\n';
                output << "local entries: ";
                for (int i = begin; i < end - 1; ++i)
                {
                    output << '(' << local_entries[i][0] << ", " << local_entries[i][1] << ')' << ", ";
                }
                output << '(' << local_entries[end - 1][0] << ", " << local_entries[end - 1][1] << ")\n";
            }
            output << "finished rank " << rank << '\n';
        }
        IBTK_MPI::barrier();
    }
}

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    const int rank = IBTK_MPI::getRank();

    // clear the output file, since we will append to it from multiple processes below:
    if (rank == 0)
    {
        std::ofstream output("output");
        output << "LData\n";
    }

    // test 1: set up an owning PETSc vector with depth 1.
    {
        std::vector<int> nonlocal_nodes;
        if (IBTK_MPI::getNodes() == 4)
        {
            switch (rank)
            {
            case 0:
                nonlocal_nodes = { 10, 20, 30 };
                break;
            case 1:
                nonlocal_nodes = { 1, 21, 31 };
                break;
            case 2:
                nonlocal_nodes = { 2, 12, 32 };
                break;
            case 3:
                nonlocal_nodes = { 3, 13 };
                break;
            default:
                TBOX_ERROR("only implemented for less than four processes");
            }
        }
        // use fewer entries on the last process
        LData l_data_1("l_data_1", rank == 3 ? 6 : 10, 1, nonlocal_nodes);

        for (int i = 0; i < IBTK_MPI::getNodes(); ++i)
        {
            boost::multi_array_ref<double, 1> entries = *l_data_1.getLocalFormArray();
            double offset = 0.0;
            for (double& entry : entries)
            {
                entry = offset + 100 * rank;
                ++offset;
            }
        }
        l_data_1.beginGhostUpdate();
        l_data_1.endGhostUpdate();
        write_ldata_info(l_data_1);
    }

    // test 2: set up an owning PETSc vector with depth 2.
    {
        std::vector<int> nonlocal_nodes;
        if (IBTK_MPI::getNodes() == 4)
        {
            switch (rank)
            {
            case 0:
                nonlocal_nodes = { 10, 11, 20, 21, 30, 31 };
                break;
            case 1:
                nonlocal_nodes = { 1, 2, 21, 22, 31, 32 };
                break;
            case 2:
                nonlocal_nodes = { 2, 3, 12, 13, 32, 33 };
                break;
            case 3:
                nonlocal_nodes = { 3, 4, 13, 14 };
                break;
            default:
                TBOX_ERROR("only implemented for less than four processes");
            }
        }
        // use fewer entries on the last process
        LData l_data_2("l_data_2", rank == 3 ? 6 : 10, 2, nonlocal_nodes);

        for (int i = 0; i < IBTK_MPI::getNodes(); ++i)
        {
            boost::multi_array_ref<double, 2> entries = *l_data_2.getLocalFormVecArray();
            for (unsigned int row_n = 0; row_n < entries.shape()[0]; ++row_n)
            {
                entries[row_n][0] = row_n + 100 * rank;
                entries[row_n][1] = row_n + 1000 * rank;
            }
        }
        l_data_2.beginGhostUpdate();
        l_data_2.endGhostUpdate();
        write_ldata_info(l_data_2);
    }
} // main
