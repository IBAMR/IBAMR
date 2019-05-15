// Copyright (c) 2019, Boyce Griffith
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

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

#include <boost/multi_array.hpp>

// Do some basic verification tests for LData to guarantee that we no longer
// leak memory.

void
write_ldata_info(LData& l_data)
{
    if (SAMRAI_MPI::getRank() == 0)
    {
        std::ofstream output("output", std::ios_base::app);
        output << "object name: " << l_data.getName() << '\n';
    }

    for (int rank = 0; rank < SAMRAI_MPI::getNodes(); ++rank)
    {
        if (SAMRAI_MPI::getRank() == rank)
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
                const int begin = SAMRAI_MPI::getRank() * 10;
                const int end = begin + (SAMRAI_MPI::getRank() == 3 ? 6 : 10);
                output << "local range: " << '[' << begin << ", " << end << ")\n";
                output << "local entries: ";
                for (int i = begin; i < end - 1; ++i) output << local_entries[i] << ", ";
                output << local_entries[end - 1] << '\n';
            }

            if (l_data.getDepth() == 2)
            {
                output << "with 2D vector:\n" << std::flush;
                boost::multi_array_ref<double, 2> local_entries = *l_data.getVecArray();
                const int begin = SAMRAI_MPI::getRank() * 10;
                const int end = begin + (SAMRAI_MPI::getRank() == 3 ? 3 : 5);
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
        SAMRAI_MPI::barrier();
    }
}

int
main(int argc, char** argv)
{
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    const int rank = SAMRAI_MPI::getRank();

    // clear the output file, since we will append to it from multiple processes below:
    if (rank == 0)
    {
        std::ofstream output("output");
        output << "LData\n";
    }

    // test 1: set up an owning PETSc vector with depth 1.
    {
        std::vector<int> nonlocal_nodes;
        if (SAMRAI_MPI::getNodes() == 4)
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

        for (int i = 0; i < SAMRAI_MPI::getNodes(); ++i)
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
        if (SAMRAI_MPI::getNodes() == 4)
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

        for (int i = 0; i < SAMRAI_MPI::getNodes(); ++i)
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

    SAMRAIManager::shutdown();
} // main
