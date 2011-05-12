// Filename: LagM3DDataWriter.C
// Created on 26 Apr 2005 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "LagM3DDataWriter.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/LagMarker.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <IndexData.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <limits>
#include <numeric>

// HDF5 INCLUDES
#include <hdf5.h>
#if (H5_VERS_MINOR == 6)
#include <H5LT.h>
#define H5Dcreate1 H5Dcreate
#define H5Gcreate1 H5Gcreate
#endif
#if (H5_VERS_MINOR == 8)
#include <hdf5_hl.h>
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// The rank of the root MPI process and the MPI tag number.
static const int M3D_MPI_ROOT = 0;

// Maximum number of fibers per group.
static const int M3D_NFG_MAX = 999;

// Buffer sizes.
static const int M3D_BUFSIZE = 128;

// The compression level for the local HDF5 files.
//
// NOTE: The deflation level may be any integer in the range [0,9], with lower
// values indicating lower (faster) compression levels and higher values
// indicating higher (slower) compression levels.
static const int M3D_DEFLATE_LEVEL = 4;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

LagM3DDataWriter::LagM3DDataWriter(
    const std::string& object_name,
    const std::string& dump_directory_name,
    const std::string& experiment_name,
    const int& experiment_number)
    : d_object_name(object_name),
      d_dump_directory_name(dump_directory_name),
      d_experiment_name(experiment_name),
      d_experiment_number(experiment_number),
      d_file_prefix(""),
      d_time_step_number(-1),
      d_hierarchy(),
      d_coarsest_ln(0),
      d_finest_ln(0),
      d_mark_idx(-1),
      d_nclouds(0),
      d_cloud_names(),
      d_cloud_nmarks(),
      d_cloud_first_mark_idx(),
      d_nblocks(d_finest_ln+1,0),
      d_block_names(d_finest_ln+1),
      d_block_nelems(d_finest_ln+1),
      d_block_periodic(d_finest_ln+1),
      d_block_nfibers(d_finest_ln+1),
      d_block_ngroups(d_finest_ln+1),
      d_block_first_lag_idx(d_finest_ln+1),
      d_coords_data(d_finest_ln+1,Pointer<LData>(NULL)),
      d_ao(d_finest_ln+1),
      d_build_vec_scatters(d_finest_ln+1),
      d_src_vec(d_finest_ln+1),
      d_dst_vec(d_finest_ln+1),
      d_vec_scatter(d_finest_ln+1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_experiment_name.size() == 3);
    TBOX_ASSERT((0 <= d_experiment_number) && (d_experiment_number <= 9999));
#endif
    std::ostringstream stream;
    stream << experiment_name << std::setfill('0') << std::setw(4) << d_experiment_number;
    d_file_prefix = stream.str();
    return;
}// LagM3DDataWriter

LagM3DDataWriter::~LagM3DDataWriter()
{
    // Destroy any remaining PETSc objects.
    int ierr;
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        for (std::map<int,Vec>::iterator it = d_dst_vec[ln].begin();
             it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = it->second;
            if (v)
            {
                ierr = VecDestroy(v);
                IBTK_CHKERRQ(ierr);
            }
        }
        for (std::map<int,VecScatter>::iterator it = d_vec_scatter[ln].begin();
             it != d_vec_scatter[ln].end(); ++it)
        {
            VecScatter& vs = it->second;
            if (vs)
            {
                ierr = VecScatterDestroy(vs);
                IBTK_CHKERRQ(ierr);
            }
        }
    }
    return;
}// ~LagM3DDataWriter

void
LagM3DDataWriter::setPatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT(hierarchy->getFinestLevelNumber() >= d_finest_ln);
#endif
    // Reset the hierarchy.
    d_hierarchy = hierarchy;
    return;
}// setPatchHierarchy

void
LagM3DDataWriter::resetLevels(
    const int coarsest_ln,
    const int finest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((coarsest_ln >= 0) && (finest_ln >= coarsest_ln));
    if (!d_hierarchy.isNull())
    {
        TBOX_ASSERT(finest_ln <= d_hierarchy->getFinestLevelNumber());
    }
#endif
    // Destroy any unneeded PETSc objects.
    int ierr;
    for (int ln = std::max(d_coarsest_ln,0); (ln <= d_finest_ln) && (ln < coarsest_ln); ++ln)
    {
        for (std::map<int,Vec>::iterator it = d_dst_vec[ln].begin();
             it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = it->second;
            if (v != static_cast<Vec>(NULL))
            {
                ierr = VecDestroy(v);  IBTK_CHKERRQ(ierr);
            }
        }
        for (std::map<int,VecScatter>::iterator it = d_vec_scatter[ln].begin();
             it != d_vec_scatter[ln].end(); ++it)
        {
            VecScatter& vs = it->second;
            if (vs != static_cast<VecScatter>(NULL))
            {
                ierr = VecScatterDestroy(vs);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    for (int ln = finest_ln+1; ln <= d_finest_ln; ++ln)
    {
        for (std::map<int,Vec>::iterator it = d_dst_vec[ln].begin();
             it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = it->second;
            if (v != static_cast<Vec>(NULL))
            {
                ierr = VecDestroy(v);  IBTK_CHKERRQ(ierr);
            }
        }
        for (std::map<int,VecScatter>::iterator it = d_vec_scatter[ln].begin();
             it != d_vec_scatter[ln].end(); ++it)
        {
            VecScatter& vs = it->second;
            if (vs != static_cast<VecScatter>(NULL))
            {
                ierr = VecScatterDestroy(vs);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    // Reset the level numbers.
    d_coarsest_ln = coarsest_ln;
    d_finest_ln   = finest_ln;

    // Resize some arrays.
    d_nblocks            .resize(d_finest_ln+1,0);
    d_block_names        .resize(d_finest_ln+1);
    d_block_nelems       .resize(d_finest_ln+1);
    d_block_periodic     .resize(d_finest_ln+1);
    d_block_nfibers      .resize(d_finest_ln+1);
    d_block_ngroups      .resize(d_finest_ln+1);
    d_block_first_lag_idx.resize(d_finest_ln+1);

    d_coords_data.resize(d_finest_ln+1,NULL);

    d_ao                .resize(d_finest_ln+1);
    d_build_vec_scatters.resize(d_finest_ln+1);
    d_src_vec           .resize(d_finest_ln+1);
    d_dst_vec           .resize(d_finest_ln+1);
    d_vec_scatter       .resize(d_finest_ln+1);
    return;
}// resetLevels

void
LagM3DDataWriter::registerLagMarkerPatchDataIndex(
    const int mark_idx)
{
    d_mark_idx = mark_idx;
    return;
}// registerLagMarkerPatchDataIndex

void
LagM3DDataWriter::registerMarkerCloud(
    const std::string& name,
    const int nmarks,
    const int first_mark_idx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(nmarks > 0);
#endif

    // Check to see if the cloud name has already been registered.
    if (find(d_cloud_names.begin(), d_cloud_names.end(), name) != d_cloud_names.end())
    {
        TBOX_ERROR(d_object_name << "::registerMarkerCloud()\n"
                   << "  marker clouds must have unique names.\n"
                   << "  a marker cloud named ``" << name << "'' has already been registered.\n");
    }

    // Record the layout of the marker cloud.
    ++d_nclouds;
    d_cloud_names         .push_back(name);
    d_cloud_nmarks        .push_back(nmarks);
    d_cloud_first_mark_idx.push_back(first_mark_idx);
    return;
}// registerMarkerCloud

void
LagM3DDataWriter::registerLogicallyCartesianBlock(
    const std::string& name,
    const IntVector<NDIM>& nelem,
    const IntVector<NDIM>& periodic,
    const int first_lag_idx,
    const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(std::min(level_number,d_coarsest_ln),std::max(level_number,d_finest_ln));
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(nelem[d] > 0);
        TBOX_ASSERT(periodic(d) == 0 || periodic(d) == 1);
    }
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif

    // Check to see if the block name has already been registered.
    if (find(d_block_names[level_number].begin(), d_block_names[level_number].end(),
             name) != d_block_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                   << "  Cartesian blocks must have unique names.\n"
                   << "  a Cartesian block named ``" << name << "'' has already been registered.\n");
    }

    // Record the layout of the logically Cartesian block.
    int nfibers, ngroups;
    if ((nelem[0] == 1 && nelem[1] == 1 && nelem[2] == 1))
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                   << "  attempting to register a single point as a logically Cartesian block.\n");
    }
    else if ((nelem[0] == 1 && nelem[1] == 1) ||
             (nelem[0] == 1 && nelem[2] == 1) ||
             (nelem[1] == 1 && nelem[2] == 1))
    {
        nfibers = 1;
        ngroups = 1;
    }
    else if ((nelem[0] == 1) || (nelem[1] == 1) || (nelem[2] == 1))
    {
        nfibers = ((nelem[0] == 1 ? 0 : nelem[0]) +
                   (nelem[1] == 1 ? 0 : nelem[1]) +
                   (nelem[2] == 1 ? 0 : nelem[2]));
        ngroups = 2;
    }
    else
    {
        nfibers = nelem[0]*nelem[1] + nelem[1]*nelem[2] + nelem[2]*nelem[0];
        ngroups = 3;
    }

    ++d_nblocks[level_number];
    d_block_names        [level_number].push_back(name);
    d_block_nelems       [level_number].push_back(nelem);
    d_block_periodic     [level_number].push_back(periodic);
    d_block_nfibers      [level_number].push_back(nfibers);
    d_block_ngroups      [level_number].push_back(ngroups);
    d_block_first_lag_idx[level_number].push_back(first_lag_idx);
    return;
}// registerLogicallyCartesianBlock

void
LagM3DDataWriter::registerCoordsData(
    Pointer<LData> coords_data,
    const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(std::min(level_number,d_coarsest_ln),std::max(level_number,d_finest_ln));
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!coords_data.isNull());
    TBOX_ASSERT(coords_data->getDepth() == NDIM);
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    d_coords_data[level_number] = coords_data;
    return;
}// registerCoordsData

void
LagM3DDataWriter::registerLagrangianAO(
    AO& ao,
    const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(std::min(level_number,d_coarsest_ln),std::max(level_number,d_finest_ln));
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= level_number &&
                d_finest_ln   >= level_number);
#endif
    d_ao[level_number] = ao;
    d_build_vec_scatters[level_number] = true;
    return;
}// registerLagrangianAO

void
LagM3DDataWriter::registerLagrangianAO(
    std::vector<AO>& ao,
    const int coarsest_ln,
    const int finest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln <= finest_ln);
#endif

    if (coarsest_ln < d_coarsest_ln || finest_ln > d_finest_ln)
    {
        resetLevels(std::min(coarsest_ln,d_coarsest_ln),std::max(finest_ln,d_finest_ln));
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_coarsest_ln <= coarsest_ln && finest_ln <= d_finest_ln);
#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        registerLagrangianAO(ao[ln], ln);
    }
    return;
}// registerLagrangianAO

void
LagM3DDataWriter::writePlotData(
    const int time_step_number,
    const double simulation_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(time_step_number >= 0);
    TBOX_ASSERT(!d_dump_directory_name.empty());
#endif

    if (time_step_number <= d_time_step_number)
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                   << "  data writer with name " << d_object_name << "\n"
                   << "  time step number: " << time_step_number
                   << " is <= last time step number: " << d_time_step_number
                   << std::endl);
    }
    d_time_step_number = time_step_number;

    if (d_dump_directory_name.empty())
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                   << "  data writer with name " << d_object_name << "\n"
                   << "  dump directory name is empty" << std::endl);
    }

    int ierr;
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int mpi_size = SAMRAI_MPI::getNodes();

    // Create the dump directory.
    Utilities::recursiveMkdir(d_dump_directory_name);

    // Determine the local file names.
    std::ostringstream stream;
    stream << std::setfill('0') << std::setw(4) << mpi_rank;
    const std::string mpi_rank_string = stream.str();

    const std::string marker_file_name = getMarkerFileName(time_step_number);
    const std::string local_marker_file_name = marker_file_name + "." + mpi_rank_string;
    const std::string marker_header_file_name = marker_file_name + ".hdr";

    const std::string fiber_file_name = getFiberFileName(time_step_number);
    const std::string local_fiber_file_name = fiber_file_name + "." + mpi_rank_string;
    const std::string fiber_header_file_name = fiber_file_name + ".hdr";

    // Determine the marker cloud and fiber offsets.
    int num_local_marker_nodes = std::accumulate(d_cloud_nmarks.begin(), d_cloud_nmarks.end(), 0);
    std::vector<int> num_marker_nodes_proc(mpi_size,0);
    SAMRAI_MPI::allGather(num_local_marker_nodes, &num_marker_nodes_proc[0]);
    const int marker_node_offset = std::accumulate(
        num_marker_nodes_proc.begin(), num_marker_nodes_proc.begin()+mpi_rank, 0);
    const int num_marker_nodes = std::accumulate(
        num_marker_nodes_proc.begin()+mpi_rank, num_marker_nodes_proc.end(), marker_node_offset);

    const int num_local_marker_clouds = d_nclouds;
    std::vector<int> num_marker_clouds_proc(mpi_size,0);
    SAMRAI_MPI::allGather(num_local_marker_clouds, &num_marker_clouds_proc[0]);
    const int marker_cloud_offset = std::accumulate(
        num_marker_clouds_proc.begin(), num_marker_clouds_proc.begin()+mpi_rank, 0);
    const int num_marker_clouds = std::accumulate(
        num_marker_clouds_proc.begin()+mpi_rank, num_marker_clouds_proc.end(), marker_cloud_offset);

    int local_marker_node_counter, local_marker_cloud_counter;

    int num_local_fibers = 0;
    for (std::vector<std::vector<int> >::const_iterator cit = d_block_nfibers.begin();
         cit != d_block_nfibers.end(); ++cit)
    {
        num_local_fibers += std::accumulate(cit->begin(), cit->end(), 0);
    }
    std::vector<int> num_fibers_proc(mpi_size,0);
    SAMRAI_MPI::allGather(num_local_fibers, &num_fibers_proc[0]);
    const int fiber_offset = std::accumulate(
        num_fibers_proc.begin(), num_fibers_proc.begin()+mpi_rank, 0);
    const int num_fibers = std::accumulate(
        num_fibers_proc.begin()+mpi_rank, num_fibers_proc.end(), fiber_offset);

    int num_local_groups = 0;
    for (std::vector<std::vector<int> >::const_iterator cit = d_block_ngroups.begin();
         cit != d_block_ngroups.end(); ++cit)
    {
        num_local_groups += std::accumulate(cit->begin(), cit->end(), 0);
    }
    std::vector<int> num_groups_proc(mpi_size,0);
    SAMRAI_MPI::allGather(num_local_groups, &num_groups_proc[0]);
    const int group_offset = std::accumulate(
        num_groups_proc.begin(), num_groups_proc.begin()+mpi_rank, 0);
    const int num_groups = std::accumulate(
        num_groups_proc.begin()+mpi_rank, num_groups_proc.end(), group_offset);

    const int num_local_layers = std::accumulate(d_nblocks.begin(), d_nblocks.end(), 0);
    std::vector<int> num_layers_proc(mpi_size,0);
    SAMRAI_MPI::allGather(num_local_layers, &num_layers_proc[0]);
    const int layer_offset = std::accumulate(
        num_layers_proc.begin(), num_layers_proc.begin()+mpi_rank, 0);
    const int num_layers = std::accumulate(
        num_layers_proc.begin()+mpi_rank, num_layers_proc.end(), layer_offset);

    int local_fiber_counter, local_group_counter, local_layer_counter;

    // Determine the length of the computational domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const grid_xLower = grid_geom->getXLower();
    const double* const grid_xUpper = grid_geom->getXUpper();
    const double L = std::min(std::min(grid_xUpper[0]-grid_xLower[0], grid_xUpper[1]-grid_xLower[1]), grid_xUpper[2]-grid_xLower[2]);

    // Create the header files.
    local_marker_node_counter = 0;
    local_marker_cloud_counter = 0;

    local_fiber_counter = 0;
    local_group_counter = 0;
    local_layer_counter = 0;
    for (int rank = 0; rank < mpi_size; ++rank)
    {
        if (rank == mpi_rank)
        {
            std::ofstream marker_header_stream(std::string(d_dump_directory_name + "/" + marker_header_file_name).c_str(), rank == 0 ? std::ios::out : std::ios::app);
            std::ofstream fiber_header_stream(std::string(d_dump_directory_name + "/" + fiber_header_file_name).c_str(), rank == 0 ? std::ios::out : std::ios::app);

            if (rank == 0)
            {
                marker_header_stream << "V3   = FORMAT VERSION\n"
                                     << std::setw(4) << int(L) << " = NG\n"
                                     << std::setw(4) << num_marker_clouds << " = NUMBER OF CLOUDS\n";

                fiber_header_stream << " " << std::setw(4) << int(L) << "           = NG\n"
                                    << "C" << std::setw(4) << num_layers << std::setw(4) << 1 << "       = MAXIMUM-LAYER-NUMBER, STARTING WITH 1\n";
            }

            for (int cloud = 0; cloud < d_nclouds; ++cloud)
            {
                const int nmarks = d_cloud_nmarks[cloud];
                marker_header_stream << std::setw(4) << nmarks << " = NUMBER OF POINTS IN CLOUD" << " " << std::setw(2) << local_marker_cloud_counter+marker_cloud_offset+1 << "\n";
                local_marker_node_counter += nmarks;
                local_marker_cloud_counter += 1;
            }

            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                for (int block = 0; block < d_nblocks[ln]; ++block)
                {
                    const int layer_number = layer_offset+local_layer_counter+1;
                    const IntVector<NDIM>& nelem = d_block_nelems[ln][block];
                    const IntVector<NDIM>& periodic = d_block_periodic[ln][block];
                    if ((nelem[0] == 1 && nelem[1] == 1) || (nelem[0] == 1 && nelem[2] == 1) || (nelem[1] == 1 && nelem[2] == 1))
                    {
                        // Output a single fiber.
                        fiber_header_stream << "C" << std::setw(4) << d_block_ngroups[ln][block] << std::setw(4) << layer_number << "       = NUMBER-OF-GROUPS LAYER-NUMBER\n"
                                            << "C" << std::setw(4) << group_offset+local_group_counter+1 << std::setw(4) << 1 << std::setw(4) << nelem.getProduct() << "   = GRP NFG NPF\n";
                    }
                    else if ((nelem[0] == 1) || (nelem[1] == 1) || (nelem[2] == 1))
                    {
                        // Output a 2D sheet of fibers.
                        fiber_header_stream << "C" << std::setw(4) << d_block_ngroups[ln][block] << std::setw(4) << layer_number << "       = NUMBER-OF-GROUPS LAYER-NUMBER\n";

                        for (int d0 = 0; d0 < NDIM; ++d0)
                        {
                            if (nelem[d0] > 1)
                            {
                                // Find the other nontrivial dimension.
                                for (int d1 = 0; d1 < NDIM; ++d1)
                                {
                                    if (d1 != d0 && nelem[d1] > 1)
                                    {
                                        fiber_header_stream << "C" << std::setw(4) << group_offset+local_group_counter+1 << std::setw(4) << nelem[d0] << std::setw(4) << nelem[d1] + (periodic[d1] ? 1 : 0) << "   = GRP NFG NPF\n";
                                        local_group_counter += 1;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        // Output a 3D volume of fibers.
                        int ngroups = 0;
                        for (int d0 = 0; d0 < NDIM; ++d0)
                        {
                            const int d1 = (d0+1)%NDIM;
                            const int nfibers_per_group = nelem[d0]*nelem[d1];
                            int split_factor = 1;
                            while(nfibers_per_group/split_factor > M3D_NFG_MAX)
                            {
                                split_factor *= 2;
                            }
                            ngroups += split_factor;
                        }

                        fiber_header_stream << "C" << std::setw(4) << ngroups << std::setw(4) << layer_number << "       = NUMBER-OF-GROUPS LAYER-NUMBER\n";
                        for (int d0 = 0; d0 < NDIM; ++d0)
                        {
                            const int d1 = (d0+1)%NDIM;
                            const int d2 = (d1+1)%NDIM;

                            // myocardial3D cannot cleanly handle headers
                            // specifying groups with more than 999 fibers.
                            const int nfibers_group = nelem[d0]*nelem[d1];
                            int split_factor = 1;
                            while(nfibers_group/split_factor > M3D_NFG_MAX)
                            {
                                split_factor *= 2;
                            }

                            const int split_nfibers_group = split_factor == 1 ? 0 : nfibers_group/split_factor;

                            for (int s = 0; s < split_factor-1; ++s)
                            {
                                fiber_header_stream << "C" << std::setw(4) << group_offset+local_group_counter+1 << std::setw(4) << split_nfibers_group << std::setw(4) << nelem[d2] + (periodic[d2] ? 1 : 0) << "   = GRP NFG NPF\n";
                            }

                            fiber_header_stream << "C" << std::setw(4) << group_offset+local_group_counter+1 << std::setw(4) << nfibers_group - (split_factor-1)*split_nfibers_group << std::setw(4) << nelem[d2] + (periodic[d2] ? 1 : 0) << "   = GRP NFG NPF\n";
                            local_group_counter += 1;
                        }
                    }
                    local_fiber_counter += d_block_nfibers[ln][block];
                    local_layer_counter += 1;
                }
            }

            if (rank == mpi_size-1)
            {
                marker_header_stream << std::setw(5) << time_step_number << "           = KLOK " << getMarkerFileName(time_step_number) << "\n"
                                     << std::setw(15) << std::fixed << std::setprecision(9) << simulation_time << std::resetiosflags(std::ios::scientific) << " = TIME\n";

                fiber_header_stream << "C*\n"
                                    << std::setw(5) << time_step_number << "           = KLOK " << getFiberFileName(time_step_number) << "\n"
                                    << std::setw(15) << std::fixed << std::setprecision(9) << simulation_time << std::resetiosflags(std::ios::scientific) << " = TIME\n"
                                    << std::setw(5) << num_fibers << "           = NUMBER OF FIBERS IN THIS DATA FILE\n";
            }
        }
        SAMRAI_MPI::barrier();
    }

    // Create the menu.text file.
    local_marker_node_counter = 0;
    local_marker_cloud_counter = 0;

    local_fiber_counter = 0;
    local_group_counter = 0;
    local_layer_counter = 0;
    if (time_step_number == 0)
    {
        const std::string menu_file_name = d_dump_directory_name + "/" + getMenuFileName();
        for (int rank = 0; rank < mpi_size; ++rank)
        {
            if (rank == mpi_rank)
            {
                std::ofstream menu_file_stream(menu_file_name.c_str(), (rank == M3D_MPI_ROOT ? std::ios::out : std::ios::app));
                if (rank == M3D_MPI_ROOT)
                {
                    menu_file_stream << "#\n"
                                     << "# layer names\n"
                                     << "#\n";
                }
                for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
                {
                    for (int block = 0; block < d_nblocks[ln]; ++block)
                    {
                        const int layer_number = layer_offset+local_layer_counter+1;
                        menu_file_stream << std::setw(32) << std::string("\"" + d_block_names[ln][block] + "\"") << "     " << std::setw(3) << layer_number << " = layer_name\n";
                        local_fiber_counter += d_block_nfibers[ln][block];
                        local_group_counter += d_block_ngroups[ln][block];
                        local_layer_counter += 1;
                    }
                }
            }
            SAMRAI_MPI::barrier();
        }

        for (int rank = 0; rank < mpi_size; ++rank)
        {
            if (rank == mpi_rank)
            {
                std::ofstream menu_file_stream(menu_file_name.c_str(), std::ios::app);
                if (rank == M3D_MPI_ROOT)
                {
                    menu_file_stream << "#\n"
                                     << "# marker cloud names\n"
                                     << "#\n";
                }
                for (int cloud = 0; cloud < d_nclouds; ++cloud)
                {
                    const int nmarks = d_cloud_nmarks[cloud];
                    menu_file_stream << std::setw(32) << std::string("\"" + d_cloud_names[cloud] + "\"") << "     " << std::setw(3) << local_marker_cloud_counter+marker_cloud_offset+1 << " = marker_name\n";
                    local_marker_node_counter += nmarks;
                    local_marker_cloud_counter += 1;
                }
            }
            SAMRAI_MPI::barrier();
        }

        if (mpi_rank == M3D_MPI_ROOT)
        {
            std::ofstream menu_file_stream(menu_file_name.c_str(), std::ios::app);
            menu_file_stream << "#\n"
                             << "# vectors names\n"
                             << "#\n"
                             << "#    intentionally blank\n"
                             << "#\n"
                             << "# activations names\n"
                             << "#\n"
                             << "#    intentionally blank\n"
                             << "#\n"
                             << "# flowmeters names\n"
                             << "#\n"
                             << "#    intentionally blank\n"
                             << "#\n"
                             << "# pressure taps names\n"
                             << "#\n"
                             << "#    intentionally blank\n"
                             << "#\n"
                             << "# colors\n"
                             << "#\n"
                             << "#    intentionally blank\n"
                             << "#\n"
                             << "# layer subdivisions\n"
                             << "#\n"
                             << "# marker cloud mass names\n"
                             << "#\n"
                             << "#    intentionally blank\n"
                             << "#\n"
                             << "# marker cloud mass numbers\n"
                             << "#\n"
                             << "#    intentionally blank\n"
                             << "#\n"
                             << "# vector cluster names\n"
                             << "#\n"
                             << "#    intentionally blank\n"
                             << "#\n"
                             << "# vector cluster numbers\n"
                             << "#\n"
                             << "#    intentionally blank\n"
                             << "#\n"
                             << "# =+=+=+=+=+=+=+=+=+= END OF DISPLAY PROGRAM DIRECTIVES =+=+=+=+=+=+=+=+=+=+=+=\n";
        }
        SAMRAI_MPI::barrier();
    }

    // Create list and cat files.
    if (mpi_rank == M3D_MPI_ROOT)
    {
        // Create or update the list file on the root MPI process.
        const std::string list_file_name = d_dump_directory_name + "/" + getListFileName();
        std::ofstream list_file_stream(list_file_name.c_str(), (time_step_number == 0 ? std::ios::out : std::ios::app));
        if (time_step_number > 0)
        {
            list_file_stream << "\n";  // we do not provide flow meter files
        }
        list_file_stream << marker_file_name << "\n" << fiber_file_name << "\n";

        // Create or update the cat script on the root MPI process.
        const std::string cat_file_name = d_dump_directory_name + "/" + getCatScriptFileName();
        std::ofstream cat_file_stream(cat_file_name.c_str(), (time_step_number == 0 ? std::ios::out : std::ios::app));
        if (time_step_number == 0)
        {
            cat_file_stream << "#!/bin/sh\n"
                            << "unset noclobber\n"
                            << "echo \"About to concatenate files in directory $PWD\"\n"
                            << "echo\n"
                            << "echo \"WARNING: this script *CAN* clobber existing files in the working directory\"\n"
                            << "echo \"         but DOES NOT modify the constituent source data files\"\n"
                            << "echo\n"
                            << "read -p \"Press <Enter> to continue...\"\n"
                            << "echo \"Concatenating files, please wait...\"\n";
        }

        // Setup commands to create a unified marker file.
        cat_file_stream << "\n";
        for (int rank = 0; rank < mpi_size; ++rank)
        {
            std::ostringstream stream;
            stream << std::setfill('0') << std::setw(4) << rank;
            const std::string rank_string = stream.str();
            const std::string local_marker_file_name = marker_file_name + "." + rank_string;
            cat_file_stream << "m3D_hdf5_marker_converter " << local_marker_file_name + ".h5" << " " << local_marker_file_name << "\n";
        }
        cat_file_stream << "\n";
        cat_file_stream << "cat " << marker_header_file_name;
        for (int rank = 0; rank < mpi_size; ++rank)
        {
            std::ostringstream stream;
            stream << std::setfill('0') << std::setw(4) << rank;
            const std::string rank_string = stream.str();
            cat_file_stream << " " << marker_file_name + "." + rank_string;
        }
        cat_file_stream << " > " << marker_file_name << "\n";
        cat_file_stream << "\n";
        for (int rank = 0; rank < mpi_size; ++rank)
        {
            std::ostringstream stream;
            stream << std::setfill('0') << std::setw(4) << rank;
            const std::string rank_string = stream.str();
            const std::string local_marker_file_name = marker_file_name + "." + rank_string;
            cat_file_stream << "rm -f " << local_marker_file_name << "\n";
        }

        // Setup commands to create a unified fiber file.
        cat_file_stream << "\n";
        for (int rank = 0; rank < mpi_size; ++rank)
        {
            std::ostringstream stream;
            stream << std::setfill('0') << std::setw(4) << rank;
            const std::string rank_string = stream.str();
            const std::string local_fiber_file_name = fiber_file_name + "." + rank_string;
            cat_file_stream << "m3D_hdf5_fiber_converter " << local_fiber_file_name + ".h5" << " " << local_fiber_file_name << "\n";
        }

        cat_file_stream << "\n";
        cat_file_stream << "cat " << fiber_header_file_name;
        for (int rank = 0; rank < mpi_size; ++rank)
        {
            std::ostringstream stream;
            stream << std::setfill('0') << std::setw(4) << rank;
            const std::string rank_string = stream.str();
            cat_file_stream << " " << fiber_file_name + "." + rank_string;
        }
        cat_file_stream << " > " << fiber_file_name << "\n";
        cat_file_stream << "\n";
        for (int rank = 0; rank < mpi_size; ++rank)
        {
            std::ostringstream stream;
            stream << std::setfill('0') << std::setw(4) << rank;
            const std::string rank_string = stream.str();
            const std::string local_fiber_file_name = fiber_file_name + "." + rank_string;
            cat_file_stream << "rm -f " << local_fiber_file_name << "\n";
        }
    }

    // Construct the VecScatter objects required to write the plot data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (d_build_vec_scatters[ln])
        {
            buildVecScatters(d_ao[ln], ln);
        }
        d_build_vec_scatters[ln] = false;
    }

    // Gather marker data.
    int max_marker_idx = 0;
    for (int cloud = 0; cloud < d_nclouds; ++cloud)
    {
        max_marker_idx = std::max(max_marker_idx,d_cloud_first_mark_idx[cloud]+d_cloud_nmarks[cloud]-1);
    }
    SAMRAI_MPI::maxReduction(&max_marker_idx);

    const int total_num_marks = max_marker_idx+1;
    std::vector<double> X_marker(NDIM*total_num_marks,std::numeric_limits<double>::max());
    if (d_mark_idx != -1)
    {
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<IndexData<NDIM,LagMarker,CellGeometry<NDIM> > > mark_data = patch->getPatchData(d_mark_idx);
                for (IndexData<NDIM,LagMarker,CellGeometry<NDIM> >::Iterator it(*mark_data); it; it++)
                {
                    const Index<NDIM>& i = it.getIndex();
                    if (patch_box.contains(i))
                    {
                        const LagMarker& mark = it();
                        const std::vector<double>& X = mark.getPositions();
                        const std::vector<int>& idx = mark.getIndices();
                        for (size_t k = 0; k < X.size()/NDIM; ++k)
                        {
                            for (int d = 0; d < NDIM; ++d)
                            {
                                /*!
                                 * \todo Add index error checking!
                                 */
                                X_marker[NDIM*idx[k]+d] = X[NDIM*k+d];
                            }
                        }
                    }
                }
            }
        }
    }
    /*!
     * \todo Make this operation more efficient!
     */
    SAMRAI_MPI::minReduction(&X_marker[0], NDIM*total_num_marks);

    // Create the HDF5 files.
    const std::string hdf5_marker_file_name = d_dump_directory_name + "/" + local_marker_file_name + ".h5";
    hid_t marker_file_id = H5Fcreate(hdf5_marker_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (marker_file_id < 0)
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                   << "  could not create marker file: " << marker_file_name << std::endl);
    }

    const std::string hdf5_fiber_file_name = d_dump_directory_name + "/" + local_fiber_file_name + ".h5";
    hid_t fiber_file_id = H5Fcreate(hdf5_fiber_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (fiber_file_id < 0)
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                   << "  could not create fiber file: " << fiber_file_name << std::endl);
    }

    // Store information about the data layout in the local HDF5 files.
    herr_t status;

    // Add the local clouds to the local marker file.
    hid_t marker_group_id = H5Gcreate1(marker_file_id, "/markers", 0);

    status = H5LTset_attribute_int(marker_file_id, "/markers", "num_local_marker_nodes", &num_local_marker_nodes, 1);
    status = H5LTset_attribute_int(marker_file_id, "/markers", "marker_node_offset", &marker_node_offset, 1);
    status = H5LTset_attribute_int(marker_file_id, "/markers", "num_marker_nodes", &num_marker_nodes, 1);

    status = H5LTset_attribute_int(marker_file_id, "/markers", "num_local_marker_clouds", &num_local_marker_clouds, 1);
    status = H5LTset_attribute_int(marker_file_id, "/markers", "marker_cloud_offset", &marker_cloud_offset, 1);
    status = H5LTset_attribute_int(marker_file_id, "/markers", "num_marker_clouds", &num_marker_clouds, 1);

    local_marker_node_counter = 0;
    local_marker_cloud_counter = 0;
    for (int cloud = 0; cloud < d_nclouds; ++cloud)
    {
        std::ostringstream dset_name_stream;
        dset_name_stream << "/markers/cloud_" << std::setw(4) << std::setfill('0') << marker_cloud_offset+local_marker_cloud_counter;
        const std::string dset_name = dset_name_stream.str();

        const int nmarks = d_cloud_nmarks[cloud];
        const int first_mark_idx = d_cloud_first_mark_idx[cloud];
        const double* const X = &X_marker[NDIM*first_mark_idx];
        const std::vector<float> buffer(X,X+NDIM*nmarks);
        const int node_offset = marker_node_offset+local_marker_node_counter;
        const int cloud_number = marker_cloud_offset+local_marker_cloud_counter;

        // Create the dataspace for the dataset.
        static const int rank = 2;
        hsize_t dims[rank] = { NDIM , nmarks };
        hid_t dataspace_id = H5Screate_simple(rank, dims, NULL);

        // Create the dataset with data compression enabled.
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t cdims[rank] = { NDIM , nmarks };
        status = H5Pset_chunk(plist_id, rank, cdims);
        status = H5Pset_deflate(plist_id, M3D_DEFLATE_LEVEL);
        hid_t dataset_id = H5Dcreate1(marker_file_id, dset_name.c_str(), H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT);

        // Write the data and related attributes to the dataset.
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer[0]);
        status = H5LTset_attribute_int(marker_file_id, dset_name.c_str(), "nmarks", &nmarks, 1);
        status = H5LTset_attribute_int(marker_file_id, dset_name.c_str(), "node_offset", &node_offset, 1);
        status = H5LTset_attribute_int(marker_file_id, dset_name.c_str(), "cloud_number", &cloud_number, 1);
        status = H5LTset_attribute_string(marker_file_id, dset_name.c_str(), "cloud_name", d_cloud_names[cloud].c_str());

        // Cleanup HDF5 data.
        status = H5Dclose(dataset_id);
        status = H5Pclose(plist_id);
        status = H5Sclose(dataspace_id);

        // Advance the counters.
        local_marker_node_counter += nmarks;
        local_marker_cloud_counter += 1;
    }
    status = H5Gclose(marker_group_id);
    status = H5Fclose(marker_file_id);

    // Add the local fibers to the local fiber file.
    hid_t fiber_group_id = H5Gcreate1(fiber_file_id, "/fibers", 0);

    status = H5LTset_attribute_int(fiber_file_id, "/fibers", "num_local_fibers", &num_local_fibers, 1);
    status = H5LTset_attribute_int(fiber_file_id, "/fibers", "fiber_offset", &fiber_offset, 1);
    status = H5LTset_attribute_int(fiber_file_id, "/fibers", "num_fibers", &num_fibers, 1);

    status = H5LTset_attribute_int(fiber_file_id, "/fibers", "num_local_groups", &num_local_groups, 1);
    status = H5LTset_attribute_int(fiber_file_id, "/fibers", "group_offset", &group_offset, 1);
    status = H5LTset_attribute_int(fiber_file_id, "/fibers", "num_groups", &num_groups, 1);

    status = H5LTset_attribute_int(fiber_file_id, "/fibers", "num_local_layers", &num_local_layers, 1);
    status = H5LTset_attribute_int(fiber_file_id, "/fibers", "layer_offset", &layer_offset, 1);
    status = H5LTset_attribute_int(fiber_file_id, "/fibers", "num_layers", &num_layers, 1);

    local_fiber_counter = 0;
    local_group_counter = 0;
    local_layer_counter = 0;
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (!d_coords_data[ln].isNull())
        {
            // Scatter the data from "global" to "local" form.
            Vec local_X_vec;
            ierr = VecDuplicate(d_dst_vec[ln][NDIM], &local_X_vec);
            IBTK_CHKERRQ(ierr);

            Vec global_X_vec = d_coords_data[ln]->getVec();
            ierr = VecScatterBegin(d_vec_scatter[ln][NDIM], global_X_vec, local_X_vec,
                                   INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterEnd(d_vec_scatter[ln][NDIM], global_X_vec, local_X_vec,
                                 INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);

            double* local_X_arr;
            ierr = VecGetArray(local_X_vec, &local_X_arr);
            IBTK_CHKERRQ(ierr);

            // Keep track of the current offset in the local Vec data.
            int offset = 0;

            // Add the local blocks to the local fiber file.
            for (int block = 0; block < d_nblocks[ln]; ++block)
            {
                std::ostringstream dset_name_stream;
                dset_name_stream << "/fibers/layer_" << std::setw(4) << std::setfill('0') << layer_offset+local_layer_counter;
                const std::string dset_name = dset_name_stream.str();

                const IntVector<NDIM>& nelem    = d_block_nelems  [ln][block];
                const IntVector<NDIM>& periodic = d_block_periodic[ln][block];
                const int ntot = nelem.getProduct();
                static const int rank = 4;
                hsize_t dims[rank] = { NDIM , nelem[0] , nelem[1] , nelem[2] };

                const double* const X = local_X_arr + NDIM*offset;
                const std::vector<float> buffer(X,X+NDIM*ntot);
                const int fiber_number = fiber_offset+local_fiber_counter;
                const int group_number = group_offset+local_group_counter;
                const int layer_number = layer_offset+local_layer_counter+1;

                // Create the dataspace for the dataset.
                hid_t dataspace_id = H5Screate_simple(rank, dims, NULL);

                // Create the dataset with data compression enabled.
                hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
                hsize_t cdims[rank] = { NDIM , nelem[0] , nelem[1] , nelem[2] };
                status = H5Pset_chunk(plist_id, rank, cdims);
                status = H5Pset_deflate(plist_id, M3D_DEFLATE_LEVEL);
                hid_t dataset_id = H5Dcreate1(fiber_file_id, dset_name.c_str(), H5T_NATIVE_FLOAT, dataspace_id, plist_id);

                // Write the data and related attributes to the dataset.
                status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer[0]);
                status = H5LTset_attribute_int(fiber_file_id, dset_name.c_str(), "nelem", &nelem[0], NDIM);
                status = H5LTset_attribute_int(fiber_file_id, dset_name.c_str(), "periodic", &periodic[0], NDIM);
                status = H5LTset_attribute_int(fiber_file_id, dset_name.c_str(), "fiber_number", &fiber_number, 1);
                status = H5LTset_attribute_int(fiber_file_id, dset_name.c_str(), "group_number", &group_number, 1);
                status = H5LTset_attribute_int(fiber_file_id, dset_name.c_str(), "layer_number", &layer_number, 1);
                status = H5LTset_attribute_int(fiber_file_id, dset_name.c_str(), "nfibers", &d_block_nfibers[ln][block], 1);
                status = H5LTset_attribute_int(fiber_file_id, dset_name.c_str(), "ngroups", &d_block_ngroups[ln][block], 1);
                status = H5LTset_attribute_string(fiber_file_id, dset_name.c_str(), "layer_name", d_block_names[ln][block].c_str());

                // Cleanup HDF5 data.
                status = H5Dclose(dataset_id);
                status = H5Pclose(plist_id);
                status = H5Sclose(dataspace_id);

                // Advance the counters.
                offset += ntot;
                local_fiber_counter += d_block_nfibers[ln][block];
                local_group_counter += d_block_ngroups[ln][block];
                local_layer_counter += 1;
            }

            // Clean up allocated data.
            ierr = VecRestoreArray(local_X_vec, &local_X_arr);  IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(local_X_vec);  IBTK_CHKERRQ(ierr);
        }
    }
    status = H5Gclose(fiber_group_id);
    status = H5Fclose(fiber_file_id);
    return;
}// writePlotData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LagM3DDataWriter::buildVecScatters(
    AO& ao,
    const int level_number)
{
    if (d_coords_data[level_number].isNull()) return;

    int ierr;

    // Setup the IS data used to generate the VecScatters that redistribute the
    // distributed data into local marker clouds, local logically Cartesian
    // blocks, and local UCD meshes.
    std::vector<int> ref_is_idxs;
    for (int block = 0; block < d_nblocks[level_number]; ++block)
    {
        const IntVector<NDIM>& nelem = d_block_nelems[level_number][block];
        const int ntot = nelem.getProduct();
        const int first_lag_idx = d_block_first_lag_idx[level_number][block];
        ref_is_idxs.reserve(ref_is_idxs.size()+ntot);

        for (int idx = first_lag_idx; idx < first_lag_idx + ntot; ++idx)
        {
            ref_is_idxs.push_back(idx);
        }
    }

    // Map Lagrangian indices to PETSc indices.
    std::vector<int> ao_dummy(1,-1);
    ierr = AOApplicationToPetsc(
        ao,
        (!ref_is_idxs.empty() ? int(ref_is_idxs.size()) : int(ao_dummy.size())),
        (!ref_is_idxs.empty() ? &ref_is_idxs[0]         : &ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    // Setup IS indices for all necessary data depths.
    std::map<int,std::vector<int> > src_is_idxs;

    src_is_idxs[NDIM].resize(ref_is_idxs.size());
    std::transform(ref_is_idxs.begin(), ref_is_idxs.end(),
                   src_is_idxs[NDIM].begin(),
                   std::bind2nd(std::multiplies<int>(),NDIM));
    d_src_vec[level_number][NDIM] = d_coords_data[level_number]->getVec();

    // Create the VecScatters to scatter data from the global PETSc Vec to
    // contiguous local subgrids.  VecScatter objects are individually created
    // for data depths as necessary.
    for (std::map<int,std::vector<int> >::iterator it = src_is_idxs.begin();
         it != src_is_idxs.end(); ++it)
    {
        const int depth = it->first;
        const std::vector<int>& idxs = it->second;

        IS src_is;
        ierr = ISCreateBlock(PETSC_COMM_WORLD, depth, idxs.size(),
                             &idxs[0], &src_is);
        IBTK_CHKERRQ(ierr);

        Vec& src_vec = d_src_vec[level_number][depth];
        Vec& dst_vec = d_dst_vec[level_number][depth];
        if (dst_vec)
        {
            ierr = VecDestroy(dst_vec);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecCreateMPI(PETSC_COMM_WORLD, depth*idxs.size(),
                            PETSC_DETERMINE, &dst_vec);
        IBTK_CHKERRQ(ierr);

        ierr = VecSetBlockSize(dst_vec, depth);
        IBTK_CHKERRQ(ierr);

        VecScatter& vec_scatter = d_vec_scatter[level_number][depth];
        if (vec_scatter)
        {
            ierr = VecScatterDestroy(vec_scatter);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecScatterCreate(src_vec, src_is, dst_vec, PETSC_NULL,
                                &vec_scatter);
        IBTK_CHKERRQ(ierr);

        ierr = ISDestroy(src_is);  IBTK_CHKERRQ(ierr);
    }
    return;
}// buildVecScatters

std::string
LagM3DDataWriter::getMarkerFileName(
    const int& timestep_number) const
{
    std::ostringstream stream;
    if (timestep_number < 999999)
    {
        stream << "_" << std::setfill('0') << std::setw(6) << timestep_number;
    }
    else
    {
        stream << std::setfill('0') << std::setw(7) << timestep_number;
    }
    return d_file_prefix + stream.str() + ".mk";
}// getMarkerFileName

std::string
LagM3DDataWriter::getFiberFileName(
    const int& timestep_number) const
{
    std::ostringstream stream;
    if (timestep_number < 999999)
    {
        stream << "_" << std::setfill('0') << std::setw(6) << timestep_number;
    }
    else
    {
        stream << std::setfill('0') << std::setw(7) << timestep_number;
    }
    return d_file_prefix + stream.str() + ".xf";
}// getFiberFileName

std::string
LagM3DDataWriter::getMenuFileName() const
{
    return "menu." + d_file_prefix;
}// getMenuFileName

std::string
LagM3DDataWriter::getListFileName() const
{
    return "list." + d_file_prefix;
}// getListFileName

std::string
LagM3DDataWriter::getCatScriptFileName() const
{
    return "cat." + d_file_prefix + ".sh";
}// getCatScriptFileName

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::LagM3DDataWriter>;

//////////////////////////////////////////////////////////////////////////////
