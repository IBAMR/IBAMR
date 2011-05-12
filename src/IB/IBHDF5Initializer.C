// Filename: IBHDF5Initializer.C
// Created on 26 Sep 2006 by Boyce Griffith
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

#include "IBHDF5Initializer.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/IBBeamForceSpec.h>
#include <ibamr/IBInstrumentationSpec.h>
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IndexUtilities.h>
#include <ibtk/LNodeIndexData.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <Index.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <fstream>
#include <iostream>

// HDF5 INCLUDES
#include <hdf5.h>
#if (H5_VERS_MINOR == 6)
#include <H5LT.h>
#define H5Dopen1 H5Dopen
#endif
#if (H5_VERS_MINOR == 8)
#include <hdf5_hl.h>
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int BUFFER_SIZE = 2097152;  // = 16 * 2^20 / 8 ~ 16 MB of double precision values
static const int STRING_BUFFER_SIZE = 256;

static const int NUM_GROUPS = 6;
static const int VERTEX_GROUP          = 0;
static const int SPRING_GROUP          = 1;
static const int BEAM_GROUP            = 2;
static const int TARGET_POINT_GROUP    = 3;
static const int MASS_POINT_GROUP      = 4;
static const int INSTRUMENTATION_GROUP = 5;

herr_t
file_info(
    hid_t loc_id,
    const char* name,
    void* op_data)
{
    std::vector<bool>& has_group = *(static_cast<std::vector<bool>*>(op_data));

    H5G_stat_t statbuf;
    H5Gget_objinfo(loc_id, name, 0, &statbuf);
    if (statbuf.type == H5G_GROUP)
    {
        if (strcmp(name,"vertex") == 0)
        {
            has_group[VERTEX_GROUP] = true;
        }
        else if (strcmp(name,"spring") == 0)
        {
            has_group[SPRING_GROUP] = true;
        }
        else if (strcmp(name,"beam") == 0)
        {
            has_group[BEAM_GROUP] = true;
        }
        else if (strcmp(name,"target_point") == 0)
        {
            has_group[TARGET_POINT_GROUP] = true;
        }
        else if (strcmp(name,"mass_point") == 0)
        {
            has_group[MASS_POINT_GROUP] = true;
        }
        else if (strcmp(name,"instrumentation") == 0)
        {
            has_group[INSTRUMENTATION_GROUP] = true;
        }
    }
    return 0;
}// file_info
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBHDF5Initializer::IBHDF5Initializer(
    const std::string& object_name,
    Pointer<Database> input_db)
    : d_object_name(object_name),
      d_use_file_batons(true),
      d_max_levels(-1),
      d_level_is_initialized(),
      d_filenames(),
      d_enable_springs(),
      d_using_uniform_spring_stiffness(),
      d_uniform_spring_stiffness(),
      d_using_uniform_spring_rest_length(),
      d_uniform_spring_rest_length(),
      d_using_uniform_spring_force_fcn_idx(),
      d_uniform_spring_force_fcn_idx(),
      d_enable_beams(),
      d_using_uniform_beam_bend_rigidity(),
      d_uniform_beam_bend_rigidity(),
      d_enable_target_points(),
      d_using_uniform_target_stiffness(),
      d_uniform_target_stiffness(),
      d_using_uniform_target_damping(),
      d_uniform_target_damping(),
      d_enable_instrumentation(),
      d_instrument_names(),
      d_cache_level_number(-1),
      d_level_num_vertex(),
      d_level_num_local_vertex(),
      d_level_vertex_offset(),
      d_level_posns(),
      d_level_vertex_idxs(),
      d_level_cell_idxs(),
      d_level_patch_nums(),
      d_level_reset_specs_set(),
      d_level_num_spring(),
      d_level_num_local_spring(),
      d_level_spring_data_map(),
      d_level_num_beam(),
      d_level_num_local_beam(),
      d_level_beam_data_map(),
      d_level_num_target_point(),
      d_level_num_local_target_point(),
      d_level_target_point_data_map(),
      d_level_num_inst_point(),
      d_level_num_local_inst_point(),
      d_level_inst_point_data_map()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
#endif
    // Register the specification objects with the StreamableManager class.
    IBSpringForceSpec::registerWithStreamableManager();
    IBBeamForceSpec::registerWithStreamableManager();
    IBTargetPointForceSpec::registerWithStreamableManager();
    IBInstrumentationSpec::registerWithStreamableManager();

    // Initialize object with data read from the input database.
    getFromInput(input_db);
    return;
}// IBHDF5Initializer

IBHDF5Initializer::~IBHDF5Initializer()
{
    // intentionally blank
    return;
}// ~IBHDF5Initializer

bool
IBHDF5Initializer::getLevelHasLagrangianData(
    const int level_number,
    const bool can_be_refined) const
{
    return !d_filenames[level_number].empty();
}// getLevelHasLagrangianData

int
IBHDF5Initializer::computeLocalNodeCountOnPatchLevel(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    if (!can_be_refined && level_number != d_max_levels-1)
    {
        TBOX_WARNING(d_object_name << "::computeLocalNodeCountOnPatchLevel():\n  Input database key `max_levels' = " << d_max_levels << " but finest level in the hierarchy is " << level_number << "\n"
                     << "  some Lagrangian data may not be properly initialized.\n");
    }

    buildLevelDataCache(hierarchy, level_number, init_data_time, can_be_refined, initial_time);
    return std::accumulate(d_level_num_local_vertex.begin(),d_level_num_local_vertex.end(),0);
}// computeLocalNodeCountOnPatchLevel

int
IBHDF5Initializer::initializeDataOnPatchLevel(
    const int lag_node_index_idx,
    const int global_index_offset,
    const int local_index_offset,
    Pointer<LMeshData>& X_data,
    Pointer<LMeshData>& U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    if (!can_be_refined && level_number != d_max_levels-1)
    {
        TBOX_WARNING(d_object_name << "::initializeDataOnPatchLevel():\n  Input database key `max_levels' = " << d_max_levels << " but finest level in the hierarchy is " << level_number << "\n"
                     << "  some Lagrangian data may not be properly initialized.\n");
    }

    buildLevelDataCache(hierarchy, level_number, init_data_time, can_be_refined, initial_time);

    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const gridXLower = grid_geom->getXLower();
    const double* const gridXUpper = grid_geom->getXUpper();

    // Loop over all vertices in the specified level and initialize the data in
    // the appropriate Cartesian grid patches.
    blitz::Array<double,2>& X_array = *X_data->getGhostedLocalFormVecArray();
    blitz::Array<double,2>& U_array = *U_data->getGhostedLocalFormVecArray();
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const int num_filenames = d_filenames[level_number].size();
    int local_idx = -1;
    int local_node_count = 0;
    for (int j = 0; j < num_filenames; ++j)
    {
        for (int k = 0; k < d_level_num_local_vertex[j]; ++k, ++local_node_count)
        {
            const std::vector<double>& X = d_level_posns[j][k];
            const std::pair<int,int>& vertex_idx = d_level_vertex_idxs[j][k];
            const Index<NDIM>& i = d_level_cell_idxs[j][k];
            const int patch_num = d_level_patch_nums[j][k];

            // Compute the index information for the present vertex.
            const int lagrangian_idx = getCanonicalLagrangianIndex(vertex_idx,global_index_offset);
            const int local_petsc_idx = ++local_idx + local_index_offset;
            const int global_petsc_idx = local_petsc_idx+global_index_offset;

            // Ensure the point lies within the physical domain.
            for (int d = 0; d < NDIM; ++d)
            {
                if (X[d] <= gridXLower[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                               << "  encountered node below lower physical boundary\n"
                               << "  level number = " << level_number << "\n"
                               << "  file name = " << d_filenames[level_number][j] << "\n"
                               << "  vertex index = " << k << "\n"
                               << "  please ensure that all nodes are within the computational domain.\n");
                }

                if (X[d] >= gridXUpper[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                               << "  encountered node above upper physical boundary\n"
                               << "  level number = " << level_number << "\n"
                               << "  file name = " << d_filenames[level_number][j] << "\n"
                               << "  vertex index = " << k << "\n"
                               << "  please ensure that all nodes are within the computational domain.\n");
                }
            }

            // Ensure the point lies with the present grid patch.
            Pointer<Patch<NDIM> > patch = level->getPatch(patch_num);
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patchXLower = patch_geom->getXLower();
            const double* const patchXUpper = patch_geom->getXUpper();
            for (int d = 0; d < NDIM; ++d)
            {
                if (X[d] < patchXLower[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                               << "  encountered node below lower patch boundary\n"
                               << "  level number = " << level_number << "\n"
                               << "  file name = " << d_filenames[level_number][j] << "\n"
                               << "  vertex index = " << k << "\n");
                }
                else if (X[d] >= patchXUpper[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                               << "  encountered node above upper patch boundary\n"
                               << "  level number = " << level_number << "\n"
                               << "  file name = " << d_filenames[level_number][j] << "\n"
                               << "  vertex index = " << k << "\n");
                }
            }

            // Initialize the position of the present vertex.
            std::copy(X.begin(),X.end(),&X_array(local_petsc_idx,0));

            // Initialize the velocity of the present vertex.
            std::fill(&U_array(local_petsc_idx,0),&U_array(local_petsc_idx,0)+NDIM,0.0);

            // Initialize the specification objects associated with the present
            // vertex.
            std::vector<Pointer<Streamable> > vertex_specs = initializeSpecs(std::make_pair(j,k), vertex_idx, global_index_offset);

            // Initialize the LNodeIndex data.
            const Box<NDIM>& patch_box = patch->getBox();
            if (!patch_box.contains(i))
            {
                TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                           << "  encountered node assigned to incorrect patch\n"
                           << "  level number = " << level_number << "\n"
                           << "  file name = " << d_filenames[level_number][j] << "\n"
                           << "  vertex index = " << k << "\n"
                           << "  vertex cell index = " << i << "\n"
                           << "  assigned patch number = " << patch_num << "\n"
                           << "  assigned patch box = " << patch_box << "\n");
            }
            Pointer<LNodeIndexData> index_data = patch->getPatchData(lag_node_index_idx);
            if (!index_data->isElement(i))
            {
                index_data->appendItemPointer(i, new LNodeIndexSet());
            }
            LNodeIndexSet* const node_set = index_data->getItem(i);
            const IntVector<NDIM> periodic_offset(0);
            const std::vector<double> periodic_displacement(NDIM,0.0);
            node_set->push_back(new LNodeIndex(lagrangian_idx, global_petsc_idx, local_petsc_idx, &X_array(local_petsc_idx,0), periodic_offset, periodic_displacement, vertex_specs));
        }
    }
    X_data->restoreArrays();
    U_data->restoreArrays();

    // Sanity check.
    const int expected_local_node_count = std::accumulate(d_level_num_local_vertex.begin(),d_level_num_local_vertex.end(),0);
    if (local_node_count != expected_local_node_count)
    {
        TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                   << "  expected local node count = " << expected_local_node_count << "\n"
                   << "  actual   local node count = " << local_node_count << "\n");
    }

    // Setup the instrument names.
    std::vector<std::string> all_instrument_names;
    for (int ln = 0; ln <= level_number; ++ln)
    {
        for (int j = 0; j < int(d_instrument_names[ln].size()); ++j)
        {
            all_instrument_names.insert(
                all_instrument_names.end(),
                d_instrument_names[ln][j].begin(), d_instrument_names[ln][j].end());
        }
    }
    IBInstrumentationSpec::setInstrumentNames(all_instrument_names);

    // Indicate that the level is initialized.
    d_level_is_initialized[level_number] = true;
    return local_node_count;
}// initializeDataOnPatchLevel

int
IBHDF5Initializer::initializeMassDataOnPatchLevel(
    const int global_index_offset,
    const int local_index_offset,
    Pointer<LMeshData>& M_data,
    Pointer<LMeshData>& K_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    TBOX_ERROR(d_object_name << "::initializeMassDataOnPatchLevel():\n  Not implemented.\n");
    int local_node_count = 0;
    return local_node_count;
}// initializeMassOnPatchLevel

void
IBHDF5Initializer::tagCellsForInitialRefinement(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    Pointer<PatchLevel<NDIM> > tag_level = hierarchy->getPatchLevel(level_number);

    // Tag cells for refinement whenever there are vertices whose initial
    // locations lie within the index space of the given patch, but on the finer
    // levels of the AMR patch hierarchy.
    std::vector<Index<NDIM> > cell_idxs;
    std::vector<int> patch_nums;
    for (int ln = level_number+1; ln < d_max_levels; ++ln)
    {
        std::vector<Index<NDIM> > level_cell_idxs;
        std::vector<int> level_patch_nums;
        findLocalPatchIndices(level_cell_idxs, level_patch_nums, ln, tag_level);
        cell_idxs.insert(cell_idxs.end(),level_cell_idxs.begin(),level_cell_idxs.end());
        patch_nums.insert(patch_nums.end(),level_patch_nums.begin(),level_patch_nums.end());
    }

    const int num_vertex = cell_idxs.size();
    for (int k = 0; k < num_vertex; ++k)
    {
        const Index<NDIM>& i = cell_idxs[k];
        const int patch_num = patch_nums[k];
        Pointer<CellData<NDIM,int> > tag_data = tag_level->getPatch(patch_num)->getPatchData(tag_index);
        const Box<NDIM>& patch_box = tag_data->getBox();
        if (!patch_box.contains(i))
        {
            TBOX_ERROR(d_object_name << "::tagCellsForInitialRefinement():\n"
                       << "  encountered node assigned to incorrect patch\n"
                       << "  level number = " << level_number << "\n"
                       << "  vertex cell index = " << i << "\n"
                       << "  assigned patch number = " << patch_num << "\n"
                       << "  assigned patch box = " << patch_box << "\n");
        }
        (*tag_data)(i) = 1;
    }
    return;
}// tagCellsForInitialRefinement

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBHDF5Initializer::findLocalPatchIndices(
    std::vector<Index<NDIM> >& cell_idxs,
    std::vector<int>& patch_nums,
    const int level_number,
    const Pointer<PatchLevel<NDIM> > level) const
{
    cell_idxs.clear();
    patch_nums.clear();

    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    const int num_filenames = d_filenames[level_number].size();
    for (int j = 0; j < num_filenames; ++j)
    {
        // Wait for the previous MPI process to finish reading the current file.
        if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

        std::string filename = d_filenames[level_number][j];
        const std::string postfix = ".h5";
        size_t pos = filename.find(postfix);
        if (pos != std::string::npos)
        {
            filename.swap(filename.erase(pos,postfix.length()));
        }
        const std::string base_group_name = "/" + filename;
        filename += postfix;

        hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Unable to open input file " << filename << "\n");
        }
        else
        {
            plog << d_object_name << ":  "
                 << "processing vertex data from input filename " << filename << "\n"
                 << "  on MPI process " << SAMRAI_MPI::getRank() << "\n";
        }

        std::vector<Index<NDIM> > file_cell_idxs;
        std::vector<int> file_patch_nums;
        findLocalPatchIndicesFromHDF5(
            file_cell_idxs, file_patch_nums,
            file_id, base_group_name, level, filename);
        cell_idxs.insert(cell_idxs.end(),file_cell_idxs.begin(),file_cell_idxs.end());
        patch_nums.insert(patch_nums.end(),file_patch_nums.begin(),file_patch_nums.end());

        H5Fclose(file_id);

        // Free the next MPI process to start reading the current file.
        if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
}// findLocalPatchIndices

void
IBHDF5Initializer::findLocalPatchIndicesFromHDF5(
    std::vector<Index<NDIM> >& cell_idxs,
    std::vector<int>& patch_nums,
    const hid_t file_id,
    const std::string& base_group_name,
    const Pointer<PatchLevel<NDIM> > level,
    const std::string& filename) const
{
    cell_idxs.clear();
    patch_nums.clear();

    // Check to see if the group for the dataset(s) exists.
    const std::string vertex_group_name = base_group_name + "/vertex";
    if (H5Gget_objinfo(file_id, vertex_group_name.c_str(), 0, NULL) >= 0)
    {
        // Open the dataset.
        const std::string posn_dset_name = vertex_group_name + "/posn";
        hid_t posn_dset = H5Dopen1(file_id, posn_dset_name.c_str());
        if (posn_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required vertex dataset in input file " << filename << "\n");
        }

        // Read in the dimensions of the dataset.
        int rank;
        hsize_t dims[2];
        H5T_class_t class_id;
        size_t type_size;

        H5LTget_dataset_ndims(file_id, posn_dset_name.c_str(), &rank);
        if (rank != 2)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid vertex dataset rank in input file " << filename << "\n");
        }

        H5LTget_dataset_info(file_id, posn_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0 || dims[1] != NDIM)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid vertex dataset dimension in input file " << filename << "\n");
        }
        const int num_vertex = int(dims[0]);

        // Define the file dataspace.
        static const int rankf = 2;
        hsize_t dimsf[rankf] = { num_vertex , NDIM };
        hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

        // Define the memory dataspace.
        static const int rankm = 2;
        hsize_t dimsm[rankm] = { BUFFER_SIZE , NDIM };
        hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

        // Read in the vertex data one block at a time.
        std::vector<double> posn_buf(NDIM*BUFFER_SIZE);
        const int num_blocks = num_vertex/BUFFER_SIZE + (num_vertex%BUFFER_SIZE == 0 ? 0 : 1);
        for (int block = 0; block < num_blocks; ++block)
        {
            // Determine whether we are reading in the last block in the file.
            const bool last_block = (block == num_blocks-1);

            // Determine the number of items to read (always BUFFER_SIZE except
            // for the final block in the file).
            const int num_vertex_block = (last_block ? num_vertex - block*BUFFER_SIZE : BUFFER_SIZE);
            TBOX_ASSERT(num_vertex_block > 0 && num_vertex_block <= BUFFER_SIZE);

            // Define the file hyperslab.
            hsize_t offsetf[rankf] = { block*BUFFER_SIZE , 0 };
            hsize_t countf[rankf] = { num_vertex_block , NDIM };
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

            // Define the memory hyperslab.
            hsize_t offsetm[rankm] = { 0 , 0 };
            hsize_t countm[rankm] = { num_vertex_block , NDIM };
            H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

            // Read data from hyperslab in the file into the hyperslab in
            // memory.
            H5Dread(posn_dset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &posn_buf[0]);

            // Setup cell indices for any local vertices in the hyperslab.
            for (int k = 0; k < num_vertex_block; ++k)
            {
                const double* const X = &posn_buf[NDIM*k];
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
                        patch->getPatchGeometry();
                    const double* const xLower = patch_geom->getXLower();
                    const double* const xUpper = patch_geom->getXUpper();
                    const bool patch_owns_node =
                        ((  xLower[0] <= X[0])&&(X[0] < xUpper[0]))
#if (NDIM > 1)
                        &&((xLower[1] <= X[1])&&(X[1] < xUpper[1]))
#if (NDIM > 2)
                        &&((xLower[2] <= X[2])&&(X[2] < xUpper[2]))
#endif
#endif
                        ;
                    if (patch_owns_node)
                    {
                        const double* const dx = patch_geom->getDx();
                        const Box<NDIM>& patch_box = patch->getBox();
                        const Index<NDIM>& patch_lower = patch_box.lower();
                        const Index<NDIM>& patch_upper = patch_box.upper();
                        const Index<NDIM> i = IndexUtilities::getCellIndex(
                            X, xLower, xUpper, dx, patch_lower, patch_upper);
                        cell_idxs.push_back(i);
                        patch_nums.push_back(p());
                        break;
                    }
                }
            }
        }

        // Cleanup HDF5 data structures.
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Dclose(posn_dset);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":\n  Cannot find vertex group in input file " << filename << "\n"
                   << "       base group name = " << base_group_name << "\n"
                   << "       vertex group name = " << vertex_group_name << "\n");
    }
    return;
}// findLocalPatchIndicesFromHDF5

void
IBHDF5Initializer::buildLevelDataCache(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    if (d_cache_level_number != level_number)
    {
        clearLevelDataCache();
    }
    else if (d_cache_level_number == level_number)
    {
        return;
    }

    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    const int num_filenames = d_filenames[level_number].size();

    d_instrument_names[level_number].resize(num_filenames);

    d_cache_level_number = level_number;
    d_level_num_vertex            .resize(num_filenames,0);
    d_level_num_local_vertex      .resize(num_filenames,0);
    d_level_vertex_offset         .resize(num_filenames,0);
    d_level_posns                 .resize(num_filenames);
    d_level_vertex_idxs           .resize(num_filenames);
    d_level_cell_idxs             .resize(num_filenames);
    d_level_patch_nums            .resize(num_filenames);
    d_level_reset_specs_set       .resize(num_filenames);
    d_level_num_spring            .resize(num_filenames,0);
    d_level_num_local_spring      .resize(num_filenames,0);
    d_level_spring_data_map       .resize(num_filenames);
    d_level_num_beam              .resize(num_filenames,0);
    d_level_num_local_beam        .resize(num_filenames,0);
    d_level_beam_data_map         .resize(num_filenames);
    d_level_num_target_point      .resize(num_filenames,0);
    d_level_num_local_target_point.resize(num_filenames,0);
    d_level_target_point_data_map .resize(num_filenames);
    d_level_num_inst_point        .resize(num_filenames,0);
    d_level_num_local_inst_point  .resize(num_filenames,0);
    d_level_inst_point_data_map   .resize(num_filenames);

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (int j = 0; j < num_filenames; ++j)
    {
        // Wait for the previous MPI process to finish reading the current file.
        if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

        std::string filename = d_filenames[level_number][j];
        const std::string postfix = ".h5";
        size_t pos = filename.find(postfix);
        if (pos != std::string::npos)
        {
            filename.swap(filename.erase(pos,postfix.length()));
        }
        const std::string base_group_name = "/" + filename;
        filename += postfix;

        hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Unable to open input file " << filename << "\n");
        }
        else
        {
            plog << d_object_name << ":  "
                 << "processing vertex data from input filename " << filename << "\n"
                 << "  on MPI process " << SAMRAI_MPI::getRank() << "\n";
        }

        // Check the file contents.
        std::vector<bool> has_group(NUM_GROUPS,false);
        H5Giterate(file_id, base_group_name.c_str(), NULL, file_info, static_cast<void*>(&has_group));

        std::set<int> local_vertex_idx_set;

        if (has_group[VERTEX_GROUP])
        {
            int& num_vertex       = d_level_num_vertex      [j];
            int& num_local_vertex = d_level_num_local_vertex[j];

            std::vector<std::vector<double> >&       posns       = d_level_posns      [j];
            std::vector<std::pair<int,int> >&        vertex_idxs = d_level_vertex_idxs[j];
            std::vector<Index<NDIM> >& cell_idxs   = d_level_cell_idxs  [j];
            std::vector<int>&                        patch_nums  = d_level_patch_nums [j];

            if (j == 0)
            {
                d_level_vertex_offset[j] = 0;
            }
            else
            {
                d_level_vertex_offset[j] = d_level_vertex_offset[j-1]+d_level_num_vertex[j-1];
            }

            buildLevelVertexDataCacheFromHDF5(
                num_vertex, num_local_vertex, posns, vertex_idxs, cell_idxs, patch_nums, local_vertex_idx_set,
                file_id, base_group_name, level, filename, j, num_filenames);
        }
        else
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find vertex group in input file " << filename << "\n");
        }

        if (d_enable_springs[level_number][j] && has_group[SPRING_GROUP])
        {
            int& num_spring       = d_level_num_spring      [j];
            int& num_local_spring = d_level_num_local_spring[j];
            std::map<int,Pointer<IBSpringForceSpec> >& spring_data_map = d_level_spring_data_map[j];

            buildLevelSpringDataCacheFromHDF5(
                num_spring, num_local_spring, spring_data_map,
                local_vertex_idx_set, file_id, base_group_name, level, filename, j, num_filenames);
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n  Cannot find spring group in input file " << filename << "\n");
        }

        if (d_enable_beams[level_number][j] && has_group[BEAM_GROUP])
        {
            int& num_beam       = d_level_num_beam      [j];
            int& num_local_beam = d_level_num_local_beam[j];
            std::map<int,Pointer<IBBeamForceSpec> >& beam_data_map = d_level_beam_data_map[j];

            buildLevelBeamDataCacheFromHDF5(
                num_beam, num_local_beam, beam_data_map,
                local_vertex_idx_set, file_id, base_group_name, level, filename, j, num_filenames);
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n  Cannot find beam group in input file " << filename << "\n");
        }

        if (d_enable_target_points[level_number][j] && has_group[TARGET_POINT_GROUP])
        {
            int& num_target_point       = d_level_num_target_point      [j];
            int& num_local_target_point = d_level_num_local_target_point[j];
            std::map<int,Pointer<IBTargetPointForceSpec> >& target_point_data_map = d_level_target_point_data_map[j];

            buildLevelTargetPointDataCacheFromHDF5(
                num_target_point, num_local_target_point, target_point_data_map,
                local_vertex_idx_set, file_id, base_group_name, level, filename, j, num_filenames);
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n  Cannot find target point group in input file " << filename << "\n");
        }

        if (d_enable_instrumentation[level_number][j] && has_group[INSTRUMENTATION_GROUP])
        {
            std::vector<std::string>& instrument_names = d_instrument_names[level_number][j];
            int& num_inst_point       = d_level_num_inst_point      [j];
            int& num_local_inst_point = d_level_num_local_inst_point[j];
            std::map<int,Pointer<IBInstrumentationSpec> >& inst_point_data_map = d_level_inst_point_data_map[j];

            buildLevelInstrumentationDataCacheFromHDF5(
                instrument_names, num_inst_point, num_local_inst_point, inst_point_data_map,
                local_vertex_idx_set, file_id, base_group_name, level, filename, j, num_filenames);
        }
        else
        {
            TBOX_WARNING(d_object_name << ":\n  Cannot find instrumentation group in input file " << filename << "\n");
        }

        H5Fclose(file_id);

        // Free the next MPI process to start reading the current file.
        if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
    }

    // Sanity checks.
    const int level_num_vertex_sum = std::accumulate(d_level_num_vertex.begin(),d_level_num_vertex.end(),0);
    const int level_num_local_vertex_sum = std::accumulate(d_level_num_local_vertex.begin(),d_level_num_local_vertex.end(),0);
    TBOX_ASSERT(level_num_vertex_sum == SAMRAI_MPI::sumReduction(level_num_local_vertex_sum));

    const int level_num_spring_sum = std::accumulate(d_level_num_spring.begin(),d_level_num_spring.end(),0);
    const int level_num_local_spring_sum = std::accumulate(d_level_num_local_spring.begin(),d_level_num_local_spring.end(),0);
    TBOX_ASSERT(level_num_spring_sum == SAMRAI_MPI::sumReduction(level_num_local_spring_sum));

    const int level_num_beam_sum = std::accumulate(d_level_num_beam.begin(),d_level_num_beam.end(),0);
    const int level_num_local_beam_sum = std::accumulate(d_level_num_local_beam.begin(),d_level_num_local_beam.end(),0);
    TBOX_ASSERT(level_num_beam_sum == SAMRAI_MPI::sumReduction(level_num_local_beam_sum));

    const int level_num_target_point_sum = std::accumulate(d_level_num_target_point.begin(),d_level_num_target_point.end(),0);
    const int level_num_local_target_point_sum = std::accumulate(d_level_num_local_target_point.begin(),d_level_num_local_target_point.end(),0);
    TBOX_ASSERT(level_num_target_point_sum == SAMRAI_MPI::sumReduction(level_num_local_target_point_sum));

    const int level_num_inst_point_sum = std::accumulate(d_level_num_inst_point.begin(),d_level_num_inst_point.end(),0);
    const int level_num_local_inst_point_sum = std::accumulate(d_level_num_local_inst_point.begin(),d_level_num_local_inst_point.end(),0);
    TBOX_ASSERT(level_num_inst_point_sum == SAMRAI_MPI::sumReduction(level_num_local_inst_point_sum));
    return;
}// buildLevelDataCache

void
IBHDF5Initializer::buildLevelVertexDataCacheFromHDF5(
    int& num_vertex,
    int& num_local_vertex,
    std::vector<std::vector<double> >& posns,
    std::vector<std::pair<int,int> >& vertex_idxs,
    std::vector<Index<NDIM> >& cell_idxs,
    std::vector<int>& patch_nums,
    std::set<int>& local_vertex_idx_set,
    const hid_t file_id,
    const std::string& base_group_name,
    const Pointer<PatchLevel<NDIM> > level,
    const std::string& filename,
    const int file_number,
    const int num_files) const
{
    num_vertex       = 0;
    num_local_vertex = 0;
    posns     .clear();
    cell_idxs .clear();
    patch_nums.clear();

    // Check to see if the group for the dataset(s) exists.
    const std::string vertex_group_name = base_group_name + "/vertex";
    if (H5Gget_objinfo(file_id, vertex_group_name.c_str(), 0, NULL) >= 0)
    {
        // Open the dataset.
        const std::string posn_dset_name = vertex_group_name + "/posn";
        hid_t posn_dset = H5Dopen1(file_id, posn_dset_name.c_str());
        if (posn_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required vertex dataset in input file " << filename << "\n");
        }

        // Read in the dimensions of the datasets.
        int rank;
        hsize_t dims[2];
        H5T_class_t class_id;
        size_t type_size;

        H5LTget_dataset_ndims(file_id, posn_dset_name.c_str(), &rank);
        if (rank != 2)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid vertex dataset rank in input file " << filename << "\n");
        }

        H5LTget_dataset_info(file_id, posn_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0 || dims[1] != NDIM)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid vertex dataset dimension in input file " << filename << "\n");
        }
        num_vertex = int(dims[0]);

        // Define the file dataspace.
        static const int rankf = 2;
        hsize_t dimsf[rankf] = { num_vertex , NDIM };
        hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

        // Define the memory dataspace.
        static const int rankm = 2;
        hsize_t dimsm[rankm] = { BUFFER_SIZE , NDIM };
        hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

        // Read in the vertex data one block at a time.
        std::vector<double> posn_buf(NDIM*BUFFER_SIZE);
        const int num_blocks = num_vertex/BUFFER_SIZE + (num_vertex%BUFFER_SIZE == 0 ? 0 : 1);
        for (int block = 0; block < num_blocks; ++block)
        {
            // Determine whether we are reading in the last block in the file.
            const bool last_block = (block == num_blocks-1);

            // Determine the number of items to read (always BUFFER_SIZE except
            // for the final block in the file).
            const int num_vertex_block = (last_block ? num_vertex - block*BUFFER_SIZE : BUFFER_SIZE);
            TBOX_ASSERT(num_vertex_block > 0 && num_vertex_block <= BUFFER_SIZE);

            // Define the file hyperslab.
            hsize_t offsetf[rankf] = { block*BUFFER_SIZE , 0 };
            hsize_t countf[rankf] = { num_vertex_block , NDIM };
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

            // Define the memory hyperslab.
            hsize_t offsetm[rankm] = { 0 , 0 };
            hsize_t countm[rankm] = { num_vertex_block , NDIM };
            H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

            // Read data from hyperslab in the file into the hyperslab in memory.
            H5Dread(posn_dset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &posn_buf[0]);

            // Setup data for all local vertices in the hyperslab.
            const int index_offset = block*BUFFER_SIZE;
            for (int k = 0; k < num_vertex_block; ++k)
            {
                const double* const X = &posn_buf[NDIM*k];
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
                        patch->getPatchGeometry();
                    const double* const xLower = patch_geom->getXLower();
                    const double* const xUpper = patch_geom->getXUpper();
                    const bool patch_owns_node =
                        ((  xLower[0] <= X[0])&&(X[0] < xUpper[0]))
#if (NDIM > 1)
                        &&((xLower[1] <= X[1])&&(X[1] < xUpper[1]))
#if (NDIM > 2)
                        &&((xLower[2] <= X[2])&&(X[2] < xUpper[2]))
#endif
#endif
                        ;
                    if (patch_owns_node)
                    {
                        ++num_local_vertex;

                        const double* const dx = patch_geom->getDx();
                        const Box<NDIM>& patch_box = patch->getBox();
                        const Index<NDIM>& patch_lower = patch_box.lower();
                        const Index<NDIM>& patch_upper = patch_box.upper();
                        const Index<NDIM> i = IndexUtilities::getCellIndex(
                            X, xLower, xUpper, dx, patch_lower, patch_upper);

                        const int index = k+index_offset;
                        posns.push_back(std::vector<double>(X,X+NDIM));
                        vertex_idxs.push_back(std::make_pair(file_number,index));
                        cell_idxs.push_back(i);
                        patch_nums.push_back(p());

                        local_vertex_idx_set.insert(local_vertex_idx_set.end(),index);

                        break;
                    }
                }
            }
        }

        // Cleanup HDF5 data structures.
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Dclose(posn_dset);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":\n  Cannot find vertex group in input file " << filename << "\n"
                   << "       base group name = " << base_group_name << "\n"
                   << "       vertex group name = " << vertex_group_name << "\n");
    }
    return;
}// buildLevelVertexDataCacheFromHDF5

void
IBHDF5Initializer::buildLevelSpringDataCacheFromHDF5(
    int& num_spring,
    int& num_local_spring,
    std::map<int,Pointer<IBSpringForceSpec> >& spring_data_map,
    const std::set<int>& local_vertex_idx_set,
    const hid_t file_id,
    const std::string& base_group_name,
    const Pointer<PatchLevel<NDIM> > level,
    const std::string& filename,
    const int file_number,
    const int num_files) const
{
    num_spring       = 0;
    num_local_spring = 0;
    spring_data_map.clear();

    // Check to see if the group for the dataset(s) exists.
    const std::string spring_group_name = base_group_name + "/spring";
    if (H5Gget_objinfo(file_id, spring_group_name.c_str(), 0, NULL) >= 0)
    {
        // Open the datasets.
        const std::string node1_idx_dset_name     = spring_group_name + "/node1_idx"    ;
        const std::string node2_idx_dset_name     = spring_group_name + "/node2_idx"    ;
        const std::string force_fcn_idx_dset_name = spring_group_name + "/force_fcn_idx";
        const std::string stiffness_dset_name     = spring_group_name + "/stiffness"    ;
        const std::string rest_length_dset_name   = spring_group_name + "/rest_length"  ;

        hid_t node1_idx_dset = H5Dopen1(file_id, node1_idx_dset_name.c_str());
        if (node1_idx_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required spring dataset in input file " << filename << "\n");
        }

        hid_t node2_idx_dset = H5Dopen1(file_id, node2_idx_dset_name.c_str());
        if (node2_idx_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required spring dataset in input file " << filename << "\n");
        }

        hid_t force_fcn_idx_dset = H5Dopen1(file_id, force_fcn_idx_dset_name.c_str());
        if (force_fcn_idx_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required spring dataset in input file " << filename << "\n");
        }

        hid_t stiffness_dset = H5Dopen1(file_id, stiffness_dset_name.c_str());
        if (stiffness_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required spring dataset in input file " << filename << "\n");
        }

        hid_t rest_length_dset = H5Dopen1(file_id, rest_length_dset_name.c_str());
        if (rest_length_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required spring dataset in input file " << filename << "\n");
        }

        // Read in the dimensions of the datasets.
        int rank;
        hsize_t dims[1];
        H5T_class_t class_id;
        size_t type_size;

        H5LTget_dataset_ndims(file_id, node1_idx_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, node1_idx_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset dimension in input file " << filename << "\n");
        }
        const int node1_idx_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, node2_idx_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, node2_idx_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset dimension in input file " << filename << "\n");
        }
        const int node2_idx_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, force_fcn_idx_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, force_fcn_idx_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset dimension in input file " << filename << "\n");
        }
        const int force_fcn_idx_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, stiffness_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, stiffness_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset dimension in input file " << filename << "\n");
        }
        const int stiffness_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, rest_length_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, rest_length_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset dimension in input file " << filename << "\n");
        }
        const int rest_length_size = int(dims[0]);

        if ((node1_idx_size != node2_idx_size    ) ||
            (node1_idx_size != force_fcn_idx_size) ||
            (node1_idx_size != stiffness_size    ) ||
            (node1_idx_size != rest_length_size  ))
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid spring dataset dimension in input file " << filename << "\n");
        }

        num_spring = node1_idx_size;

        // Define the file dataspace.
        static const int rankf = 1;
        hsize_t dimsf[rankf] = { num_spring };
        hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

        // Define the memory dataspace.
        static const int rankm = 1;
        hsize_t dimsm[rankm] = { BUFFER_SIZE };
        hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

        // Read in the spring data one block at a time.
        std::vector<int> node1_idx_buf(BUFFER_SIZE), node2_idx_buf(BUFFER_SIZE), force_fcn_idx_buf(BUFFER_SIZE);
        std::vector<double> stiffness_buf(BUFFER_SIZE), rest_length_buf(BUFFER_SIZE);
        const int num_blocks = num_spring/BUFFER_SIZE + (num_spring%BUFFER_SIZE == 0 ? 0 : 1);
        for (int block = 0; block < num_blocks; ++block)
        {
            // Determine whether we are reading in the last block in the file.
            const bool last_block = (block == num_blocks-1);

            // Determine the number of items to read (always BUFFER_SIZE except
            // for the final block in the file).
            const int num_spring_block = (last_block ? num_spring - block*BUFFER_SIZE : BUFFER_SIZE);
            TBOX_ASSERT(num_spring_block > 0 && num_spring_block <= BUFFER_SIZE);

            // Define the file hyperslab.
            hsize_t offsetf[rankf] = { block*BUFFER_SIZE };
            hsize_t countf[rankf] = { num_spring_block };
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

            // Define the memory hyperslab.
            hsize_t offsetm[rankm] = { 0 };
            hsize_t countm[rankm] = { num_spring_block };
            H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

            // Read data from hyperslab in the file into the hyperslab in memory.
            H5Dread(node1_idx_dset    , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node1_idx_buf    [0]);
            H5Dread(node2_idx_dset    , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node2_idx_buf    [0]);
            H5Dread(force_fcn_idx_dset, H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &force_fcn_idx_buf[0]);
            H5Dread(stiffness_dset    , H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &stiffness_buf    [0]);
            H5Dread(rest_length_dset  , H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &rest_length_buf  [0]);

            // Setup data for all local springs in the hyperslab.
            //
            // Note that in the spring map, each spring is associated with only
            // the "master" (first) vertex in the spring.  Consequently, spring
            // specifications are only constructed for and assigned to this
            // "master" vertex.  Appropriate forces are automatically applied to
            // the "slave" (second) vertex in
            // IBSpringForceGen::computeLagrangianForce().
            for (int k = 0; k < num_spring_block; ++k)
            {
                const int& node1_idx = node1_idx_buf[k];
                const int& node2_idx = node2_idx_buf[k];
                const int& force_fcn_idx = force_fcn_idx_buf[k];
                const double& stiffness = stiffness_buf[k];
                const double& rest_length = rest_length_buf[k];
                if (local_vertex_idx_set.find(node1_idx) != local_vertex_idx_set.end())
                {
                    ++num_local_spring;

                    std::map<int,Pointer<IBSpringForceSpec> >::iterator it = spring_data_map.find(node1_idx);
                    if (it == spring_data_map.end())
                    {
                        it = spring_data_map.insert(spring_data_map.end(), std::make_pair(node1_idx,new IBSpringForceSpec()));
                    }

                    Pointer<IBSpringForceSpec> old_spring_data = (*it).second;
                    TBOX_ASSERT(old_spring_data->getMasterNodeIndex() == -1 ||
                                old_spring_data->getMasterNodeIndex() == node1_idx);

                    std::vector<int>&    slave_idxs     = old_spring_data->getSlaveNodeIndices();
                    std::vector<int>&    force_fcn_idxs = old_spring_data->getForceFunctionIndices();
                    std::vector<double>& stiffnesses    = old_spring_data->getStiffnesses();
                    std::vector<double>& rest_lengths   = old_spring_data->getRestingLengths();

                    slave_idxs    .push_back(node2_idx    );
                    force_fcn_idxs.push_back(force_fcn_idx);
                    stiffnesses   .push_back(stiffness    );
                    rest_lengths  .push_back(rest_length  );

                    Pointer<IBSpringForceSpec> new_spring_data = new IBSpringForceSpec(
                        node1_idx, slave_idxs, force_fcn_idxs, stiffnesses, rest_lengths);

                    (*it).second = new_spring_data;
                }
            }
        }

        // Cleanup HDF5 data structures.
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Dclose(node1_idx_dset);
        H5Dclose(node2_idx_dset);
        H5Dclose(force_fcn_idx_dset);
        H5Dclose(stiffness_dset);
        H5Dclose(rest_length_dset);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":\n  Cannot find spring group in input file " << filename << "\n"
                   << "       base group name = " << base_group_name << "\n"
                   << "       spring group name = " << spring_group_name << "\n");
    }
    return;
}// buildLevelSpringDataCacheFromHDF5

void
IBHDF5Initializer::buildLevelBeamDataCacheFromHDF5(
    int& num_beam,
    int& num_local_beam,
    std::map<int,Pointer<IBBeamForceSpec> >& beam_data_map,
    const std::set<int>& local_vertex_idx_set,
    const hid_t file_id,
    const std::string& base_group_name,
    const Pointer<PatchLevel<NDIM> > level,
    const std::string& filename,
    const int file_number,
    const int num_files) const
{
    num_beam       = 0;
    num_local_beam = 0;
    beam_data_map.clear();

    // Check to see if the group for the dataset(s) exists.
    const std::string beam_group_name = base_group_name + "/beam";
    if (H5Gget_objinfo(file_id, beam_group_name.c_str(), 0, NULL) >= 0)
    {
        // Open the datasets.
        const std::string node1_idx_dset_name     = beam_group_name + "/node1_idx"    ;
        const std::string node2_idx_dset_name     = beam_group_name + "/node2_idx"    ;
        const std::string node3_idx_dset_name     = beam_group_name + "/node3_idx"    ;
        const std::string bend_rigidity_dset_name = beam_group_name + "/bend_rigidity";
        std::string rest_curvature_dset_name[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "_" << d;
            rest_curvature_dset_name[d] = beam_group_name + "/rest_curvature" + os.str();
        }

        hid_t node1_idx_dset = H5Dopen1(file_id, node1_idx_dset_name.c_str());
        if (node1_idx_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required beam dataset in input file " << filename << "\n");
        }

        hid_t node2_idx_dset = H5Dopen1(file_id, node2_idx_dset_name.c_str());
        if (node2_idx_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required beam dataset in input file " << filename << "\n");
        }

        hid_t node3_idx_dset = H5Dopen1(file_id, node3_idx_dset_name.c_str());
        if (node3_idx_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required beam dataset in input file " << filename << "\n");
        }

        hid_t bend_rigidity_dset = H5Dopen1(file_id, bend_rigidity_dset_name.c_str());
        if (bend_rigidity_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required beam dataset in input file " << filename << "\n");
        }

        hid_t rest_curvature_dset[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            rest_curvature_dset[d] = H5Dopen1(file_id, rest_curvature_dset_name[d].c_str());
            if (rest_curvature_dset < 0)
            {
                TBOX_ERROR(d_object_name << ":\n  Cannot find required beam dataset in input file " << filename << "\n");
            }
        }

        // Read in the dimensions of the datasets.
        int rank;
        hsize_t dims[1];
        H5T_class_t class_id;
        size_t type_size;

        H5LTget_dataset_ndims(file_id, node1_idx_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, node1_idx_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset dimension in input file " << filename << "\n");
        }
        const int node1_idx_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, node2_idx_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, node2_idx_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset dimension in input file " << filename << "\n");
        }
        const int node2_idx_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, node3_idx_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, node3_idx_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset dimension in input file " << filename << "\n");
        }
        const int node3_idx_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, bend_rigidity_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, bend_rigidity_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset dimension in input file " << filename << "\n");
        }
        const int bend_rigidity_size = int(dims[0]);

        int rest_curvature_size[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            H5LTget_dataset_ndims(file_id, rest_curvature_dset_name[d].c_str(), &rank);
            if (rank != 1)
            {
                TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset rank in input file " << filename << "\n");
            }
            H5LTget_dataset_info(file_id, rest_curvature_dset_name[d].c_str(), dims, &class_id, &type_size);
            if (dims[0] <= 0)
            {
                TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset dimension in input file " << filename << "\n");
            }
            rest_curvature_size[d] = int(dims[0]);
        }

        if ((node1_idx_size != node2_idx_size    ) ||
            (node1_idx_size != node3_idx_size    ) ||
            (node1_idx_size != bend_rigidity_size))
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset dimension in input file " << filename << "\n");
        }
        for (int d = 0; d < NDIM; ++d)
        {
            if (node1_idx_size != rest_curvature_size[d])
            {
                TBOX_ERROR(d_object_name << ":\n  Invalid beam dataset dimension in input file " << filename << "\n");
            }
        }

        num_beam = node1_idx_size;

        // Define the file dataspace.
        static const int rankf = 1;
        hsize_t dimsf[rankf] = { num_beam };
        hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

        // Define the memory dataspace.
        static const int rankm = 1;
        hsize_t dimsm[rankm] = { BUFFER_SIZE };
        hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

        // Read in the beam data one block at a time.
        std::vector<int> node1_idx_buf(BUFFER_SIZE), node2_idx_buf(BUFFER_SIZE), node3_idx_buf(BUFFER_SIZE);
        std::vector<double> bend_rigidity_buf(BUFFER_SIZE);
        std::vector<std::vector<double> > rest_curvature_buf(NDIM,std::vector<double>(BUFFER_SIZE));
        const int num_blocks = num_beam/BUFFER_SIZE + (num_beam%BUFFER_SIZE == 0 ? 0 : 1);
        for (int block = 0; block < num_blocks; ++block)
        {
            // Determine whether we are reading in the last block in the file.
            const bool last_block = (block == num_blocks-1);

            // Determine the number of items to read (always BUFFER_SIZE except
            // for the final block in the file).
            const int num_beam_block = (last_block ? num_beam - block*BUFFER_SIZE : BUFFER_SIZE);
            TBOX_ASSERT(num_beam_block > 0 && num_beam_block <= BUFFER_SIZE);

            // Define the file hyperslab.
            hsize_t offsetf[rankf] = { block*BUFFER_SIZE };
            hsize_t countf[rankf] = { num_beam_block };
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

            // Define the memory hyperslab.
            hsize_t offsetm[rankm] = { 0 };
            hsize_t countm[rankm] = { num_beam_block };
            H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

            // Read data from hyperslab in the file into the hyperslab in memory.
            H5Dread(node1_idx_dset    , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node1_idx_buf    [0]);
            H5Dread(node2_idx_dset    , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node2_idx_buf    [0]);
            H5Dread(node3_idx_dset    , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node3_idx_buf    [0]);
            H5Dread(bend_rigidity_dset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &bend_rigidity_buf[0]);
            for (int d = 0; d < NDIM; ++d)
            {
                H5Dread(rest_curvature_dset[d], H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &rest_curvature_buf[d][0]);
            }

            // Setup data for all local beams in the hyperslab.
            //
            // Note that in the beam map, each beam is associated with only the
            // "current" vertex in the beam.  Consequently, beam specifications
            // are only constructed for and assigned to this "current" vertex.
            // Appropriate forces are automatically applied to the "previous"
            // and "next" vertices in IBBeamForceGen::computeLagrangianForce().
            for (int k = 0; k < num_beam_block; ++k)
            {
                const int& node1_idx = node1_idx_buf[k];
                const int& node2_idx = node2_idx_buf[k];
                const int& node3_idx = node3_idx_buf[k];
                const double& bend_rigidity = bend_rigidity_buf[k];
                std::vector<double> rest_curvature(NDIM);
                for (int d = 0; d < NDIM; ++d)
                {
                    rest_curvature[d] = rest_curvature_buf[d][k];
                }
                if (local_vertex_idx_set.find(node2_idx) != local_vertex_idx_set.end())
                {
                    ++num_local_beam;

                    std::map<int,Pointer<IBBeamForceSpec> >::iterator it = beam_data_map.find(node2_idx);
                    if (it == beam_data_map.end())
                    {
                        it = beam_data_map.insert(beam_data_map.end(), std::make_pair(node2_idx,new IBBeamForceSpec()));
                    }

                    Pointer<IBBeamForceSpec> old_beam_data = (*it).second;
                    TBOX_ASSERT(old_beam_data->getMasterNodeIndex() == -1 ||
                                old_beam_data->getMasterNodeIndex() == node2_idx);

                    std::vector<std::pair<int,int> >&  neighbor_idxs   = old_beam_data->getNeighborNodeIndices();
                    std::vector<double>&               bend_rigidities = old_beam_data->getBendingRigidities();
                    std::vector<std::vector<double> >& rest_curvatures = old_beam_data->getMeshDependentCurvatures();

                    neighbor_idxs  .push_back(std::make_pair(node1_idx,node3_idx));
                    bend_rigidities.push_back(bend_rigidity                      );
                    rest_curvatures.push_back(rest_curvature                     );

                    Pointer<IBBeamForceSpec> new_beam_data = new IBBeamForceSpec(
                        node2_idx, neighbor_idxs, bend_rigidities, rest_curvatures);

                    (*it).second = new_beam_data;
                }
            }
        }

        // Cleanup HDF5 data structures.
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Dclose(node1_idx_dset);
        H5Dclose(node2_idx_dset);
        H5Dclose(node3_idx_dset);
        H5Dclose(bend_rigidity_dset);
        for (int d = 0; d < NDIM; ++d)
        {
            H5Dclose(rest_curvature_dset[d]);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << ":\n  Cannot find beam group in input file " << filename << "\n"
                   << "       base group name = " << base_group_name << "\n"
                   << "       beam group name = " << beam_group_name << "\n");
    }
    return;
}// buildLevelBeamDataCacheFromHDF5

void
IBHDF5Initializer::buildLevelTargetPointDataCacheFromHDF5(
    int& num_target_point,
    int& num_local_target_point,
    std::map<int,Pointer<IBTargetPointForceSpec> >& target_point_data_map,
    const std::set<int>& local_vertex_idx_set,
    const hid_t file_id,
    const std::string& base_group_name,
    const Pointer<PatchLevel<NDIM> > level,
    const std::string& filename,
    const int file_number,
    const int num_files) const
{
    num_target_point       = 0;
    num_local_target_point = 0;
    target_point_data_map.clear();

    // Check to see if the group for the dataset(s) exists.
    const std::string target_point_group_name = base_group_name + "/target_point";
    if (H5Gget_objinfo(file_id, target_point_group_name.c_str(), 0, NULL) >= 0)
    {
        // Open the datasets.
        const std::string node_idx_dset_name  = target_point_group_name + "/node_idx" ;
        const std::string stiffness_dset_name = target_point_group_name + "/stiffness";
        const std::string damping_dset_name   = target_point_group_name + "/damping";

        hid_t node_idx_dset = H5Dopen1(file_id, node_idx_dset_name.c_str());
        if (node_idx_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required target point dataset in input file " << filename << "\n");
        }

        hid_t stiffness_dset = H5Dopen1(file_id, stiffness_dset_name.c_str());
        if (stiffness_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required target point dataset in input file " << filename << "\n");
        }

        hid_t damping_dset = H5Dopen1(file_id, damping_dset_name.c_str());
        if (damping_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required target point dataset in input file " << filename << "\n");
        }

        // Read in the dimensions of the datasets.
        int rank;
        hsize_t dims[1];
        H5T_class_t class_id;
        size_t type_size;

        H5LTget_dataset_ndims(file_id, node_idx_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid target point dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, node_idx_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid target point dataset dimension in input file " << filename << "\n");
        }
        const int node_idx_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, stiffness_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid target point dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, stiffness_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid target point dataset dimension in input file " << filename << "\n");
        }
        const int stiffness_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, damping_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid target point dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, damping_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid target point dataset dimension in input file " << filename << "\n");
        }
        const int damping_size = int(dims[0]);

        if ((node_idx_size != stiffness_size) ||
            (node_idx_size != damping_size  ) )
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid target point dataset dimension in input file " << filename << "\n");
        }

        num_target_point = node_idx_size;

        // Define the file dataspace.
        static const int rankf = 1;
        hsize_t dimsf[rankf] = { num_target_point };
        hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

        // Define the memory dataspace.
        static const int rankm = 1;
        hsize_t dimsm[rankm] = { BUFFER_SIZE };
        hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

        // Read in the target point data one block at a time.
        std::vector<int> node_idx_buf(BUFFER_SIZE);
        std::vector<double> stiffness_buf(BUFFER_SIZE), damping_buf(BUFFER_SIZE);
        const int num_blocks = num_target_point/BUFFER_SIZE + (num_target_point%BUFFER_SIZE == 0 ? 0 : 1);
        for (int block = 0; block < num_blocks; ++block)
        {
            // Determine whether we are reading in the last block in the file.
            const bool last_block = (block == num_blocks-1);

            // Determine the number of items to read (always BUFFER_SIZE except
            // for the final block in the file).
            const int num_target_point_block = (last_block ? num_target_point - block*BUFFER_SIZE : BUFFER_SIZE);
            TBOX_ASSERT(num_target_point_block > 0 && num_target_point_block <= BUFFER_SIZE);

            // Define the file hyperslab.
            hsize_t offsetf[rankf] = { block*BUFFER_SIZE };
            hsize_t countf[rankf] = { num_target_point_block };
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

            // Define the memory hyperslab.
            hsize_t offsetm[rankm] = { 0 };
            hsize_t countm[rankm] = { num_target_point_block };
            H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

            // Read data from hyperslab in the file into the hyperslab in memory.
            H5Dread(node_idx_dset , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node_idx_buf [0]);
            H5Dread(stiffness_dset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &stiffness_buf[0]);
            H5Dread(damping_dset  , H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &damping_buf  [0]);

            // Setup data for all local target_points in the hyperslab.
            for (int k = 0; k < num_target_point_block; ++k)
            {
                const int& node_idx = node_idx_buf[k];
                const double& stiffness = stiffness_buf[k];
                const double& damping = damping_buf[k];
                if (local_vertex_idx_set.find(node_idx) != local_vertex_idx_set.end())
                {
                    ++num_local_target_point;

                    Pointer<IBTargetPointForceSpec> target_point_data;
                    std::map<int,Pointer<IBTargetPointForceSpec> >::iterator it = target_point_data_map.find(node_idx);
                    if (it == target_point_data_map.end())
                    {
                        target_point_data = (*target_point_data_map.insert(target_point_data_map.end(), std::make_pair(node_idx,new IBTargetPointForceSpec()))).second;
                        TBOX_ASSERT(target_point_data->getMasterNodeIndex() == -1);
                    }
                    else
                    {
                        target_point_data = (*it).second;
                        TBOX_ASSERT(target_point_data->getMasterNodeIndex() == node_idx);
                    }
                    target_point_data->getMasterNodeIndex() = node_idx;
                    target_point_data->getStiffness() = stiffness;
                    target_point_data->getDamping() = damping;
                }
            }
        }

        // Cleanup HDF5 data structures.
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Dclose(node_idx_dset);
        H5Dclose(stiffness_dset);
        H5Dclose(damping_dset);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":\n  Cannot find target point group in input file " << filename << "\n"
                   << "       base group name = " << base_group_name << "\n"
                   << "       target point group name = " << target_point_group_name << "\n");
    }
    return;
}// buildLevelTargetPointDataCacheFromHDF5

void
IBHDF5Initializer::buildLevelInstrumentationDataCacheFromHDF5(
    std::vector<std::string>& instrument_names,
    int& num_inst_point,
    int& num_local_inst_point,
    std::map<int,Pointer<IBInstrumentationSpec> >& inst_point_data_map,
    const std::set<int>& local_vertex_idx_set,
    const hid_t file_id,
    const std::string& base_group_name,
    const Pointer<PatchLevel<NDIM> > level,
    const std::string& filename,
    const int file_number,
    const int num_files) const
{
    num_inst_point       = 0;
    num_local_inst_point = 0;
    inst_point_data_map.clear();

    // Check to see if the group for the dataset(s) exists.
    const std::string instrumentation_group_name = base_group_name + "/instrumentation";
    if (H5Gget_objinfo(file_id, instrumentation_group_name.c_str(), 0, NULL) >= 0)
    {
        // Read the instrument names.
        const std::string num_inst_dset_name = instrumentation_group_name + "/num_inst";
        hid_t num_inst_dset = H5Dopen1(file_id, num_inst_dset_name.c_str());
        if (num_inst_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required instrumentation dataset in input file " << filename << "\n");
        }
        H5Dclose(num_inst_dset);

        int num_inst = -1;
        H5LTread_dataset_int(file_id, num_inst_dset_name.c_str(), &num_inst);
        if (num_inst <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid number of instruments in input file " << filename << "\n");
        }
        for (int k = 0; k < num_inst; ++k)
        {
            std::ostringstream num_stream;
            num_stream << k;
            const std::string name_dset_name = instrumentation_group_name + "/name_" + num_stream.str();
            char buffer[STRING_BUFFER_SIZE];
            H5LTread_dataset_string(file_id, name_dset_name.c_str(), buffer);
            instrument_names.push_back(std::string(buffer));
        }

        // Open the datasets.
        const std::string node_idx_dset_name       = instrumentation_group_name + "/node_idx"      ;
        const std::string meter_idx_dset_name      = instrumentation_group_name + "/meter_idx"     ;
        const std::string meter_node_idx_dset_name = instrumentation_group_name + "/meter_node_idx";

        hid_t node_idx_dset = H5Dopen1(file_id, node_idx_dset_name.c_str());
        if (node_idx_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required instrumentation dataset in input file " << filename << "\n");
        }

        hid_t meter_idx_dset = H5Dopen1(file_id, meter_idx_dset_name.c_str());
        if (meter_idx_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required instrumentation dataset in input file " << filename << "\n");
        }

        hid_t meter_node_idx_dset = H5Dopen1(file_id, meter_node_idx_dset_name.c_str());
        if (meter_node_idx_dset < 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Cannot find required instrumentation dataset in input file " << filename << "\n");
        }

        // Read in the dimensions of the datasets.
        int rank;
        hsize_t dims[1];
        H5T_class_t class_id;
        size_t type_size;

        H5LTget_dataset_ndims(file_id, node_idx_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid instrumentation dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, node_idx_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid instrumentation dataset dimension in input file " << filename << "\n");
        }
        const int node_idx_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, meter_idx_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid instrumentation dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, meter_idx_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid instrumentation dataset dimension in input file " << filename << "\n");
        }
        const int meter_idx_size = int(dims[0]);

        H5LTget_dataset_ndims(file_id, meter_node_idx_dset_name.c_str(), &rank);
        if (rank != 1)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid instrumentation dataset rank in input file " << filename << "\n");
        }
        H5LTget_dataset_info(file_id, meter_node_idx_dset_name.c_str(), dims, &class_id, &type_size);
        if (dims[0] <= 0)
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid instrumentation dataset dimension in input file " << filename << "\n");
        }
        const int meter_node_idx_size = int(dims[0]);

        if ((node_idx_size != meter_idx_size) || (node_idx_size != meter_node_idx_size))
        {
            TBOX_ERROR(d_object_name << ":\n  Invalid instrumentation dataset dimension in input file " << filename << "\n");
        }

        num_inst_point = node_idx_size;

        // Define the file dataspace.
        static const int rankf = 1;
        hsize_t dimsf[rankf] = { num_inst_point };
        hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

        // Define the memory dataspace.
        static const int rankm = 1;
        hsize_t dimsm[rankm] = { BUFFER_SIZE };
        hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

        // Read in the instrumentation data one block at a time.
        std::vector<int> node_idx_buf(BUFFER_SIZE), meter_idx_buf(BUFFER_SIZE), meter_node_idx_buf(BUFFER_SIZE);
        const int num_blocks = num_inst_point/BUFFER_SIZE + (num_inst_point%BUFFER_SIZE == 0 ? 0 : 1);
        for (int block = 0; block < num_blocks; ++block)
        {
            // Determine whether we are reading in the last block in the file.
            const bool last_block = (block == num_blocks-1);

            // Determine the number of items to read (always BUFFER_SIZE except
            // for the final block in the file).
            const int num_inst_point_block = (last_block ? num_inst_point - block*BUFFER_SIZE : BUFFER_SIZE);
            TBOX_ASSERT(num_inst_point_block > 0 && num_inst_point_block <= BUFFER_SIZE);

            // Define the file hyperslab.
            hsize_t offsetf[rankf] = { block*BUFFER_SIZE };
            hsize_t countf[rankf] = { num_inst_point_block };
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

            // Define the memory hyperslab.
            hsize_t offsetm[rankm] = { 0 };
            hsize_t countm[rankm] = { num_inst_point_block };
            H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

            // Read data from hyperslab in the file into the hyperslab in memory.
            H5Dread(node_idx_dset      , H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, &node_idx_buf      [0]);
            H5Dread(meter_idx_dset     , H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, &meter_idx_buf     [0]);
            H5Dread(meter_node_idx_dset, H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, &meter_node_idx_buf[0]);

            // Setup data for all local instrumented points in the hyperslab.
            for (int k = 0; k < num_inst_point_block; ++k)
            {
                const int& node_idx = node_idx_buf[k];
                const int& meter_idx = meter_idx_buf[k];
                const int& meter_node_idx = meter_node_idx_buf[k];
                if (local_vertex_idx_set.find(node_idx) != local_vertex_idx_set.end())
                {
                    ++num_local_inst_point;

                    Pointer<IBInstrumentationSpec> inst_point_data;
                    std::map<int,Pointer<IBInstrumentationSpec> >::iterator it = inst_point_data_map.find(node_idx);
                    if (it == inst_point_data_map.end())
                    {
                        inst_point_data = (*inst_point_data_map.insert(inst_point_data_map.end(), std::make_pair(node_idx,new IBInstrumentationSpec()))).second;
                        TBOX_ASSERT(inst_point_data->getMasterNodeIndex() == -1);
                    }
                    else
                    {
                        inst_point_data = (*it).second;
                        TBOX_ASSERT(inst_point_data->getMasterNodeIndex() == node_idx);
                    }
                    inst_point_data->getMasterNodeIndex() = node_idx;
                    inst_point_data->getMeterIndex() = meter_idx;
                    inst_point_data->getNodeIndex() = meter_node_idx;
                }
            }
        }

        // Cleanup HDF5 data structures.
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Dclose(node_idx_dset);
        H5Dclose(meter_idx_dset);
        H5Dclose(meter_node_idx_dset);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":\n  Cannot find instrumentation group in input file " << filename << "\n"
                   << "       base_group name = " << base_group_name << "\n"
                   << "       instrumentation group name = " << instrumentation_group_name << "\n");
    }
    return;
}// buildLevelInstrumentationDataCacheFromHDF5

int
IBHDF5Initializer::getCanonicalLagrangianIndex(
    const std::pair<int,int>& global_vertex_idx,
    const int global_index_offset) const
{
    return d_level_vertex_offset[global_vertex_idx.first]+global_vertex_idx.second+global_index_offset;
}// getCanonicalLagrangianIndex

void
IBHDF5Initializer::clearLevelDataCache()
{
    d_cache_level_number = -1;
    d_level_num_vertex            .clear();
    d_level_num_local_vertex      .clear();
    d_level_vertex_offset         .clear();
    d_level_posns                 .clear();
    d_level_vertex_idxs           .clear();
    d_level_cell_idxs             .clear();
    d_level_patch_nums            .clear();
    d_level_reset_specs_set       .clear();
    d_level_num_spring            .clear();
    d_level_num_local_spring      .clear();
    d_level_spring_data_map       .clear();
    d_level_num_beam              .clear();
    d_level_num_local_beam        .clear();
    d_level_beam_data_map         .clear();
    d_level_num_target_point      .clear();
    d_level_num_local_target_point.clear();
    d_level_target_point_data_map .clear();
    d_level_num_inst_point        .clear();
    d_level_num_local_inst_point  .clear();
    d_level_inst_point_data_map   .clear();
    return;
}// clearLevelDataCache

std::vector<Pointer<Streamable> >
IBHDF5Initializer::initializeSpecs(
    const std::pair<int,int>& local_vertex_idx,
    const std::pair<int,int>& global_vertex_idx,
    const int global_index_offset)
{
    std::vector<Pointer<Streamable> > vertex_specs;

    const int& ln = d_cache_level_number;
    const int& j_local = local_vertex_idx.first;
    const int& k_local = local_vertex_idx.second;
    const int& j = global_vertex_idx.first;
    const int& k = global_vertex_idx.second;
    const int vertex_offset = d_level_vertex_offset[j]+global_index_offset;

    const bool specs_need_reset = d_level_reset_specs_set[j].find(k) == d_level_reset_specs_set[j].end();
    TBOX_ASSERT(specs_need_reset);

    // Initialize any spring specifications associated with the present vertex.
    if (d_enable_springs[ln][j])
    {
        std::map<int,Pointer<IBSpringForceSpec> >::const_iterator cit =
            d_level_spring_data_map[j].find(k);
        if (cit != d_level_spring_data_map[j].end())
        {
            Pointer<IBSpringForceSpec> spring_data = (*cit).second;
            if (specs_need_reset)
            {
                if (d_using_uniform_spring_force_fcn_idx[ln][j])
                {
                    std::vector<int>& force_fcn_idxs = spring_data->getForceFunctionIndices();
                    std::fill(force_fcn_idxs.begin(), force_fcn_idxs.end(), d_uniform_spring_force_fcn_idx[ln][j]);
                }
                if (d_using_uniform_spring_stiffness[ln][j])
                {
                    std::vector<double>& stiffnesses = spring_data->getStiffnesses();
                    std::fill(stiffnesses.begin(), stiffnesses.end(), d_uniform_spring_stiffness[ln][j]);
                }
                if (d_using_uniform_spring_rest_length[ln][j])
                {
                    std::vector<double>& rest_lengths = spring_data->getRestingLengths();
                    std::fill(rest_lengths.begin(), rest_lengths.end(), d_uniform_spring_rest_length[ln][j]);
                }
                spring_data->getMasterNodeIndex() += vertex_offset;
                std::vector<int>& slave_idxs = spring_data->getSlaveNodeIndices();
                for (std::vector<int>::iterator it = slave_idxs.begin();
                     it != slave_idxs.end(); ++it)
                {
                    (*it) += vertex_offset;
                }
            }
            vertex_specs.push_back(spring_data);
        }
    }

    // Initialize any beam specifications associated with the present vertex.
    if (d_enable_beams[ln][j])
    {
        std::map<int,Pointer<IBBeamForceSpec> >::const_iterator cit =
            d_level_beam_data_map[j].find(k);
        if (cit != d_level_beam_data_map[j].end())
        {
            Pointer<IBBeamForceSpec> beam_data = (*cit).second;
            if (specs_need_reset)
            {
                if (d_using_uniform_beam_bend_rigidity[ln][j])
                {
                    std::vector<double>& bend_rigidities = beam_data->getBendingRigidities();
                    std::fill(bend_rigidities.begin(), bend_rigidities.end(), d_uniform_beam_bend_rigidity[ln][j]);
                }
                beam_data->getMasterNodeIndex() += vertex_offset;
                std::vector<std::pair<int,int> >& neighbor_idxs = beam_data->getNeighborNodeIndices();
                for (std::vector<std::pair<int,int> >::iterator it = neighbor_idxs.begin();
                     it != neighbor_idxs.end(); ++it)
                {
                    (*it).first  += vertex_offset;
                    (*it).second += vertex_offset;
                }
            }
            vertex_specs.push_back(beam_data);
        }
    }

    // Initialize any target point specifications associated with the present
    // vertex.
    if (d_enable_target_points[ln][j])
    {
        if (d_using_uniform_target_stiffness[ln][j] && d_using_uniform_target_damping[ln][j])
        {
            Pointer<IBTargetPointForceSpec> target_point_data = new IBTargetPointForceSpec(
                k+vertex_offset, d_uniform_target_stiffness[ln][j], d_uniform_target_damping[ln][j], d_level_posns[j_local][k_local]);
            vertex_specs.push_back(target_point_data);
        }
        else if (d_using_uniform_target_stiffness[ln][j] || d_using_uniform_target_damping[ln][j])
        {
            TBOX_ERROR(d_object_name << ":\n  Uniform properties for target point springs incorrectly specified\n");
        }
        else
        {
            std::map<int,Pointer<IBTargetPointForceSpec> >::const_iterator cit =
                d_level_target_point_data_map[j].find(k);
            if (cit != d_level_target_point_data_map[j].end())
            {
                Pointer<IBTargetPointForceSpec> target_point_data = (*cit).second;
                if (specs_need_reset)
                {
                    target_point_data->getMasterNodeIndex() += vertex_offset;
                    target_point_data->getTargetPointPosition() = d_level_posns[j_local][k_local];
                }
                vertex_specs.push_back(target_point_data);
            }
        }
    }

    // Initialize any instrumentation specifications associated with the present
    // vertex.
    int instrument_index_offset = 0;
    for (int coarser_ln = 0; coarser_ln < ln; ++coarser_ln)
    {
        for (int file_number = 0; file_number < int(d_instrument_names[coarser_ln].size());
             ++file_number)
        {
            instrument_index_offset += d_instrument_names[coarser_ln][file_number].size();
        }
    }
    for (int file_number = 0; file_number < j; ++file_number)
    {
        instrument_index_offset += d_instrument_names[ln][file_number].size();
    }

    if (d_enable_instrumentation[ln][j])
    {
        std::map<int,Pointer<IBInstrumentationSpec> >::const_iterator cit =
            d_level_inst_point_data_map[j].find(k);
        if (cit != d_level_inst_point_data_map[j].end())
        {
            Pointer<IBInstrumentationSpec> inst_point_data = (*cit).second;
            if (specs_need_reset)
            {
                inst_point_data->getMasterNodeIndex() += vertex_offset;
                inst_point_data->getMeterIndex() += instrument_index_offset;
            }
            vertex_specs.push_back(inst_point_data);
        }
    }

    // Indicate that the spec objects associated with the present vertex have
    // been reset.
    d_level_reset_specs_set[j].insert(d_level_reset_specs_set[j].end(),k);
    return vertex_specs;
}// initializeSpecs

void
IBHDF5Initializer::getFromInput(
    Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    // Determine whether to use "batons" to prevent multiple MPI processes from
    // reading the same file at once.
    d_use_file_batons = db->getBoolWithDefault("use_file_batons",d_use_file_batons);

    // Determine the (maximum) number of levels in the locally refined grid.
    // Note that each piece of the Lagrangian structure must be assigned to a
    // particular level of the grid.
    if (db->keyExists("max_levels"))
    {
        d_max_levels = db->getInteger("max_levels");
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Key data `max_levels' not found in input.");
    }

    if (d_max_levels < 1)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Key data `max_levels' found in input is < 1.");
    }

    // Resize the vectors that are indexed by the level number.
    d_level_is_initialized.resize(d_max_levels,false);

    d_filenames.resize(d_max_levels);

    d_enable_springs.resize(d_max_levels);
    d_using_uniform_spring_stiffness.resize(d_max_levels);
    d_uniform_spring_stiffness.resize(d_max_levels);
    d_using_uniform_spring_rest_length.resize(d_max_levels);
    d_uniform_spring_rest_length.resize(d_max_levels);
    d_using_uniform_spring_force_fcn_idx.resize(d_max_levels);
    d_uniform_spring_force_fcn_idx.resize(d_max_levels);

    d_enable_beams.resize(d_max_levels);
    d_using_uniform_beam_bend_rigidity.resize(d_max_levels);
    d_uniform_beam_bend_rigidity.resize(d_max_levels);

    d_enable_target_points.resize(d_max_levels);
    d_using_uniform_target_stiffness.resize(d_max_levels);
    d_uniform_target_stiffness.resize(d_max_levels);
    d_using_uniform_target_damping.resize(d_max_levels);
    d_uniform_target_damping.resize(d_max_levels);

    d_enable_instrumentation.resize(d_max_levels);
    d_instrument_names.resize(d_max_levels);

    // Determine the various input file names.
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        std::ostringstream db_key_name_stream;
        db_key_name_stream << "filenames_" << ln;
        const std::string db_key_name = db_key_name_stream.str();
        if (db->keyExists(db_key_name))
        {
            const int n_files = db->getArraySize(db_key_name);
            d_filenames[ln].resize(n_files);
            db->getStringArray(db_key_name, &d_filenames[ln][0], n_files);
        }
        else
        {
            TBOX_WARNING(d_object_name << ":  "
                         << "Key data `" + db_key_name + "' not found in input.");
        }
    }

    // Read in any sub-databases associated with the input file names.
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_filenames = d_filenames[ln].size();

        d_enable_springs[ln].resize(num_filenames,true);

        d_using_uniform_spring_stiffness[ln].resize(num_filenames,false);
        d_uniform_spring_stiffness[ln].resize(num_filenames,-1.0);

        d_using_uniform_spring_rest_length[ln].resize(num_filenames,false);
        d_uniform_spring_rest_length[ln].resize(num_filenames,-1.0);

        d_using_uniform_spring_force_fcn_idx[ln].resize(num_filenames,false);
        d_uniform_spring_force_fcn_idx[ln].resize(num_filenames,-1);

        d_enable_beams[ln].resize(num_filenames,true);

        d_using_uniform_beam_bend_rigidity[ln].resize(num_filenames,false);
        d_uniform_beam_bend_rigidity[ln].resize(num_filenames,-1.0);

        d_enable_target_points[ln].resize(num_filenames,true);

        d_using_uniform_target_stiffness[ln].resize(num_filenames,false);
        d_uniform_target_stiffness[ln].resize(num_filenames,-1.0);

        d_using_uniform_target_damping[ln].resize(num_filenames,false);
        d_uniform_target_damping[ln].resize(num_filenames,-1.0);

        d_enable_instrumentation[ln].resize(num_filenames,true);

        for (int j = 0; j < num_filenames; ++j)
        {
            const std::string& filename = d_filenames[ln][j];
            if (db->isDatabase(filename))
            {
                Pointer<Database> sub_db = db->getDatabase(filename);

                // Determine whether to enable or disable any particular
                // features.
                if (sub_db->keyExists("enable_springs"))
                {
                    d_enable_springs[ln][j] = sub_db->getBool("enable_springs");
                }
                if (sub_db->keyExists("enable_beams"))
                {
                    d_enable_beams[ln][j] = sub_db->getBool("enable_beams");
                }
                if (sub_db->keyExists("enable_target_points"))
                {
                    d_enable_target_points[ln][j] = sub_db->getBool("enable_target_points");
                }
                if (sub_db->keyExists("enable_instrumentation"))
                {
                    d_enable_instrumentation[ln][j] = sub_db->getBool("enable_instrumentation");
                }

                // Determine whether to use uniform values for any particular
                // structure attributes.
                if (sub_db->keyExists("uniform_spring_stiffness"))
                {
                    d_using_uniform_spring_stiffness[ln][j] = true;
                    d_uniform_spring_stiffness[ln][j] = sub_db->getDouble("uniform_spring_stiffness");

                    if (d_uniform_spring_stiffness[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_spring_stiffness' in database " << filename << "\n"
                                   << "  spring constant is negative\n");
                    }
                }
                if (sub_db->keyExists("uniform_spring_rest_length"))
                {
                    d_using_uniform_spring_rest_length[ln][j] = true;
                    d_uniform_spring_rest_length[ln][j] = sub_db->getDouble("uniform_spring_rest_length");

                    if (d_uniform_spring_rest_length[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_spring_rest_length' in database " << filename << "\n"
                                   << "  spring resting length is negative\n");
                    }
                }
                if (sub_db->keyExists("uniform_spring_force_fcn_idx"))
                {
                    d_using_uniform_spring_force_fcn_idx[ln][j] = true;
                    d_uniform_spring_force_fcn_idx[ln][j] = sub_db->getInteger("uniform_spring_force_fcn_idx");
                }

                if (sub_db->keyExists("uniform_beam_bend_rigidity"))
                {
                    d_using_uniform_beam_bend_rigidity[ln][j] = true;
                    d_uniform_beam_bend_rigidity[ln][j] = sub_db->getDouble("uniform_beam_bend_rigidity");

                    if (d_uniform_beam_bend_rigidity[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_beam_bend_rigidity' in database " << filename << "\n"
                                   << "  beam bending rigidity is negative\n");
                    }
                }

                if (sub_db->keyExists("uniform_target_stiffness"))
                {
                    d_using_uniform_target_stiffness[ln][j] = true;
                    d_uniform_target_stiffness[ln][j] = sub_db->getDouble("uniform_target_stiffness");

                    if (d_uniform_target_stiffness[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_target_stiffness' in database " << filename << "\n"
                                   << "  target point spring constant is negative\n");
                    }
                }

                if (sub_db->keyExists("uniform_target_damping"))
                {
                    d_using_uniform_target_damping[ln][j] = true;
                    d_uniform_target_damping[ln][j] = sub_db->getDouble("uniform_target_damping");

                    if (d_uniform_target_damping[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_target_damping' in database " << filename << "\n"
                                   << "  target point spring constant is negative\n");
                    }
                }
            }
        }
    }

    // Output the names of the input files to be read along with additional
    // debugging information.
    pout << d_object_name << ":  Reading from input files: \n";
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_filenames = d_filenames[ln].size();
        for (int j = 0; j < num_filenames; ++j)
        {
            const std::string& filename = d_filenames[ln][j];
            pout << "  filename: " << filename << "\n"
                 << "  assigned to level " << ln << " of the Cartesian grid patch hierarchy\n";
            if (!d_enable_springs[ln][j])
            {
                pout << "  NOTE: spring forces are DISABLED for " << filename << "\n";
            }
            else
            {
                if (d_using_uniform_spring_stiffness[ln][j])
                {
                    pout << "  NOTE: uniform spring stiffnesses are being employed for the structure named " << filename << "\n"
                         << "        any stiffness information in input file " << filename << " will be IGNORED\n";
                }
                if (d_using_uniform_spring_rest_length[ln][j])
                {
                    pout << "  NOTE: uniform spring resting lengths are being employed for the structure named " << filename << "\n"
                         << "        any resting length information in input file " << filename << " will be IGNORED\n";
                }
                if (d_using_uniform_spring_force_fcn_idx[ln][j])
                {
                    pout << "  NOTE: uniform spring force functions are being employed for the structure named " << filename << "\n"
                         << "        any force function index information in input file " << filename << " will be IGNORED\n";
                }
            }

            if (!d_enable_beams[ln][j])
            {
                pout << "  NOTE: beam forces are DISABLED for " << filename << "\n";
            }
            else
            {
                if (d_using_uniform_beam_bend_rigidity[ln][j])
                {
                    pout << "  NOTE: uniform beam bending rigidities are being employed for the structure named " << filename << "\n"
                         << "        any stiffness information in input file " << filename << " will be IGNORED\n";
                }
            }

            if (!d_enable_target_points[ln][j])
            {
                pout << "  NOTE: target point penalty forces are DISABLED for " << filename << "\n";
            }
            else
            {
                if (d_using_uniform_target_stiffness[ln][j])
                {
                    pout << "  NOTE: uniform target point stiffnesses are being employed for the structure named " << filename << "\n"
                         << "        any target point stiffness information in input file " << filename << " will be IGNORED\n";
                }
                if (d_using_uniform_target_damping[ln][j])
                {
                    pout << "  NOTE: uniform target point damping factors are being employed for the structure named " << filename << "\n"
                         << "        any target point damping factor information in input file " << filename << " will be IGNORED\n";
                }
            }

            if (!d_enable_instrumentation[ln][j])
            {
                pout << "  NOTE: instrumentation is DISABLED for " << filename << "\n";
            }

            pout << "\n";
        }
    }
    return;
}// getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBHDF5Initializer>;

//////////////////////////////////////////////////////////////////////////////
