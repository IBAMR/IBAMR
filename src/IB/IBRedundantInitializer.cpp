// Filename: IBRedundantInitializer.cpp
// Created on 24 May 2018 by Boyce Griffith
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

#include <algorithm>
#include <cmath>
#include <ios>
#include <iosfwd>
#include <istream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <stddef.h>
#include <string>
#include <utility>
#include <vector>

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "boost/array.hpp"
#include "boost/math/special_functions/round.hpp"
#include "boost/multi_array.hpp"
#include "ibamr/IBAnchorPointSpec.h"
#include "ibamr/IBBeamForceSpec.h"
#include "ibamr/IBInstrumentationSpec.h"
#include "ibamr/IBRedundantInitializer.h"
#include "ibamr/IBRodForceSpec.h"
#include "ibamr/IBSourceSpec.h"
#include "ibamr/IBSpringForceSpec.h"
#include "ibamr/IBStandardSourceGen.h"
#include "ibamr/IBTargetPointForceSpec.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IndexUtilities.h"
#include "ibtk/LData.h"
#include "ibtk/LIndexSetData.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeSet.h"
#include "ibtk/LNodeSetData.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/Streamable.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"
#include <../contrib/eigen/Eigen/src/Core/products/GeneralBlockPanelKernel.h>

namespace IBTK
{
class LDataManager;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

IBRedundantInitializer::IBRedundantInitializer(const std::string& object_name, Pointer<Database> input_db)
    : d_object_name(object_name),
      d_max_levels(-1),
      d_level_is_initialized(),
      d_silo_writer(NULL),
      d_base_filename(),
      d_length_scale_factor(1.0),
      d_posn_shift(Vector::Zero()),
      d_num_vertex(),
      d_vertex_offset(),
      d_vertex_posn(),
      d_spring_edge_map(),
      d_spring_spec_data(),
      d_xspring_edge_map(),
      d_xspring_spec_data(),
      d_beam_spec_data(),
      d_rod_edge_map(),
      d_rod_spec_data(),
      d_target_spec_data(),
      d_anchor_spec_data(),
      d_bdry_mass_spec_data(),
      d_directors(),
      d_instrument_idx(),
      d_source_idx(),
      d_global_index_offset(),
      d_init_structure_on_level_fcn(NULL),
      d_init_spring_on_level_fcn(NULL),
      d_init_beam_on_level_fcn(NULL),
      d_init_director_and_rod_on_level_fcn(NULL),
      d_init_boundary_mass_on_level_fcn(NULL),
      d_init_target_pt_on_level_fcn(NULL),
      d_data_processed(false)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
#endif

    // Register the specification objects with the StreamableManager class.
    IBAnchorPointSpec::registerWithStreamableManager();
    IBBeamForceSpec::registerWithStreamableManager();
    IBInstrumentationSpec::registerWithStreamableManager();
    IBRodForceSpec::registerWithStreamableManager();
    IBSourceSpec::registerWithStreamableManager();
    IBSpringForceSpec::registerWithStreamableManager();
    IBTargetPointForceSpec::registerWithStreamableManager();

    getFromInput(input_db);
    return;
} // IBRedundantInitializer

IBRedundantInitializer::~IBRedundantInitializer()
{
    pout << d_object_name << ":  Deallocating initialization data.\n";
    return;
} // ~IBRedundantInitializer

void
IBRedundantInitializer::registerLSiloDataWriter(Pointer<LSiloDataWriter> silo_writer)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(silo_writer);
#endif

    // Cache a pointer to the data writer.
    d_silo_writer = silo_writer;

    // Check to see if we are starting from a restart file.
    RestartManager* restart_manager = RestartManager::getManager();
    const bool is_from_restart = restart_manager->isFromRestart();

    // Initialize the Silo data writer only if we are not starting from a
    // restart file.
    if (!is_from_restart)
    {
        // Check if data has been processed.
        init();

        for (int ln = 0; ln < d_max_levels; ++ln)
        {
            if (d_level_is_initialized[ln])
            {
                initializeLSiloDataWriter(ln);
            }
        }
    }
    return;
} // registerLSiloDataWriter

void
IBRedundantInitializer::init()
{
    if (d_data_processed)
    {
        return;
    }
    else
    {
        // Process structure information.
        initializeStructurePosition();
        initializeSprings();
        initializeBeams();
        initializeDirectorAndRods();
        initializeBoundaryMass();
        initializeTargetPts();
    }

    // Indicate that we have processed data.
    d_data_processed = true;

    return;
}

bool
IBRedundantInitializer::getLevelHasLagrangianData(const int level_number, const bool /*can_be_refined*/) const
{
    return !d_num_vertex[level_number].empty();
} // getLevelHasLagrangianData

unsigned int
IBRedundantInitializer::computeGlobalNodeCountOnPatchLevel(const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                           const int level_number,
                                                           const double /*init_data_time*/,
                                                           const bool /*can_be_refined*/,
                                                           const bool /*initial_time*/)
{
    // Check if data has been processed.
    init();

    return std::accumulate(d_num_vertex[level_number].begin(), d_num_vertex[level_number].end(), 0);
}

unsigned int
IBRedundantInitializer::computeLocalNodeCountOnPatchLevel(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                          const int level_number,
                                                          const double /*init_data_time*/,
                                                          const bool /*can_be_refined*/,
                                                          const bool /*initial_time*/)
{
    // Check if data has been processed.
    init();
    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();

    // Loop over all patches in the specified level of the patch level and count
    // the number of local vertices.
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        // Count the number of vertices whose initial locations will be within
        // the given patch.
        std::vector<std::pair<int, int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, hierarchy);
        local_node_count += patch_vertices.size();
    }
    return local_node_count;
} // computeLocalNodeCountOnPatchLevel

void
IBRedundantInitializer::initializeStructureIndexingOnPatchLevel(
    std::map<int, std::string>& strct_id_to_strct_name_map,
    std::map<int, std::pair<int, int> >& strct_id_to_lag_idx_range_map,
    const int level_number,
    const double /*init_data_time*/,
    const bool /*can_be_refined*/,
    const bool /*initial_time*/,
    LDataManager* const /*l_data_manager*/)
{
    // Check if data has been processed.
    init();

    int offset = 0;
    for (int j = 0; j < static_cast<int>(d_base_filename[level_number].size()); ++j)
    {
        strct_id_to_strct_name_map[j] = d_base_filename[level_number][j];
        strct_id_to_lag_idx_range_map[j] = std::make_pair(offset, offset + d_num_vertex[level_number][j]);
        offset += d_num_vertex[level_number][j];
    }
    return;
} // initializeStructureIndexingOnPatchLevel

void
IBRedundantInitializer::registerInitStructureFunction(InitStructureOnLevel fcn)
{
    d_init_structure_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitSpringDataFunction(InitSpringDataOnLevel fcn)
{
    d_init_spring_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitBeamDataFunction(InitBeamDataOnLevel fcn)
{
    d_init_beam_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitDirectorAndRodFunction(InitDirectorAndRodOnLevel fcn)
{
    d_init_director_and_rod_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitBoundaryMassFunction(InitBoundaryMassOnLevel fcn)
{
    d_init_boundary_mass_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitTargetPtFunction(InitTargetPtOnLevel fcn)
{
    d_init_target_pt_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitAnchorPtFunction(InitAnchorPtOnLevel fcn)
{
    d_init_anchor_pt_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::initializeStructurePosition()
{
    if (!d_init_structure_on_level_fcn)
    {
        TBOX_ERROR("IBRedundantInitializer::initializeStructurePosition()\n"
                   << " no function registered to initialize structure.\n");
    }
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_num_vertex[ln].resize(num_base_filename, 0);
        d_vertex_offset[ln].resize(num_base_filename, std::numeric_limits<int>::max());
        d_vertex_posn[ln].resize(num_base_filename);
        for (unsigned int j = 0; j < num_base_filename; ++j)
        {
            if (j == 0)
            {
                d_vertex_offset[ln][j] = 0;
            }
            else
            {
                d_vertex_offset[ln][j] = d_vertex_offset[ln][j - 1] + d_num_vertex[ln][j - 1];
            }

            d_init_structure_on_level_fcn(j, ln, d_num_vertex[ln][j], d_vertex_posn[ln][j]);
#if !defined(NDEBUG)
            if (d_num_vertex[ln][j] <= 0)
            {
                TBOX_ERROR(d_object_name << ":\n Invalid number of vertices " << d_num_vertex[ln][j] << " of structure "
                                         << j << " on level " << ln << std::endl);
            }
#endif

            // Shift and scale the position of structures
            for (int k = 0; k < d_num_vertex[ln][j]; ++k)
            {
                Point& X = d_vertex_posn[ln][j][k];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X[d] = d_length_scale_factor * (X[d] + d_posn_shift[d]);
                }
            }
        }
    }

    return;
} // initializeStructurePosition

void
IBRedundantInitializer::initializeSprings()
{
    if (!d_init_spring_on_level_fcn) return;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_spring_edge_map[ln].resize(num_base_filename);
        d_spring_spec_data[ln].resize(num_base_filename);
        for (unsigned int j = 0; j < num_base_filename; ++j)
        {
            d_init_spring_on_level_fcn(j, ln, d_spring_edge_map[ln][j], d_spring_spec_data[ln][j]);
        }
    }
    return;
} // initializeSprings

void
IBRedundantInitializer::initializeBeams()
{
    if (!d_init_beam_on_level_fcn) return;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_beam_spec_data[ln].resize(num_base_filename);
        for (unsigned int j = 0; j < num_base_filename; ++j)
        {
            d_init_beam_on_level_fcn(j, ln, d_beam_spec_data[ln][j]);
        }
    }
    return;
} // initializeBeams

void
IBRedundantInitializer::initializeTargetPts()
{
    if (!d_init_target_pt_on_level_fcn) return;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_target_spec_data[ln].resize(num_base_filename);
        for (unsigned int j = 0; j < num_base_filename; ++j)
        {
            d_init_target_pt_on_level_fcn(j, ln, d_target_spec_data[ln][j]);
        }
    }
    return;
} // initializeTargetPts

void
IBRedundantInitializer::initializeDirectorAndRods()
{
    if (!d_init_director_and_rod_on_level_fcn) return;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_directors[ln].resize(num_base_filename);
        d_rod_edge_map[ln].resize(num_base_filename);
        d_rod_spec_data[ln].resize(num_base_filename);
        for (unsigned int j = 0; j < num_base_filename; ++j)
        {
            d_init_director_and_rod_on_level_fcn(
                j, ln, d_directors[ln][j], d_rod_edge_map[ln][j], d_rod_spec_data[ln][j]);
        }
    }
    return;
} // initializeDirectorAndRods

void
IBRedundantInitializer::initializeBoundaryMass()
{
    if (!d_init_boundary_mass_on_level_fcn) return;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_bdry_mass_spec_data[ln].resize(num_base_filename);
        for (unsigned int j = 0; j < num_base_filename; ++j)
        {
            d_init_boundary_mass_on_level_fcn(j, ln, d_bdry_mass_spec_data[ln][j]);
        }
    }
    return;
}

void
IBRedundantInitializer::initializeAnchorPts()
{
    if (!d_init_anchor_pt_on_level_fcn) return;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_anchor_spec_data[ln].resize(num_base_filename);
        for (unsigned int j = 0; j < num_base_filename; ++j)
        {
            d_init_anchor_pt_on_level_fcn(j, ln, d_anchor_spec_data[ln][j]);
        }
    }
    return;
}

unsigned int
IBRedundantInitializer::initializeDataOnPatchLevel(const int lag_node_index_idx,
                                                   const unsigned int global_index_offset,
                                                   const unsigned int local_index_offset,
                                                   Pointer<LData> X_data,
                                                   Pointer<LData> U_data,
                                                   const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                   const int level_number,
                                                   const double /*init_data_time*/,
                                                   const bool /*can_be_refined*/,
                                                   const bool /*initial_time*/,
                                                   LDataManager* const /*l_data_manager*/)
{
    // Check if data has been processed.
    init();

    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
    Vector domain_length;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        domain_length[d] = domain_x_upper[d] - domain_x_lower[d];
    }

    // Set the global index offset.  This is equal to the number of Lagrangian
    // indices that have already been initialized on the specified level.
    d_global_index_offset[level_number] = global_index_offset;

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    boost::multi_array_ref<double, 2>& X_array = *X_data->getLocalFormVecArray();
    boost::multi_array_ref<double, 2>& U_array = *U_data->getLocalFormVecArray();
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(ratio);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();

        Pointer<LNodeSetData> index_data = patch->getPatchData(lag_node_index_idx);

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int, int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, hierarchy);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin(); it != patch_vertices.end();
             ++it)
        {
            const std::pair<int, int>& point_idx = (*it);
            const int lagrangian_idx = getCanonicalLagrangianIndex(point_idx, level_number) + global_index_offset;
            const int local_petsc_idx = ++local_idx + local_index_offset;
            const int global_petsc_idx = local_petsc_idx + global_index_offset;

            // Get the coordinates and periodic shifters of the present vertex.
            Point X_real = getVertexPosn(point_idx, level_number);
            Point X = getShiftedVertexPosn(point_idx, level_number, domain_x_lower, domain_x_upper, periodic_shift);
            Vector periodic_displacement = X_real - X;
            IntVector<NDIM> periodic_offset;
            for (int d = 0; d < NDIM; ++d)
            {
                periodic_offset[d] = boost::math::round(periodic_displacement[d] / patch_dx[d]);
            }

            // Ensure that all points are initially within the computational
            // domain.
            for (int d = 0; d < NDIM; ++d)
            {
                if (!periodic_shift[d] && X[d] < domain_x_lower[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                                             << "  encountered node below lower physical boundary\n"
                                             << "  please ensure that all nodes are within the "
                                                "computational domain."
                                             << std::endl);
                }

                if (!periodic_shift[d] && X[d] >= domain_x_upper[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                                             << "  encountered node above upper physical boundary\n"
                                             << "  please ensure that all nodes are within the "
                                                "computational domain."
                                             << std::endl);
                }
            }

            // Set X_array.
            for (int d = 0; d < NDIM; ++d)
            {
                X_array[local_petsc_idx][d] = X[d];
            }

            // Get the index of the cell in which the present vertex is
            // initially located.
            const CellIndex<NDIM> idx = IndexUtilities::getCellIndex(X, grid_geom, ratio);

            // Initialize the specification objects associated with the present
            // vertex.
            std::vector<Pointer<Streamable> > node_data =
                initializeNodeData(point_idx, global_index_offset, level_number);
            for (std::vector<Pointer<Streamable> >::iterator it = node_data.begin(); it != node_data.end(); ++it)
            {
                (*it)->registerPeriodicShift(periodic_offset, periodic_displacement);
            }

            // Create or retrieve a pointer to the LNodeSet associated with the
            // current Cartesian grid cell.
            if (!index_data->isElement(idx))
            {
                index_data->appendItemPointer(idx, new LNodeSet());
            }
            LNodeSet* const node_set = index_data->getItem(idx);
            node_set->push_back(new LNode(lagrangian_idx,
                                          global_petsc_idx,
                                          local_petsc_idx,
                                          /*initial*/ periodic_offset,
                                          /*current*/ periodic_offset,
                                          /*initial*/ periodic_displacement,
                                          /*current*/ periodic_displacement,
                                          node_data));

            // Initialize the velocity of the present vertex.
            std::fill(&U_array[local_petsc_idx][0], &U_array[local_petsc_idx][0] + NDIM, 0.0);
        }
    }
    X_data->restoreArrays();
    U_data->restoreArrays();

    d_level_is_initialized[level_number] = true;

    // If a Lagrangian Silo data writer is registered with the initializer,
    // setup the visualization data corresponding to the present level of the
    // locally refined grid.
    if (d_silo_writer)
    {
        initializeLSiloDataWriter(level_number);
    }
    return local_node_count;
} // initializeDataOnPatchLevel

unsigned int
IBRedundantInitializer::initializeMassDataOnPatchLevel(const unsigned int /*global_index_offset*/,
                                                       const unsigned int local_index_offset,
                                                       Pointer<LData> M_data,
                                                       Pointer<LData> K_data,
                                                       const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                       const int level_number,
                                                       const double /*init_data_time*/,
                                                       const bool /*can_be_refined*/,
                                                       const bool /*initial_time*/,
                                                       LDataManager* const /*l_data_manager*/)
{
    // Check if data has been processed.
    init();

    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    boost::multi_array_ref<double, 1>& M_array = *M_data->getLocalFormArray();
    boost::multi_array_ref<double, 1>& K_array = *K_data->getLocalFormArray();
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int, int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, hierarchy);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin(); it != patch_vertices.end();
             ++it)
        {
            const std::pair<int, int>& point_idx = (*it);
            const int local_petsc_idx = ++local_idx + local_index_offset;

            // Initialize the mass and penalty stiffness coefficient
            // corresponding to the present vertex.
            const BdryMassSpec& spec = getVertexBdryMassSpec(point_idx, level_number);
            const double M = spec.bdry_mass;
            const double K = spec.stiffness;

            // Avoid division by zero at massless nodes.
            if (MathUtilities<double>::equalEps(M, 0.0))
            {
                M_array[local_petsc_idx] = std::numeric_limits<double>::epsilon();
                K_array[local_petsc_idx] = 0.0;
            }
            else
            {
                M_array[local_petsc_idx] = M;
                K_array[local_petsc_idx] = K;
            }
        }
    }
    M_data->restoreArrays();
    K_data->restoreArrays();
    return local_node_count;
} // initializeMassOnPatchLevel

unsigned int
IBRedundantInitializer::initializeDirectorDataOnPatchLevel(const unsigned int /*global_index_offset*/,
                                                           const unsigned int local_index_offset,
                                                           Pointer<LData> D_data,
                                                           const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                           const int level_number,
                                                           const double /*init_data_time*/,
                                                           const bool /*can_be_refined*/,
                                                           const bool /*initial_time*/,
                                                           LDataManager* const /*l_data_manager*/)
{
    // Check if data has been processed.
    init();

    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    boost::multi_array_ref<double, 2>& D_array = *D_data->getLocalFormVecArray();
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int, int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, hierarchy);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin(); it != patch_vertices.end();
             ++it)
        {
            const std::pair<int, int>& point_idx = (*it);
            const int local_petsc_idx = ++local_idx + local_index_offset;

            // Initialize the director corresponding to the present vertex.
            const std::vector<double>& D = getVertexDirectors(point_idx, level_number);
            for (int d = 0; d < 3 * 3; ++d)
            {
                D_array[local_petsc_idx][d] = D[d];
            }
        }
    }
    D_data->restoreArrays();
    return local_node_count;
} // initializeDirectorOnPatchLevel

void
IBRedundantInitializer::tagCellsForInitialRefinement(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                     const int level_number,
                                                     const double /*error_data_time*/,
                                                     const int tag_index)
{
    // Check if data has been processed.
    init();

    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();

    // Loop over all patches in the specified level of the patch level and tag
    // cells for refinement wherever there are vertices assigned to a finer
    // level of the Cartesian grid.
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(ratio);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const Box<NDIM>& patch_box = patch->getBox();

        Pointer<CellData<NDIM, int> > tag_data = patch->getPatchData(tag_index);

        // Tag cells for refinement whenever there are vertices whose initial
        // locations will be within the index space of the given patch, but on
        // the finer levels of the AMR patch hierarchy.
        for (int ln = level_number + 1; ln < d_max_levels; ++ln)
        {
            std::vector<std::pair<int, int> > patch_vertices;
            getPatchVerticesAtLevel(patch_vertices, patch, hierarchy, ln);
            for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin();
                 it != patch_vertices.end();
                 ++it)
            {
                const std::pair<int, int>& point_idx = (*it);

                // Get the coordinates of the present vertex.
                const Point& X = getShiftedVertexPosn(point_idx, ln, domain_x_lower, domain_x_upper, periodic_shift);

                // Get the index of the cell in which the present vertex is
                // initially located.
                const CellIndex<NDIM> i = IndexUtilities::getCellIndex(X, grid_geom, ratio);

                // Tag the cell for refinement.
                if (patch_box.contains(i)) (*tag_data)(i) = 1;
            }
        }
    }
    return;
} // tagCellsForInitialRefinement

void
IBRedundantInitializer::setStructureNamesOnLevel(const int& level_num, const std::vector<std::string>& strct_names)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level_num >= 0);
    TBOX_ASSERT(level_num < d_max_levels);
    TBOX_ASSERT(d_level_is_initialized[level_num]);
#endif
    d_base_filename[level_num] = strct_names;
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBRedundantInitializer::initializeLSiloDataWriter(const int level_number)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level_number >= 0);
    TBOX_ASSERT(level_number < d_max_levels);
    TBOX_ASSERT(d_level_is_initialized[level_number]);
#endif

    // WARNING: This code does not work if the global node offset is nonzero on
    // any of the levels of the locally refined Cartesian grid.
    if (d_global_index_offset[level_number] != 0)
    {
        TBOX_ERROR("This is broken --- please submit a bug report if you encounter this error.\n");
    }

    // WARNING: For now, we just register the visualization data on MPI process
    // 0.  This will fail if the structure is too large to be stored in the
    // memory available to a single MPI process.
    if (SAMRAI_MPI::getRank() == 0)
    {
        for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_num_vertex[level_number][j] > 0)
            {
                const std::string postfix = "_vertices";
                d_silo_writer->registerMarkerCloud(d_base_filename[level_number][j] + postfix,
                                                   d_num_vertex[level_number][j],
                                                   d_vertex_offset[level_number][j],
                                                   level_number);
            }
        }

        bool registered_spring_edge_map = false;
        for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_spring_edge_map[level_number][j].size() > 0)
            {
                registered_spring_edge_map = true;
                const std::string postfix = "_mesh";
                d_silo_writer->registerUnstructuredMesh(
                    d_base_filename[level_number][j] + postfix, d_spring_edge_map[level_number][j], level_number);
            }
        }

        for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_xspring_edge_map[level_number][j].size() > 0)
            {
                const std::string postfix = "_xmesh";
                d_silo_writer->registerUnstructuredMesh(
                    d_base_filename[level_number][j] + postfix, d_xspring_edge_map[level_number][j], level_number);
            }
        }

        for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_rod_edge_map[level_number][j].size() > 0)
            {
                const std::string postfix = registered_spring_edge_map ? "_rod_mesh" : "_mesh";
                d_silo_writer->registerUnstructuredMesh(
                    d_base_filename[level_number][j] + postfix, d_rod_edge_map[level_number][j], level_number);
            }
        }
    }
    return;
} // initializeLSiloDataWriter

void
IBRedundantInitializer::getPatchVertices(std::vector<std::pair<int, int> >& patch_vertices,
                                         const Pointer<Patch<NDIM> > patch,
                                         const Pointer<PatchHierarchy<NDIM> > hierarchy) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(patch->inHierarchy());
#endif
    const int level_number = patch->getPatchLevelNumber();
    getPatchVerticesAtLevel(patch_vertices, patch, hierarchy, level_number);
    return;
} // getPatchVertices

void
IBRedundantInitializer::getPatchVerticesAtLevel(std::vector<std::pair<int, int> >& patch_vertices,
                                                const Pointer<Patch<NDIM> > patch,
                                                const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                const int vertex_level_number) const
{
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
#if !defined(NDEBUG)
    TBOX_ASSERT(patch->inHierarchy());
#endif
    const int level_number = patch->getPatchLevelNumber();
    const Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(ratio);

    // Loop over all of the vertices to determine the indices of those vertices
    // within the present patch.
    //
    // NOTE: This is clearly not the best way to do this, but it will work for
    // now.
    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    for (unsigned int j = 0; j < d_num_vertex[vertex_level_number].size(); ++j)
    {
        for (int k = 0; k < d_num_vertex[vertex_level_number][j]; ++k)
        {
            std::pair<int, int> point_index(j, k);
            const Point& X =
                getShiftedVertexPosn(point_index, vertex_level_number, domain_x_lower, domain_x_upper, periodic_shift);
            const CellIndex<NDIM> idx = IndexUtilities::getCellIndex(X, grid_geom, ratio);
            if (patch_box.contains(idx)) patch_vertices.push_back(point_index);
        }
    }
    return;
} // getPatchVerticesAtLevel

int
IBRedundantInitializer::getCanonicalLagrangianIndex(const std::pair<int, int>& point_index,
                                                    const int level_number) const
{
    return d_vertex_offset[level_number][point_index.first] + point_index.second;
} // getCanonicalLagrangianIndex

Point
IBRedundantInitializer::getVertexPosn(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_vertex_posn[level_number][point_index.first][point_index.second];
} // getVertexPosn

Point
IBRedundantInitializer::getShiftedVertexPosn(const std::pair<int, int>& point_index,
                                             const int level_number,
                                             const double* const domain_x_lower,
                                             const double* const domain_x_upper,
                                             const IntVector<NDIM>& periodic_shift) const
{
    Point X = getVertexPosn(point_index, level_number);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (periodic_shift[d])
        {
            double domain_length = domain_x_upper[d] - domain_x_lower[d];
            while (X[d] < domain_x_lower[d]) X[d] += domain_length;
            while (X[d] >= domain_x_upper[d]) X[d] -= domain_length;
            TBOX_ASSERT(X[d] >= domain_x_lower[d] && X[d] < domain_x_upper[d]);
            X[d] = std::max(X[d], domain_x_lower[d]);
            X[d] = std::min(X[d], domain_x_upper[d] - std::numeric_limits<double>::epsilon());
        }
    }
    return X;
} // getShiftedVertexPosn

const IBRedundantInitializer::TargetSpec&
IBRedundantInitializer::getVertexTargetSpec(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_target_spec_data[level_number][point_index.first][point_index.second];
} // getVertexTargetSpec

const IBRedundantInitializer::AnchorSpec&
IBRedundantInitializer::getVertexAnchorSpec(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_anchor_spec_data[level_number][point_index.first][point_index.second];
} // getVertexAnchorSpec

const IBRedundantInitializer::BdryMassSpec&
IBRedundantInitializer::getVertexBdryMassSpec(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_bdry_mass_spec_data[level_number][point_index.first][point_index.second];
} // getVertexBdryMassSpec

const std::vector<double>&
IBRedundantInitializer::getVertexDirectors(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_directors[level_number][point_index.first][point_index.second];
} // getVertexDirectors

std::pair<int, int>
IBRedundantInitializer::getVertexInstrumentationIndices(const std::pair<int, int>& point_index,
                                                        const int level_number) const
{
    std::map<int, std::pair<int, int> >::const_iterator it =
        d_instrument_idx[level_number][point_index.first].find(point_index.second);
    if (it != d_instrument_idx[level_number][point_index.first].end())
    {
        return it->second;
    }
    else
    {
        return std::make_pair(-1, -1);
    }
} // getVertexInstrumentationIndices

int
IBRedundantInitializer::getVertexSourceIndices(const std::pair<int, int>& point_index, const int level_number) const
{
    std::map<int, int>::const_iterator it = d_source_idx[level_number][point_index.first].find(point_index.second);
    if (it != d_source_idx[level_number][point_index.first].end())
    {
        return it->second;
    }
    else
    {
        return -1;
    }
} // getVertexSourceIndices

std::vector<Pointer<Streamable> >
IBRedundantInitializer::initializeNodeData(const std::pair<int, int>& point_index,
                                           const unsigned int global_index_offset,
                                           const int level_number) const
{
    std::vector<Pointer<Streamable> > node_data;

    const int j = point_index.first;
    const int mastr_idx = getCanonicalLagrangianIndex(point_index, level_number);

    // Initialize any spring specifications associated with the present vertex.
    {
        std::vector<int> slave_idxs, force_fcn_idxs;
        std::vector<std::vector<double> > parameters;
        for (std::multimap<int, Edge>::const_iterator it = d_spring_edge_map[level_number][j].lower_bound(mastr_idx);
             it != d_spring_edge_map[level_number][j].upper_bound(mastr_idx);
             ++it)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(mastr_idx == it->first);
#endif
            // The connectivity information.
            const Edge& e = it->second;
            if (e.first == mastr_idx)
            {
                slave_idxs.push_back(e.second + global_index_offset);
            }
            else
            {
                slave_idxs.push_back(e.first + global_index_offset);
            }

            // The material properties.
            const SpringSpec& spec_data = d_spring_spec_data[level_number][j].find(e)->second;
            parameters.push_back(spec_data.parameters);
            force_fcn_idxs.push_back(spec_data.force_fcn_idx);
        }
        const size_t num_base_filename = d_base_filename[level_number].size();
        for (unsigned int j = 0; j < num_base_filename; ++j)
        {
            for (std::multimap<int, Edge>::const_iterator it =
                     d_xspring_edge_map[level_number][j].lower_bound(mastr_idx);
                 it != d_xspring_edge_map[level_number][j].upper_bound(mastr_idx);
                 ++it)
            {
#if !defined(NDEBUG)
                TBOX_ASSERT(mastr_idx == it->first);
#endif
                // The connectivity information.
                const Edge& e = it->second;
                if (e.first == mastr_idx)
                {
                    slave_idxs.push_back(e.second + global_index_offset);
                }
                else
                {
                    slave_idxs.push_back(e.first + global_index_offset);
                }

                // The material properties.
                const XSpringSpec& spec_data = d_xspring_spec_data[level_number][j].find(e)->second;
                parameters.push_back(spec_data.parameters);
                force_fcn_idxs.push_back(spec_data.force_fcn_idx);
            }
        }
        if (slave_idxs.size() > 0)
        {
            node_data.push_back(new IBSpringForceSpec(mastr_idx, slave_idxs, force_fcn_idxs, parameters));
        }
    }

    // Initialize any beam specifications associated with the present vertex.
    {
        std::vector<std::pair<int, int> > beam_neighbor_idxs;
        std::vector<double> beam_bend_rigidity;
        std::vector<Vector> beam_mesh_dependent_curvature;
        for (std::multimap<int, BeamSpec>::const_iterator it = d_beam_spec_data[level_number][j].lower_bound(mastr_idx);
             it != d_beam_spec_data[level_number][j].upper_bound(mastr_idx);
             ++it)
        {
            const BeamSpec& spec_data = it->second;
            beam_neighbor_idxs.push_back(spec_data.neighbor_idxs);
            beam_bend_rigidity.push_back(spec_data.bend_rigidity);
            beam_mesh_dependent_curvature.push_back(spec_data.curvature);
        }
        if (!beam_neighbor_idxs.empty())
        {
            node_data.push_back(
                new IBBeamForceSpec(mastr_idx, beam_neighbor_idxs, beam_bend_rigidity, beam_mesh_dependent_curvature));
        }
    }

    // Initialize any rod specifications associated with the present vertex.
    {
        std::vector<int> rod_next_idxs;
        std::vector<boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> > rod_material_params;
        for (std::multimap<int, Edge>::const_iterator it = d_rod_edge_map[level_number][j].lower_bound(mastr_idx);
             it != d_rod_edge_map[level_number][j].upper_bound(mastr_idx);
             ++it)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(mastr_idx == it->first);
#endif
            // The connectivity information.
            const Edge& e = it->second;
            if (e.first == mastr_idx)
            {
                rod_next_idxs.push_back(e.second + global_index_offset);
            }
            else
            {
                rod_next_idxs.push_back(e.first + global_index_offset);
            }

            // The material properties.
            const RodSpec& spec_data = d_rod_spec_data[level_number][j].find(e)->second;
            rod_material_params.push_back(spec_data.properties);
        }
        if (!rod_next_idxs.empty())
        {
            node_data.push_back(new IBRodForceSpec(mastr_idx, rod_next_idxs, rod_material_params));
        }
    }

    // Initialize any target point specifications associated with the present
    // vertex.
    {
        const TargetSpec& spec_data = getVertexTargetSpec(point_index, level_number);
        const double kappa_target = spec_data.stiffness;
        const double eta_target = spec_data.damping;
        const Point& X_target = getVertexPosn(point_index, level_number);
        node_data.push_back(new IBTargetPointForceSpec(mastr_idx, kappa_target, eta_target, X_target));
    }

    // Initialize any anchor point specifications associated with the present
    // vertex.
    {
        const AnchorSpec& spec_data = getVertexAnchorSpec(point_index, level_number);
        const bool is_anchor_point = spec_data.is_anchor_point;
        if (is_anchor_point)
        {
            node_data.push_back(new IBAnchorPointSpec(mastr_idx));
        }
    }

    // Initialize any instrumentation specifications associated with the present
    // vertex.
    {
        const std::pair<int, int> inst_idx = getVertexInstrumentationIndices(point_index, level_number);
        if (inst_idx.first != -1 && inst_idx.second != -1)
        {
            node_data.push_back(new IBInstrumentationSpec(mastr_idx, inst_idx.first, inst_idx.second));
        }
    }

    // Initialize any source specifications associated with the present
    // vertex.
    {
        const int source_idx = getVertexSourceIndices(point_index, level_number);
        if (source_idx != -1)
        {
            node_data.push_back(new IBSourceSpec(mastr_idx, source_idx));
        }
    }

    return node_data;
} // initializeNodeData

void
IBRedundantInitializer::getFromInput(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
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

    d_level_is_initialized.resize(d_max_levels, false);
    d_base_filename.resize(d_max_levels);
    d_num_vertex.resize(d_max_levels);
    d_spring_edge_map.resize(d_max_levels);
    d_spring_spec_data.resize(d_max_levels);
    d_xspring_edge_map.resize(d_max_levels);
    d_xspring_spec_data.resize(d_max_levels);
    d_beam_spec_data.resize(d_max_levels);
    d_rod_edge_map.resize(d_max_levels);
    d_rod_spec_data.resize(d_max_levels);
    d_target_spec_data.resize(d_max_levels);
    d_anchor_spec_data.resize(d_max_levels);
    d_bdry_mass_spec_data.resize(d_max_levels);
    d_directors.resize(d_max_levels);

    // Determine the various input file names.
    //
    // Prefer to use the new ``structure_names'' key, but revert to the
    // level-by-level ``base_filenames'' keys if necessary.
    if (db->keyExists("structure_names"))
    {
        const int num_strcts = db->getArraySize("structure_names");
        std::vector<std::string> structure_names(num_strcts);
        db->getStringArray("structure_names", &structure_names[0], num_strcts);
        for (int n = 0; n < num_strcts; ++n)
        {
            const std::string& strct_name = structure_names[n];
            if (db->keyExists(strct_name))
            {
                Pointer<Database> sub_db = db->getDatabase((strct_name));
                if (sub_db->keyExists("level_number"))
                {
                    const int ln = sub_db->getInteger("level_number");
                    if (ln < 0)
                    {
                        TBOX_ERROR(d_object_name << ":  "
                                                 << "Key data `level_number' associated with structure `" << strct_name
                                                 << "' is negative.");
                    }
                    else if (ln > d_max_levels)
                    {
                        TBOX_ERROR(d_object_name << ":  "
                                                 << "Key data `level_number' associated with structure `" << strct_name
                                                 << "' is greater than the expected maximum level number "
                                                 << d_max_levels << ".");
                    }
                    d_base_filename[ln].push_back(strct_name);
                }
                else
                {
                    TBOX_ERROR(d_object_name << ":  "
                                             << "Key data `level_number' not found in structure `" << strct_name
                                             << "' input.");
                }
            }
            else
            {
                TBOX_ERROR(d_object_name << ":  "
                                         << "Key data `" << strct_name << "' not found in input.");
            }
        }
    }
    else if (db->keyExists("structure_levels"))
    {
        const int num_levels = db->getArraySize("structure_levels");
        std::vector<int> strct_levels(num_levels);
        db->getIntegerArray("structure_levels", &strct_levels[0], num_levels);
        std::vector<int> num_strcts(num_levels);
        db->getIntegerArray("num_structures_levels", &num_strcts[0], num_levels);

        for (int k = 0; k < num_levels; ++k)
        {
            const int ln = strct_levels[k];
            const int num_strcts_on_ln = num_strcts[k];
            d_base_filename[ln].resize(num_strcts_on_ln);

            for (int s = 0; s < num_strcts_on_ln; ++s)
            {
                std::ostringstream strct_name_stream;
                strct_name_stream << "body_" << s << "_ln_" << ln;
                d_base_filename[ln][s] = strct_name_stream.str();
            }
        }
    }
    else
    {
        for (int ln = 0; ln < d_max_levels; ++ln)
        {
            std::ostringstream db_key_name_stream;
            db_key_name_stream << "base_filenames_" << ln;
            const std::string db_key_name = db_key_name_stream.str();
            if (db->keyExists(db_key_name))
            {
                const int num_files = db->getArraySize(db_key_name);
                d_base_filename[ln].resize(num_files);
                db->getStringArray(db_key_name, &d_base_filename[ln][0], num_files);
            }
            else
            {
                TBOX_WARNING(d_object_name << ":  "
                                           << "Key data `" << db_key_name << "' not found in input.");
            }
        }
    }

    return;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
