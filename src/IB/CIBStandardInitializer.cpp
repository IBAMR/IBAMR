//CIBStandardInitilizer.cpp for suspension of strucrutes modified by Baky and Floren Balboa Aug 2015

// Filename: CIBStandardInitializer.cpp
// Created on 22 Nov 2006 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Geometry>

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
#include "boost/multi_array.hpp"
#include "ibamr/IBAnchorPointSpec.h"
#include "ibamr/IBBeamForceSpec.h"
#include "ibamr/IBInstrumentationSpec.h"
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


#include "ibamr/CIBStandardInitializer.h"

namespace IBTK
{
  class LDataManager;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////



namespace
  {
    inline std::string discard_comments(const std::string& input_string)
    {
      // Create a copy of the input string, but without any text following a '!',
      // '#', or '%' character.
      std::string output_string = input_string;
      std::istringstream string_stream;

      // Discard any text following a '!' character.
      string_stream.str(output_string);
      std::getline(string_stream, output_string, '!');
      string_stream.clear();

      // Discard any text following a '#' character.
      string_stream.str(output_string);
      std::getline(string_stream, output_string, '#');
      string_stream.clear();

      // Discard any text following a '%' character.
      string_stream.str(output_string);
      std::getline(string_stream, output_string, '%');
      string_stream.clear();
      return output_string;
    } // discard_comments

    inline int round(double x)
    {
      return floor(x + 0.5);
    } // round
  }

/////////////////////////////// PUBLIC ///////////////////////////////////////

  CIBStandardInitializer::CIBStandardInitializer(const std::string& object_name, Pointer<Database> input_db)
    : d_object_name(object_name), d_use_file_batons(true), d_max_levels(-1), d_level_is_initialized(),
      d_silo_writer(NULL), d_base_filename(), d_length_scale_factor(1.0), d_posn_shift(Vector::Zero()), d_num_vertex(),
      d_vertex_offset(), d_vertex_posn(), d_enable_springs(), d_spring_edge_map(), d_spring_spec_data(),
      d_using_uniform_spring_stiffness(), d_uniform_spring_stiffness(), d_using_uniform_spring_rest_length(),
      d_uniform_spring_rest_length(), d_using_uniform_spring_force_fcn_idx(), d_uniform_spring_force_fcn_idx(),
      d_enable_xsprings(), d_xspring_edge_map(), d_xspring_spec_data(), d_using_uniform_xspring_stiffness(),
      d_uniform_xspring_stiffness(), d_using_uniform_xspring_rest_length(), d_uniform_xspring_rest_length(),
      d_using_uniform_xspring_force_fcn_idx(), d_uniform_xspring_force_fcn_idx(), d_enable_beams(), d_beam_spec_data(),
      d_using_uniform_beam_bend_rigidity(), d_uniform_beam_bend_rigidity(), d_using_uniform_beam_curvature(),
      d_uniform_beam_curvature(), d_enable_rods(), d_rod_edge_map(), d_rod_spec_data(),
      d_using_uniform_rod_properties(), d_uniform_rod_properties(), d_enable_target_points(), d_target_spec_data(),
      d_using_uniform_target_stiffness(), d_uniform_target_stiffness(), d_using_uniform_target_damping(),
      d_uniform_target_damping(), d_enable_anchor_points(), d_anchor_spec_data(), d_enable_bdry_mass(),
      d_bdry_mass_spec_data(), d_using_uniform_bdry_mass(), d_uniform_bdry_mass(),
      d_using_uniform_bdry_mass_stiffness(), d_uniform_bdry_mass_stiffness(), d_directors(), d_enable_instrumentation(),
      d_instrument_idx(), d_enable_sources(), d_source_idx(), d_global_index_offset()
  {
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
#endif

    //Register the specification objects with the StreamableManager class.
    IBAnchorPointSpec::registerWithStreamableManager();
    IBBeamForceSpec::registerWithStreamableManager();
    IBInstrumentationSpec::registerWithStreamableManager();
    IBRodForceSpec::registerWithStreamableManager();
    IBSourceSpec::registerWithStreamableManager();
    IBSpringForceSpec::registerWithStreamableManager();
    IBTargetPointForceSpec::registerWithStreamableManager();

    // Initialize object with data read from the input database.
    getFromInput(input_db);
    
    // Check to see if we are starting from a restart file.
    RestartManager* restart_manager = RestartManager::getManager();
    const bool is_from_restart = restart_manager->isFromRestart();
    
    // Process the input files only if we are not starting from a restart file.
    if (!is_from_restart)
    {
      // Process the vertex information.
      readVertexFiles(".vertex");

      // Process the spring information.
      readSpringFiles(".spring", /*input_uses_global_idxs*/ false);

      // Process the crosslink spring ("x-spring") information.
      readXSpringFiles(".xspring", /*input_uses_global_idxs*/ true);

      // Process the beam information.
      readBeamFiles(".beam", /*input_uses_global_idxs*/ false);

      // Process the rod information.
      readRodFiles(".rod", /*input_uses_global_idxs*/ false);

      // Process the target point information.
      readTargetPointFiles(".target");

      // Process the anchor point information.
      readAnchorPointFiles(".anchor");

      // Process the mass information.
      readBoundaryMassFiles(".mass");

      // Process the directors information.
      readDirectorFiles(".director");

      // Process the instrumentation information.
      readInstrumentationFiles(".inst");

      // Process the source information.
      readSourceFiles(".source");
    }
    return;
  } // CIBStandardInitializer

  CIBStandardInitializer::~CIBStandardInitializer()
  {
    pout << d_object_name << ":  Deallocating initialization data.\n";
    return;
  } // ~CIBStandardInitializer

  void CIBStandardInitializer::registerLSiloDataWriter(Pointer<LSiloDataWriter> silo_writer)
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

  bool CIBStandardInitializer::getLevelHasLagrangianData(const int level_number, const bool /*can_be_refined*/) const
  {
    return !d_num_vertex[level_number].empty();
  } // getLevelHasLagrangianData

  unsigned int
  CIBStandardInitializer::computeGlobalNodeCountOnPatchLevel(const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
							    const int level_number,
							    const double /*init_data_time*/,
							    const bool /*can_be_refined*/,
							    const bool /*initial_time*/)
  {
    return std::accumulate(d_num_vertex[level_number].begin(), d_num_vertex[level_number].end(), 0);
  }

  unsigned int CIBStandardInitializer::computeLocalNodeCountOnPatchLevel(const Pointer<PatchHierarchy<NDIM> > hierarchy,
									const int level_number,
									const double /*init_data_time*/,
									const bool can_be_refined,
									const bool /*initial_time*/)
  {
    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();

    // Loop over all patches in the specified level of the patch level and count
    // the number of local vertices.
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
      Pointer<Patch<NDIM> > patch = level->getPatch(p());

      // Count the number of vertices whose initial locations will be within
      // the given patch.
      std::vector<std::pair<int, int> > patch_vertices;
      getPatchVertices(
	patch_vertices, patch, level_number, can_be_refined, domain_x_lower, domain_x_upper, periodic_shift);
      local_node_count += patch_vertices.size();
    }
    return local_node_count;
  } // computeLocalNodeCountOnPatchLevel

  void CIBStandardInitializer::initializeStructureIndexingOnPatchLevel(
    std::map<int, std::string>& strct_id_to_strct_name_map,
    std::map<int, std::pair<int, int> >& strct_id_to_lag_idx_range_map,
    const int level_number,
    const double /*init_data_time*/,
    const bool /*can_be_refined*/,
    const bool /*initial_time*/,
    LDataManager* const /*l_data_manager*/)
  {
    int offset = 0;
    for (int j = 0; j < static_cast<int>(d_base_filename[level_number].size()); ++j)
    {
      strct_id_to_strct_name_map[j] = d_base_filename[level_number][j];
      strct_id_to_lag_idx_range_map[j] = std::make_pair(offset, offset + d_num_vertex[level_number][j]);
      offset += d_num_vertex[level_number][j];
    }
    return;
  } // initializeStructureIndexingOnPatchLevel

  unsigned int CIBStandardInitializer::initializeDataOnPatchLevel(const int lag_node_index_idx,
								 const unsigned int global_index_offset,
								 const unsigned int local_index_offset,
								 Pointer<LData> X_data,
								 Pointer<LData> U_data,
								 const Pointer<PatchHierarchy<NDIM> > hierarchy,
								 const int level_number,
								 const double /*init_data_time*/,
								 const bool can_be_refined,
								 const bool /*initial_time*/,
								 LDataManager* const /*l_data_manager*/)
  {
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
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
      Pointer<Patch<NDIM> > patch = level->getPatch(p());
      const Box<NDIM>& patch_box = patch->getBox();
      const CellIndex<NDIM>& patch_lower = patch_box.lower();
      const CellIndex<NDIM>& patch_upper = patch_box.upper();
      const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
      const double* const patch_x_lower = patch_geom->getXLower();
      const double* const patch_x_upper = patch_geom->getXUpper();
      const double* const patch_dx = patch_geom->getDx();

      Pointer<LNodeSetData> index_data = patch->getPatchData(lag_node_index_idx);

      // Initialize the vertices whose initial locations will be within the
      // given patch.
      std::vector<std::pair<int, int> > patch_vertices;
      getPatchVertices(
	patch_vertices, patch, level_number, can_be_refined, domain_x_lower, domain_x_upper, periodic_shift);
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
	  periodic_offset[d] = round(periodic_displacement[d] / patch_dx[d]);
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
		       "computational domain." << std::endl);
	  }

	  if (!periodic_shift[d] && X[d] >= domain_x_upper[d])
	  {
	    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
		       << "  encountered node above upper physical boundary\n"
		       << "  please ensure that all nodes are within the "
		       "computational domain." << std::endl);
	  }
	}

	// Set X_array.
	for (int d = 0; d < NDIM; ++d)
	{
	  X_array[local_petsc_idx][d] = X[d];
	}

	// Get the index of the cell in which the present vertex is
	// initially located.
	const CellIndex<NDIM> idx =
	  IndexUtilities::getCellIndex(X, patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper);

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
	node_set->push_back(new LNode(
			      lagrangian_idx, global_petsc_idx, local_petsc_idx, periodic_offset, periodic_displacement, node_data));

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

  unsigned int CIBStandardInitializer::initializeMassDataOnPatchLevel(const unsigned int /*global_index_offset*/,
								     const unsigned int local_index_offset,
								     Pointer<LData> M_data,
								     Pointer<LData> K_data,
								     const Pointer<PatchHierarchy<NDIM> > hierarchy,
								     const int level_number,
								     const double /*init_data_time*/,
								     const bool can_be_refined,
								     const bool /*initial_time*/,
								     LDataManager* const /*l_data_manager*/)
  {
    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    boost::multi_array_ref<double, 1>& M_array = *M_data->getLocalFormArray();
    boost::multi_array_ref<double, 1>& K_array = *K_data->getLocalFormArray();
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
      Pointer<Patch<NDIM> > patch = level->getPatch(p());

      // Initialize the vertices whose initial locations will be within the
      // given patch.
      std::vector<std::pair<int, int> > patch_vertices;
      getPatchVertices(
	patch_vertices, patch, level_number, can_be_refined, domain_x_lower, domain_x_upper, periodic_shift);
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

  unsigned int CIBStandardInitializer::initializeDirectorDataOnPatchLevel(const unsigned int /*global_index_offset*/,
									 const unsigned int local_index_offset,
									 Pointer<LData> D_data,
									 const Pointer<PatchHierarchy<NDIM> > hierarchy,
									 const int level_number,
									 const double /*init_data_time*/,
									 const bool can_be_refined,
									 const bool /*initial_time*/,
									 LDataManager* const /*l_data_manager*/)
  {
    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    boost::multi_array_ref<double, 2>& D_array = *D_data->getLocalFormVecArray();
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
      Pointer<Patch<NDIM> > patch = level->getPatch(p());

      // Initialize the vertices whose initial locations will be within the
      // given patch.
      std::vector<std::pair<int, int> > patch_vertices;
      getPatchVertices(
	patch_vertices, patch, level_number, can_be_refined, domain_x_lower, domain_x_upper, periodic_shift);
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

  void CIBStandardInitializer::tagCellsForInitialRefinement(const Pointer<PatchHierarchy<NDIM> > hierarchy,
							   const int level_number,
							   const double /*error_data_time*/,
							   const int tag_index)
  {
    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();

    // Loop over all patches in the specified level of the patch level and tag
    // cells for refinement wherever there are vertices assigned to a finer
    // level of the Cartesian grid.
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(level->getRatio());
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
      Pointer<Patch<NDIM> > patch = level->getPatch(p());
      const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
      const Box<NDIM>& patch_box = patch->getBox();
      const CellIndex<NDIM>& patch_lower = patch_box.lower();
      const CellIndex<NDIM>& patch_upper = patch_box.upper();
      const double* const x_lower = patch_geom->getXLower();
      const double* const x_upper = patch_geom->getXUpper();
      const double* const dx = patch_geom->getDx();

      Pointer<CellData<NDIM, int> > tag_data = patch->getPatchData(tag_index);

      // Tag cells for refinement whenever there are vertices whose initial
      // locations will be within the index space of the given patch, but on
      // the finer levels of the AMR patch hierarchy.
      const bool can_be_refined = level_number + 2 < d_max_levels;
      for (int ln = level_number + 1; ln < d_max_levels; ++ln)
      {
	std::vector<std::pair<int, int> > patch_vertices;
	getPatchVertices(patch_vertices, patch, ln, can_be_refined, domain_x_lower, domain_x_upper, periodic_shift);
	for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin();
	     it != patch_vertices.end();
	     ++it)
	{
	  const std::pair<int, int>& point_idx = (*it);

	  // Get the coordinates of the present vertex.
	  const Point& X = getShiftedVertexPosn(point_idx, ln, domain_x_lower, domain_x_upper, periodic_shift);

	  // Get the index of the cell in which the present vertex is
	  // initially located.
	  const CellIndex<NDIM> i =
	    IndexUtilities::getCellIndex(X, x_lower, x_upper, dx, patch_lower, patch_upper);

	  // Tag the cell for refinement.
	  if (patch_box.contains(i)) (*tag_data)(i) = 1;
	}
      }
    }
    return;
  } // tagCellsForInitialRefinement

 //Baky modifications
Point CIBStandardInitializer::getPrototypeVertexPosn(const int level_number, const int structID, const int nodeID) const
{
    int prototypeID = d_struct_prototype_map[level_number][structID];
    return d_vertex_posn[level_number][prototypeID][nodeID];
} // getVertexPosn

Point CIBStandardInitializer::getInitialVertexPosn(const int level_number, const int structID, const int nodeID) const
{
    std::pair<int, int> point_index(structID,nodeID);
    return getVertexPosn(point_index, level_number);
}

Eigen::Vector3d CIBStandardInitializer::getInitialCOMStructure(const int level_number, const int structID) const
{
    return  d_X_com[level_number][structID];
};

Eigen::Quaterniond* CIBStandardInitializer::getStructureQuaternion(const int level_number, const int structID)
{
    return &d_quaternion[level_number][structID];
} // getVertexPosn

void CIBStandardInitializer::getClonesParameters(int& num_structs_types, 
						 std::vector<int>& structs_clones_num)
{
    num_structs_types = d_num_structs_types;
    structs_clones_num = d_structs_clones_num;
};
/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

  void CIBStandardInitializer::initializeLSiloDataWriter(const int level_number)
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



void CIBStandardInitializer::readClonesCenterPositions(const std::string& fileNameCenters,
						       const unsigned num_clones,
						       double *coords,
						       double *direction)
{
    std::ifstream file_stream;
    file_stream.open(fileNameCenters.c_str(), std::ios::in);
    if (file_stream.is_open())
    {
	std::string line_string;
	std::ifstream fileCenters(fileNameCenters.c_str());
	for(unsigned i=0; i<num_clones; ++i)
	{
	    getline(fileCenters,line_string);
	    std::istringstream line_stream(line_string);
	    //read positions
	    for (unsigned int d = 0; d < NDIM; ++d)
	    {
		if (!(line_stream >> coords[i*NDIM+d]))
		{
		    TBOX_ERROR("CIBStandardInitializer:: Invalid entry in structure center position file encountered on line "
			       << i + 2 << " of file " << fileNameCenters << std::endl);
		}
	    }
	    //read orientation directions
	    for (unsigned int d = 0; d < NDIM; ++d)
	    {
		if (!(line_stream >> direction[i*NDIM+d]))
		{
		    //pout<<"assuming the defaulf orientation vector (1,0,0) of the structure clone encountered on line"
		    //<< i + 2 << " of file " << fileNameCenters << std::endl;
		    direction[i*NDIM] = 1.0;
		    for (unsigned int dd = 1; dd < NDIM; ++dd) direction[i*NDIM+dd]=0;
		    break;
		}
	    }
	}
    }
    else
    {
	if (num_clones > 1)
	{
	    TBOX_ERROR("CIBStandardInitializer:: Cannot find required structure center position file: " <<  fileNameCenters << std::endl);
	}else
	{
	    //use a default position and a default orientation for a single structure
	    for (unsigned int d = 0; d < NDIM; ++d)
	    {
		coords[d] = 0;
		direction[d] = 0;
	    }
	    direction[0] = 1.0;
	}

    }
    file_stream.close();
}

void CIBStandardInitializer::readVertexFiles(const std::string& extension)
{
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      const size_t num_base_filename = d_base_filename[ln].size();
      d_num_vertex[ln].resize(num_base_filename, 0);
      d_vertex_offset[ln].resize(num_base_filename, std::numeric_limits<int>::max());
      d_vertex_posn[ln].resize(num_base_filename);
      d_quaternion[ln].resize(num_base_filename);
      d_X_com[ln].resize(num_base_filename);
      d_struct_prototype_map[ln].resize(num_base_filename);
      
      unsigned struct_begin_idx=0;
      for (unsigned itype = 0; itype < d_num_structs_types; ++itype)
      {
	  d_struct_prototype_map[ln][struct_begin_idx] = struct_begin_idx;
	  const unsigned num_clones =  d_structs_clones_num[itype];
	  double* shift =new double[NDIM*num_clones];
	  double* dircs =new double[NDIM*num_clones];
	  //center of mass of the structure prototype
	  Eigen::Vector3d X_com = Eigen::Vector3d::Zero();

	  //Read initial positions of all clones for the given structure
	  std::string fileNameCenters = d_base_filename[ln][struct_begin_idx];

	  if (fileNameCenters.size()) fileNameCenters.erase(fileNameCenters.size() - 1); //erase number zero
	  fileNameCenters +=".clones";


	  //unsigned int j = 0; //Baky - reading file only once

	  // Wait for the previous MPI process to finish reading the current file.
	  if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, struct_begin_idx);

	  readClonesCenterPositions(fileNameCenters,num_clones, shift, dircs);

	  if (itype == 0)
	  {
	      d_vertex_offset[ln][0] = 0;
	  }
	  else
	  {
	      d_vertex_offset[ln][struct_begin_idx] = d_vertex_offset[ln][struct_begin_idx - 1] + d_num_vertex[ln][struct_begin_idx - 1];
	  }

	  //BAKY - Read vertex file
	  std::string vertex_filename = d_base_filename[ln][struct_begin_idx];
	  if (vertex_filename.size()) vertex_filename.erase(vertex_filename.size() - 1); //erase number zero
	  vertex_filename +=extension;

	  std::ifstream file_stream;
	  file_stream.open(vertex_filename.c_str(), std::ios::in);
	  if (file_stream.is_open())
	  {
	      plog << d_object_name << ":  "
		   << "processing vertex data from ASCII input file named " << vertex_filename << std::endl
		   << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	      
	      // The first entry in the file is the number of vertices.
	      if (!std::getline(file_stream, line_string))
	      {
		  TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
			     "before line 1 of file " << vertex_filename << std::endl);
	      }
	      else
	      {
		  line_string = discard_comments(line_string);
		  std::istringstream line_stream(line_string);
		  if (!(line_stream >> d_num_vertex[ln][struct_begin_idx]))
		  {
		      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
				 "encountered on line 1 of file " << vertex_filename << std::endl);
		  }
	      }

	      if (d_num_vertex[ln][struct_begin_idx] <= 0)
	      {
		  TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
			     << vertex_filename << std::endl);
	      }
	      
	      // Each successive line provides the initial position of each
	      // vertex in the input file.
	      d_vertex_posn[ln][struct_begin_idx].resize(d_num_vertex[ln][struct_begin_idx]);

	      for (int k = 0; k < d_num_vertex[ln][struct_begin_idx]; ++k)
	      {
		  Point& X = d_vertex_posn[ln][struct_begin_idx][k];
		  if (!std::getline(file_stream, line_string))
		  {
		      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k + 2
				 << " of file " << vertex_filename << std::endl);
		  }
		  else
		  {
		      line_string = discard_comments(line_string);
		      std::istringstream line_stream(line_string);
		      for (unsigned int d = 0; d < NDIM; ++d)
		      {
			  if (!(line_stream >> X[d]))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
					 << k + 2 << " of file " << vertex_filename << std::endl);
			  }
		      }
		  }
		  for (unsigned int d = 0; d < NDIM; ++d) X_com[d] +=X[d];
	      }

	      X_com /=d_num_vertex[ln][struct_begin_idx];

	      //put protoypes CM in the origin
	      for (int k = 0; k < d_num_vertex[ln][struct_begin_idx]; ++k)
	      {
		  Point& X = d_vertex_posn[ln][struct_begin_idx][k];
		  for (unsigned int d = 0; d < NDIM; ++d) 
		  {
		      X[d] -=X_com[d];
		      X[d] *=d_length_scale_factor;
		  }
	      }

	      // Close the input file.
	      file_stream.close();

	      plog << d_object_name << ":  "
		   << "read " << d_num_vertex[ln][struct_begin_idx] << " vertices from ASCII input file named " << vertex_filename
		   << std::endl << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	  }
	  else
	  {
	      TBOX_ERROR(d_object_name << ":\n  Cannot find required vertex file: " << vertex_filename << std::endl);
	  }

	  // Free the next MPI process to start reading the current file.
	  if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, struct_begin_idx);
      
        
	  //Baky distribute other structures  
	  for (unsigned j = 1; j < num_clones; ++j)
	  {
	      d_struct_prototype_map[ln][struct_begin_idx+j]=struct_begin_idx;
	      d_vertex_offset[ln][struct_begin_idx+j] = d_vertex_offset[ln][struct_begin_idx+j - 1] + d_num_vertex[ln][struct_begin_idx+j-1];	      
	      d_num_vertex[ln][struct_begin_idx+j]=d_num_vertex[ln][struct_begin_idx];
	      
	      // Each successive line provides the initial position of each
	      // vertex in the input file.

	      //create rotation according to input
	      Eigen::Vector3d x_axis = Eigen::Vector3d::Zero();
	      x_axis[0]=1.0;
	      Eigen::Vector3d dir_vec = Eigen::Vector3d::Zero();
	      for (unsigned int d = 0; d < NDIM; ++d) dir_vec[d] = dircs[j*NDIM+d];
	      double vecnorm=dir_vec.norm();
	      if (!vecnorm)  TBOX_ERROR("CIBStandardInitializer:: Direction vector cannot be zero!" << std::endl);
	      dir_vec /=vecnorm;
	      Eigen::Quaterniond q;
	      q.setFromTwoVectors(x_axis, dir_vec);
	      Eigen::Matrix3d rot_mat = (q.normalized()).toRotationMatrix();
	      d_quaternion[ln][struct_begin_idx+j]=q.normalized();

	      for (unsigned int d = 0; d < NDIM; ++d) d_X_com[ln][struct_begin_idx+j][d] = d_length_scale_factor * (shift[j*NDIM+d] + d_posn_shift[d]);

	  }//j


	  // correct zero structure coords
	  //create rotation according to input
	  Eigen::Vector3d x_axis = Eigen::Vector3d::Zero();
	  x_axis[0]=1.0;
	  Eigen::Vector3d dir_vec = Eigen::Vector3d::Zero();
	  for (unsigned int d = 0; d < NDIM; ++d) dir_vec[d] = dircs[d];
	  double vecnorm=dir_vec.norm();
	  if (!vecnorm)  TBOX_ERROR("CIBStandardInitializer:: Direction vector cannot be zero!" << std::endl);
	  dir_vec /=vecnorm;
	  Eigen::Quaterniond q;
	  q.setFromTwoVectors(x_axis, dir_vec);
	  Eigen::Matrix3d rot_mat = (q.normalized()).toRotationMatrix();
	  d_quaternion[ln][struct_begin_idx]=q.normalized();

	  // for (int k = 0; k < d_num_vertex[ln][struct_begin_idx]; ++k)
	  // {
	  //     Point& X0 = d_vertex_posn[ln][struct_begin_idx][k];
	  //     Eigen::Vector3d x0 = Eigen::Vector3d::Zero();
	  //     Eigen::Vector3d x  = Eigen::Vector3d::Zero();
	  //     for (unsigned int d = 0; d < NDIM; ++d) x0[d] = X0[d]-X_com[d]; 
	  //     x = rot_mat*x0;

	  //     for (unsigned int d = 0; d < NDIM; ++d)
	  //     {
	  // 	  X0[d] = d_length_scale_factor * (x[d] + shift[d]+ d_posn_shift[d]);
	  //     }
	  // }
	  for (unsigned int d = 0; d < NDIM; ++d) d_X_com[ln][struct_begin_idx][d] = d_length_scale_factor * (shift[d]+ d_posn_shift[d]);

	  struct_begin_idx +=num_clones;
	  delete [] shift;
	  delete [] dircs;
      }//itype
    }//ln
    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();

    return;
} // readVertexFiles
  

  void CIBStandardInitializer::readSpringFiles(const std::string& extension, const bool input_uses_global_idxs)
  {
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      const size_t num_base_filename = d_base_filename[ln].size();
      d_spring_edge_map[ln].resize(num_base_filename);
      d_spring_spec_data[ln].resize(num_base_filename);

      unsigned struct_begin_idx=0;
      unsigned j=0;
      for(unsigned itype=0;itype< d_num_structs_types;itype++)
      {
	  for(int istruct=0;istruct < d_structs_clones_num[itype];istruct++)
	  {
	      bool warned = false;
	      
	      // Determine min/max index ranges.
	      const int min_idx = 0;
	      const int max_idx =
		  (input_uses_global_idxs ? std::accumulate(d_num_vertex[ln].begin(), d_num_vertex[ln].end(), 0) :
		   d_num_vertex[ln][j]);
	      
	      // Wait for the previous MPI process to finish reading the current file.
	      if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);
	      
	      // Ensure that the file exists.
	      std::string spring_filename = d_base_filename[ln][struct_begin_idx];
	      if (spring_filename.size()) spring_filename.erase(spring_filename.size() - 1); //erase number zero
	      spring_filename +=extension;

	      std::ifstream file_stream;
	      file_stream.open(spring_filename.c_str(), std::ios::in);
	      if (file_stream.is_open())
	      {
		  plog << d_object_name << ":  "
		       << "processing spring data from ASCII input file named " << spring_filename << std::endl
		       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

		  // The first line in the file indicates the number of edges in the input
		  // file.
		  int num_edges = -1;
		  if (!std::getline(file_stream, line_string))
		  {
		      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
				 "before line 1 of file " << spring_filename << std::endl);
		  }
		  else
		  {
		      line_string = discard_comments(line_string);
		      std::istringstream line_stream(line_string);
		      if (!(line_stream >> num_edges))
		      {
			  TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
				     "encountered on line 1 of file " << spring_filename << std::endl);
		      }
		  }

		  if (num_edges <= 0)
		  {
		      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
				 << spring_filename << std::endl);
		  }

		  // Each successive line provides the connectivity and material parameter
		  // information for each spring in the structure.
		  for (int k = 0; k < num_edges; ++k)
		  {
		      Edge e;
		      std::vector<double> parameters(2);
		      int force_fcn_idx = 0;
		      if (!std::getline(file_stream, line_string))
		      {
			  TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k + 2
				     << " of file " << spring_filename << std::endl);
		      }
		      else
		      {
			  line_string = discard_comments(line_string);
			  std::istringstream line_stream(line_string);
			  if (!(line_stream >> e.first))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << spring_filename << std::endl);
			  }
			  else if ((e.first < min_idx) || (e.first >= max_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << spring_filename << std::endl << "  vertex index "
					 << e.first << " is out of range" << std::endl);
			  }

			  if (!(line_stream >> e.second))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << spring_filename << std::endl);
			  }
			  else if ((e.second < min_idx) || (e.second >= max_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << spring_filename << std::endl << "  vertex index "
					 << e.second << " is out of range" << std::endl);
			  }

			  if (!(line_stream >> parameters[0]))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << spring_filename << std::endl);
			  }
			  else if (parameters[0] < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << spring_filename << std::endl
					 << "  spring constant is negative" << std::endl);
			  }

			  if (!(line_stream >> parameters[1]))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << spring_filename << std::endl);
			  }
			  else if (parameters[1] < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << spring_filename << std::endl
					 << "  spring resting length is negative" << std::endl);
			  }
			  parameters[1] *= d_length_scale_factor;

			  if (!(line_stream >> force_fcn_idx))
			  {
			      force_fcn_idx = 0; // default force function specification.
			  }

			  double param;
			  while (line_stream >> param)
			  {
			      parameters.push_back(param);
			  }
		      }

		      // Modify kappa and length according to whether uniform
		      // values are to be employed for this particular structure.
		      if (d_using_uniform_spring_stiffness[ln][j])
		      {
			  parameters[0] = d_uniform_spring_stiffness[ln][j];
		      }
		      if (d_using_uniform_spring_rest_length[ln][j])
		      {
			  parameters[1] = d_uniform_spring_rest_length[ln][j];
		      }
		      if (d_using_uniform_spring_force_fcn_idx[ln][j])
		      {
			  force_fcn_idx = d_uniform_spring_force_fcn_idx[ln][j];
		      }

		      // Check to see if the spring constant is zero and, if so,
		      // emit a warning.
		      if (!warned && d_enable_springs[ln][j] &&
			  (parameters[0] == 0.0 || MathUtilities<double>::equalEps(parameters[0], 0.0)))
		      {
			  TBOX_WARNING(d_object_name << ":\n  Spring with zero spring constant "
				       "encountered in ASCII input file named " << spring_filename << "."
				       << std::endl);
			  warned = true;
		      }

		      // Correct the edge numbers to be in the global Lagrangian indexing
		      // scheme.
		      if (!input_uses_global_idxs)
		      {
			  e.first += d_vertex_offset[ln][j];
			  e.second += d_vertex_offset[ln][j];
		      }

		      // Initialize the map data corresponding to the present edge.
		      //
		      // Note that in the edge map, each edge is associated with only the
		      // first vertex.
		      if (e.first > e.second)
		      {
			  std::swap<int>(e.first, e.second);
		      }
		      d_spring_edge_map[ln][j].insert(std::make_pair(e.first, e));
		      SpringSpec spec_data;
		      spec_data.parameters = parameters;
		      spec_data.force_fcn_idx = force_fcn_idx;
		      d_spring_spec_data[ln][j].insert(std::make_pair(e, spec_data));
		  }

		  // Close the input file.
		  file_stream.close();

		  plog << d_object_name << ":  "
		       << "read " << num_edges << " edges from ASCII input file named " << spring_filename << std::endl
		       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	      }

	      // Free the next MPI process to start reading the current file.
	      if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);

	      j++;
	  }
	  struct_begin_idx += d_structs_clones_num[itype];
      }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
  } // readSpringFiles

  void CIBStandardInitializer::readXSpringFiles(const std::string& extension, const bool input_uses_global_idxs)
  {
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      const size_t num_base_filename = d_base_filename[ln].size();
      d_xspring_edge_map[ln].resize(num_base_filename);
      d_xspring_spec_data[ln].resize(num_base_filename);

      unsigned struct_begin_idx=0;

      unsigned j = 0;
      for(unsigned itype=0;itype< d_num_structs_types;itype++)
      {
	  for(int istruct=0;istruct < d_structs_clones_num[itype];istruct++)
	  {
	      bool warned = false;

	      // Determine min/max index ranges.
	      const int min_idx = 0;
	      const int max_idx =
		  (input_uses_global_idxs ? std::accumulate(d_num_vertex[ln].begin(), d_num_vertex[ln].end(), 0) :
		   d_num_vertex[ln][j]);

	      // Wait for the previous MPI process to finish reading the current file.
	      if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

	      // Ensure that the file exists.
	      std::string xspring_filename = d_base_filename[ln][struct_begin_idx];
	      if (xspring_filename.size()) xspring_filename.erase(xspring_filename.size() - 1); 
	      xspring_filename +=extension;

	      std::ifstream file_stream;
	      file_stream.open(xspring_filename.c_str(), std::ios::in);
	      if (file_stream.is_open())
	      {
		  plog << d_object_name << ":  "
		       << "processing crosslink spring data from ASCII input file named " << xspring_filename << std::endl
		       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

		  // The first line in the file indicates the number of edges in the input
		  // file.
		  int num_edges = -1;
		  if (!std::getline(file_stream, line_string))
		  {
		      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
				 "before line 1 of file " << xspring_filename << std::endl);
		  }
		  else
		  {
		      line_string = discard_comments(line_string);
		      std::istringstream line_stream(line_string);
		      if (!(line_stream >> num_edges))
		      {
			  TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
				     "encountered on line 1 of file " << xspring_filename << std::endl);
		      }
		  }

		  if (num_edges <= 0)
		  {
		      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
				 << xspring_filename << std::endl);
		  }

		  // Each successive line provides the connectivity and material parameter
		  // information for each crosslink spring in the structure.
		  for (int k = 0; k < num_edges; ++k)
		  {
		      Edge e;
		      std::vector<double> parameters(2);
		      int force_fcn_idx = 0;
		      if (!std::getline(file_stream, line_string))
		      {
			  TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k + 2
				     << " of file " << xspring_filename << std::endl);
		      }
		      else
		      {
			  line_string = discard_comments(line_string);
			  std::istringstream line_stream(line_string);
			  if (!(line_stream >> e.first))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << xspring_filename << std::endl);
			  }
			  else if ((e.first < min_idx) || (e.first >= max_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << xspring_filename << std::endl
					 << "  vertex index " << e.first << " is out of range"
					 << std::endl);
			  }

			  if (!(line_stream >> e.second))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << xspring_filename << std::endl);
			  }
			  else if ((e.second < min_idx) || (e.second >= max_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << xspring_filename << std::endl
					 << "  vertex index " << e.second << " is out of range"
					 << std::endl);
			  }

			  if (!(line_stream >> parameters[0]))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << xspring_filename << std::endl);
			  }
			  else if (parameters[0] < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << xspring_filename << std::endl
					 << "  spring constant is negative" << std::endl);
			  }

			  if (!(line_stream >> parameters[1]))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << xspring_filename << std::endl);
			  }
			  else if (parameters[1] < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << xspring_filename << std::endl
					 << "  spring resting length is negative" << std::endl);
			  }
			  parameters[1] *= d_length_scale_factor;

			  if (!(line_stream >> force_fcn_idx))
			  {
			      force_fcn_idx = 0; // default force function specification.
			  }

			  double param;
			  while (line_stream >> param)
			  {
			      parameters.push_back(param);
			  }
		      }

		      // Modify kappa and length according to whether uniform
		      // values are to be employed for this particular structure.
		      if (d_using_uniform_xspring_stiffness[ln][j])
		      {
			  parameters[0] = d_uniform_xspring_stiffness[ln][j];
		      }
		      if (d_using_uniform_xspring_rest_length[ln][j])
		      {
			  parameters[1] = d_uniform_xspring_rest_length[ln][j];
		      }
		      if (d_using_uniform_xspring_force_fcn_idx[ln][j])
		      {
			  force_fcn_idx = d_uniform_xspring_force_fcn_idx[ln][j];
		      }

		      // Check to see if the spring constant is zero and, if so,
		      // emit a warning.
		      if (!warned && d_enable_xsprings[ln][j] &&
			  (parameters[0] == 0.0 || MathUtilities<double>::equalEps(parameters[0], 0.0)))
		      {
			  TBOX_WARNING(d_object_name << ":\n  Crosslink spring with zero spring "
				       "constant encountered in ASCII input file "
				       "named " << xspring_filename << "." << std::endl);
			  warned = true;
		      }

		      // Correct the edge numbers to be in the global Lagrangian indexing
		      // scheme.
		      if (!input_uses_global_idxs)
		      {
			  e.first += d_vertex_offset[ln][j];
			  e.second += d_vertex_offset[ln][j];
		      }

		      // Initialize the map data corresponding to the present edge.
		      //
		      // Note that in the edge map, each edge is associated with only the
		      // first vertex.
		      if (e.first > e.second)
		      {
			  std::swap<int>(e.first, e.second);
		      }
		      d_xspring_edge_map[ln][j].insert(std::make_pair(e.first, e));
		      XSpringSpec spec_data;
		      spec_data.parameters = parameters;
		      spec_data.force_fcn_idx = force_fcn_idx;
		      d_xspring_spec_data[ln][j].insert(std::make_pair(e, spec_data));
		  }

		  // Close the input file.
		  file_stream.close();

		  plog << d_object_name << ":  "
		       << "read " << num_edges << " edges from ASCII input file named " << xspring_filename << std::endl
		       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	      }

	      // Free the next MPI process to start reading the current file.
	      if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);
	  
	      j++;
	  }
	  struct_begin_idx += d_structs_clones_num[itype];
      }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
  } // readXSpringFiles

  void CIBStandardInitializer::readBeamFiles(const std::string& extension, const bool input_uses_global_idxs)
  {
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      const size_t num_base_filename = d_base_filename[ln].size();
      d_beam_spec_data[ln].resize(num_base_filename);

      unsigned struct_begin_idx=0;

      unsigned j = 0;
      for(unsigned itype=0;itype< d_num_structs_types;itype++)
      {
	  for(int istruct=0;istruct < d_structs_clones_num[itype];istruct++)
	  {
	      bool warned = false;

	      // Determine min/max index ranges.
	      const int min_idx = 0;
	      const int max_idx =
		  (input_uses_global_idxs ? std::accumulate(d_num_vertex[ln].begin(), d_num_vertex[ln].end(), 0) :
		   d_num_vertex[ln][j]);

	      // Wait for the previous MPI process to finish reading the current file.
	      if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

	      std::string beam_filename = d_base_filename[ln][struct_begin_idx];
	      if (beam_filename.size()) beam_filename.erase(beam_filename.size() - 1); 
	      beam_filename +=extension;

	      std::ifstream file_stream;
	      file_stream.open(beam_filename.c_str(), std::ios::in);
	      if (file_stream.is_open())
	      {
		  plog << d_object_name << ":  "
		       << "processing beam data from ASCII input file named " << beam_filename << std::endl
		       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

		  // The first line in the file indicates the number of beams in
		  // the input file.
		  int num_beams = -1;
		  if (!std::getline(file_stream, line_string))
		  {
		      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
				 "before line 1 of file " << beam_filename << std::endl);
		  }
		  else
		  {
		      line_string = discard_comments(line_string);
		      std::istringstream line_stream(line_string);
		      if (!(line_stream >> num_beams))
		      {
			  TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
				     "encountered on line 1 of file " << beam_filename << std::endl);
		      }
		  }

		  if (num_beams <= 0)
		  {
		      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
				 << beam_filename << std::endl);
		  }

		  // Each successive line provides the connectivity and material
		  // parameter information for each beam in the structure.
		  for (int k = 0; k < num_beams; ++k)
		  {
		      int prev_idx = std::numeric_limits<int>::max(), curr_idx = std::numeric_limits<int>::max(),
			  next_idx = std::numeric_limits<int>::max();
		      double bend = 0.0;
		      Vector curv(Vector::Zero());
		      if (!std::getline(file_stream, line_string))
		      {
			  TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k + 2
				     << " of file " << beam_filename << std::endl);
		      }
		      else
		      {
			  line_string = discard_comments(line_string);
			  std::istringstream line_stream(line_string);
			  if (!(line_stream >> prev_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << beam_filename << std::endl);
			  }
			  else if ((prev_idx < min_idx) || (prev_idx >= max_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << beam_filename << std::endl << "  vertex index "
					 << prev_idx << " is out of range" << std::endl);
			  }

			  if (!(line_stream >> curr_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << beam_filename << std::endl);
			  }
			  else if ((curr_idx < min_idx) || (curr_idx >= max_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << beam_filename << std::endl << "  vertex index "
					 << curr_idx << " is out of range" << std::endl);
			  }

			  if (!(line_stream >> next_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << beam_filename << std::endl);
			  }
			  else if ((next_idx < min_idx) || (next_idx >= max_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << beam_filename << std::endl << "  vertex index "
					 << next_idx << " is out of range" << std::endl);
			  }

			  if (!(line_stream >> bend))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << beam_filename << std::endl);
			  }
			  else if (bend < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << beam_filename << std::endl
					 << "  beam constant is negative" << std::endl);
			  }

			  bool curv_found_in_input = false;
			  for (unsigned int d = 0; d < NDIM; ++d)
			  {
			      double c;
			      if (!(line_stream >> c))
			      {
				  if (curv_found_in_input)
				  {
				      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
						 "encountered on line " << k + 2 << " of file "
						 << beam_filename << std::endl
						 << "  incomplete beam curvature specification"
						 << std::endl);
				  }
			      }
			      else
			      {
				  curv_found_in_input = true;
				  curv[d] = c;
			      }
			  }
		      }

		      // Modify bend and curvature according to whether uniform
		      // values are to be employed for this particular structure.
		      if (d_using_uniform_beam_bend_rigidity[ln][j])
		      {
			  bend = d_uniform_beam_bend_rigidity[ln][j];
		      }
		      if (d_using_uniform_beam_curvature[ln][j])
		      {
			  curv = d_uniform_beam_curvature[ln][j];
		      }

		      // Check to see if the bending rigidity is zero and, if so,
		      // emit a warning.
		      if (!warned && d_enable_beams[ln][j] && (bend == 0.0 || MathUtilities<double>::equalEps(bend, 0.0)))
		      {
			  TBOX_WARNING(d_object_name << ":\n  Beam with zero bending rigidity "
				       "encountered in ASCII input file named " << beam_filename << "."
				       << std::endl);
			  warned = true;
		      }

		      // Correct the node numbers to be in the global Lagrangian
		      // indexing scheme.
		      if (!input_uses_global_idxs)
		      {
			  prev_idx += d_vertex_offset[ln][j];
			  curr_idx += d_vertex_offset[ln][j];
			  next_idx += d_vertex_offset[ln][j];
		      }

		      // Initialize the map data corresponding to the present
		      // beam.
		      //
		      // Note that in the beam property map, each edge is
		      // associated with only the "current" vertex.
		      BeamSpec spec_data;
		      spec_data.neighbor_idxs = std::make_pair(next_idx, prev_idx);
		      spec_data.bend_rigidity = bend;
		      spec_data.curvature = curv;
		      d_beam_spec_data[ln][j].insert(std::make_pair(curr_idx, spec_data));
		  }

		  // Close the input file.
		  file_stream.close();

		  plog << d_object_name << ":  "
		       << "read " << num_beams << " beams from ASCII input file named " << beam_filename << std::endl
		       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	      }

	      // Free the next MPI process to start reading the current file.
	      if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);

	      j++;
	  }
	  struct_begin_idx += d_structs_clones_num[itype];
      }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
  } // readBeamFiles

  void CIBStandardInitializer::readRodFiles(const std::string& extension, const bool input_uses_global_idxs)
  {
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      const size_t num_base_filename = d_base_filename[ln].size();
      d_rod_edge_map[ln].resize(num_base_filename);
      d_rod_spec_data[ln].resize(num_base_filename);

      unsigned struct_begin_idx=0;

      unsigned j = 0;
      for(unsigned itype=0;itype< d_num_structs_types;itype++)
      {
	  for(int istruct=0;istruct < d_structs_clones_num[itype];istruct++)
	  {
	      // Determine min/max index ranges.
	      const int min_idx = 0;
	      const int max_idx =
		  (input_uses_global_idxs ? std::accumulate(d_num_vertex[ln].begin(), d_num_vertex[ln].end(), 0) :
		   d_num_vertex[ln][j]);

	      // Wait for the previous MPI process to finish reading the current file.
	      if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

	      std::string rod_filename = d_base_filename[ln][struct_begin_idx];
	      if (rod_filename.size()) rod_filename.erase(rod_filename.size() - 1-std::max(0,(int)log10(1e-16+struct_begin_idx))); //erase number
	      rod_filename +=extension;

	      std::ifstream file_stream;
	      file_stream.open(rod_filename.c_str(), std::ios::in);
	      if (file_stream.is_open())
	      {
		  plog << d_object_name << ":  "
		       << "processing rod data from ASCII input file named " << rod_filename << std::endl
		       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

		  // The first line in the file indicates the number of rods in
		  // the input file.
		  int num_rods = -1;
		  if (!std::getline(file_stream, line_string))
		  {
		      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
				 "before line 1 of file " << rod_filename << std::endl);
		  }
		  else
		  {
		      line_string = discard_comments(line_string);
		      std::istringstream line_stream(line_string);
		      if (!(line_stream >> num_rods))
		      {
			  TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
				     "encountered on line 1 of file " << rod_filename << std::endl);
		      }
		  }

		  if (num_rods <= 0)
		  {
		      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
				 << rod_filename << std::endl);
		  }

		  // Each successive line provides the connectivity and material
		  // parameter information for each rod in the structure.
		  for (int k = 0; k < num_rods; ++k)
		  {
		      int curr_idx = std::numeric_limits<int>::max(), next_idx = std::numeric_limits<int>::max();
		      boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> properties;
		      double& ds = properties[0];
		      double& a1 = properties[1];
		      double& a2 = properties[2];
		      double& a3 = properties[3];
		      double& b1 = properties[4];
		      double& b2 = properties[5];
		      double& b3 = properties[6];
		      double& kappa1 = properties[7];
		      double& kappa2 = properties[8];
		      double& tau = properties[9];

		      if (!std::getline(file_stream, line_string))
		      {
			  TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k + 2
				     << " of file " << rod_filename << std::endl);
		      }
		      else
		      {
			  line_string = discard_comments(line_string);
			  std::istringstream line_stream(line_string);

			  if (!(line_stream >> curr_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl);
			  }
			  else if ((curr_idx < min_idx) || (curr_idx >= max_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl << "  vertex index "
					 << curr_idx << " is out of range" << std::endl);
			  }

			  if (!(line_stream >> next_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl);
			  }
			  else if ((next_idx < min_idx) || (next_idx >= max_idx))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl << "  vertex index "
					 << next_idx << " is out of range" << std::endl);
			  }

			  if (!(line_stream >> ds))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl);
			  }
			  else if (ds < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl
					 << "  rod material constant ds is negative" << std::endl);
			  }

			  if (!(line_stream >> a1))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl);
			  }
			  else if (a1 < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl
					 << "  rod material constant a1 is negative" << std::endl);
			  }

			  if (!(line_stream >> a2))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl);
			  }
			  else if (a2 < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl
					 << "  rod material constant a2 is negative" << std::endl);
			  }

			  if (!(line_stream >> a3))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl);
			  }
			  else if (a3 < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl
					 << "  rod material constant a3 is negative" << std::endl);
			  }

			  if (!(line_stream >> b1))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl);
			  }
			  else if (b1 < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl
					 << "  rod material constant b1 is negative" << std::endl);
			  }

			  if (!(line_stream >> b2))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl);
			  }
			  else if (b2 < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl
					 << "  rod material constant b2 is negative" << std::endl);
			  }

			  if (!(line_stream >> b3))
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl);
			  }
			  else if (b3 < 0.0)
			  {
			      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
					 << " of file " << rod_filename << std::endl
					 << "  rod material constant b3 is negative" << std::endl);
			  }

			  bool curvature_data_found_in_input = false;

			  if (!(line_stream >> kappa1))
			  {
			      kappa1 = 0.0;
			  }
			  else
			  {
			      curvature_data_found_in_input = true;
			  }

			  if (!(line_stream >> kappa2))
			  {
			      kappa2 = 0.0;
			      if (curvature_data_found_in_input)
			      {
				  TBOX_WARNING(d_object_name << ":\n  Potentially invalid entry in input file "
					       "encountered on line " << k + 2 << " of file "
					       << rod_filename << std::endl
					       << "  intrinsic curvature kappa1 was specified "
					       "but kappa2 was not" << std::endl);
			      }
			  }
			  else
			  {
			      curvature_data_found_in_input = true;
			  }

			  if (!(line_stream >> tau))
			  {
			      tau = 0.0;
			      if (curvature_data_found_in_input)
			      {
				  TBOX_WARNING(d_object_name << ":\n  Potentially invalid entry in input file "
					       "encountered on line " << k + 2 << " of file "
					       << rod_filename << std::endl
					       << "  intrinsic curvatures kappa1 and kappa2 "
					       "were specified but intrinsic twist tau was "
					       "not" << std::endl);
			      }
			  }
			  else
			  {
			      curvature_data_found_in_input = true;
			  }
		      }

		      // Modify properties according to whether uniform values are
		      // to be employed for this particular structure.
		      if (d_using_uniform_rod_properties[ln][j])
		      {
			  properties = d_uniform_rod_properties[ln][j];
		      }

		      // Correct the node numbers to be in the global Lagrangian
		      // indexing scheme.
		      if (!input_uses_global_idxs)
		      {
			  curr_idx += d_vertex_offset[ln][j];
			  next_idx += d_vertex_offset[ln][j];
		      }

		      // Initialize the map data corresponding to the present rod.
		      //
		      // Note that in the rod property map, each edge is
		      // associated with only the "current" vertex.
		      Edge e;
		      e.first = curr_idx;
		      e.second = next_idx;
		      d_rod_edge_map[ln][j].insert(std::make_pair(e.first, e));
		      RodSpec rod_spec;
		      rod_spec.properties = properties;
		      d_rod_spec_data[ln][j].insert(std::make_pair(e, rod_spec));
		  }

		  // Close the input file.
		  file_stream.close();

		  plog << d_object_name << ":  "
		       << "read " << num_rods << " rods from ASCII input file named " << rod_filename << std::endl
		       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	      }

	      // Free the next MPI process to start reading the current file.
	      if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);
	      
	      j++;
	  }
	  struct_begin_idx += d_structs_clones_num[itype];
      }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
  } // readRodFiles

  void CIBStandardInitializer::readTargetPointFiles(const std::string& extension)
  {
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      const size_t num_base_filename = d_base_filename[ln].size();
      d_target_spec_data[ln].resize(num_base_filename);
      for (unsigned int j = 0; j < num_base_filename; ++j)
      {
	bool warned = false;
	  
	// Determine min/max index ranges.
	 const int min_idx = 0;
	const int max_idx = d_num_vertex[ln][j];

	// Wait for the previous MPI process to finish reading the current file.
	if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

	TargetSpec default_spec;
	default_spec.stiffness = 0.0;
	default_spec.damping = 0.0;
	d_target_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);

	const std::string target_point_stiffness_filename = d_base_filename[ln][j] + extension;
	std::ifstream file_stream;
	file_stream.open(target_point_stiffness_filename.c_str(), std::ios::in);
	if (file_stream.is_open())
	{
	  plog << d_object_name << ":  "
	       << "processing target point data from ASCII input file named " << target_point_stiffness_filename
	       << std::endl << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

	  // The first line in the file indicates the number of target
	  // point specifications in the input file.
	  int num_target_points = -1;
	  if (!std::getline(file_stream, line_string))
	  {
	    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
		       "before line 1 of file " << target_point_stiffness_filename
		       << std::endl);
	  }
	  else
	  {
	    line_string = discard_comments(line_string);
	    std::istringstream line_stream(line_string);
	    if (!(line_stream >> num_target_points))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
			 "encountered on line 1 of file " << target_point_stiffness_filename
			 << std::endl);
	    }
	  }

	  if (num_target_points <= 0)
	  {
	    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
		       << target_point_stiffness_filename << std::endl);
	  }

	  // Each successive line indicates the vertex number and spring
	  // constant associated with any target points.
	  for (int k = 0; k < num_target_points; ++k)
	  {
	    int n = std::numeric_limits<int>::max();
	    if (!std::getline(file_stream, line_string))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k + 2
			 << " of file " << target_point_stiffness_filename << std::endl);
	    }
	    else
	    {
	      line_string = discard_comments(line_string);
	      std::istringstream line_stream(line_string);
	      if (!(line_stream >> n))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << target_point_stiffness_filename << std::endl);
	      }
	      else if ((n < min_idx) || (n >= max_idx))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << target_point_stiffness_filename << std::endl
			   << "  vertex index " << n << " is out of range" << std::endl);
	      }

	      if (!(line_stream >> d_target_spec_data[ln][j][n].stiffness))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << target_point_stiffness_filename << std::endl);
	      }
	      else if (d_target_spec_data[ln][j][n].stiffness < 0.0)
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << target_point_stiffness_filename << std::endl
			   << "  target point spring constant is negative" << std::endl);
	      }

	      if (!(line_stream >> d_target_spec_data[ln][j][n].damping))
	      {
		d_target_spec_data[ln][j][n].damping = 0.0;
	      }
	      else if (d_target_spec_data[ln][j][n].damping < 0.0)
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << target_point_stiffness_filename << std::endl
			   << "  target point damping coefficient is negative" << std::endl);
	      }
	    }

	    // Check to see if the penalty spring constant is zero and,
	    // if so, emit a warning.
	    const double kappa = d_target_spec_data[ln][j][n].stiffness;
	    if (!warned && d_enable_target_points[ln][j] &&
		(kappa == 0.0 || MathUtilities<double>::equalEps(kappa, 0.0)))
	    {
	      TBOX_WARNING(d_object_name << ":\n  Target point with zero penalty spring "
			   "constant encountered in ASCII input file "
			   "named " << target_point_stiffness_filename << "." << std::endl);
	      warned = true;
	    }
	  }

	  // Close the input file.
	  file_stream.close();

	  plog << d_object_name << ":  "
	       << "read " << num_target_points << " target points from ASCII input file named "
	       << target_point_stiffness_filename << std::endl << "  on MPI process " << SAMRAI_MPI::getRank()
	       << std::endl;
	}

	// Modify the target point stiffness constants according to whether
	// target point penalty forces are enabled, or whether uniform
	// values are to be employed, for this particular structure.
	if (!d_enable_target_points[ln][j])
	{
	  for (int k = 0; k < d_num_vertex[ln][j]; ++k)
	  {
	    d_target_spec_data[ln][j][k].stiffness = 0.0;
	    d_target_spec_data[ln][j][k].damping = 0.0;
	  }
	}
	else
	{
	  if (d_using_uniform_target_stiffness[ln][j])
	  {
	    for (int k = 0; k < d_num_vertex[ln][j]; ++k)
	    {
	      d_target_spec_data[ln][j][k].stiffness = d_uniform_target_stiffness[ln][j];
	    }
	  }
	  if (d_using_uniform_target_damping[ln][j])
	  {
	    for (int k = 0; k < d_num_vertex[ln][j]; ++k)
	    {
	      d_target_spec_data[ln][j][k].damping = d_uniform_target_damping[ln][j];
	    }
	  }
	}

	// Free the next MPI process to start reading the current file.
	if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);
      }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
  } // readTargetPointFiles

  void CIBStandardInitializer::readAnchorPointFiles(const std::string& extension)
  {
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      const size_t num_base_filename = d_base_filename[ln].size();
      d_anchor_spec_data[ln].resize(num_base_filename);
      for (unsigned int j = 0; j < num_base_filename; ++j)
      {
	// Determine min/max index ranges.
	const int min_idx = 0;
	const int max_idx = d_num_vertex[ln][j];

	// Wait for the previous MPI process to finish reading the current file.
	if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

	AnchorSpec default_spec;
	default_spec.is_anchor_point = false;
	d_anchor_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);

	const std::string anchor_point_filename = d_base_filename[ln][j] + extension;
	std::ifstream file_stream;
	file_stream.open(anchor_point_filename.c_str(), std::ios::in);
	if (file_stream.is_open())
	{
	  plog << d_object_name << ":  "
	       << "processing anchor point data from ASCII input file named " << anchor_point_filename
	       << std::endl << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

	  // The first line in the file indicates the number of anchor
	  // points in the input file.
	  int num_anchor_pts = -1;
	  if (!std::getline(file_stream, line_string))
	  {
	    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
		       "before line 1 of file " << anchor_point_filename << std::endl);
	  }
	  else
	  {
	    line_string = discard_comments(line_string);
	    std::istringstream line_stream(line_string);
	    if (!(line_stream >> num_anchor_pts))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
			 "encountered on line 1 of file " << anchor_point_filename
			 << std::endl);
	    }
	  }

	  if (num_anchor_pts <= 0)
	  {
	    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
		       << anchor_point_filename << std::endl);
	  }

	  // Each successive line indicates the vertex number of the
	  // anchor points in the input file.
	  for (int k = 0; k < num_anchor_pts; ++k)
	  {
	    int n;
	    if (!std::getline(file_stream, line_string))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k + 2
			 << " of file " << anchor_point_filename << std::endl);
	    }
	    else
	    {
	      line_string = discard_comments(line_string);
	      std::istringstream line_stream(line_string);
	      if (!(line_stream >> n))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << anchor_point_filename << std::endl);
	      }
	      else if ((n < min_idx) || (n >= max_idx))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << anchor_point_filename << std::endl
			   << "  vertex index " << n << " is out of range" << std::endl);
	      }

	      d_anchor_spec_data[ln][j][n].is_anchor_point = true;
	    }
	  }

	  // Close the input file.
	  file_stream.close();

	  plog << d_object_name << ":  "
	       << "read " << num_anchor_pts << " anchor points from ASCII input file named "
	       << anchor_point_filename << std::endl << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	}

	// Free the next MPI process to start reading the current file.
	if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);
      }
    }
    return;
  } // readAnchorPointFiles

  void CIBStandardInitializer::readBoundaryMassFiles(const std::string& extension)
  {
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      const size_t num_base_filename = d_base_filename[ln].size();
      d_bdry_mass_spec_data[ln].resize(num_base_filename);
      for (unsigned int j = 0; j < num_base_filename; ++j)
      {
	// Determine min/max index ranges.
	const int min_idx = 0;
	const int max_idx = d_num_vertex[ln][j];

	// Wait for the previous MPI process to finish reading the current file.
	if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

	BdryMassSpec default_spec;
	default_spec.bdry_mass = 0.0;
	default_spec.stiffness = 0.0;
	d_bdry_mass_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);

	const std::string bdry_mass_filename = d_base_filename[ln][j] + extension;
	std::ifstream file_stream;
	file_stream.open(bdry_mass_filename.c_str(), std::ios::in);
	if (file_stream.is_open())
	{
	  plog << d_object_name << ":  "
	       << "processing boundary mass data from ASCII input file named " << bdry_mass_filename << std::endl
	       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

	  // The first line in the file indicates the number of massive IB
	  // points in the input file.
	  int num_bdry_mass_pts = -1;
	  if (!std::getline(file_stream, line_string))
	  {
	    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
		       "before line 1 of file " << bdry_mass_filename << std::endl);
	  }
	  else
	  {
	    line_string = discard_comments(line_string);
	    std::istringstream line_stream(line_string);
	    if (!(line_stream >> num_bdry_mass_pts))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
			 "encountered on line 1 of file " << bdry_mass_filename
			 << std::endl);
	    }
	  }

	  if (num_bdry_mass_pts <= 0)
	  {
	    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
		       << bdry_mass_filename << std::endl);
	  }

	  // Each successive line indicates the vertex number, mass, and
	  // penalty spring constant associated with any massive IB
	  // points.
	  for (int k = 0; k < num_bdry_mass_pts; ++k)
	  {
	    int n;
	    if (!std::getline(file_stream, line_string))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k + 2
			 << " of file " << bdry_mass_filename << std::endl);
	    }
	    else
	    {
	      line_string = discard_comments(line_string);
	      std::istringstream line_stream(line_string);
	      if (!(line_stream >> n))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << bdry_mass_filename << std::endl);
	      }
	      else if ((n < min_idx) || (n >= max_idx))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << bdry_mass_filename << std::endl
			   << "  vertex index " << n << " is out of range" << std::endl);
	      }

	      if (!(line_stream >> d_bdry_mass_spec_data[ln][j][n].bdry_mass))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << bdry_mass_filename << std::endl);
	      }
	      else if (d_bdry_mass_spec_data[ln][j][n].bdry_mass < 0.0)
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << bdry_mass_filename << std::endl
			   << "  boundary mass is negative" << std::endl);
	      }

	      if (!(line_stream >> d_bdry_mass_spec_data[ln][j][n].stiffness))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << bdry_mass_filename << std::endl);
	      }
	      else if (d_bdry_mass_spec_data[ln][j][n].stiffness < 0.0)
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k + 2
			   << " of file " << bdry_mass_filename << std::endl
			   << "  boundary mass spring constant is negative" << std::endl);
	      }
	    }
	  }

	  // Close the input file.
	  file_stream.close();

	  plog << d_object_name << ":  "
	       << "read " << num_bdry_mass_pts << " boundary mass points from ASCII input file named "
	       << bdry_mass_filename << std::endl << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	}

	// Modify the boundary masses and boundary mass stiffness constants
	// according to whether boundary mass is enabled, or whether uniform
	// values are to be employed, for this particular structure.
	if (!d_enable_bdry_mass[ln][j])
	{
	  for (int k = 0; k < d_num_vertex[ln][j]; ++k)
	  {
	    d_bdry_mass_spec_data[ln][j][k].bdry_mass = 0.0;
	    d_bdry_mass_spec_data[ln][j][k].stiffness = 0.0;
	  }
	}
	else
	{
	  if (d_using_uniform_bdry_mass[ln][j])
	  {
	    for (int k = 0; k < d_num_vertex[ln][j]; ++k)
	    {
	      d_bdry_mass_spec_data[ln][j][k].bdry_mass = d_uniform_bdry_mass[ln][j];
	    }
	  }
	  if (d_using_uniform_bdry_mass_stiffness[ln][j])
	  {
	    for (int k = 0; k < d_num_vertex[ln][j]; ++k)
	    {
	      d_bdry_mass_spec_data[ln][j][k].stiffness = d_uniform_bdry_mass_stiffness[ln][j];
	    }
	  }
	}

	// Free the next MPI process to start reading the current file.
	if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);
      }
    }
    return;
  } // readBoundaryMassFiles

  void CIBStandardInitializer::readDirectorFiles(const std::string& extension)
  {
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      const size_t num_base_filename = d_base_filename[ln].size();
      d_directors[ln].resize(num_base_filename);
      for (unsigned int j = 0; j < num_base_filename; ++j)
      {
	// Wait for the previous MPI process to finish reading the current file.
	if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

	d_directors[ln][j].resize(d_num_vertex[ln][j], std::vector<double>(3 * 3, 0.0));

	const std::string directors_filename = d_base_filename[ln][j] + extension;
	std::ifstream file_stream;
	file_stream.open(directors_filename.c_str(), std::ios::in);
	if (file_stream.is_open())
	{
	  plog << d_object_name << ":  "
	       << "processing director data from ASCII input file named " << directors_filename << std::endl
	       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

	  // The first line in the file indicates the number of sets of
	  // directors in the input file.
	  int num_directors_pts = -1;
	  if (!std::getline(file_stream, line_string))
	  {
	    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
		       "before line 1 of file " << directors_filename << std::endl);
	  }
	  else
	  {
	    line_string = discard_comments(line_string);
	    std::istringstream line_stream(line_string);
	    if (!(line_stream >> num_directors_pts))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
			 "encountered on line 1 of file " << directors_filename
			 << std::endl);
	    }
	  }

	  if (num_directors_pts != d_num_vertex[ln][j])
	  {
	    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
		       << directors_filename << std::endl);
	  }

	  // Each successive set of three lines indicates the initial
	  // configuration of a triad.
	  for (int k = 0; k < num_directors_pts; ++k)
	  {
	    for (int n = 0; n < 3; ++n)
	    {
	      if (!std::getline(file_stream, line_string))
	      {
		TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line "
			   << 3 * k + n + 2 << " of file " << directors_filename
			   << std::endl);
	      }
	      else
	      {
		line_string = discard_comments(line_string);
		std::istringstream line_stream(line_string);
		double D_norm_squared = 0.0;
		for (int d = 0; d < 3; ++d)
		{
		  if (!(line_stream >> d_directors[ln][j][k][3 * n + d]))
		  {
		    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input "
			       "file encountered on line " << 3 * k + n + 2
			       << " of file " << directors_filename << std::endl);
		  }
		  D_norm_squared += d_directors[ln][j][k][3 * n + d] * d_directors[ln][j][k][3 * n + d];
		}
		const double D_norm = sqrt(D_norm_squared);
		if (!MathUtilities<double>::equalEps(D_norm, 1.0))
		{
		  TBOX_WARNING(d_object_name << ":\n  Director vector on line " << 3 * k + n + 2
			       << " of file " << directors_filename
			       << " is not normalized; norm = " << D_norm << std::endl);
		  for (int d = 0; d < 3; ++d)
		  {
		    d_directors[ln][j][k][3 * n + d] /= D_norm;
		  }
		}
	      }
	    }
	  }

	  // Close the input file.
	  file_stream.close();

	  plog << d_object_name << ":  "
	       << "read " << num_directors_pts << " director triads from ASCII input file named "
	       << directors_filename << std::endl << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	}

	// Free the next MPI process to start reading the current file.
	if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);
      }
    }
    return;
  } // readDirectorFiles

  void CIBStandardInitializer::readInstrumentationFiles(const std::string& extension)
  {
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    int instrument_offset = 0;
    std::vector<std::string> instrument_names;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      const size_t num_base_filename = d_base_filename[ln].size();
      d_instrument_idx[ln].resize(num_base_filename);
      for (unsigned int j = 0; j < num_base_filename; ++j)
      {
	// Determine min/max index ranges.
	const int min_idx = 0;
	const int max_idx = d_num_vertex[ln][j];

	// Wait for the previous MPI process to finish reading the current file.
	if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

	const std::string inst_filename = d_base_filename[ln][j] + extension;
	std::ifstream file_stream;
	file_stream.open(inst_filename.c_str(), std::ios::in);
	if (file_stream.is_open() && d_enable_instrumentation[ln][j])
	{
	  plog << d_object_name << ":  "
	       << "processing instrumentation data from ASCII input file named " << inst_filename << std::endl
	       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

	  // The first line in the file indicates the number of
	  // instruments in the input file.
	  int num_inst = -1;
	  if (!std::getline(file_stream, line_string))
	  {
	    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
		       "before line 1 of file " << inst_filename << std::endl);
	  }
	  else
	  {
	    line_string = discard_comments(line_string);
	    std::istringstream line_stream(line_string);
	    if (!(line_stream >> num_inst))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
			 "encountered on line 1 of file " << inst_filename << std::endl);
	    }
	  }

	  if (num_inst <= 0)
	  {
	    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
		       << inst_filename << std::endl);
	  }

	  // The next several lines in the file indicate the names of the
	  // instruments in the input file.
	  for (int m = 0; m < num_inst; ++m)
	  {
	    if (!std::getline(file_stream, line_string))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << m + 2
			 << " of file " << inst_filename << std::endl);
	    }
	    else
	    {
	      line_string = discard_comments(line_string);

	      // trim leading whitespace
	      std::string::size_type notwhite = line_string.find_first_not_of(" \t\n");
	      line_string.erase(0, notwhite);

	      // trim trailing whitespace
	      notwhite = line_string.find_last_not_of(" \t\n");
	      line_string.erase(notwhite + 1);

	      instrument_names.push_back(line_string);
	    }
	  }

	  // The next line in the file indicates the number of
	  // instrumented IB points in the input file.
	  int num_inst_pts = -1;
	  if (!std::getline(file_stream, line_string))
	  {
	    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line "
		       << num_inst + 2 << " of file " << inst_filename << std::endl);
	  }
	  else
	  {
	    line_string = discard_comments(line_string);
	    std::istringstream line_stream(line_string);
	    if (!(line_stream >> num_inst_pts))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			 << num_inst + 2 << " of file " << inst_filename << std::endl);
	    }
	  }

	  if (num_inst_pts <= 0)
	  {
	    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << num_inst + 2
		       << " of file " << inst_filename << std::endl);
	  }

	  // Each successive line indicates the vertex number, meter
	  // number, and meter node indices of each of the instrumented IB
	  // points in the input file.
	  std::vector<bool> encountered_instrument_idx;
	  std::map<size_t, std::vector<bool> > encountered_node_idx;
	  for (int k = 0; k < num_inst_pts; ++k)
	  {
	    int n;
	    if (!std::getline(file_stream, line_string))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line "
			 << num_inst + k + 3 << " of file " << inst_filename << std::endl);
	    }
	    else
	    {
	      line_string = discard_comments(line_string);
	      std::istringstream line_stream(line_string);
	      if (!(line_stream >> n))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << num_inst + k + 3 << " of file " << inst_filename << std::endl);
	      }
	      else if ((n < min_idx) || (n >= max_idx))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << num_inst + k + 3 << " of file " << inst_filename << std::endl
			   << "  vertex index " << n << " is out of range" << std::endl);
	      }

	      std::pair<int, int>& idx = d_instrument_idx[ln][j][n];

	      if (!(line_stream >> idx.first))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << num_inst + k + 3 << " of file " << inst_filename << std::endl);
	      }
	      else if (idx.first < 0 || idx.first >= num_inst)
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << num_inst + k + 3 << " of file " << inst_filename << std::endl
			   << "  meter index " << idx.first << " is out of range"
			   << std::endl);
	      }

	      if (idx.first >= static_cast<int>(encountered_instrument_idx.size()))
	      {
		encountered_instrument_idx.resize(idx.first + 1, false);
	      }
	      encountered_instrument_idx[idx.first] = true;

	      if (!(line_stream >> idx.second))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << num_inst + k + 3 << " of file " << inst_filename << std::endl);
	      }
	      else if (idx.second < 0)
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << num_inst + k + 3 << " of file " << inst_filename << std::endl
			   << "  meter node index is negative" << std::endl);
	      }

	      if (idx.second >= static_cast<int>(encountered_node_idx[idx.first].size()))
	      {
		encountered_node_idx[idx.first].resize(idx.second + 1, false);
	      }
	      encountered_node_idx[idx.first][idx.second] = true;

	      // Correct the instrument index to account for
	      // instrument indices from earlier files.
	      idx.first += instrument_offset;
	    }
	  }

	  // Ensure that a complete range of instrument indices were found
	  // in the input file.
	  for (std::vector<bool>::iterator meter_it = encountered_instrument_idx.begin();
	       meter_it != encountered_instrument_idx.end();
	       ++meter_it)
	  {
	    const size_t meter_idx = std::distance(encountered_instrument_idx.begin(), meter_it);
	    if ((*meter_it) == false)
	    {
	      TBOX_ERROR(d_object_name << ":\n  "
			 << "  Instrument index " << meter_idx << " not found in input file "
			 << inst_filename << std::endl);
	    }

	    std::vector<bool>& meter_node_idxs = encountered_node_idx[meter_idx];
	    for (std::vector<bool>::iterator node_it = meter_node_idxs.begin();
		 node_it != meter_node_idxs.end();
		 ++node_it)
	    {
	      const size_t node_idx = std::distance(meter_node_idxs.begin(), node_it);
	      if ((*node_it) == false)
	      {
		TBOX_ERROR(d_object_name << ":\n  "
			   << "  Node index " << node_idx << " associated with meter index "
			   << meter_idx << " not found in input file " << inst_filename
			   << std::endl);
	      }
	    }
	  }

	  if (static_cast<int>(encountered_instrument_idx.size()) != num_inst)
	  {
	    TBOX_ERROR(d_object_name << ":\n  "
		       << "  Not all anticipated instrument indices were found in input file "
		       << inst_filename << "  Expected to find " << num_inst
		       << " distinct meter indices in input file" << std::endl);
	  }

	  // Increment the meter offset.
	  instrument_offset += encountered_instrument_idx.size();

	  // Close the input file.
	  file_stream.close();

	  plog << d_object_name << ":  "
	       << "read " << num_inst_pts << " instrumentation points from ASCII input file named "
	       << inst_filename << std::endl << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	}

	// Free the next MPI process to start reading the current file.
	if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);
      }
    }
    IBInstrumentationSpec::setInstrumentNames(instrument_names);
    return;
  } // readInstrumentationFiles

  void CIBStandardInitializer::readSourceFiles(const std::string& extension)
  {
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
      int source_offset = 0;
      std::vector<std::string> source_names;
      std::vector<double> source_radii;
      const size_t num_base_filename = d_base_filename[ln].size();
      d_source_idx[ln].resize(num_base_filename);
      for (unsigned int j = 0; j < num_base_filename; ++j)
      {
	// Determine min/max index ranges.
	const int min_idx = 0;
	const int max_idx = d_num_vertex[ln][j];

	// Wait for the previous MPI process to finish reading the current file.
	if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

	const std::string source_filename = d_base_filename[ln][j] + extension;
	std::ifstream file_stream;
	file_stream.open(source_filename.c_str(), std::ios::in);
	if (file_stream.is_open() && d_enable_sources[ln][j])
	{
	  plog << d_object_name << ":  "
	       << "processing source data from ASCII input file named " << source_filename << std::endl
	       << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

	  // The first line in the file indicates the number of sources in
	  // the input file.
	  int num_source = -1;
	  if (!std::getline(file_stream, line_string))
	  {
	    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
		       "before line 1 of file " << source_filename << std::endl);
	  }
	  else
	  {
	    line_string = discard_comments(line_string);
	    std::istringstream line_stream(line_string);
	    if (!(line_stream >> num_source))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
			 "encountered on line 1 of file " << source_filename << std::endl);
	    }
	  }

	  if (num_source <= 0)
	  {
	    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
		       << source_filename << std::endl);
	  }

	  // The next several lines in the file indicate the names of the
	  // sources in the input file.
	  for (int m = 0; m < num_source; ++m)
	  {
	    if (!std::getline(file_stream, line_string))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << m + 2
			 << " of file " << source_filename << std::endl);
	    }
	    else
	    {
	      line_string = discard_comments(line_string);

	      // trim leading whitespace
	      std::string::size_type notwhite = line_string.find_first_not_of(" \t\n");
	      line_string.erase(0, notwhite);

	      // trim trailing whitespace
	      notwhite = line_string.find_last_not_of(" \t\n");
	      line_string.erase(notwhite + 1);

	      source_names.push_back(line_string);
	    }
	  }

	  // The next several lines in the file indicate the sizes of the
	  // sources in the input file.
	  for (int m = 0; m < num_source; ++m)
	  {
	    if (!std::getline(file_stream, line_string))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << m + 2
			 << " of file " << source_filename << std::endl);
	    }
	    else
	    {
	      line_string = discard_comments(line_string);
	      std::istringstream line_stream(line_string);
	      double r;
	      if (!(line_stream >> r) || r <= 0.0)
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << num_source + m + 2 << " of file " << source_filename
			   << std::endl);
	      }
	      source_radii.push_back(r);
	    }
	  }

	  // The next line in the file indicates the number of source
	  // points in the input file.
	  int num_source_pts = -1;
	  if (!std::getline(file_stream, line_string))
	  {
	    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line "
		       << 2 * num_source + 2 << " of file " << source_filename << std::endl);
	  }
	  else
	  {
	    line_string = discard_comments(line_string);
	    std::istringstream line_stream(line_string);
	    if (!(line_stream >> num_source_pts) || (num_source_pts <= 0))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			 << 2 * num_source + 2 << " of file " << source_filename << std::endl);
	    }
	  }

	  // Each successive line indicates the vertex number and source
	  // number of the IB points in the input file.
	  for (int k = 0; k < num_source_pts; ++k)
	  {
	    int n;
	    if (!std::getline(file_stream, line_string))
	    {
	      TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line "
			 << 2 * num_source + k + 3 << " of file " << source_filename
			 << std::endl);
	    }
	    else
	    {
	      line_string = discard_comments(line_string);
	      std::istringstream line_stream(line_string);
	      if (!(line_stream >> n))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << 2 * num_source + k + 3 << " of file " << source_filename
			   << std::endl);
	      }
	      else if ((n < min_idx) || (n >= max_idx))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << 2 * num_source + k + 3 << " of file " << source_filename
			   << std::endl << "  vertex index " << n << " is out of range"
			   << std::endl);
	      }

	      int& source_idx = d_source_idx[ln][j][n];

	      if (!(line_stream >> source_idx))
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << 2 * num_source + k + 3 << " of file " << source_filename
			   << std::endl);
	      }
	      else if (source_idx < 0 || source_idx >= num_source)
	      {
		TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
			   << 2 * num_source + k + 3 << " of file " << source_filename
			   << std::endl << "  meter index " << source_idx
			   << " is out of range" << std::endl);
	      }

	      // Correct the source index to account for source
	      // indices from earlier files.
	      source_idx += source_offset;
	    }
	  }

	  // Increment the meter offset.
	  source_offset += num_source;

	  // Close the input file.
	  file_stream.close();

	  plog << d_object_name << ":  "
	       << "read " << num_source_pts << " source points from ASCII input file named " << source_filename
	       << std::endl << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
	}

	// Free the next MPI process to start reading the current file.
	if (d_use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);
      }
      IBStandardSourceGen::setNumSources(ln, source_offset);
      IBStandardSourceGen::setSourceNames(ln, source_names);
      IBStandardSourceGen::setSourceRadii(ln, source_radii);
    }
    return;
  } // readSourceFiles

  void CIBStandardInitializer::getPatchVertices(std::vector<std::pair<int, int> >& patch_vertices,
					       const Pointer<Patch<NDIM> > patch,
					       const int level_number,
					       const bool /*can_be_refined*/,
					       const double* const domain_x_lower,
					       const double* const domain_x_upper,
					       const IntVector<NDIM>& periodic_shift) const
  {
    // Loop over all of the vertices to determine the indices of those vertices
    // within the present patch.
    //
    // NOTE: This is clearly not the best way to do this, but it will work for
    // now.
    const Box<NDIM>& patch_box = patch->getBox();
    const CellIndex<NDIM>& patch_lower = patch_box.lower();
    const CellIndex<NDIM>& patch_upper = patch_box.upper();
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const patch_x_lower = patch_geom->getXLower();
    const double* const patch_x_upper = patch_geom->getXUpper();
    const double* const patch_dx = patch_geom->getDx();
    for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
    {
      for (int k = 0; k < d_num_vertex[level_number][j]; ++k)
      {
	std::pair<int, int> point_index(j, k);
	const Point& X =
	  getShiftedVertexPosn(point_index, level_number, domain_x_lower, domain_x_upper, periodic_shift);
	const CellIndex<NDIM> idx =
	  IndexUtilities::getCellIndex(X, patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper);
	if (patch_box.contains(idx)) patch_vertices.push_back(point_index);
      }
    }
    return;
  } // getPatchVertices

  int CIBStandardInitializer::getCanonicalLagrangianIndex(const std::pair<int, int>& point_index,
							 const int level_number) const
  {
    return d_vertex_offset[level_number][point_index.first] + point_index.second;
  } // getCanonicalLagrangianIndex

  Point CIBStandardInitializer::getVertexPosn(const std::pair<int, int>& point_index, const int level_number) const
  {
      int prototypeID = d_struct_prototype_map[level_number][point_index.first];

      Point X;
      Eigen::Vector3d X_com = d_X_com[level_number][point_index.first];
      Eigen::Matrix3d rot_mat = (d_quaternion[level_number][point_index.first]).toRotationMatrix();
	      
      const Point X0 = d_vertex_posn[level_number][prototypeID][point_index.second];
      Eigen::Vector3d x0;
      for (unsigned int d = 0; d < NDIM; ++d) x0[d] = X0[d]; 

      Eigen::Vector3d x = rot_mat*x0;
      for (unsigned int d = 0; d < NDIM; ++d) X[d] = X_com[d] + x[d];
      
      return X;      
  } // getVertexPosn

  Point CIBStandardInitializer::getShiftedVertexPosn(const std::pair<int, int>& point_index,
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

  const CIBStandardInitializer::TargetSpec&
  CIBStandardInitializer::getVertexTargetSpec(const std::pair<int, int>& point_index, const int level_number) const
  {
    return d_target_spec_data[level_number][point_index.first][point_index.second];
  } // getVertexTargetSpec

  const CIBStandardInitializer::AnchorSpec&
  CIBStandardInitializer::getVertexAnchorSpec(const std::pair<int, int>& point_index, const int level_number) const
  {
    return d_anchor_spec_data[level_number][point_index.first][point_index.second];
  } // getVertexAnchorSpec

  const CIBStandardInitializer::BdryMassSpec&
  CIBStandardInitializer::getVertexBdryMassSpec(const std::pair<int, int>& point_index, const int level_number) const
  {
    return d_bdry_mass_spec_data[level_number][point_index.first][point_index.second];
  } // getVertexBdryMassSpec

  const std::vector<double>& CIBStandardInitializer::getVertexDirectors(const std::pair<int, int>& point_index,
								       const int level_number) const
  {
    return d_directors[level_number][point_index.first][point_index.second];
  } // getVertexDirectors

  std::pair<int, int> CIBStandardInitializer::getVertexInstrumentationIndices(const std::pair<int, int>& point_index,
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

  int CIBStandardInitializer::getVertexSourceIndices(const std::pair<int, int>& point_index, const int level_number) const
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

  std::vector<Pointer<Streamable> > CIBStandardInitializer::initializeNodeData(const std::pair<int, int>& point_index,
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
      if (d_enable_springs[level_number][j])
      {
	for (std::multimap<int, Edge>::const_iterator it =
	       d_spring_edge_map[level_number][j].lower_bound(mastr_idx);
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
      }
      const size_t num_base_filename = d_base_filename[level_number].size();
      for (unsigned int j = 0; j < num_base_filename; ++j)
      {
	if (!d_enable_xsprings[level_number][j]) continue;
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
    if (d_enable_beams[level_number][j])
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
    if (d_enable_rods[level_number][j])
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
    if (d_enable_target_points[level_number][j])
    {
      const TargetSpec& spec_data = getVertexTargetSpec(point_index, level_number);
      const double kappa_target = spec_data.stiffness;
      const double eta_target = spec_data.damping;
      const Point& X_target = getVertexPosn(point_index, level_number);
      node_data.push_back(new IBTargetPointForceSpec(mastr_idx, kappa_target, eta_target, X_target));
    }

    // Initialize any anchor point specifications associated with the present
    // vertex.
    if (d_enable_anchor_points[level_number][j])
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
    if (d_enable_instrumentation[level_number][j])
    {
      const std::pair<int, int> inst_idx = getVertexInstrumentationIndices(point_index, level_number);
      if (inst_idx.first != -1 && inst_idx.second != -1)
      {
	node_data.push_back(new IBInstrumentationSpec(mastr_idx, inst_idx.first, inst_idx.second));
      }
    }

    // Initialize any source specifications associated with the present
    // vertex.
    if (d_enable_sources[level_number][j])
    {
      const int source_idx = getVertexSourceIndices(point_index, level_number);
      if (source_idx != -1)
      {
	node_data.push_back(new IBSourceSpec(mastr_idx, source_idx));
      }
    }
    return node_data;
  } // initializeNodeData

void CIBStandardInitializer::getFromInput(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif

    // Determine whether to use "batons" to prevent multiple MPI processes from
    // reading the same file at once.
    if (db->keyExists("use_file_batons")) d_use_file_batons = db->getBool("use_file_batons");

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
    d_level_is_initialized.resize(d_max_levels, false);

    d_base_filename.resize(d_max_levels);

    d_num_vertex.resize(d_max_levels);
    d_vertex_offset.resize(d_max_levels);
    d_vertex_posn.resize(d_max_levels);
    d_quaternion.resize(d_max_levels);
    d_X_com.resize(d_max_levels);
    d_struct_prototype_map.resize(d_max_levels);
    d_enable_springs.resize(d_max_levels);
    d_spring_edge_map.resize(d_max_levels);
    d_spring_spec_data.resize(d_max_levels);
    d_using_uniform_spring_stiffness.resize(d_max_levels);
    d_uniform_spring_stiffness.resize(d_max_levels);
    d_using_uniform_spring_rest_length.resize(d_max_levels);
    d_uniform_spring_rest_length.resize(d_max_levels);
    d_using_uniform_spring_force_fcn_idx.resize(d_max_levels);
    d_uniform_spring_force_fcn_idx.resize(d_max_levels);

    d_enable_xsprings.resize(d_max_levels);
    d_xspring_edge_map.resize(d_max_levels);
    d_xspring_spec_data.resize(d_max_levels);
    d_using_uniform_xspring_stiffness.resize(d_max_levels);
    d_uniform_xspring_stiffness.resize(d_max_levels);
    d_using_uniform_xspring_rest_length.resize(d_max_levels);
    d_uniform_xspring_rest_length.resize(d_max_levels);
    d_using_uniform_xspring_force_fcn_idx.resize(d_max_levels);
    d_uniform_xspring_force_fcn_idx.resize(d_max_levels);

    d_enable_beams.resize(d_max_levels);
    d_beam_spec_data.resize(d_max_levels);
    d_using_uniform_beam_bend_rigidity.resize(d_max_levels);
    d_uniform_beam_bend_rigidity.resize(d_max_levels);
    d_using_uniform_beam_curvature.resize(d_max_levels);
    d_uniform_beam_curvature.resize(d_max_levels);

    d_enable_rods.resize(d_max_levels);
    d_rod_edge_map.resize(d_max_levels);
    d_rod_spec_data.resize(d_max_levels);
    d_using_uniform_rod_properties.resize(d_max_levels);
    d_uniform_rod_properties.resize(d_max_levels);

    d_enable_target_points.resize(d_max_levels);
    d_target_spec_data.resize(d_max_levels);
    d_using_uniform_target_stiffness.resize(d_max_levels);
    d_uniform_target_stiffness.resize(d_max_levels);
    d_using_uniform_target_damping.resize(d_max_levels);
    d_uniform_target_damping.resize(d_max_levels);

    d_enable_anchor_points.resize(d_max_levels);
    d_anchor_spec_data.resize(d_max_levels);

    d_enable_bdry_mass.resize(d_max_levels);
    d_bdry_mass_spec_data.resize(d_max_levels);
    d_using_uniform_bdry_mass.resize(d_max_levels);
    d_uniform_bdry_mass.resize(d_max_levels);
    d_using_uniform_bdry_mass_stiffness.resize(d_max_levels);
    d_uniform_bdry_mass_stiffness.resize(d_max_levels);

    d_directors.resize(d_max_levels);

    d_enable_instrumentation.resize(d_max_levels);
    d_instrument_idx.resize(d_max_levels);

    d_enable_sources.resize(d_max_levels);
    d_source_idx.resize(d_max_levels);

    d_global_index_offset.resize(d_max_levels);

    // Determine the various input file names.
    //
    // Prefer to use the new ``structure_names'' key, but revert to the
    // level-by-level ``base_filenames'' keys if necessary.
    if (db->keyExists("structure_names"))
    {
	//Baky code changing here 
	//we are reading only one file and number of structures setting standard names as "name"+id

	d_num_structs_types = db->getArraySize("structure_names");
	d_structs_clones_num.resize(d_num_structs_types);

	std::vector<std::string> structure_type_names(d_num_structs_types);

	db->getStringArray("structure_names", &structure_type_names[0], d_num_structs_types);
	db->getIntegerArray("structs_clones_num", &d_structs_clones_num[0], d_num_structs_types);
	
	int num_strcts=0; 
	for(unsigned itype=0;itype< d_num_structs_types;itype++) num_strcts +=d_structs_clones_num[itype];
      
	for(unsigned itype=0;itype< d_num_structs_types;itype++)
	{
	    std::string struct_name = structure_type_names[itype];
	  
	    if (db->keyExists(struct_name))
	    {
		Pointer<Database> sub_db = db->getDatabase((struct_name));
		if (sub_db->keyExists("level_number"))
		{
		    const int ln = sub_db->getInteger("level_number");
		    if (ln < 0)
		    {
			TBOX_ERROR(d_object_name << ":  "
				   << "Key data `level_number' associated with structure `" << struct_name
				   << "' is negative.");
		    }
		    else if (ln > d_max_levels)
		    {
			TBOX_ERROR(d_object_name << ":  "
				   << "Key data `level_number' associated with structure `" << struct_name
				   << "' is greater than the expected maximum level number "
				   << d_max_levels << ".");
		    }
		    for(int istruct=0;istruct < d_structs_clones_num[itype];istruct++)
		    {
			std::stringstream convert;
			convert<<istruct;
			std::string structure_name = struct_name + convert.str();
			d_base_filename[ln].push_back(structure_name);
		    }
		}
		else
		{
		    TBOX_ERROR(d_object_name << ":  "
			       << "Key data `level_number' not found in structure `" << struct_name
			       << "' input.");
		}
	    }
	    else
	    {
		TBOX_ERROR(d_object_name << ":  "
			   << "Key data `" << struct_name << "' not found in input.");
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
    // Read in any shift and scale factors.
    if (db->keyExists("length_scale_factor"))
    {
	d_length_scale_factor = db->getDouble("length_scale_factor");
    }

    if (db->keyExists("posn_shift"))
    {
	db->getDoubleArray("posn_shift", &d_posn_shift[0], NDIM);
    }

    // Read in any sub-databases associated with the input file names.
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
	const size_t num_base_filename = d_base_filename[ln].size();
	d_enable_springs[ln].resize(num_base_filename, true);
	d_using_uniform_spring_stiffness[ln].resize(num_base_filename, false);
	d_uniform_spring_stiffness[ln].resize(num_base_filename, -1.0);
	d_using_uniform_spring_rest_length[ln].resize(num_base_filename, false);
	d_uniform_spring_rest_length[ln].resize(num_base_filename, -1.0);
	d_using_uniform_spring_force_fcn_idx[ln].resize(num_base_filename, false);
	d_uniform_spring_force_fcn_idx[ln].resize(num_base_filename, -1);

	d_enable_xsprings[ln].resize(num_base_filename, true);
	d_using_uniform_xspring_stiffness[ln].resize(num_base_filename, false);
	d_uniform_xspring_stiffness[ln].resize(num_base_filename, -1.0);
	d_using_uniform_xspring_rest_length[ln].resize(num_base_filename, false);
	d_uniform_xspring_rest_length[ln].resize(num_base_filename, -1.0);
	d_using_uniform_xspring_force_fcn_idx[ln].resize(num_base_filename, false);
	d_uniform_xspring_force_fcn_idx[ln].resize(num_base_filename, -1);

	d_enable_beams[ln].resize(num_base_filename, true);
	d_using_uniform_beam_bend_rigidity[ln].resize(num_base_filename, false);
	d_uniform_beam_bend_rigidity[ln].resize(num_base_filename, -1.0);
	d_using_uniform_beam_curvature[ln].resize(num_base_filename, false);
	d_uniform_beam_curvature[ln].resize(num_base_filename, Vector::Zero());

	d_enable_rods[ln].resize(num_base_filename, true);
	d_using_uniform_rod_properties[ln].resize(num_base_filename, false);
	d_uniform_rod_properties[ln].resize(num_base_filename,
					    array_constant<double, IBRodForceSpec::NUM_MATERIAL_PARAMS>(0.0));

	d_enable_target_points[ln].resize(num_base_filename, true);
	d_using_uniform_target_stiffness[ln].resize(num_base_filename, false);
	d_uniform_target_stiffness[ln].resize(num_base_filename, -1.0);
	d_using_uniform_target_damping[ln].resize(num_base_filename, false);
	d_uniform_target_damping[ln].resize(num_base_filename, -1.0);

	d_enable_anchor_points[ln].resize(num_base_filename, true);

	d_enable_bdry_mass[ln].resize(num_base_filename, true);
	d_using_uniform_bdry_mass[ln].resize(num_base_filename, false);
	d_uniform_bdry_mass[ln].resize(num_base_filename, -1.0);
	d_using_uniform_bdry_mass_stiffness[ln].resize(num_base_filename, false);
	d_uniform_bdry_mass_stiffness[ln].resize(num_base_filename, -1.0);

	d_enable_instrumentation[ln].resize(num_base_filename, true);

	d_enable_sources[ln].resize(num_base_filename, true);
	unsigned j=0;
	for(unsigned itype=0;itype< d_num_structs_types;itype++)
	{
	    for(int istruct=0;istruct < d_structs_clones_num[itype];istruct++)
	    {
		std::string base_filename = d_base_filename[ln][j];
		if (base_filename.size()) base_filename.erase(base_filename.size() - 1-std::max(0,(int)log10(1e-16+istruct))); //erase number
		
		if (db->isDatabase(base_filename))
		{
		    Pointer<Database> sub_db = db->getDatabase(base_filename);
		    
		    // Determine whether to enable or disable any particular
		    // features.
		    if (sub_db->keyExists("enable_springs"))
		    {
			d_enable_springs[ln][j] = sub_db->getBool("enable_springs");
		    }
		    if (sub_db->keyExists("enable_xsprings"))
		    {
			d_enable_xsprings[ln][j] = sub_db->getBool("enable_xsprings");
		    }
		    if (sub_db->keyExists("enable_beams"))
		    {
			d_enable_beams[ln][j] = sub_db->getBool("enable_beams");
		    }
		    if (sub_db->keyExists("enable_rods"))
		    {
			d_enable_rods[ln][j] = sub_db->getBool("enable_rods");
		    }
		    if (sub_db->keyExists("enable_target_points"))
		    {
			d_enable_target_points[ln][j] = sub_db->getBool("enable_target_points");
		    }
		    if (sub_db->keyExists("enable_anchor_points"))
		    {
			d_enable_anchor_points[ln][j] = sub_db->getBool("enable_anchor_points");
		    }
		    if (sub_db->keyExists("enable_bdry_mass"))
		    {
			d_enable_bdry_mass[ln][j] = sub_db->getBool("enable_bdry_mass");
		    }
		    if (sub_db->keyExists("enable_instrumentation"))
		    {
			d_enable_instrumentation[ln][j] = sub_db->getBool("enable_instrumentation");
		    }
		    if (sub_db->keyExists("enable_sources"))
		    {
			d_enable_sources[ln][j] = sub_db->getBool("enable_sources");
		    }
		    
		    // Determine whether to use uniform values for any particular
		    // structure attributes.
		    if (sub_db->keyExists("uniform_spring_stiffness"))
		    {
			d_using_uniform_spring_stiffness[ln][j] = true;
			d_uniform_spring_stiffness[ln][j] = sub_db->getDouble("uniform_spring_stiffness");
			if (d_uniform_spring_stiffness[ln][j] < 0.0)
			{
			    TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_spring_stiffness' "
				       "in database " << base_filename << std::endl
				       << "  spring constant is negative" << std::endl);
			}
		    }
		    if (sub_db->keyExists("uniform_spring_rest_length"))
		    {
			d_using_uniform_spring_rest_length[ln][j] = true;
			d_uniform_spring_rest_length[ln][j] = sub_db->getDouble("uniform_spring_rest_length");
			if (d_uniform_spring_rest_length[ln][j] < 0.0)
			{
			    TBOX_ERROR(d_object_name << ":\n  Invalid entry for key "
				       "`uniform_spring_rest_length' in database " << base_filename
				       << std::endl << "  spring resting length is negative" << std::endl);
			}
		    }
		    if (sub_db->keyExists("uniform_spring_force_fcn_idx"))
		    {
			d_using_uniform_spring_force_fcn_idx[ln][j] = true;
			d_uniform_spring_force_fcn_idx[ln][j] = sub_db->getInteger("uniform_spring_force_fcn_idx");
		    }
		    
		    if (sub_db->keyExists("uniform_xspring_stiffness"))
		    {
			d_using_uniform_xspring_stiffness[ln][j] = true;
			d_uniform_xspring_stiffness[ln][j] = sub_db->getDouble("uniform_xspring_stiffness");
			if (d_uniform_xspring_stiffness[ln][j] < 0.0)
			{
			    TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_xspring_stiffness' "
				       "in database " << base_filename << std::endl
				       << "  spring constant is negative" << std::endl);
			}
		    }
		    if (sub_db->keyExists("uniform_xspring_rest_length"))
		    {
			d_using_uniform_xspring_rest_length[ln][j] = true;
			d_uniform_xspring_rest_length[ln][j] = sub_db->getDouble("uniform_xspring_rest_length");
			if (d_uniform_xspring_rest_length[ln][j] < 0.0)
			{
			    TBOX_ERROR(d_object_name << ":\n  Invalid entry for key "
				       "`uniform_xspring_rest_length' in database " << base_filename
				       << std::endl << "  spring resting length is negative" << std::endl);
			}
		    }
		    if (sub_db->keyExists("uniform_xspring_force_fcn_idx"))
		    {
			d_using_uniform_xspring_force_fcn_idx[ln][j] = true;
			d_uniform_xspring_force_fcn_idx[ln][j] = sub_db->getInteger("uniform_xspring_force_fcn_idx");
		    }

		    if (sub_db->keyExists("uniform_beam_bend_rigidity"))
		    {
			d_using_uniform_beam_bend_rigidity[ln][j] = true;
			d_uniform_beam_bend_rigidity[ln][j] = sub_db->getDouble("uniform_beam_bend_rigidity");
			if (d_uniform_beam_bend_rigidity[ln][j] < 0.0)
			{
			    TBOX_ERROR(d_object_name << ":\n  Invalid entry for key "
				       "`uniform_beam_bend_rigidity' in database " << base_filename
				       << std::endl << "  beam bending rigidity is negative" << std::endl);
			}
		    }
		    if (sub_db->keyExists("uniform_beam_curvature"))
		    {
			d_using_uniform_beam_curvature[ln][j] = true;
			sub_db->getDoubleArray("uniform_beam_curvature", d_uniform_beam_curvature[ln][j].data(), NDIM);
		    }

		    if (sub_db->keyExists("uniform_rod_properties"))
		    {
			d_using_uniform_rod_properties[ln][j] = true;
			sub_db->getDoubleArray("uniform_rod_properties",
					       &d_uniform_rod_properties[ln][j][0],
					       IBRodForceSpec::NUM_MATERIAL_PARAMS);
		    }

		    if (sub_db->keyExists("uniform_target_stiffness"))
		    {
			d_using_uniform_target_stiffness[ln][j] = true;
			d_uniform_target_stiffness[ln][j] = sub_db->getDouble("uniform_target_stiffness");
			if (d_uniform_target_stiffness[ln][j] < 0.0)
			{
			    TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_target_stiffness' "
				       "in database " << base_filename << std::endl
				       << "  target point spring constant is negative" << std::endl);
			}
		    }
		    if (sub_db->keyExists("uniform_target_damping"))
		    {
			d_using_uniform_target_damping[ln][j] = true;
			d_uniform_target_damping[ln][j] = sub_db->getDouble("uniform_target_damping");
			if (d_uniform_target_damping[ln][j] < 0.0)
			{
			    TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_target_damping' in "
				       "database " << base_filename << std::endl
				       << "  target point spring constant is negative" << std::endl);
			}
		    }

		    if (sub_db->keyExists("uniform_bdry_mass"))
		    {
			d_using_uniform_bdry_mass[ln][j] = true;
			d_uniform_bdry_mass[ln][j] = sub_db->getDouble("uniform_bdry_mass");

			if (d_uniform_bdry_mass[ln][j] < 0.0)
			{
			    TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_bdry_mass' in database "
				       << base_filename << std::endl << "  boundary mass is negative"
				       << std::endl);
			}
		    }
		    if (sub_db->keyExists("uniform_bdry_mass_stiffness"))
		    {
			d_using_uniform_bdry_mass_stiffness[ln][j] = true;
			d_uniform_bdry_mass_stiffness[ln][j] = sub_db->getDouble("uniform_bdry_mass_stiffness");

			if (d_uniform_bdry_mass_stiffness[ln][j] < 0.0)
			{
			    TBOX_ERROR(d_object_name << ":\n  Invalid entry for key "
				       "`uniform_bdry_mass_stiffness' in database " << base_filename
				       << std::endl << "  boundary mass spring constant is negative"
				       << std::endl);
			}
		    }
		}

		j++;
	    }
	}
    }
    
    // Output the names of the input files to be read along with additional
    // debugging information.
    pout << d_object_name << ":  Reading from input files.\n";
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
	const size_t num_base_filename = d_base_filename[ln].size();

	unsigned itype=0;
	for (unsigned int j = 0; j < num_base_filename; j += d_structs_clones_num[itype++])/// ++j)
	{
	    std::string base_filename = d_base_filename[ln][j];
	    if (base_filename.size()) base_filename.erase(base_filename.size() - 1); //erase number zero
	    pout << "  base filename: " << base_filename << "\n"
		 << "  assigned to level " << ln << " of the Cartesian grid patch hierarchy\n";
	    if (!d_enable_springs[ln][j])
	    {
		pout << "  NOTE: spring forces are DISABLED for the structure named " << base_filename << "\n";
	    }
	    else
	    {
		if (d_using_uniform_spring_stiffness[ln][j])
		{
		    pout << "  NOTE: UNIFORM spring stiffnesses are being employed for the "
			"structure named " << base_filename << "\n";
		}
		if (d_using_uniform_spring_rest_length[ln][j])
		{
		    pout << "  NOTE: UNIFORM spring resting lengths are being employed for the "
			"structure named " << base_filename << "\n";
		}
		if (d_using_uniform_spring_force_fcn_idx[ln][j])
		{
		    pout << "  NOTE: UNIFORM spring force functions are being employed for the "
			"structure named " << base_filename << "\n";
		}
	    }

	    if (!d_enable_xsprings[ln][j])
	    {
		pout << "  NOTE: crosslink spring forces are DISABLED for the structure named " << base_filename
		     << "\n";
	    }
	    else
	    {
		if (d_using_uniform_xspring_stiffness[ln][j])
		{
		    pout << "  NOTE: UNIFORM crosslink spring stiffnesses are being employed "
			"for "
			"the structure named " << base_filename << "\n";
		}
		if (d_using_uniform_xspring_rest_length[ln][j])
		{
		    pout << "  NOTE: UNIFORM crosslink spring resting lengths are being "
			"employed "
			"for the structure named " << base_filename << "\n";
		}
		if (d_using_uniform_xspring_force_fcn_idx[ln][j])
		{
		    pout << "  NOTE: UNIFORM crosslink spring force functions are being "
			"employed "
			"for the structure named " << base_filename << "\n";
		}
	    }
	    
	    if (!d_enable_beams[ln][j])
	    {
		pout << "  NOTE: beam forces are DISABLED for the structure named " << base_filename << "\n";
	    }
	    else
	    {
		if (d_using_uniform_beam_bend_rigidity[ln][j])
		{
		    pout << "  NOTE: UNIFORM beam bending rigidities are being employed for "
			"the "
			"structure named " << base_filename << "\n";
		}
		if (d_using_uniform_beam_curvature[ln][j])
		{
		    pout << "  NOTE: UNIFORM beam curvatures are being employed for the "
			"structure "
			"named " << base_filename << "\n";
		}
	    }

	    if (!d_enable_rods[ln][j])
	    {
		pout << "  NOTE: rod forces are DISABLED for the structure named " << base_filename << "\n";
	    }
	    else
	    {
		if (d_using_uniform_rod_properties[ln][j])
		{
		    pout << "  NOTE: UNIFORM rod material properties are being employed for "
			"the "
			"structure named " << base_filename << "\n";
		}
	    }

	    if (!d_enable_target_points[ln][j])
	    {
		pout << "  NOTE: target point penalty forces are DISABLED for the structure "
		    "named " << base_filename << "\n";
	    }
	    else
	    {
		if (d_using_uniform_target_stiffness[ln][j])
		{
		    pout << "  NOTE: UNIFORM target point stiffnesses are being employed for "
			"the "
			"structure named " << base_filename << "\n";
		}
		if (d_using_uniform_target_damping[ln][j])
		{
		    pout << "  NOTE: UNIFORM target point damping factors are being employed "
			"for "
			"the structure named " << base_filename << "\n";
		}
	    }

	    if (!d_enable_anchor_points[ln][j])
	    {
		pout << "  NOTE: anchor points are DISABLED for the structure named " << base_filename << "\n";
	    }

	    if (!d_enable_bdry_mass[ln][j])
	    {
		pout << "  NOTE: massive boundary points are DISABLED for the structure named " << base_filename
		     << "\n";
	    }
	    else
	    {
		if (d_using_uniform_bdry_mass[ln][j])
		{
		    pout << "  NOTE: UNIFORM boundary point masses are being employed for the "
			"structure named " << base_filename << "\n";
		}
		if (d_using_uniform_bdry_mass_stiffness[ln][j])
		{
		    pout << "  NOTE: UNIFORM massive boundary point stiffnesses are being "
			"employed "
			"for the structure named " << base_filename << "\n";
		}
	    }

	    if (!d_enable_instrumentation[ln][j])
	    {
		pout << "  NOTE: instrumentation is DISABLED for the structure named " << base_filename << "\n";
	    }

	    if (!d_enable_sources[ln][j])
	    {
		pout << "  NOTE: sources/sinks are DISABLED for the structure named " << base_filename << "\n";
	    }

	    pout << "\n";
	}
    }
    return;
} // getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
