// Filename: IBHydrodynamicSurfaceForceEvaluator.cpp
// Created on 22 Oct 2016 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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

#include "HydroForceEval.h"

#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CellIndex.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/Array.h>
#include <tbox/RestartManager.h>
#include <tbox/Utilities.h>
#include <ibtk/IndexUtilities.h>
#include <fstream> 
#include <set>
#include <sstream>


/////////////////////////////// NAMESPACE ////////////////////////////////////

using namespace SAMRAI::geom;
using namespace SAMRAI::hier;
using namespace SAMRAI::pdat;
using namespace SAMRAI::tbox;
using namespace IBTK;

namespace {
	
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
	
} // namespace

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

HydroForceEval::HydroForceEval(const std::string& object_name,
							   Pointer<Database> input_db,
							   Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                               LDataManager* const l_data_manager)
 : d_object_name(object_name),
   d_mu(0.0), 
   d_num_structs(0),
   d_print_interval(1),
   d_surf_print_interval(1),
   d_dir_name("HydroForceDump"), 
   d_finest_ln(0)
{
    getGridData(patch_hierarchy);
	getFromInput(input_db, l_data_manager);
	readVertexFiles(".ghost"); 
	readElemFiles(".elem");
	
	// Do printing operation for processor 0 only.
	bool from_restart = RestartManager::getManager()->isFromRestart();
	if (!from_restart)
	{
		Utilities::recursiveMkdir(d_dir_name);
        
        const int nodes = SAMRAI_MPI::getNodes();
        for (int n = 0; n < nodes; ++n)
            Utilities::recursiveMkdir(d_dir_name + "/node" + std::to_string(n));
	}
	
	if (!SAMRAI_MPI::getRank())
	{
		d_force_stream.resize(d_num_structs);
		for (int struct_no = 0; struct_no < d_num_structs; ++struct_no)
		{
			std::ostringstream force_stream; 
			force_stream << d_dir_name + '/' + d_struct_names[struct_no] + "_HF.txt";
			if (from_restart) d_force_stream[struct_no] = 
				new std::ofstream(force_stream.str().c_str(), std::fstream::app);
			else d_force_stream[struct_no] = 
				new std::ofstream(force_stream.str().c_str(), std::fstream::out);
		}
	}
	
	return;
} // HydroForceEval

HydroForceEval::~HydroForceEval()
{
	// Delete printing streams. 
	if (!SAMRAI_MPI::getRank())
	{
		for (int struct_no = 0; struct_no < d_num_structs; ++struct_no)
			delete (d_force_stream[struct_no]);
	}
	
	return;
} // ~HydroForceEval


void HydroForceEval::calcHydroForce(const int u_idx, const int p_idx, const int /*f_idx*/,
									Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                    LDataManager* const l_data_manager, const double time,
                                    const int iteration_num)
{
	// WARNING: The following calculation is only implemented for 2D, not 3D. 
    
	bool print = !(iteration_num % d_print_interval); 
	bool surf_print = !(iteration_num % d_surf_print_interval); 
	if (!print && !surf_print) return; 
	
	// Determine the grid extents.
	Pointer<PatchLevel<NDIM> > patch_level = patch_hierarchy->getPatchLevel(d_finest_ln);
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_level->getGridGeometry();
    const IntVector<NDIM>& ratio = patch_level->getRatio();
	
	std::vector<std::pair<int, Box<NDIM> > > patch_boxes;
	for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
	{
		Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
		const int patch_num = patch->getPatchNumber();
		const Box<NDIM>& box = patch->getBox();
		patch_boxes.push_back(std::make_pair(patch_num, box));
	}
    
    // Calculate structure center of mass position.
    calcStructCOM(l_data_manager);
    
    // Local elements.
    std::vector<std::set<Elem, elem_cmp> > elem_set(d_num_structs);
    
	// Initialize net pressure and viscous forces and torques to zero. 
    std::vector<Vector> p_net_force(d_num_structs, Vector::Zero()), v_net_force(d_num_structs, Vector::Zero());
    std::vector<double> p_net_torque(d_num_structs, 0.0), v_net_torque(d_num_structs, 0.0);
#if !defined(NDEBUG)
	std::vector<int> local_patch_found(d_num_structs, 0);
#endif
	for (int struct_no = 0; struct_no < d_num_structs; ++struct_no)
	{
		for (int k = 0; k < d_num_elem[struct_no]; ++k)
		{
			const std::pair<int, int> mstr_slv_pair = d_elem_conn[struct_no][k];
			
			const Point X_mstr = d_vertex_posn[struct_no][mstr_slv_pair.first];
			const Point X_slv = d_vertex_posn[struct_no][mstr_slv_pair.second];
            Vector tangent = X_mstr - X_slv;
			Vector normal; 
			normal[0] = -tangent[1];
			normal[1] = tangent[0];
			const Point X_eval = 0.5 * (X_mstr + X_slv); 
            
			const CellIndex<NDIM> cell_idx = IndexUtilities::getCellIndex(&X_eval[0], grid_geom, ratio);
			bool found_local_patch = false;
			int patch_num = -1;
			for (std::vector<std::pair<int, Box<NDIM> > >::const_iterator it = patch_boxes.begin(); 
				 it != patch_boxes.end() && !found_local_patch; ++it)
			{
				const Box<NDIM>& box = it->second;
				found_local_patch = box.contains(cell_idx);
				if (found_local_patch) 
				{
#if !defined(NDEBUG)
					local_patch_found[struct_no] += k + 1; 
#endif
                    Elem elem;
                    elem.posn_idx = k;
                    elem.X = X_eval;
                    
					// Get patch containing X_eval. 
					patch_num = it->first;
					Pointer<Patch<NDIM> > patch = patch_level->getPatch(patch_num);
					
					// Get velocity and pressure data. 
					Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(p_idx);
					Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
               
					// Force due to pressure.
                    elem.pressure = (*p_data)(cell_idx);
                    Vector& p_force = elem.pressure_force;
                    p_force = -(*p_data)(cell_idx) * normal;
                    
					// Force due to viscosity. 
					Matrix vel_grad_tensor = Matrix::Zero();
					for (unsigned int d = 0; d < NDIM; ++d)
					{
						vel_grad_tensor(d, d) = ((*u_data)(SideIndex<NDIM>(cell_idx, d, SideIndex<NDIM>::Upper)) - 
												(*u_data)(SideIndex<NDIM>(cell_idx, d, SideIndex<NDIM>::Lower))) / d_mesh_width[d];
						for (unsigned int axis = 0; axis < NDIM; ++axis)
						{
							if (axis == d) continue; 
							
							CellIndex<NDIM> cell_upper_nbr_idx = cell_idx, cell_lower_nbr_idx = cell_idx;
							cell_upper_nbr_idx(axis) += 1;
							cell_lower_nbr_idx(axis) -= 1;
							
							vel_grad_tensor(d, axis) = 0.25 * (((*u_data)(SideIndex<NDIM>(cell_upper_nbr_idx, d, SideIndex<NDIM>::Lower)) +
														(*u_data)(SideIndex<NDIM>(cell_upper_nbr_idx, d, SideIndex<NDIM>::Upper))) - 
														((*u_data)(SideIndex<NDIM>(cell_lower_nbr_idx, d, SideIndex<NDIM>::Lower)) +
														(*u_data)(SideIndex<NDIM>(cell_lower_nbr_idx, d, SideIndex<NDIM>::Upper)))) / d_mesh_width[axis];
						}
					}
					
					Matrix def_grad_tensor = vel_grad_tensor;
					for (unsigned int d = 0; d < NDIM; ++d)
					{
						for (unsigned int axis = 0; axis < NDIM; ++axis)
						{
							if (axis == d) continue; 
							
							def_grad_tensor(d, axis) = 0.5 * (vel_grad_tensor(d, axis) + vel_grad_tensor(axis, d));
							def_grad_tensor(axis, d) = def_grad_tensor(d, axis);
						}
					}
                    elem.vorticity = vel_grad_tensor(1,0) - vel_grad_tensor(0,1);
                    Vector& v_force = elem.viscous_force;
					v_force = 2.0 * d_mu * def_grad_tensor * normal;
                    
                    // Torque w.r.t. c.o.m.
                    const int struct_id = d_struct_map[struct_no];
                    const Point& X_com = d_X_com[struct_id];
					Vector R = X_eval - X_com;
                    elem.pressure_torque = R[0] * p_force[1] - R[1] * p_force[0];
                    elem.viscous_torque = R[0] * v_force[1] - R[1] * v_force[0];
                    
					// Increment p & v force and torque vectors.
					for (unsigned int d = 0; d < NDIM; ++d)
					{
						p_net_force[struct_no][d] += p_force[d];
						v_net_force[struct_no][d] += v_force[d];
					}
                    p_net_torque[struct_no] += elem.pressure_torque;
                    v_net_torque[struct_no] += elem.viscous_torque;
                    
                    // Store element information for printing.
                    elem_set[struct_no].insert(elem);
                }
			}
		}
		// Sum across all processors. 
		SAMRAI_MPI::sumReduction(&p_net_force[struct_no][0], NDIM);
		SAMRAI_MPI::sumReduction(&v_net_force[struct_no][0], NDIM);
	}
    SAMRAI_MPI::sumReduction(&p_net_torque[0], d_num_structs);
    SAMRAI_MPI::sumReduction(&v_net_torque[0], d_num_structs);
	
#if !defined(NDEBUG)	
	// Check all surface points found. 
	SAMRAI_MPI::sumReduction(&local_patch_found[0], d_num_structs);
	for (int struct_no = 0; struct_no < d_num_structs; ++struct_no)
	{
		TBOX_ASSERT(local_patch_found[struct_no] == 
					0.5 * d_num_elem[struct_no] * (d_num_elem[struct_no] + 1)); 
	}
#endif 

	// Print net force and torque data. 
	if (print && !SAMRAI_MPI::getRank())
	{
		for (int struct_no = 0; struct_no < d_num_structs; ++struct_no)
		{
			*d_force_stream[struct_no] << time << '\t';
			for (int d = 0; d < NDIM; ++d) *d_force_stream[struct_no] << p_net_force[struct_no][d] << '\t';
			for (int d = 0; d < NDIM; ++d) *d_force_stream[struct_no] << v_net_force[struct_no][d] << '\t';
            *d_force_stream[struct_no] << p_net_torque[struct_no] << '\t';
            *d_force_stream[struct_no] << v_net_torque[struct_no] << '\t';
			*d_force_stream[struct_no] << std::endl;
		}
	}
    
    // Print element data.
    if (surf_print) printData(elem_set, time, iteration_num);

	return;
}
/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void HydroForceEval::getFromInput(Pointer<Database> db, LDataManager* const l_data_manager)
{
	// Get printing intervals. 
	d_print_interval = db->getIntegerWithDefault("print_interval", d_print_interval);
	d_surf_print_interval = db->getIntegerWithDefault("surf_print_interval", d_surf_print_interval);
	
	// Get the fluid viscosity. 
	if (db->keyExists("mu_fluid"))
    {
        d_mu = db->getDouble("mu_fluid");
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                                 << "Key data `mu_fluid' not found in input.");
    }
	
	if (db->keyExists("lifted_struct_names"))
    {
		d_num_structs = db->getArraySize("lifted_struct_names");
		std::vector<std::string> struct_names(d_num_structs); 
        db->getStringArray("lifted_struct_names", &struct_names[0], d_num_structs);
        for (int n = 0; n < d_num_structs; ++n)
        {
            d_struct_names.push_back(struct_names[n]);
        }
    }
    else
    {
         TBOX_ERROR(d_object_name << ":  "
                                  << "Key data `lifted_struct_names' not found in input.");
    }
    
    if (db->keyExists("struct_names"))
    {
        const int num_structs = db->getArraySize("struct_names");
        TBOX_ASSERT(num_structs == d_num_structs);
        std::vector<std::string> struct_names(d_num_structs);
        db->getStringArray("struct_names", &struct_names[0], d_num_structs);
        for (int n = 0; n < d_num_structs; ++n)
        {
            const std::string& struct_name = struct_names[n];
            const int struct_id = l_data_manager->getLagrangianStructureID(struct_name, d_finest_ln);
            if (struct_id == -1)
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "Structure `"
                           << struct_name
                           << "' not found on finest level.");
            }
            else
            {
                d_struct_map.insert(std::make_pair(n, struct_id));
                d_lag_struct_id.insert(struct_id);
                Point X_com = l_data_manager->computeLagrangianStructureCenterOfMass(struct_id, d_finest_ln);
                d_X_com.insert(std::make_pair(struct_id, X_com));
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Key data `struct_names' not found in input.");
    }
    
    if (db->keyExists("is_stationary"))
    {
        const int num_structs = db->getArraySize("is_stationary");
        TBOX_ASSERT(num_structs == d_num_structs);
        Array<bool> is_stationary = db->getBoolArray("is_stationary");
        for (int n = 0; n < d_num_structs; ++n)
        {
            const int struct_id = d_struct_map[n];
            d_is_stationary.insert(std::make_pair(struct_id, is_stationary[n]));
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Key data `is_stationary' not found in input.");
    }
    
    d_dir_name = db->getStringWithDefault("output_dirname", d_dir_name) + "/";
	
	return; 
} // getFromInput

void HydroForceEval::readVertexFiles(const std::string& extension)
{
	std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;
	const bool use_file_batons = true; 

	d_num_vertex.resize(d_num_structs, 0);
	d_vertex_posn.resize(d_num_structs);
	for (int j = 0; j < d_num_structs; ++j)
	{
		// Wait for the previous MPI process to finish reading the current file.
		if (use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

		// Ensure that the file exists.
		const std::string vertex_filename = d_struct_names[j] + extension;
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
											"before line 1 of file "
											<< vertex_filename
											<< std::endl);
			}
			else
			{
				line_string = discard_comments(line_string);
				std::istringstream line_stream(line_string);
				if (!(line_stream >> d_num_vertex[j]))
				{
					TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
												"encountered on line 1 of file "
												<< vertex_filename
												<< std::endl);
				}
			}

			if (d_num_vertex[j] <= 0)
			{
				TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
											<< vertex_filename
											<< std::endl);
			}

			// Each successive line provides the initial position of each
			// vertex in the input file.
			d_vertex_posn[j].resize(d_num_vertex[j]);
			for (int k = 0; k < d_num_vertex[j]; ++k)
			{
				Point& X = d_vertex_posn[j][k];
				if (!std::getline(file_stream, line_string))
				{
					TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k + 2
												<< " of file "
												<< vertex_filename
												<< std::endl);
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
														<< k + 2
														<< " of file "
														<< vertex_filename
														<< std::endl);
						}
					}
				}
			}

			// Close the input file.
			file_stream.close();

			plog << d_object_name << ":  "
					<< "read " << d_num_vertex[j] << " vertices from ASCII input file named " << vertex_filename
					<< std::endl
					<< "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
		}
		else
		{
			TBOX_ERROR(d_object_name << ":\n  Cannot find required vertex file: " << vertex_filename << std::endl);
		}

		// Free the next MPI process to start reading the current file.
		if (use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);
	}
        
    // Synchronize the processes.
    if (use_file_batons) SAMRAI_MPI::barrier();
	
    return;
} // readVertexFiles

void HydroForceEval::readElemFiles(const std::string& extension)
{
	std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;
	const bool use_file_batons = true; 

	d_num_elem.resize(d_num_structs); 
	d_elem_conn.resize(d_num_structs);
	for (int j = 0; j < d_num_structs; ++j)
	{
		// Wait for the previous MPI process to finish reading the current file.
		if (use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank - 1, false, j);

		// Ensure that the file exists.
		const std::string elem_filename = d_struct_names[j] + extension;
		std::ifstream file_stream;
		file_stream.open(elem_filename.c_str(), std::ios::in);
		if (file_stream.is_open())
		{
			plog << d_object_name << ":  "
				 << "processing vertex data from ASCII input file named " << elem_filename << std::endl
				 << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

			// The first entry in the file is the number of vertices.
			if (!std::getline(file_stream, line_string))
			{
				TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered "
											"before line 1 of file "
											<< elem_filename
											<< std::endl);
			}
			else
			{
				line_string = discard_comments(line_string);
				std::istringstream line_stream(line_string);
				if (!(line_stream >> d_num_elem[j]))
				{
					TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file "
												"encountered on line 1 of file "
												<< elem_filename
												<< std::endl);
				}
			}

			if (d_num_elem[j] <= 0)
			{
				TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file "
											<< elem_filename
											<< std::endl);
			}

			// Each successive line provides the initial position of each
			// vertex in the input file.
			d_elem_conn[j].resize(d_num_elem[j]);
			for (int k = 0; k < d_num_elem[j]; ++k)
			{
				std::pair<int, int>& elem = d_elem_conn[j][k];
				if (!std::getline(file_stream, line_string))
				{
					TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k + 2
												<< " of file "
												<< elem_filename
												<< std::endl);
				}
				else
				{
					line_string = discard_comments(line_string);
					std::istringstream line_stream(line_string);
					if (!(line_stream >> elem.first))
					{
						TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
													<< k + 2
													<< " of file "
													<< elem_filename
													<< std::endl);
					}
					else if (elem.first < 0 || elem.first >= d_num_elem[j])
					{
						TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
													<< k + 2
													<< " of file "
													<< elem_filename
													<< std::endl);
					}
					if (!(line_stream >> elem.second))
					{
						TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
													<< k + 2
													<< " of file "
													<< elem_filename
													<< std::endl);
					}
					else if (elem.second < 0 || elem.second >= d_num_elem[j])
					{
						TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line "
													<< k + 2
													<< " of file "
													<< elem_filename
													<< std::endl);
					}
				}
			}

			// Close the input file.
			file_stream.close();

			plog << d_object_name << ":  "
					<< "read " << d_num_elem[j] << " vertices from ASCII input file named " << elem_filename
					<< std::endl
					<< "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
		}
		else
		{
			TBOX_ERROR(d_object_name << ":\n  Cannot find required vertex file: " << elem_filename << std::endl);
		}

		// Free the next MPI process to start reading the current file.
		if (use_file_batons && rank != nodes - 1) SAMRAI_MPI::send(&flag, sz, rank + 1, false, j);
	}
        
    // Synchronize the processes.
    if (use_file_batons) SAMRAI_MPI::barrier();
	
	return; 
} // readElemFiles
    
void HydroForceEval::calcStructCOM(LDataManager* const l_data_manager)
{
    // Get center of mass of structure.
    for (std::set<int>::const_iterator it = d_lag_struct_id.begin(); it != d_lag_struct_id.end(); ++it)
    {
        const int struct_id = *it;
        if (d_is_stationary[struct_id]) continue;
        Point X_com = l_data_manager->computeLagrangianStructureCenterOfMass(struct_id, d_finest_ln);
        d_X_com[struct_id] = X_com;
    }
    
    return;
} // calcStructCOM

void HydroForceEval::getGridData(Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
{
	// Get the finest level number. 
	d_finest_ln = patch_hierarchy->getFinestLevelNumber();
	
	// Determine the grid extents.
	Pointer<PatchLevel<NDIM> > patch_level = patch_hierarchy->getPatchLevel(d_finest_ln);
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_level->getGridGeometry();
    const double* const dx0 = grid_geom->getDx();
    const IntVector<NDIM>& ratio = patch_level->getRatio();
	for (unsigned int d = 0; d < NDIM; ++d) 
		d_mesh_width[d] = dx0[d] / static_cast<double>(ratio(d));
	
	return; 
} // getGridData
    
    void HydroForceEval::printData(const std::vector<std::set<Elem, elem_cmp> > elem_set,
                                   const double time, const int iteration_num)
{
    const int rank = SAMRAI_MPI::getRank();
    
    std::ofstream hf_stream;
    for (int struct_no = 0; struct_no < d_num_structs; ++struct_no)
    {
		if (!elem_set[struct_no].size()) continue; 
		hf_stream.open(d_dir_name + "/node" + std::to_string(rank) + "/" + d_struct_names[struct_no] + 
			"_" + std::to_string(iteration_num), std::fstream::out);
		hf_stream << time << '\n';
		hf_stream << elem_set[struct_no].size() << '\n';
        for (std::set<Elem, elem_cmp>::const_iterator it = elem_set[struct_no].begin(); it != elem_set[struct_no].end(); ++it)
        {
            const Elem& elem = *it;
            hf_stream << elem.posn_idx << '\t';
            for (unsigned int d = 0; d < NDIM; ++d) hf_stream << elem.X(d) << '\t';
            hf_stream << elem.pressure << '\t';
            hf_stream << elem.vorticity << '\t';

            for (unsigned int d = 0; d < NDIM; ++d) hf_stream << elem.pressure_force(d) << '\t';
            for (unsigned int d = 0; d < NDIM; ++d) hf_stream << elem.viscous_force(d) << '\t';
            hf_stream << elem.pressure_torque << '\t';
            hf_stream << elem.viscous_torque << '\t';
            hf_stream << '\n';
        }
		hf_stream.close();
    }
    
    return;
} // printData

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
