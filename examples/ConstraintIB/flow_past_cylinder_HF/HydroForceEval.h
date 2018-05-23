// Filename: HydroForceEval.h
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

#ifndef included_HydroForceEval
#define included_HydroForceEval

#include <ibtk/ibtk_utilities.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <PatchHierarchy.h> 
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/SAMRAI_MPI.h>

#include <map>
#include <string>
#include <vector>
#include <utility>

namespace IBAMR
{
/////////////////////////////// CLASS DEFINITION /////////////////////////////
/*!
 * \brief Class IBHydrodynamicSurfaceForceEvaluator computes hydrodynamic force and
 * torque on the surface of the immersed bodies. See 
 * <A HREF="https://arxiv.org/pdf/1610.04398.pdf">Computing the Force Distribution 
 * on the Surface of Complex Deforming Geometries using Vortex Methods and Brinkman 
 * Penalization</A> by Verma et al.
 *
 * \note A file with line element connectivity in 2D and triangle element connectivity 
 * in 3D with an extension .elem must be supplied. 
 */
class HydroForceEval
{
/////////////////////////////// PUBLIC ///////////////////////////////////////
public:
    /*!
     * \brief Default constructor.
     */
    HydroForceEval(const std::string& object_name,
		SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
		SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
        IBTK::LDataManager* const l_data_manager);
	
	/*!
     * \brief Destructor.
     */
    ~HydroForceEval();
	
	/*!
	 * \brief Calculate surface forces and torques. 
	 */
	void calcHydroForce(const int u_idx, const int p_idx, const int f_idx,
						SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                        IBTK::LDataManager* const l_data_manager, const double time,
                        const int iteration_num);
	
/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////
private:	
	/*!
     * \brief Get values from input file.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
                      IBTK::LDataManager* const l_data_manager);
	
	/*!
     * \brief Read the surface vertex data from one or more input files.
     */
    void readVertexFiles(const std::string& extension);
	
	/*!
     * \brief Read the surface element data from one or more input files.
     */
    void readElemFiles(const std::string& extension);
    
    /*!
     * \brief Calculate the center of mass of the body which is encompassed by the lifted surface.
     */
    void calcStructCOM(IBTK::LDataManager* const l_data_manager);
	
	/*!
     * \brief Get the finest level number and the finest mesh width. 
     */
    void getGridData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);
    
    /*!
     * \brief Structure to encapsulate the element data.
     */
    struct Elem {
        int posn_idx;
        IBTK::Point X;
        double pressure;
        double vorticity;
        IBTK::Vector pressure_force;
        IBTK::Vector viscous_force;
        double pressure_torque;
        double viscous_torque;
    }; // Elem
    struct elem_cmp {
        bool operator() (const Elem& lhs, const Elem& rhs) const{
            return lhs.posn_idx < rhs.posn_idx;
        }
    };
    
    /*!
     * \brief Print element data. 
     */
    void printData(const std::vector<std::set<Elem, elem_cmp> > surf_elem_set,
                   const double time, const int iteration_num);
	
	/*!
     * \brief Object name.
     */
    std::string d_object_name;
	
	/*!
	 * \brief Fluid viscosity. 
	 */
	double d_mu; 
	
	/*!
	 * \brief Number of structures. 
	 */
	int d_num_structs;
	
	/*!
     * \brief The base filenames of the lifted structures.
	 * 
	 * \note The structure must lie on the finesest grid level. 
     */
    std::vector<std::string> d_struct_names;
    
    /*!
     * \brief The base filenames of the structures encompassed by the lifted structures.
     *
     * \note The structure must lie on the finesest grid level.
     */
    std::map<int, int> d_struct_map;
    
    /*!
     * \brief Is structure stationary. 
     */
    std::map<int, int> d_is_stationary; 
    
    /*!
     * \brief Struture IDs.
     */
    std::set<int> d_lag_struct_id;
    
    /*!
     * \brief Structue C.O.M.
     */
    std::map<int, IBTK::Point> d_X_com;
	
	/*!
	 * \brief Number of verticies in a structure. 
	 */
	std::vector<int> d_num_vertex; 
	
	/*!
	 * \brief Structure vertex positions. 
	 */
	std::vector<std::vector<IBTK::Point> > d_vertex_posn; 
	
	/*!
	 * \brief Number of elements in a structure. 
	 */
	std::vector<int> d_num_elem; 
	
	/*!
	 * \brief Structure element vertex master-slave pairs. 
	 * 
	 * \note First is the master and second is slave. Order is 
	 * used to determine the direction of the surface normal. 
	 */
	std::vector<std::vector<std::pair<int, int> > > d_elem_conn;
    
	/*!
	 * \brief Printing frequency for net and surface distribution 
	 * quantities.
	 */
	int d_print_interval, d_surf_print_interval;
	
	/*!
	 * \brief Output file name string.
	 */
	std::string d_dir_name;
	
	/*!
	 * File streams associated for the output.
	 */
	std::vector<std::ofstream*> d_force_stream;
	
	/*!
	 * \brief Finest level number where the ghost surface lies. 
	 */
	int d_finest_ln;
	
	/*!
	 * \brief Mesh width on the finest level.
	 */
	IBTK::Vector d_mesh_width; 
	
}; // HydroForceEval
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_HydroForceEval
