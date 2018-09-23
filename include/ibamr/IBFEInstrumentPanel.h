// Filename: IBFEInstrumentPanel.h
// Created on 12 April 2018 by Charles Puelz
//
// Copyright (c) 2002-2018, Charles Puelz and Boyce Griffith
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

#ifndef included_IBAMR_IBFEInstrumentPanel
#define included_IBAMR_IBFEInstrumentPanel

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <fstream>

#include "boost/multi_array.hpp"
#include "ibamr/IBFEMethod.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/ibtk_utilities.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/point.h"
#include "libmesh/serial_mesh.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBFEInstrumentPanel provides support for meters to measure flow
 * and pressure.
 */
class IBFEInstrumentPanel
{
public:
    /*!
     * \brief Constructor.
     */
    IBFEInstrumentPanel(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, int part);

    /*!
     * \brief Destructor.
     */
    ~IBFEInstrumentPanel();

    /*!
     * \brief Get data from input file.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * \brief initialize data which don't depend on the Cartesian grid hierarchy.
     *
     * \note this is done only once.
     */
    void initializeHierarchyIndependentData(IBAMR::IBFEMethod* ib_method_ops);

    /*!
     * \brief initialize data which depends on the Cartesian grid hierarchy.
     *
     * \note this is done each time we compute things with the meter.
     */
    void initializeHierarchyDependentData(IBAMR::IBFEMethod* ib_method_ops,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);
    /*!
     * \brief read instrument data.
     */
    void readInstrumentData(int U_data_idx,
                            int P_data_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                            double data_time);

    /*!
     * \return The instrument dump interval
     */
    int getInstrumentDumpInterval() const;
    
    /*!
     * \return The name of the directory for output
     */
    std::string getPlotDirectoryName() const;

    /*!
     * \return The number of meter meshes
     */
    int getNumberOfMeterMeshes() const;

    /*!
     * \return A reference to the jjth meter mesh
     */
    libMesh::MeshBase& getMeterMesh(const unsigned int jj) const;

    /*!
     * \return A reference to the EquationSystems object for jjth meter mesh
     */
    libMesh::EquationSystems& getMeterMeshEquationSystems(const unsigned int jj) const;

    /*!
     * \return The name for jjth meter mesh
     */
    const std::string& getMeterMeshName(const unsigned int jj) const;

    /*!
     * \return An ordered vector of nodes for the jjth meter mesh
     */
    const std::vector<libMesh::Point>& getMeterMeshNodes(const unsigned int jj) const;

    /*!
     * \return The quadrature order for the jjth meter mesh
     */
    libMesh::Order getMeterMeshQuadOrder(const unsigned int jj) const;

private:
    /*!
     * \brief initialize data which depend on the FE equation systems for
     * the meter mesh.  this includes computing the max radius of the mesh
     * in the current configuration, and updating the mesh's system data.
     */
    void initializeSystemDependentData(IBAMR::IBFEMethod* ib_method_ops, int meter_mesh_number);

    /*!
     * \brief write out data to file.
     */
    void outputData(double data_time);

    /*!
     * \brief write out nodes.
     */
    void outputNodes();

    /*!
     * \brief vector storing radius for each meter
     */
    std::vector<double> d_meter_radii;

    /*!
     * \brief get the maximum radius of the meter in its current configuration
     */
    double getMeterRadius(int meter_mesh_number);

    /*!
     * \brief number of mesh meters.
     */
    unsigned int d_num_meters;

    /*!
     * \brief quadrature order used for the meter meshes.
     */
    std::vector<libMesh::Order> d_quad_order;

    /*!
     * \brief quadrature order from input file
     */
    libMesh::Order d_input_quad_order;

    /*!
     * \brief whether we use a grid based quadrature rule.
     * if false, then we default to using a high order Gauss
     * quadrature.
     */
    bool d_use_adaptive_quadrature;

    /*!
     * \brief quadrature type used for the meter meshes.
     */
    libMesh::QuadratureType d_quad_type;

    /*!
     * \brief part ID where the meter mesh lives, i.e. its parent mesh.
     */
    unsigned int d_part;

    /*!
     * \brief true if meter meshes and other data are built and initialized.
     */
    bool d_initialized;

    /*!
     * \brief number of nodes in the perimeter of the meter mesh.
     */
    std::vector<unsigned int> d_num_nodes;

    /*!
     * \brief vectors to store the dof indices for the velocity and displacement
     * systems in the parent mesh.  this is used to ensure the velocity
     * and displacement systems for the meter mesh have the same values as
     * in the parent mesh.
     * dimension 1 = number of meter meshes
     * dimension 2 = number of mesh nodes
     * dimension 3 = NDIM
     */
    std::vector<std::vector<std::vector<libMesh::dof_id_type> > > d_U_dof_idx;
    std::vector<std::vector<std::vector<libMesh::dof_id_type> > > d_dX_dof_idx;

    /*!
     * \brief a vector containing the nodes of each meter mesh.
     */
    std::vector<std::vector<libMesh::Point> > d_nodes;

    /*!
     * \brief a vector storing the dof indices for each meter mesh.
     */
    std::vector<std::vector<libMesh::dof_id_type> > d_node_dof_IDs;

    /*!
     * \brief contains pointers to the equation systems for the meter mesh.
     */
    std::vector<libMesh::EquationSystems*> d_meter_systems;

    /*!
     * \brief vector of meter mesh pointers.
     */
    std::vector<libMesh::SerialMesh*> d_meter_meshes;

    /*!
     * \brief names for each meter mesh.
     */
    std::vector<std::string> d_meter_mesh_names;

    /*!
     * \brief contains the nodeset IDs on the parent mesh, for the nodesets
     * from which the meter meshes are built.
     */
    SAMRAI::tbox::Array<int> d_nodeset_IDs_for_meters;

    /*!
     * \brief things for data io.
     */
    int d_instrument_dump_interval;
    std::vector<double> d_flow_values, d_mean_pressure_values;
    std::string d_plot_directory_name;
    std::ofstream d_mean_pressure_stream;
    std::ofstream d_flux_stream;

    /*!
     * \brief for ordering objects in a multimap.
     */
    struct IndexFortranOrder : public std::binary_function<SAMRAI::hier::Index<NDIM>, SAMRAI::hier::Index<NDIM>, bool>
    {
        inline bool operator()(const SAMRAI::hier::Index<NDIM>& lhs, const SAMRAI::hier::Index<NDIM>& rhs) const
        {
            return (lhs(0) < rhs(0)
#if (NDIM > 1)
                    || (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                    || (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
            );
        }
    };

    /*!
     * \brief struct for storing information about the quadrature points.
     */
    struct QuadPointStruct
    {
        int meter_num;
        IBTK::Vector normal;
        IBTK::Vector qp_xyz_current;
        double JxW;
    };

    /*!
     * \brief a multimap which associates SAMRAI indices with quadrature point structures.
     */
    typedef std::multimap<SAMRAI::hier::Index<NDIM>, QuadPointStruct, IndexFortranOrder> QuadPointMap;
    std::vector<QuadPointMap> d_quad_point_map;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBFEInstrumentPanel
