// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_IBFEInstrumentPanel
#define included_IBAMR_IBFEInstrumentPanel

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include "ibamr/IBFEMethod.h"

#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/serial_mesh.h"

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
     * \brief Default destructor.
     */
    ~IBFEInstrumentPanel() = default;

    /*!
     * \brief Initialize hierarchy-independent data.
     *
     * The data initialized by this method are assumed \em not to change during
     * the course of a simulation.
     */
    void initializeHierarchyIndependentData(IBAMR::IBFEMethod* ib_method_ops);

    /*!
     * \brief Compute the flow rates and pressures in the various distributed
     * internal flow meters and pressure gauges.
     */
    void readInstrumentData(int U_data_idx,
                            int P_data_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                            IBAMR::IBFEMethod* ib_method_ops,
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
     * \return A reference to the specified meter mesh
     */
    libMesh::MeshBase& getMeterMesh(unsigned int meter_idx) const;

    /*!
     * \return A reference to the EquationSystems object for specified meter mesh
     */
    libMesh::EquationSystems& getMeterMeshEquationSystems(unsigned int meter_idx) const;

    /*!
     * \return The name for specified meter mesh
     */
    const std::string& getMeterMeshName(unsigned int meter_idx) const;

    /*!
     * \return The quadrature rule for the specified meter mesh.
     */
    libMesh::QuadratureType getMeterMeshQuadType(unsigned int meter_idx) const;

    /*!
     * \return The quadrature order for the specified meter mesh
     */
    libMesh::Order getMeterMeshQuadOrder(unsigned int meter_idx) const;

    /*!
     * \return The flow rates for each meter mesh.
     */
    const std::vector<double>& getMeterFlowRates() const;

    /*!
     * \return The mean pressures for each meter mesh.
     */
    const std::vector<double>& getMeterMeanPressures() const;

    /*!
     * \return The centroid pressures for each meter mesh.
     */
    const std::vector<double>& getMeterCentroidPressures() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBFEInstrumentPanel() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBFEInstrumentPanel(const IBFEInstrumentPanel& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFEInstrumentPanel& operator=(const IBFEInstrumentPanel& that) = delete;

    /*!
     * \brief Get data from input file.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * \brief Initialize data that depend on the FE equation systems for
     * the meter mesh.  this includes computing the max radius of the mesh
     * in the current configuration, and updating the mesh's system data.
     */
    void resetMeterConfiguration(IBAMR::IBFEMethod* ib_method_ops, int meter_mesh_number);

    /*!
     * \brief Initialize mappings from Cartesian grid patches to meter quadrature points.
     */
    void computeMeterQuadratureData(std::vector<std::map<int, std::vector<int> > >& meter_idx_map,
                                    std::vector<std::map<int, std::vector<IBTK::Vector> > >& meter_x_map,
                                    std::vector<std::map<int, std::vector<IBTK::Vector> > >& meter_u_corr_map,
                                    std::vector<std::map<int, std::vector<IBTK::Vector> > >& meter_normal_map,
                                    std::vector<std::map<int, std::vector<double> > >& meter_JxW_map,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                    const IBFEMethod* ib_method_ops);

    /*!
     * \brief Write out data to file.
     */
    void outputData(double data_time);

    /*!
     * \brief Get the maximum radius of the meter in its current configuration.
     */
    double getMeterRadius(int meter_idx);

    /*!
     * \brief The number of meters.
     */
    unsigned int d_num_meters = 0;

    /*!
     * \brief Quadrature rule for computing surface integrals.
     */
    libMesh::QuadratureType d_quad_type = libMesh::QGRID;

    /*!
     * \brief Default quadrature order for computing surface integrals.
     */
    libMesh::Order d_default_quad_order = libMesh::FIFTH;

    /*!
     * \brief Flag tracking whether to use an adaptive quadrature order rule.
     *
     * \note This only works with QGRID.
     */
    bool d_use_adaptive_quadrature = true;

    /*!
     * \brief Interpolation type for velocity data.
     */
    std::string d_u_interp_fcn = "PIECEWISE_LINEAR";

    /*!
     * \brief Interpolation type for pressure data.
     */
    std::string d_p_interp_fcn = "PIECEWISE_LINEAR";

    /*!
     * \brief part ID where the meter mesh lives, i.e., its parent mesh.
     */
    unsigned int d_part;

    /*!
     * \brief Flag tracking if meter meshes and other data are built and initialized.
     */
    bool d_initialized = false;

    /*!
     * \brief Actual quadrature order used for each meter mesh.
     */
    std::vector<libMesh::Order> d_quad_order;

    /*!
     * \brief The number of nodes in the perimeter of the meter mesh.
     *
     * \note Each mesh also has an extra node at the centroid of the perimeter nodes.
     */
    std::vector<unsigned int> d_num_perim_nodes;

    /*!
     * \brief Position and velocity DOF indices that link the meter mesh to the parent mesh data structures.
     */
    std::vector<std::vector<std::array<libMesh::dof_id_type, NDIM> > > d_x_dof_idx, d_u_dof_idx;

    /*!
     * \brief The meter meshes.
     */
    std::vector<std::unique_ptr<libMesh::SerialMesh> > d_meter_meshes;

    /*!
     * \brief Equation systems for the meter meshes.
     */
    std::vector<std::unique_ptr<libMesh::EquationSystems> > d_meter_systems;

    /*!
     * \brief Names for each meter mesh.
     */
    std::vector<std::string> d_meter_mesh_names;

    /*!
     * \brief Nodeset IDs on the structure mesh that specify the mesh nodes that form the perimeter of each meter mesh.
     */
    SAMRAI::tbox::Array<int> d_perimeter_nodeset_ids;

    /*!
     * \brief Meter radius data.
     */
    std::vector<double> d_meter_radii;

    /*!
     * \brief Flow rate data for each meter.
     */
    std::vector<double> d_flow_rate_values;

    /*!
     * \brief Mean pressure data for each meter.
     */
    std::vector<double> d_mean_pressure_values;

    /*!
     * \brief Centroid pressure data for each meter.
     */
    std::vector<double> d_centroid_pressure_values;

    /*!
     * \brief I/O settings.
     */
    int d_instrument_dump_interval;
    std::string d_plot_directory_name;
    std::ofstream d_flux_stream, d_mean_pressure_stream, d_centroid_pressure_stream;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBAMR_HAVE_LIBMESH
#endif //#ifndef included_IBAMR_IBFEInstrumentPanel
