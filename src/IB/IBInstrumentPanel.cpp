// Filename: IBInstrumentPanel.cpp
// Created on 12 May 2007 by Boyce Griffith
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

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "BasePatchLevel.h"
#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "Eigen/Geometry" // IWYU pragma: keep
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "SideIndex.h"
#include "boost/array.hpp"
#include "boost/multi_array.hpp"
#include "ibamr/IBInstrumentPanel.h"
#include "ibamr/IBInstrumentationSpec.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/ibtk_utilities.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#if defined(IBAMR_HAVE_SILO)
#include <silo.h>
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_initialize_hierarchy_independent_data;
static Timer* t_initialize_hierarchy_dependent_data;
static Timer* t_read_instrument_data;
static Timer* t_write_plot_data;

// The rank of the root MPI process and the MPI tag number.
static const unsigned int SILO_MPI_ROOT = 0;

// The name of the Silo dumps and database filenames.
static const int SILO_NAME_BUFSIZE = 128;
static const std::string VISIT_DUMPS_FILENAME = "meter_data.visit";
static const std::string SILO_DUMP_DIR_PREFIX = "meter_data.cycle_";
static const std::string SILO_SUMMARY_FILE_PREFIX = "meter_data.cycle_";
static const std::string SILO_SUMMARY_FILE_POSTFIX = ".summary.silo";
static const std::string SILO_PROCESSOR_FILE_PREFIX = "meter_data.proc_";
static const std::string SILO_PROCESSOR_FILE_POSTFIX = ".silo";

void init_meter_elements(boost::multi_array<Point, 2>& X_web,
                         boost::multi_array<Vector, 2>& dA_web,
                         const boost::multi_array<Point, 1>& X_perimeter,
                         const Point& X_centroid)
{
#if (NDIM == 2)
    TBOX_ERROR("no support for 2D flow meters at this time!\n");
    NULL_USE(X_web);
    NULL_USE(dA_web);
    NULL_USE(X_perimeter);
    NULL_USE(X_centroid);
#endif
#if (NDIM == 3)
#if !defined(NDEBUG)
    TBOX_ASSERT(X_web.shape()[0] == X_perimeter.shape()[0]);
    TBOX_ASSERT(dA_web.shape()[0] == X_perimeter.shape()[0]);
#endif
    const int num_perimeter_nodes = static_cast<int>(X_web.shape()[0]);
    const int num_web_nodes = static_cast<int>(X_web.shape()[1]);
    for (int m = 0; m < num_perimeter_nodes; ++m)
    {
        const Point& X_perimeter0(X_perimeter[m]);
        const Vector dX0((X_centroid - X_perimeter0) / static_cast<double>(num_web_nodes));

        const Point& X_perimeter1(X_perimeter[(m + 1) % num_perimeter_nodes]);
        const Vector dX1((X_centroid - X_perimeter1) / static_cast<double>(num_web_nodes));

        // Away from the center of the web, each web patch is a planar
        // quadrilateral.  At the web centroid, the quadrilateral is degenerate,
        // i.e., it is a triangle.
        for (int n = 0; n < num_web_nodes; ++n)
        {
            // Compute the four vertices of the quadrilateral web patch.
            //
            // Note that here the vertices are placed in "standard" (i.e.,
            // "counter-clockwise") orientation.
            const Point X0(X_perimeter0 + static_cast<double>(n) * dX0);
            const Point X1(X_perimeter1 + static_cast<double>(n) * dX1);
            const Point X2(X_perimeter1 + static_cast<double>(n + 1) * dX1);
            const Point X3(X_perimeter0 + static_cast<double>(n + 1) * dX0);

            // Compute the midpoints of the edges of the quadrilateral.
            const Point X01(0.5 * (X0 + X1));
            const Point X12(0.5 * (X1 + X2));
            const Point X23(0.5 * (X2 + X3));
            const Point X30(0.5 * (X3 + X0));

            // Construct a parametric representation of the lines connecting the
            // midpoints of the edges.
            const Point& l0 = X01;
            const Vector d0 = X23 - X01;

            const Point& l1 = X12;
            const Vector d1 = X30 - X12;

            // Compute the centroid as the intersection of the lines connecting
            // the midpoints of the edges.
            const double d0d0 = d0.dot(d0);
            const double d0d1 = d0.dot(d1);
            const double d1d1 = d1.dot(d1);
            const double d0l0 = d0.dot(l0);
            const double d0l1 = d0.dot(l1);
            const double d1l0 = d1.dot(l0);
            const double d1l1 = d1.dot(l1);
            const double t = (-d0l0 * d1d1 + d0l1 * d1d1 + d0d1 * d1l0 - d0d1 * d1l1) / (-d0d1 * d0d1 + d1d1 * d0d0);
            const double s = (d1l0 * d0d0 - d0d1 * d0l0 + d0d1 * d0l1 - d1l1 * d0d0) / (-d0d1 * d0d1 + d1d1 * d0d0);
            X_web[m][n] = 0.5 * (l0 + t * d0 + l1 + s * d1);

            // Compute the area-weighted normal to the quadrilateral web patch,
            // i.e.,
            //
            //    dA = 0.5*((X2-X0) X (X3-X1))
            //
            // Note that by construction, the quadrilateral is guaranteed to lie
            // within a plane.  Also, note that if X2 == X3, the following is
            // simply the formula for the area-weighted normal to a triangle.
            dA_web[m][n] = 0.5 * (X2 - X0).cross(X3 - X1);
        }
    }
#endif
    return;
} // init_meter_elements

double
compute_flow_correction(const boost::multi_array<Vector, 1>& U_perimeter,
                        const Vector& U_centroid,
                        const boost::multi_array<Point, 1>& X_perimeter,
                        const Point& X_centroid)
{
    double U_dot_dA = 0.0;
#if (NDIM == 2)
    TBOX_ERROR("no support for 2D flow meters at this time!\n");
    NULL_USE(U_perimeter);
    NULL_USE(U_centroid);
    NULL_USE(X_perimeter);
    NULL_USE(X_centroid);
#endif
#if (NDIM == 3)
    const int num_perimeter_nodes = static_cast<int>(X_perimeter.shape()[0]);
    for (int m = 0; m < num_perimeter_nodes; ++m)
    {
        const Vector& U_perimeter0(U_perimeter[m]);
        const Point& X_perimeter0(X_perimeter[m]);

        const Vector& U_perimeter1(U_perimeter[(m + 1) % num_perimeter_nodes]);
        const Point& X_perimeter1(X_perimeter[(m + 1) % num_perimeter_nodes]);

        // Compute the linear interpolation of the velocity at the center of the
        // triangle.
        const Vector U = (U_perimeter0 + U_perimeter1 + U_centroid) / 3.0;

        // Compute the area weighted normal to the triangle.
        const Vector dA = 0.5 * (X_centroid - X_perimeter0).cross(X_centroid - X_perimeter1);

        // Compute the contribution to U_dot_dA.
        U_dot_dA += U.dot(dA);
    }
#endif
    return U_dot_dA;
} // compute_flow_correction

#if defined(IBAMR_HAVE_SILO)
/*!
 * \brief Build a local mesh database entry corresponding to a meter web.
 */
void
build_meter_web(DBfile* dbfile,
                std::string& dirname,
                const boost::multi_array<Point, 2>& X_web,
                const boost::multi_array<Vector, 2>& dA_web,
                const int timestep,
                const double simulation_time)
{
    const int npoints = static_cast<int>(X_web.num_elements());

    std::vector<float> block_X(NDIM * npoints);
    std::vector<float> block_dA(NDIM * npoints);

    for (unsigned int m = 0, i = 0; m < X_web.shape()[0]; ++m)
    {
        for (unsigned int n = 0; n < X_web.shape()[1]; ++n, ++i)
        {
            // Get the coordinate and normal vector data.
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                block_X[d * npoints + i] = float(X_web[m][n][d]);
                block_dA[d * npoints + i] = float(dA_web[m][n][d]);
            }
        }
    }

    // Set the working directory in the Silo database.
    if (DBSetDir(dbfile, dirname.c_str()) == -1)
    {
        TBOX_ERROR("IBInstrumentPanel::build_meter_web():\n"
                   << "  Could not set directory "
                   << dirname
                   << std::endl);
    }

    // Write out the variables.
    int cycle = timestep;
    float time = float(simulation_time);
    double dtime = simulation_time;

    static const int MAX_OPTS = 3;
    DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    DBAddOption(optlist, DBOPT_TIME, &time);
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

    const char* meshname = "mesh";
    std::vector<float*> coords(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        coords[d] = &block_X[d * npoints];
    }

    int ndims = NDIM;

    DBPutPointmesh(dbfile, meshname, ndims, &coords[0], npoints, DB_FLOAT, optlist);

    const char* varname = "scaled_normal";
    std::vector<float*> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        vars[d] = &block_dA[d * npoints];
    }

    DBPutPointvar(dbfile, varname, meshname, ndims, &vars[0], npoints, DB_FLOAT, optlist);

    DBFreeOptlist(optlist);

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        TBOX_ERROR("IBInstrumentPanel::build_meter_web():\n"
                   << "  Could not return to the base directory from subdirectory "
                   << dirname
                   << std::endl);
    }
    return;
} // build_meter_web
#endif

double
linear_interp(const Point& X,
              const Index<NDIM>& i_cell,
              const Point& X_cell,
              const CellData<NDIM, double>& v,
              const Index<NDIM>& /*patch_lower*/,
              const Index<NDIM>& /*patch_upper*/,
              const double* const /*x_lower*/,
              const double* const /*x_upper*/,
              const double* const dx)
{
    boost::array<bool, NDIM> is_lower;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        is_lower[d] = X[d] < X_cell[d];
    }
    double U = 0.0;
#if (NDIM == 3)
    for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
    {
#endif
        for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
        {
            for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
            {
                const Point X_center(X_cell[0] + static_cast<double>(i_shift0) * dx[0],
                                     X_cell[1] + static_cast<double>(i_shift1) * dx[1]
#if (NDIM == 3)
                                     ,
                                     X_cell[2] + static_cast<double>(i_shift2) * dx[2]
#endif
                                     );
                const double wgt =
                    (((X[0] < X_center[0] ? X[0] - (X_center[0] - dx[0]) : (X_center[0] + dx[0]) - X[0]) / dx[0]) *
                     ((X[1] < X_center[1] ? X[1] - (X_center[1] - dx[1]) : (X_center[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                     *
                     ((X[2] < X_center[2] ? X[2] - (X_center[2] - dx[2]) : (X_center[2] + dx[2]) - X[2]) / dx[2])
#endif
                         );
                const Index<NDIM> i(i_shift0 + i_cell(0),
                                    i_shift1 + i_cell(1)
#if (NDIM == 3)
                                        ,
                                    i_shift2 + i_cell(2)
#endif
                                        );
                const CellIndex<NDIM> i_c(i);
                U += v(i_c) * wgt;
            }
        }
#if (NDIM == 3)
    }
#endif
    return U;
} // linear_interp

template <int N>
Eigen::Matrix<double, N, 1>
linear_interp(const Point& X,
              const Index<NDIM>& i_cell,
              const Point& X_cell,
              const CellData<NDIM, double>& v,
              const Index<NDIM>& /*patch_lower*/,
              const Index<NDIM>& /*patch_upper*/,
              const double* const /*x_lower*/,
              const double* const /*x_upper*/,
              const double* const dx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v.getDepth() == N);
#endif
    boost::array<bool, NDIM> is_lower;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        is_lower[d] = X[d] < X_cell[d];
    }
    Eigen::Matrix<double, N, 1> U(Eigen::Matrix<double, N, 1>::Zero());
#if (NDIM == 3)
    for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
    {
#endif
        for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
        {
            for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
            {
                const Point X_center(X_cell[0] + static_cast<double>(i_shift0) * dx[0],
                                     X_cell[1] + static_cast<double>(i_shift1) * dx[1]
#if (NDIM == 3)
                                     ,
                                     X_cell[2] + static_cast<double>(i_shift2) * dx[2]
#endif
                                     );
                const double wgt =
                    (((X[0] < X_center[0] ? X[0] - (X_center[0] - dx[0]) : (X_center[0] + dx[0]) - X[0]) / dx[0]) *
                     ((X[1] < X_center[1] ? X[1] - (X_center[1] - dx[1]) : (X_center[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                     *
                     ((X[2] < X_center[2] ? X[2] - (X_center[2] - dx[2]) : (X_center[2] + dx[2]) - X[2]) / dx[2])
#endif
                         );
                const Index<NDIM> i(i_shift0 + i_cell(0),
                                    i_shift1 + i_cell(1)
#if (NDIM == 3)
                                        ,
                                    i_shift2 + i_cell(2)
#endif
                                        );
                const CellIndex<NDIM> i_c(i);
                for (int k = 0; k < N; ++k)
                {
                    U[k] += v(i_c, k) * wgt;
                }
            }
        }
#if (NDIM == 3)
    }
#endif
    return U;
} // linear_interp

Vector
linear_interp(const Point& X,
              const Index<NDIM>& i_cell,
              const Point& X_cell,
              const SideData<NDIM, double>& v,
              const Index<NDIM>& /*patch_lower*/,
              const Index<NDIM>& /*patch_upper*/,
              const double* const /*x_lower*/,
              const double* const /*x_upper*/,
              const double* const dx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v.getDepth() == 1);
#endif
    Vector U(Vector::Zero());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        boost::array<bool, NDIM> is_lower;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d == axis)
            {
                is_lower[d] = false;
            }
            else
            {
                is_lower[d] = X[d] < X_cell[d];
            }
        }
#if (NDIM == 3)
        for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
        {
#endif
            for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
            {
                for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
                {
                    const Point X_side(X_cell[0] + (static_cast<double>(i_shift0) + (axis == 0 ? -0.5 : 0.0)) * dx[0],
                                       X_cell[1] + (static_cast<double>(i_shift1) + (axis == 1 ? -0.5 : 0.0)) * dx[1]
#if (NDIM == 3)
                                       ,
                                       X_cell[2] + (static_cast<double>(i_shift2) + (axis == 2 ? -0.5 : 0.0)) * dx[2]
#endif
                                       );
                    const double wgt =
                        (((X[0] < X_side[0] ? X[0] - (X_side[0] - dx[0]) : (X_side[0] + dx[0]) - X[0]) / dx[0]) *
                         ((X[1] < X_side[1] ? X[1] - (X_side[1] - dx[1]) : (X_side[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                         *
                         ((X[2] < X_side[2] ? X[2] - (X_side[2] - dx[2]) : (X_side[2] + dx[2]) - X[2]) / dx[2])
#endif
                             );
                    const Index<NDIM> i(i_shift0 + i_cell(0),
                                        i_shift1 + i_cell(1)
#if (NDIM == 3)
                                            ,
                                        i_shift2 + i_cell(2)
#endif
                                            );
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    U[axis] += v(i_s) * wgt;
                }
            }
#if (NDIM == 3)
        }
#endif
    }
    return U;
} // linear_interp
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInstrumentPanel::IBInstrumentPanel(const std::string& object_name, Pointer<Database> input_db)
    : d_object_name(object_name),
      d_initialized(false),
      d_num_meters(0),
      d_num_perimeter_nodes(),
      d_X_centroid(),
      d_X_perimeter(),
      d_X_web(),
      d_dA_web(),
      d_instrument_read_timestep_num(-1),
      d_instrument_read_time(std::numeric_limits<double>::quiet_NaN()),
      d_max_instrument_name_len(-1),
      d_instrument_names(),
      d_flow_values(),
      d_mean_pres_values(),
      d_point_pres_values(),
      d_web_patch_map(),
      d_web_centroid_map(),
      d_plot_directory_name(NDIM == 2 ? "viz_inst2d" : "viz_inst3d"),
      d_output_log_file(false),
      d_log_file_name(NDIM == 2 ? "inst2d.log" : "inst3d.log"),
      d_log_file_stream(),
      d_flow_conv(1.0),
      d_pres_conv(1.0),
      d_flow_units(""),
      d_pres_units("")
{
#if defined(IBAMR_HAVE_SILO)
// intentionally blank
#else
    TBOX_WARNING("IBInstrumentPanel::IBInstrumentPanel(): SILO is not installed; cannot write data." << std::endl);
#endif

    // Initialize object with data read from the input database.
    if (input_db) getFromInput(input_db);

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_initialize_hierarchy_independent_data =
            TimerManager::getManager()->getTimer("IBAMR::IBInstrumentPanel::initializeHierarchyIndependentData()");
        t_initialize_hierarchy_dependent_data =
            TimerManager::getManager()->getTimer("IBAMR::IBInstrumentPanel::initializeHierarchyDependentData()");
        t_read_instrument_data = TimerManager::getManager()->getTimer("IBAMR::IBInstrumentPanel::readInstrumentData()");
        t_write_plot_data = TimerManager::getManager()->getTimer("IBAMR::IBInstrumentPanel::writePlotData()"););
    return;
} // IBInstrumentPanel

IBInstrumentPanel::~IBInstrumentPanel()
{
    // Close the log file stream.
    if (SAMRAI_MPI::getRank() == 0)
    {
        d_log_file_stream.close();
    }
    return;
} // ~IBInstrumentPanel

const std::vector<std::string>&
IBInstrumentPanel::getInstrumentNames() const
{
    return d_instrument_names;
} // getInstrumentNames

const double&
IBInstrumentPanel::getInstrumentDataReadTime() const
{
    return d_instrument_read_time;
} // getInstrumentDataReadTime

const std::vector<double>&
IBInstrumentPanel::getFlowValues() const
{
    return d_flow_values;
} // getFlowValues

const std::vector<double>&
IBInstrumentPanel::getMeanPressureValues() const
{
    return d_mean_pres_values;
} // getMeanPressureValues

const std::vector<double>&
IBInstrumentPanel::getPointwisePressureValues() const
{
    return d_point_pres_values;
} // getPointwisePressureValues

bool
IBInstrumentPanel::isInstrumented() const
{
    if (!d_initialized)
    {
        TBOX_WARNING(d_object_name << "::isInstrumented():\n"
                                   << "  instrument data has not been initialized."
                                   << std::endl);
        return false;
    }
    return (d_num_meters > 0);
} // isInstrumented

void
IBInstrumentPanel::initializeHierarchyIndependentData(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                      LDataManager* const l_data_manager)
{
    IBAMR_TIMER_START(t_initialize_hierarchy_independent_data);

    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Determine how many flow meters/pressure gauges are present in the local
    // data.
    int max_meter_index = -1;
    std::vector<int> max_node_index;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (l_data_manager->levelContainsLagrangianData(ln))
        {
            const Pointer<LMesh> mesh = l_data_manager->getLMesh(ln);
            const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
            for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
            {
                const LNode* const node_idx = *cit;
                const IBInstrumentationSpec* const spec = node_idx->getNodeDataItem<IBInstrumentationSpec>();
                if (spec)
                {
                    const int m = spec->getMeterIndex();
                    max_meter_index = std::max(m, max_meter_index);

                    const int n = spec->getNodeIndex();
                    max_node_index.resize(max_meter_index + 1, -1);
                    max_node_index[m] = std::max(n, max_node_index[m]);
                }
            }
        }
    }

    // Communicate local data to all processes.
    d_num_meters = SAMRAI_MPI::maxReduction(max_meter_index) + 1;
    max_node_index.resize(d_num_meters, -1);
    d_num_perimeter_nodes.clear();
    d_num_perimeter_nodes.resize(d_num_meters, -1);
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        d_num_perimeter_nodes[m] = max_node_index[m] + 1;
    }
    SAMRAI_MPI::maxReduction(d_num_meters > 0 ? &d_num_perimeter_nodes[0] : NULL, d_num_meters);
#if !defined(NDEBUG)
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        TBOX_ASSERT(d_num_perimeter_nodes[m] > 0);
    }
#endif

    // Resize arrays.
    d_X_centroid.resize(d_num_meters);
    d_X_perimeter.resize(d_num_meters);
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        d_X_perimeter[m].resize(boost::extents[d_num_perimeter_nodes[m]]);
    }
    d_X_web.resize(d_num_meters);
    d_dA_web.resize(d_num_meters);
    d_instrument_names = IBInstrumentationSpec::getInstrumentNames();
    if (d_instrument_names.size() != d_num_meters)
    {
        TBOX_WARNING(d_object_name << "::initializeHierarchyIndependentData():\n"
                                   << "  instrument names are not initialized\n"
                                   << "  using default names"
                                   << std::endl);
        d_instrument_names.resize(d_num_meters);
        for (unsigned int m = 0; m < d_num_meters; ++m)
        {
            std::ostringstream meter_stream;
            meter_stream << "meter_" << m;
            d_instrument_names[m] = meter_stream.str();
        }
        IBInstrumentationSpec::setInstrumentNames(d_instrument_names);
    }
    d_flow_values.resize(d_num_meters, std::numeric_limits<double>::quiet_NaN());
    d_mean_pres_values.resize(d_num_meters, std::numeric_limits<double>::quiet_NaN());
    d_point_pres_values.resize(d_num_meters, std::numeric_limits<double>::quiet_NaN());

    // Open the log file stream.
    d_max_instrument_name_len = 0;
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        d_max_instrument_name_len =
            std::max(d_max_instrument_name_len, static_cast<int>(d_instrument_names[m].length()));
    }

    if (d_output_log_file && SAMRAI_MPI::getRank() == 0 && !d_log_file_stream.is_open())
    {
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        if (from_restart)
        {
            d_log_file_stream.open(d_log_file_name.c_str(), std::ios::app);
        }
        else
        {
            d_log_file_stream.open(d_log_file_name.c_str(), std::ios::out);
            if (d_flow_units != "") d_log_file_stream << "flow units: " << d_flow_units << "    ";
            if (d_pres_units != "") d_log_file_stream << "pressure units: " << d_pres_units << "    ";
            if (d_flow_units != "" || d_pres_units != "") d_log_file_stream << "\n";
            d_log_file_stream << std::string(d_max_instrument_name_len, ' ') << "  time       "
                              << "  x_centroid "
                              << "  y_centroid " << (NDIM == 3 ? "  z_centroid " : "") << "  flow rate   "
                              << "  mean pres.  "
                              << "  point pres. "
                              << "\n";
        }
    }

    // Indicate that the hierarchy-independent data has been initialized.
    d_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_hierarchy_independent_data);
    return;
} // initializeHierarchyIndependentData

void
IBInstrumentPanel::initializeHierarchyDependentData(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                    LDataManager* const l_data_manager,
                                                    const int timestep_num,
                                                    const double data_time)
{
    if (!d_initialized)
    {
        initializeHierarchyIndependentData(hierarchy, l_data_manager);
    }
    if (d_num_meters == 0) return;

    IBAMR_TIMER_START(t_initialize_hierarchy_dependent_data);

    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Keep track of the timestep number and time at which the instrument data
    // was collected.
    d_instrument_read_timestep_num = timestep_num;
    d_instrument_read_time = data_time;

    // Loop over all local nodes to determine the positions of the local
    // perimeter nodes.
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            d_X_perimeter[m][n] = Point::Zero();
        }
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (l_data_manager->levelContainsLagrangianData(ln))
        {
            // Extract the local position array.
            Pointer<LData> X_data = l_data_manager->getLData(LDataManager::POSN_DATA_NAME, ln);
            Vec X_vec = X_data->getVec();
            double* X_arr;
            int ierr = VecGetArray(X_vec, &X_arr);
            IBTK_CHKERRQ(ierr);

            // Store the local positions of the perimeter nodes.
            const Pointer<LMesh> mesh = l_data_manager->getLMesh(ln);
            const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
            for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
            {
                const LNode* const node_idx = *cit;
                const IBInstrumentationSpec* const spec = node_idx->getNodeDataItem<IBInstrumentationSpec>();
                if (spec)
                {
                    const int& petsc_idx = node_idx->getLocalPETScIndex();
                    const double* const X = &X_arr[NDIM * petsc_idx];
                    const int m = spec->getMeterIndex();
                    const int n = spec->getNodeIndex();
                    std::copy(X, X + NDIM, d_X_perimeter[m][n].data());
                }
            }

            // Restore the local position array.
            ierr = VecGetArray(X_vec, &X_arr);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Set the positions of all perimeter nodes on all processes.
    std::vector<double> X_perimeter_flattened;
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            X_perimeter_flattened.insert(
                X_perimeter_flattened.end(), d_X_perimeter[m][n].data(), d_X_perimeter[m][n].data() + NDIM);
        }
    }
    SAMRAI_MPI::sumReduction(&X_perimeter_flattened[0], static_cast<int>(X_perimeter_flattened.size()));
    for (unsigned int m = 0, k = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n, ++k)
        {
            std::copy(&X_perimeter_flattened[NDIM * k],
                      (&X_perimeter_flattened[NDIM * k]) + NDIM,
                      d_X_perimeter[m][n].data());
        }
    }

    // Determine the centroid of each perimeter.
    std::fill(d_X_centroid.begin(), d_X_centroid.end(), Point::Zero());
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            d_X_centroid[m] += d_X_perimeter[m][n];
        }
        d_X_centroid[m] /= static_cast<double>(d_num_perimeter_nodes[m]);
    }

    // Determine the maximum distance from perimeter nodes to centroids.
    std::vector<double> r_max(d_num_meters, 0.0);
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            const Vector r(d_X_perimeter[m][n] - d_X_centroid[m]);
            r_max[m] = std::max(r_max[m], r.norm());
        }
    }

    // Determine the finest grid spacing in the Cartesian grid hierarchy.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domainXLower = grid_geom->getXLower();
    const double* const domainXUpper = grid_geom->getXUpper();
    const double* const dx_coarsest = grid_geom->getDx();
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
    const Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];

    const IntVector<NDIM>& ratio_to_level_zero = hierarchy->getPatchLevel(finest_ln)->getRatio();
    boost::array<double, NDIM> dx_finest;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d] / static_cast<double>(ratio_to_level_zero(d));
    }
    const double h_finest = *std::min_element(dx_finest.begin(), dx_finest.end());

    // Build the meter web patch centroids and area elements.
    //
    // Note that we set the number of web nodes in each meter to that the
    // spacing is approximately half a meshwidth on the finest level of the
    // Cartesian grid hierarchy.
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        const int num_web_nodes = 2 * static_cast<int>(ceil(r_max[m] / h_finest));
        d_X_web[m].resize(boost::extents[d_num_perimeter_nodes[m]][num_web_nodes]);
        d_dA_web[m].resize(boost::extents[d_num_perimeter_nodes[m]][num_web_nodes]);
        init_meter_elements(d_X_web[m], d_dA_web[m], d_X_perimeter[m], d_X_centroid[m]);
    }

    // Setup the mappings from cell indices to the web patch and web centroid
    // data.
    //
    // NOTE: Each meter web patch/centroid is assigned to precisely one
    // Cartesian grid cell in precisely one level.  In particular, each web
    // patch is assigned to whichever grid cell is the finest cell that contains
    // the region of physical space in which the centroid of the web patch is
    // located.  Similarly, each web centroid is assigned to which ever grid
    // cell is the finest cell that contains the region of physical space in
    // which the web centroid is located.
    d_web_patch_map.clear();
    d_web_patch_map.resize(finest_ln + 1);
    d_web_centroid_map.clear();
    d_web_centroid_map.resize(finest_ln + 1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = level->getRatio();
        const Box<NDIM> domain_box_level = Box<NDIM>::refine(domain_box, ratio);
        const Index<NDIM>& domain_box_level_lower = domain_box_level.lower();
        const Index<NDIM>& domain_box_level_upper = domain_box_level.upper();
        boost::array<double, NDIM> dx;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d] / static_cast<double>(ratio(d));
        }

        Pointer<PatchLevel<NDIM> > finer_level =
            (ln < finest_ln ? hierarchy->getPatchLevel(ln + 1) : Pointer<BasePatchLevel<NDIM> >(NULL));
        const IntVector<NDIM>& finer_ratio = (ln < finest_ln ? finer_level->getRatio() : IntVector<NDIM>(1));
        const Box<NDIM> finer_domain_box_level = Box<NDIM>::refine(domain_box, finer_ratio);
        const Index<NDIM>& finer_domain_box_level_lower = finer_domain_box_level.lower();
        const Index<NDIM>& finer_domain_box_level_upper = finer_domain_box_level.upper();
        boost::array<double, NDIM> finer_dx;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            finer_dx[d] = dx_coarsest[d] / static_cast<double>(finer_ratio(d));
        }

        for (unsigned int l = 0; l < d_num_meters; ++l)
        {
            // Setup the web patch mapping.
            for (unsigned int m = 0; m < d_X_web[l].shape()[0]; ++m)
            {
                for (unsigned int n = 0; n < d_X_web[l].shape()[1]; ++n)
                {
                    const Point& X = d_X_web[l][m][n];
                    const Index<NDIM> i = IndexUtilities::getCellIndex(
                        X, domainXLower, domainXUpper, dx.data(), domain_box_level_lower, domain_box_level_upper);
                    const Index<NDIM> finer_i = IndexUtilities::getCellIndex(X,
                                                                             domainXLower,
                                                                             domainXUpper,
                                                                             finer_dx.data(),
                                                                             finer_domain_box_level_lower,
                                                                             finer_domain_box_level_upper);
                    if (level->getBoxes().contains(i) &&
                        (ln == finest_ln || !finer_level->getBoxes().contains(finer_i)))
                    {
                        WebPatch p;
                        p.meter_num = l;
                        p.X = &d_X_web[l][m][n];
                        p.dA = &d_dA_web[l][m][n];
                        d_web_patch_map[ln].insert(std::make_pair(i, p));
                    }
                }
            }

            // Setup the web centroid mapping.
            const Point& X = d_X_centroid[l];
            const Index<NDIM> i = IndexUtilities::getCellIndex(
                X, domainXLower, domainXUpper, dx.data(), domain_box_level_lower, domain_box_level_upper);
            const Index<NDIM> finer_i = IndexUtilities::getCellIndex(X,
                                                                     domainXLower,
                                                                     domainXUpper,
                                                                     finer_dx.data(),
                                                                     finer_domain_box_level_lower,
                                                                     finer_domain_box_level_upper);
            if (level->getBoxes().contains(i) && (ln == finest_ln || !finer_level->getBoxes().contains(finer_i)))
            {
                WebCentroid c;
                c.meter_num = l;
                c.X = &d_X_centroid[l];
                d_web_centroid_map[ln].insert(std::make_pair(i, c));
            }
        }
    }

    IBAMR_TIMER_STOP(t_initialize_hierarchy_dependent_data);
    return;
} // initializeHierarchyDependentData

void
IBInstrumentPanel::readInstrumentData(const int U_data_idx,
                                      const int P_data_idx,
                                      const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                      LDataManager* const l_data_manager,
                                      const int timestep_num,
                                      const double data_time)
{
    if (d_num_meters == 0) return;

    IBAMR_TIMER_START(t_read_instrument_data);

    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Ensure we are collecting data at the
    if (timestep_num != d_instrument_read_timestep_num)
    {
        TBOX_ERROR(d_object_name << "::readInstrumentData():\n"
                                 << "  time step number: "
                                 << timestep_num
                                 << " is != instrumentation time step number: "
                                 << d_instrument_read_timestep_num
                                 << std::endl);
    }

    if (!MathUtilities<double>::equalEps(data_time, d_instrument_read_time))
    {
        TBOX_ERROR(d_object_name << "::readInstrumentData():\n"
                                 << "  data read time: "
                                 << data_time
                                 << " is != instrumentation data read time: "
                                 << d_instrument_read_time
                                 << std::endl);
    }

    // Reset the instrument values.
    std::fill(d_flow_values.begin(), d_flow_values.end(), 0.0);
    std::fill(d_mean_pres_values.begin(), d_mean_pres_values.end(), 0.0);
    std::fill(d_point_pres_values.begin(), d_point_pres_values.end(), 0.0);
    std::vector<double> A(d_num_meters, 0.0);

    // Compute the local contributions to the flux of U through the flow meter,
    // the average value of P in the flow meter, and the pointwise value of P at
    // the centroid of the meter.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Index<NDIM>& patch_lower = patch_box.lower();
            const Index<NDIM>& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const x_lower = pgeom->getXLower();
            const double* const x_upper = pgeom->getXUpper();
            const double* const dx = pgeom->getDx();

            Pointer<CellData<NDIM, double> > U_cc_data = patch->getPatchData(U_data_idx);
            Pointer<SideData<NDIM, double> > U_sc_data = patch->getPatchData(U_data_idx);
            Pointer<CellData<NDIM, double> > P_cc_data = patch->getPatchData(P_data_idx);

            for (Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const Index<NDIM>& i = b();
                std::pair<WebPatchMap::const_iterator, WebPatchMap::const_iterator> patch_range =
                    d_web_patch_map[ln].equal_range(i);
                if (patch_range.first != patch_range.second)
                {
                    const Point X_cell(x_lower[0] + dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5),
                                       x_lower[1] + dx[1] * (static_cast<double>(i(1) - patch_lower(1)) + 0.5)
#if (NDIM == 3)
                                           ,
                                       x_lower[2] + dx[2] * (static_cast<double>(i(2) - patch_lower(2)) + 0.5)
#endif
                                           );
                    if (U_cc_data)
                    {
                        for (WebPatchMap::const_iterator it = patch_range.first; it != patch_range.second; ++it)
                        {
                            const int& meter_num = it->second.meter_num;
                            const Point& X = *(it->second.X);
                            const Vector& dA = *(it->second.dA);
                            const Vector U = linear_interp<NDIM>(
                                X, i, X_cell, *U_cc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                            d_flow_values[meter_num] += U.dot(dA);
                        }
                    }
                    if (U_sc_data)
                    {
                        for (WebPatchMap::const_iterator it = patch_range.first; it != patch_range.second; ++it)
                        {
                            const int& meter_num = it->second.meter_num;
                            const Point& X = *(it->second.X);
                            const Vector& dA = *(it->second.dA);
                            const Vector U =
                                linear_interp(X, i, X_cell, *U_sc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                            d_flow_values[meter_num] += U.dot(dA);
                        }
                    }
                    if (P_cc_data)
                    {
                        for (WebPatchMap::const_iterator it = patch_range.first; it != patch_range.second; ++it)
                        {
                            const int& meter_num = it->second.meter_num;
                            const Point& X = *(it->second.X);
                            const Vector& dA = *(it->second.dA);
                            double P =
                                linear_interp(X, i, X_cell, *P_cc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                            d_mean_pres_values[meter_num] += P * dA.norm();
                            A[meter_num] += dA.norm();
                        }
                    }
                }

                std::pair<WebCentroidMap::const_iterator, WebCentroidMap::const_iterator> centroid_range =
                    d_web_centroid_map[ln].equal_range(i);
                if (centroid_range.first != centroid_range.second)
                {
                    const Point X_cell(x_lower[0] + dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5),
                                       x_lower[1] + dx[1] * (static_cast<double>(i(1) - patch_lower(1)) + 0.5)
#if (NDIM == 3)
                                           ,
                                       x_lower[2] + dx[2] * (static_cast<double>(i(2) - patch_lower(2)) + 0.5)
#endif
                                           );
                    if (P_cc_data)
                    {
                        for (WebCentroidMap::const_iterator it = centroid_range.first; it != centroid_range.second;
                             ++it)
                        {
                            const int& meter_num = it->second.meter_num;
                            const Point& X = *(it->second.X);
                            const double P =
                                linear_interp(X, i, X_cell, *P_cc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                            d_point_pres_values[meter_num] = P;
                        }
                    }
                }
            }
        }
    }

    // Synchronize the values across all processes.
    SAMRAI_MPI::sumReduction(&d_flow_values[0], d_num_meters);
    SAMRAI_MPI::sumReduction(&d_mean_pres_values[0], d_num_meters);
    SAMRAI_MPI::sumReduction(&d_point_pres_values[0], d_num_meters);
    SAMRAI_MPI::sumReduction(&A[0], d_num_meters);

    // Normalize the mean pressure.
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        d_mean_pres_values[m] /= A[m];
    }

    // Loop over all local nodes to determine the velocities of the local
    // perimeter nodes.
    std::vector<boost::multi_array<Vector, 1> > U_perimeter(d_num_meters);
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        U_perimeter[m].resize(boost::extents[d_num_perimeter_nodes[m]]);
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            U_perimeter[m][n] = Vector::Zero();
        }
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (l_data_manager->levelContainsLagrangianData(ln))
        {
            // Extract the local velocity array.
            Pointer<LData> U_data = l_data_manager->getLData(LDataManager::VEL_DATA_NAME, ln);
            Vec U_vec = U_data->getVec();
            double* U_arr;
            int ierr = VecGetArray(U_vec, &U_arr);
            IBTK_CHKERRQ(ierr);

            // Store the local velocities of the perimeter nodes.
            const Pointer<LMesh> mesh = l_data_manager->getLMesh(ln);
            const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
            for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
            {
                const LNode* const node_idx = *cit;
                const IBInstrumentationSpec* const spec = node_idx->getNodeDataItem<IBInstrumentationSpec>();
                if (spec)
                {
                    const int& petsc_idx = node_idx->getLocalPETScIndex();
                    const double* const U = &U_arr[NDIM * petsc_idx];
                    const int m = spec->getMeterIndex();
                    const int n = spec->getNodeIndex();
                    std::copy(U, U + NDIM, U_perimeter[m][n].data());
                }
            }

            // Restore the local position array.
            ierr = VecGetArray(U_vec, &U_arr);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Set the velocities of all perimeter nodes on all processes.
    std::vector<double> U_perimeter_flattened;
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            U_perimeter_flattened.insert(
                U_perimeter_flattened.end(), U_perimeter[m][n].data(), U_perimeter[m][n].data() + NDIM);
        }
    }
    SAMRAI_MPI::sumReduction(&U_perimeter_flattened[0], static_cast<int>(U_perimeter_flattened.size()));
    for (unsigned int m = 0, k = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n, ++k)
        {
            std::copy(
                &U_perimeter_flattened[NDIM * k], (&U_perimeter_flattened[NDIM * k]) + NDIM, U_perimeter[m][n].data());
        }
    }

    // Determine the velocity of the centroid of each perimeter.
    std::vector<Vector> U_centroid(d_num_meters, Vector::Zero());
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            U_centroid[m] += U_perimeter[m][n];
        }
        U_centroid[m] /= static_cast<double>(d_num_perimeter_nodes[m]);
    }

    // Correct for the relative motion of the flow meters.
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        d_flow_values[m] -= compute_flow_correction(U_perimeter[m], U_centroid[m], d_X_perimeter[m], d_X_centroid[m]);
    }

    // Output meter data.
    plog << std::string(d_max_instrument_name_len + 94, '*') << "\n";
    plog << d_object_name << "::readInstrumentData():\n";
    plog << std::string(d_max_instrument_name_len, ' ') << "  time       "
         << "  x_centroid "
         << "  y_centroid " << (NDIM == 3 ? "  z_centroid " : "") << "  flow rate   "
         << "  mean pres.  "
         << "  point pres. "
         << "\n";

    outputLogData(plog);
    if (d_flow_units != "") plog << "flow units: " << d_flow_units << "    ";
    if (d_pres_units != "") plog << "pressure units: " << d_pres_units;
    if (d_flow_units != "" || d_pres_units != "") plog << "\n";
    plog << std::string(d_max_instrument_name_len + 94, '*') << "\n";

    if (d_output_log_file && SAMRAI_MPI::getRank() == 0)
    {
        outputLogData(d_log_file_stream);
        d_log_file_stream.flush();
    }

    IBAMR_TIMER_STOP(t_read_instrument_data);
    return;
} // readInstrumentData

void
IBInstrumentPanel::setPlotDirectory(const std::string& plot_directory_name)
{
    d_plot_directory_name = plot_directory_name;
    return;
} // setPlotDirectory

void
IBInstrumentPanel::writePlotData(const int timestep_num, const double simulation_time)
{
    if (d_num_meters == 0) return;

    IBAMR_TIMER_START(t_write_plot_data);
#if defined(IBAMR_HAVE_SILO)
#if !defined(NDEBUG)
    TBOX_ASSERT(timestep_num >= 0);
    TBOX_ASSERT(!d_plot_directory_name.empty());
#endif

    if (timestep_num != d_instrument_read_timestep_num)
    {
        TBOX_ERROR(d_object_name << "::writePlotData():\n"
                                 << "  time step number: "
                                 << timestep_num
                                 << " is != last time step number: "
                                 << d_instrument_read_timestep_num
                                 << std::endl);
    }

    if (d_plot_directory_name.empty())
    {
        TBOX_ERROR(d_object_name << "::writePlotData():\n"
                                 << "  dump directory name is empty"
                                 << std::endl);
    }

    char temp_buf[SILO_NAME_BUFSIZE];
    std::string current_file_name;
    DBfile* dbfile;
    const unsigned int mpi_rank = SAMRAI_MPI::getRank();
    const unsigned int mpi_nodes = SAMRAI_MPI::getNodes();

    // Create the working directory.
    sprintf(temp_buf, "%06d", d_instrument_read_timestep_num);
    std::string current_dump_directory_name = SILO_DUMP_DIR_PREFIX + temp_buf;
    std::string dump_dirname = d_plot_directory_name + "/" + current_dump_directory_name;

    Utilities::recursiveMkdir(dump_dirname);

    // Create one local DBfile per MPI process.
    sprintf(temp_buf, "%04d", mpi_rank);
    current_file_name = dump_dirname + "/" + SILO_PROCESSOR_FILE_PREFIX;
    current_file_name += temp_buf;
    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

    if (!(dbfile = DBCreate(current_file_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB)))
    {
        TBOX_ERROR(d_object_name + "::writePlotData():\n"
                   << "  Could not create DBfile named "
                   << current_file_name
                   << std::endl);
    }

    // Output the web data on the available MPI processes.
    for (unsigned int meter = 0; meter < d_num_meters; ++meter)
    {
        if (meter % mpi_nodes == mpi_rank)
        {
            std::string dirname = d_instrument_names[meter];
            if (DBMkDir(dbfile, dirname.c_str()) == -1)
            {
                TBOX_ERROR(d_object_name + "::writePlotData():\n"
                           << "  Could not create directory named "
                           << dirname
                           << std::endl);
            }
            build_meter_web(dbfile, dirname, d_X_web[meter], d_dA_web[meter], timestep_num, simulation_time);
        }
    }

    DBClose(dbfile);

    if (mpi_rank == SILO_MPI_ROOT)
    {
        // Create and initialize the multimesh Silo database on the root MPI
        // process.
        sprintf(temp_buf, "%06d", d_instrument_read_timestep_num);
        std::string summary_file_name =
            dump_dirname + "/" + SILO_SUMMARY_FILE_PREFIX + temp_buf + SILO_SUMMARY_FILE_POSTFIX;
        if (!(dbfile = DBCreate(summary_file_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB)))
        {
            TBOX_ERROR(d_object_name + "::writePlotData():\n"
                       << "  Could not create DBfile named "
                       << summary_file_name
                       << std::endl);
        }

        int cycle = timestep_num;
        float time = float(simulation_time);
        double dtime = simulation_time;

        static const int MAX_OPTS = 3;
        DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
        DBAddOption(optlist, DBOPT_CYCLE, &cycle);
        DBAddOption(optlist, DBOPT_TIME, &time);
        DBAddOption(optlist, DBOPT_DTIME, &dtime);

        for (unsigned int meter = 0; meter < d_num_meters; ++meter)
        {
            const int proc = meter % mpi_nodes;
            sprintf(temp_buf, "%04d", proc);
            current_file_name = SILO_PROCESSOR_FILE_PREFIX;
            current_file_name += temp_buf;
            current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

            std::string meshname = current_file_name + ":" + d_instrument_names[meter] + "/mesh";
            char* meshname_ptr = const_cast<char*>(meshname.c_str());
            int meshtype = DB_POINTMESH;

            std::string meter_name = d_instrument_names[meter];

            DBPutMultimesh(dbfile, meter_name.c_str(), 1, &meshname_ptr, &meshtype, optlist);

            if (DBMkDir(dbfile, meter_name.c_str()) == -1)
            {
                TBOX_ERROR(d_object_name + "::writePlotData():\n"
                           << "  Could not create directory named "
                           << meter_name
                           << std::endl);
            }

            std::string varname = current_file_name + ":" + d_instrument_names[meter] + "/scaled_normal";
            char* varname_ptr = const_cast<char*>(varname.c_str());
            int vartype = DB_POINTVAR;

            std::string var_name = d_instrument_names[meter] + "_normal";

            DBPutMultivar(dbfile, var_name.c_str(), 1, &varname_ptr, &vartype, optlist);
        }

        DBClose(dbfile);

        // Create or update the dumps file on the root MPI process.
        static bool summary_file_opened = false;
        std::string path = d_plot_directory_name + "/" + VISIT_DUMPS_FILENAME;
        sprintf(temp_buf, "%06d", d_instrument_read_timestep_num);
        std::string file =
            current_dump_directory_name + "/" + SILO_SUMMARY_FILE_PREFIX + temp_buf + SILO_SUMMARY_FILE_POSTFIX;
        if (!summary_file_opened)
        {
            summary_file_opened = true;
            std::ofstream sfile(path.c_str(), std::ios::out);
            sfile << file << std::endl;
            sfile.close();
        }
        else
        {
            std::ofstream sfile(path.c_str(), std::ios::app);
            sfile << file << std::endl;
            sfile.close();
        }
    }
#else
    TBOX_WARNING("IBInstrumentPanel::writePlotData(): SILO is not installed; cannot write data." << std::endl);
#endif // if defined(IBAMR_HAVE_SILO)
    IBAMR_TIMER_STOP(t_write_plot_data);
    return;
} // writePlotData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBInstrumentPanel::getFromInput(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    if (db->keyExists("plot_directory_name")) d_plot_directory_name = db->getString("plot_directory_name");
    if (db->keyExists("output_log_file")) d_output_log_file = db->getBool("output_log_file");
    if (d_output_log_file)
    {
        if (db->keyExists("log_file_name")) d_log_file_name = db->getString("log_file_name");
    }
    if (db->keyExists("flow_conv")) d_flow_conv = db->getDouble("flow_conv");
    if (db->keyExists("pres_conv")) d_pres_conv = db->getDouble("pres_conv");
    if (db->keyExists("flow_units")) d_flow_units = db->getString("flow_units");
    if (db->keyExists("pres_units")) d_pres_units = db->getString("pres_units");
    return;
} // getFromInput

void
IBInstrumentPanel::outputLogData(std::ostream& os)
{
    for (unsigned int m = 0; m < d_num_meters; ++m)
    {
        std::string meter_name(d_instrument_names[m]);
        meter_name.resize(d_max_instrument_name_len, ' ');
        os << meter_name;

        os.setf(std::ios_base::scientific);
        os.precision(5);
        os << "  " << d_instrument_read_time;

        os.setf(std::ios_base::scientific);
        os.precision(5);
        os << "  " << d_X_centroid[m][0];
        os.setf(std::ios_base::scientific);
        os.precision(5);
        os << "  " << d_X_centroid[m][1];
#if (NDIM == 3)
        os.setf(std::ios_base::scientific);
        os.precision(5);
        os << "  " << d_X_centroid[m][2];
#endif

        os.setf(std::ios_base::scientific);
        os.setf(std::ios_base::showpos);
        os.precision(5);
        os << "  " << d_flow_conv * d_flow_values[m];

        os.setf(std::ios_base::scientific);
        os.setf(std::ios_base::showpos);
        os.precision(5);
        os << "  " << d_pres_conv * d_mean_pres_values[m];

        os.setf(std::ios_base::scientific);
        os.setf(std::ios_base::showpos);
        os.precision(5);
        os << "  " << d_pres_conv * d_point_pres_values[m];

        os << "\n";

        os.unsetf(std::ios_base::scientific);
        os.unsetf(std::ios_base::showpos);
    }
    return;
} // outputLogData

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
