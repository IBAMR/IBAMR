// Filename: IBInstrumentPanel.C
// Last modified: <13.Dec.2009 15:53:28 griffith@griffith-macbook-pro.local>
// Created on 12 May 2007 by Boyce Griffith (boyce@trasnaform2.local)

#include "IBInstrumentPanel.h"

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
#include <ibamr/IBInstrumentationSpec.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LNodeIndexData.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// BLITZ++ INCLUDES
#include <blitz/tinyvec-et.h>

// SILO INCLUDES
#if HAVE_LIBSILO
extern "C"
{
#include <silo.h>
}
#endif

// C++ STDLIB INCLUDES
#include <algorithm>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy_independent_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy_dependent_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_read_instrument_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_write_plot_data;

// The rank of the root MPI process and the MPI tag number.
static const int SILO_MPI_ROOT = 0;
static const int SILO_MPI_TAG  = 0;

// The name of the Silo dumps and database filenames.
static const int SILO_NAME_BUFSIZE = 128;
static const std::string VISIT_DUMPS_FILENAME = "meter_data.visit";
static const std::string SILO_DUMP_DIR_PREFIX = "meter_data.cycle_";
static const std::string SILO_SUMMARY_FILE_PREFIX= "meter_data.cycle_";
static const std::string SILO_SUMMARY_FILE_POSTFIX = ".summary.silo";
static const std::string SILO_PROCESSOR_FILE_PREFIX = "meter_data.proc_";
static const std::string SILO_PROCESSOR_FILE_POSTFIX = ".silo";

void
init_meter_elements(
    blitz::Array<blitz::TinyVector<double,NDIM>,2>& X_web,
    blitz::Array<blitz::TinyVector<double,NDIM>,2>& dA_web,
    const blitz::Array<blitz::TinyVector<double,NDIM>,1>& X_perimeter,
    const blitz::TinyVector<double,NDIM>& X_centroid)
{
#if (NDIM == 2)
    TBOX_ERROR("no support for 2D flow meters at this time!\n");
#endif
#if (NDIM == 3)
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(X_web.extent(0) == X_perimeter.extent(0));
    TBOX_ASSERT(dA_web.extent(0) == X_perimeter.extent(0));
#endif
    const int num_perimeter_nodes = X_web.extent(0);
    const int num_web_nodes = X_web.extent(1);
    for (int m = 0; m < num_perimeter_nodes; ++m)
    {
        const blitz::TinyVector<double,NDIM>& X_perimeter0(X_perimeter(m));
        const blitz::TinyVector<double,NDIM> dX0((X_centroid-X_perimeter0)/double(num_web_nodes));

        const blitz::TinyVector<double,NDIM>& X_perimeter1(X_perimeter((m+1)%num_perimeter_nodes));
        const blitz::TinyVector<double,NDIM> dX1((X_centroid-X_perimeter1)/double(num_web_nodes));

        // Away from the center of the web, each web patch is a planar
        // quadrilateral.  At the web centroid, the quadrilateral is degenerate,
        // i.e., it becomes a triangle.
        //
        // The only change required between these two cases is that the centroid
        // of a nondegenerate quadrilateral is different from that of a
        // degenerate quadrilateral.
        for (int n = 0; n < num_web_nodes; ++n)
        {
            // Compute the four vertices of the quadrilateral web patch.
            //
            // Note that here the vertices are placed in "standard" (i.e.,
            // "counter-clockwise") orientation.
            const blitz::TinyVector<double,NDIM> X0(X_perimeter0+double(n  )*dX0);
            const blitz::TinyVector<double,NDIM> X1(X_perimeter1+double(n  )*dX1);
            const blitz::TinyVector<double,NDIM> X2(X_perimeter1+double(n+1)*dX1);
            const blitz::TinyVector<double,NDIM> X3(X_perimeter0+double(n+1)*dX0);

            // Compute the centroid of the quadrilateral web patch.
            X_web(m,n) = ((n+1 < num_web_nodes) ? blitz::TinyVector<double,NDIM>((X0+X1+X2+X3)/4.0) : blitz::TinyVector<double,NDIM>((X0+X1+X2)/3.0));

            // Compute the area-weighted normal to the quadrilateral web patch,
            // i.e.,
            //
            //    dA = 0.5*((X2-X0) X (X3-X1))
            //
            // Note that by construction, the quadrilateral is guaranteed to lie
            // within a plane.  Also, note that if X2 == X3, the following is
            // simply the forumla for the area-weighted normal to a triangle.
            dA_web(m,n) = 0.5*blitz::cross(blitz::TinyVector<double,NDIM>(X2-X0),blitz::TinyVector<double,NDIM>(X3-X1));
        }
    }
#endif
    return;
}// init_meter_elements

double
compute_flow_correction(
    const blitz::Array<blitz::TinyVector<double,NDIM>,1>& U_perimeter,
    const blitz::TinyVector<double,NDIM>& U_centroid,
    const blitz::Array<blitz::TinyVector<double,NDIM>,1>& X_perimeter,
    const blitz::TinyVector<double,NDIM>& X_centroid)
{
    double U_dot_dA = 0.0;
#if (NDIM == 2)
    TBOX_ERROR("no support for 2D flow meters at this time!\n");
#endif
#if (NDIM == 3)
    const int num_perimeter_nodes = X_perimeter.extent(0);
    for (int m = 0; m < num_perimeter_nodes; ++m)
    {
        const blitz::TinyVector<double,NDIM>& U_perimeter0(U_perimeter(m));
        const blitz::TinyVector<double,NDIM>& X_perimeter0(X_perimeter(m));

        const blitz::TinyVector<double,NDIM>& U_perimeter1(U_perimeter((m+1)%num_perimeter_nodes));
        const blitz::TinyVector<double,NDIM>& X_perimeter1(X_perimeter((m+1)%num_perimeter_nodes));

        // Compute the linear interpolation of the velocity at the center of the
        // triangle.
        const blitz::TinyVector<double,NDIM> U = blitz::TinyVector<double,NDIM>((U_perimeter0+U_perimeter1+U_centroid)/3.0);

        // Compute the area weighted normal to the triangle.
        const blitz::TinyVector<double,NDIM> dA = 0.5*blitz::cross(blitz::TinyVector<double,NDIM>(X_centroid-X_perimeter0),blitz::TinyVector<double,NDIM>(X_centroid-X_perimeter1));

        // Compute the contribution to U_dot_dA.
        U_dot_dA += blitz::dot(U,dA);
    }
#endif
    return U_dot_dA;
}// compute_flow_correction

#if HAVE_LIBSILO
/*!
 * \brief Build a local mesh database entry corresponding to a meter web.
 */
void
build_meter_web(
    DBfile* dbfile,
    std::string& dirname,
    const blitz::Array<blitz::TinyVector<double,NDIM>,2>& X_web,
    const blitz::Array<blitz::TinyVector<double,NDIM>,2>& dA_web,
    const int timestep,
    const double simulation_time)
{
    const int npoints = X_web.numElements();

    std::vector<float> block_X (NDIM*npoints);
    std::vector<float> block_dA(NDIM*npoints);

    for (int m = 0, i = 0; m < X_web.extent(0); ++m)
    {
        for (int n = 0; n < X_web.extent(1); ++n, ++i)
        {
            // Get the coordinate and normal vector data.
            for (int d = 0; d < NDIM; ++d)
            {
                block_X [d*npoints+i] = float( X_web(m,n)(d));
                block_dA[d*npoints+i] = float(dA_web(m,n)(d));
            }
        }
    }

    // Set the working directory in the Silo database.
    if (DBSetDir(dbfile, dirname.c_str()) == -1)
    {
        TBOX_ERROR("IBInstrumentPanel::build_meter_web():\n"
                   << "  Could not set directory " << dirname << std::endl);
    }

    // Write out the variables.
    int    cycle = timestep;
    float  time  = simulation_time;
    double dtime = simulation_time;

    static const int MAX_OPTS = 3;
    DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    DBAddOption(optlist, DBOPT_TIME , &time);
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

    const char* meshname = "mesh";
    std::vector<float*> coords(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        coords[d] = &block_X[d*npoints];
    }

    int ndims = NDIM;

    DBPutPointmesh(dbfile, meshname, ndims, &coords[0], npoints, DB_FLOAT, optlist);

    const char* varname = "scaled_normal";
    std::vector<float*> vars(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        vars[d] = &block_dA[d*npoints];
    }

    DBPutPointvar(dbfile, varname, meshname, ndims, &vars[0], npoints, DB_FLOAT, optlist);

    DBFreeOptlist(optlist);

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        TBOX_ERROR("IBInstrumentPanel::build_meter_web():\n"
                   << "  Could not return to the base directory from subdirectory " << dirname << std::endl);
    }
    return;
}// build_meter_web
#endif

template <int DEPTH>
inline blitz::TinyVector<double,DEPTH>
linear_interp(
    const blitz::TinyVector<double,NDIM>& X,
    const SAMRAI::hier::Index<NDIM>& i_cell,
    const blitz::TinyVector<double,NDIM>& X_cell,
    const SAMRAI::pdat::CellData<NDIM,double>& v,
    const SAMRAI::hier::Index<NDIM>& patch_lower,
    const SAMRAI::hier::Index<NDIM>& patch_upper,
    const double* const xLower,
    const double* const xUpper,
    const double* const dx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(v.getDepth() == DEPTH);
#endif
    const blitz::TinyVector<bool,NDIM> is_lower(X < X_cell);
    blitz::TinyVector<double,DEPTH> U(0.0);
#if (NDIM == 3)
    for (int i_shift2 = (is_lower(2) ? -1 : 0); i_shift2 <= (is_lower(2) ? 0 : 1); ++i_shift2)
    {
#endif
        for (int i_shift1 = (is_lower(1) ? -1 : 0); i_shift1 <= (is_lower(1) ? 0 : 1); ++i_shift1)
        {
            for (int i_shift0 = (is_lower(0) ? -1 : 0); i_shift0 <= (is_lower(0) ? 0 : 1); ++i_shift0)
            {
                const blitz::TinyVector<double,NDIM> X_center(X_cell(0)+double(i_shift0)*dx[0],
                                                              X_cell(1)+double(i_shift1)*dx[1]
#if (NDIM == 3)
                                                              ,
                                                              X_cell(2)+double(i_shift2)*dx[2]
#endif
                                                              );
                const double wgt = (((X(0) < X_center(0) ? X(0) - (X_center(0)-dx[0]) : (X_center(0)+dx[0]) - X(0))/dx[0])*
                                    ((X(1) < X_center(1) ? X(1) - (X_center(1)-dx[1]) : (X_center(1)+dx[1]) - X(1))/dx[1])
#if (NDIM == 3)
                                    *
                                    ((X(2) < X_center(2) ? X(2) - (X_center(2)-dx[2]) : (X_center(2)+dx[2]) - X(2))/dx[2])
#endif
                                    );
                const SAMRAI::hier::Index<NDIM> i(i_shift0+i_cell(0),
                                                  i_shift1+i_cell(1)
#if (NDIM == 3)
                                                  ,
                                                  i_shift2+i_cell(2)
#endif
                                                  );
                const SAMRAI::pdat::CellIndex<NDIM> i_c(i);
                for (int d = 0; d < DEPTH; ++d)
                {
                    U(d) += v(i_c,d)*wgt;
                }
            }
        }
#if (NDIM == 3)
    }
#endif
    return U;
}// linear_interp

inline blitz::TinyVector<double,NDIM>
linear_interp(
    const blitz::TinyVector<double,NDIM>& X,
    const SAMRAI::hier::Index<NDIM>& i_cell,
    const blitz::TinyVector<double,NDIM>& X_cell,
    const SAMRAI::pdat::SideData<NDIM,double>& v,
    const SAMRAI::hier::Index<NDIM>& patch_lower,
    const SAMRAI::hier::Index<NDIM>& patch_upper,
    const double* const xLower,
    const double* const xUpper,
    const double* const dx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(v.getDepth() == 1);
#endif
    blitz::TinyVector<double,NDIM> U(0.0);
    for (int axis = 0; axis < NDIM; ++axis)
    {
        blitz::TinyVector<bool,NDIM> is_lower(X < X_cell);
        is_lower(axis) = false;
#if (NDIM == 3)
        for (int i_shift2 = (is_lower(2) ? -1 : 0); i_shift2 <= (is_lower(2) ? 0 : 1); ++i_shift2)
        {
#endif
            for (int i_shift1 = (is_lower(1) ? -1 : 0); i_shift1 <= (is_lower(1) ? 0 : 1); ++i_shift1)
            {
                for (int i_shift0 = (is_lower(0) ? -1 : 0); i_shift0 <= (is_lower(0) ? 0 : 1); ++i_shift0)
                {
                    const blitz::TinyVector<double,NDIM> X_side(X_cell(0)+(double(i_shift0) + (axis == 0 ? -0.5 : 0.0))*dx[0],
                                                                X_cell(1)+(double(i_shift1) + (axis == 1 ? -0.5 : 0.0))*dx[1]
#if (NDIM == 3)
                                                                ,
                                                                X_cell(2)+(double(i_shift2) + (axis == 2 ? -0.5 : 0.0))*dx[2]
#endif
                                                                );
                    const double wgt = (((X(0) < X_side(0) ? X(0) - (X_side(0)-dx[0]) : (X_side(0)+dx[0]) - X(0))/dx[0])*
                                        ((X(1) < X_side(1) ? X(1) - (X_side(1)-dx[1]) : (X_side(1)+dx[1]) - X(1))/dx[1])
#if (NDIM == 3)
                                        *
                                        ((X(2) < X_side(2) ? X(2) - (X_side(2)-dx[2]) : (X_side(2)+dx[2]) - X(2))/dx[2])
#endif
                                        );
                    const SAMRAI::hier::Index<NDIM> i(i_shift0+i_cell(0),
                                                      i_shift1+i_cell(1)
#if (NDIM == 3)
                                                      ,
                                                      i_shift2+i_cell(2)
#endif
                                                      );
                    const SAMRAI::pdat::SideIndex<NDIM> i_s(i, axis, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                    U(axis) += v(i_s)*wgt;
                }
            }
#if (NDIM == 3)
        }
#endif
    }
    return U;
}// linear_interp
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInstrumentPanel::IBInstrumentPanel(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
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
#if HAVE_LIBSILO
    // intentionally blank
#else
    TBOX_WARNING("IBInstrumentPanel::IBInstrumentPanel(): SILO is not installed; cannot write data." << std::endl);
#endif

    // Initialize object with data read from the input database.
    if (!input_db.isNull())
    {
        getFromInput(input_db);
    }

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_independent_data = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBInstrumentPanel::initializeHierarchyIndependentData()");
        t_initialize_hierarchy_dependent_data = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBInstrumentPanel::initializeHierarchyDependentData()");
        t_read_instrument_data = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBInstrumentPanel::readInstrumentData()");
        t_write_plot_data = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBAMR::IBInstrumentPanel::writePlotData()");
    }
    return;
}// IBInstrumentPanel

IBInstrumentPanel::~IBInstrumentPanel()
{
    // Close the log file stream.
    if (SAMRAI::tbox::SAMRAI_MPI::getRank() == 0)
    {
        d_log_file_stream.close();
    }
    return;
}// ~IBInstrumentPanel

const std::vector<std::string>&
IBInstrumentPanel::getInstrumentNames() const
{
    return d_instrument_names;
}// getInstrumentNames

const double&
IBInstrumentPanel::getInstrumentDataReadTime() const
{
    return d_instrument_read_time;
}// getInstrumentDataReadTime

const std::vector<double>&
IBInstrumentPanel::getFlowValues() const
{
    return d_flow_values;
}// getFlowValues

const std::vector<double>&
IBInstrumentPanel::getMeanPressureValues() const
{
    return d_mean_pres_values;
}// getMeanPressureValues

const std::vector<double>&
IBInstrumentPanel::getPointwisePressureValues() const
{
    return d_point_pres_values;
}// getPointwisePressureValues

bool
IBInstrumentPanel::isInstrumented() const
{
    if (!d_initialized)
    {
        TBOX_WARNING(d_object_name << "::isInstrumented():\n"
                     << "  instrument data has not been initialized." << std::endl);
        return false;
    }
    return (d_num_meters > 0);
}// isInstrumented

void
IBInstrumentPanel::initializeHierarchyIndependentData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    IBTK::LDataManager* const lag_manager)
{
    t_initialize_hierarchy_independent_data->start();

    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager-> getLNodeIndexPatchDescriptorIndex();

    // Determine how many flow meters/pressure gauges are present in the local
    // data.
    int max_meter_index = -1;
    std::vector<int> max_node_index;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (lag_manager->levelContainsLagrangianData(ln))
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
                for (IBTK::LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
                     it != idx_data->lnode_index_end(); ++it)
                {
                    const IBTK::LNodeIndex& node_idx = *it;
                    SAMRAI::tbox::Pointer<IBInstrumentationSpec> spec = node_idx.getStashData<IBInstrumentationSpec>();
                    if (!spec.isNull())
                    {
                        const int m = spec->getMeterIndex();
                        max_meter_index = std::max(m, max_meter_index);

                        const int n = spec->getNodeIndex();
                        max_node_index.resize(max_meter_index+1,-1);
                        max_node_index[m] = std::max(n, max_node_index[m]);
                    }
                }
            }
        }
    }

    // Communicate local data to all processes.
    d_num_meters = SAMRAI::tbox::SAMRAI_MPI::maxReduction(max_meter_index)+1;
    max_node_index.resize(d_num_meters,-1);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_num_meters >= 0);
#endif
    d_num_perimeter_nodes.clear();
    d_num_perimeter_nodes.resize(d_num_meters,-1);
    for (int m = 0; m < d_num_meters; ++m)
    {
        d_num_perimeter_nodes[m] = max_node_index[m]+1;
    }
    SAMRAI::tbox::SAMRAI_MPI::maxReduction(d_num_meters > 0 ? &d_num_perimeter_nodes[0] : NULL,
                                           d_num_meters);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int m = 0; m < d_num_meters; ++m)
    {
        TBOX_ASSERT(d_num_perimeter_nodes[m] > 0);
    }
#endif

    // Resize arrays.
    d_X_centroid.resize(d_num_meters);
    d_X_perimeter.resize(d_num_meters);
    for (int m = 0; m < d_num_meters; ++m)
    {
        d_X_perimeter[m].resize(d_num_perimeter_nodes[m]);
    }
    d_X_web.resize(d_num_meters);
    d_dA_web.resize(d_num_meters);
    d_instrument_names = IBInstrumentationSpec::getInstrumentNames();
    if (int(d_instrument_names.size()) != d_num_meters)
    {
        TBOX_WARNING(d_object_name << "::initializeHierarchyIndependentData():\n"
                     << "  instrument names are not initialized\n"
                     << "  using default names" << std::endl);
        d_instrument_names.resize(d_num_meters);
        for (int m = 0; m < d_num_meters; ++m)
        {
            std::ostringstream meter_stream;
            meter_stream << "meter_" << m;
            d_instrument_names[m] = meter_stream.str();
        }
        IBInstrumentationSpec::setInstrumentNames(d_instrument_names);
    }
    d_flow_values      .resize(d_num_meters,std::numeric_limits<double>::quiet_NaN());
    d_mean_pres_values .resize(d_num_meters,std::numeric_limits<double>::quiet_NaN());
    d_point_pres_values.resize(d_num_meters,std::numeric_limits<double>::quiet_NaN());

    // Open the log file stream.
    d_max_instrument_name_len = 0;
    for (int m = 0; m < d_num_meters; ++m)
    {
        d_max_instrument_name_len = std::max(d_max_instrument_name_len, int(d_instrument_names[m].length()));
    }

    if (d_output_log_file && SAMRAI::tbox::SAMRAI_MPI::getRank() == 0 && !d_log_file_stream.is_open())
    {
        const bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
        if (from_restart)
        {
            d_log_file_stream.open(d_log_file_name.c_str(),std::ios::app);
        }
        else
        {
            d_log_file_stream.open(d_log_file_name.c_str(),std::ios::out);
            if (d_flow_units != "") d_log_file_stream << "flow units: " << d_flow_units << "    ";
            if (d_pres_units != "") d_log_file_stream << "pressure units: " << d_pres_units << "    ";
            if (d_flow_units != "" || d_pres_units != "") d_log_file_stream << "\n";
            d_log_file_stream << std::string(d_max_instrument_name_len,' ')
                              << "  time       "
                              << "  x_centroid "
                              << "  y_centroid "
                              << (NDIM == 3 ? "  z_centroid " : "")
                              << "  flow rate   "
                              << "  mean pres.  "
                              << "  point pres. "
                              << "\n";
        }
    }

    // Indicate that the hierarchy-independent data has been initialized.
    d_initialized = true;

    t_initialize_hierarchy_independent_data->stop();
    return;
}// initializeHierarchyIndependentData

void
IBInstrumentPanel::initializeHierarchyDependentData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    IBTK::LDataManager* const lag_manager,
    const int timestep_num,
    const double data_time)
{
    if (!d_initialized)
    {
        initializeHierarchyIndependentData(hierarchy, lag_manager);
    }
    if (d_num_meters == 0) return;

    t_initialize_hierarchy_dependent_data->start();

    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Keep track of the timestep number and time at which the instrument data
    // was collected.
    d_instrument_read_timestep_num = timestep_num;
    d_instrument_read_time = data_time;

    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Loop over all local nodes to determine the positions of the local
    // perimeter nodes.
    for (int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            d_X_perimeter[m](n) = 0.0;
        }
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (lag_manager->levelContainsLagrangianData(ln))
        {
            // Extract the local position array.
            SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data = lag_manager->getLNodeLevelData(IBTK::LDataManager::POSN_DATA_NAME,ln);
            Vec X_vec = X_data->getGlobalVec();
            double* X_arr;
            int ierr = VecGetArray(X_vec, &X_arr);  IBTK_CHKERRQ(ierr);

            // Store the local positions of the perimeter nodes.
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
                for (IBTK::LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
                     it != idx_data->lnode_index_end(); ++it)
                {
                    const IBTK::LNodeIndex& node_idx = *it;
                    SAMRAI::tbox::Pointer<IBInstrumentationSpec> spec = node_idx.getStashData<IBInstrumentationSpec>();
                    if (!spec.isNull())
                    {
                        const int& petsc_idx = node_idx.getLocalPETScIndex();
                        const double* const X = &X_arr[NDIM*petsc_idx];
                        const int m = spec->getMeterIndex();
                        const int n = spec->getNodeIndex();
                        std::copy(X,X+NDIM,d_X_perimeter[m](n).data());
                    }
                }
            }

            // Restore the local position array.
            ierr = VecGetArray(X_vec, &X_arr);  IBTK_CHKERRQ(ierr);
        }
    }

    // Set the positions of all perimeter nodes on all processes.
    std::vector<double> X_perimeter_flattened;
    for (int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            X_perimeter_flattened.insert(X_perimeter_flattened.end(),d_X_perimeter[m](n).data(),d_X_perimeter[m](n).data()+NDIM);
        }
    }
    SAMRAI::tbox::SAMRAI_MPI::sumReduction(&X_perimeter_flattened[0],X_perimeter_flattened.size());
    for (int m = 0, k = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n, ++k)
        {
            std::copy(&X_perimeter_flattened[NDIM*k],(&X_perimeter_flattened[NDIM*k])+NDIM,d_X_perimeter[m](n).data());
        }
    }

    // Determine the centroid of each perimeter.
    std::fill(d_X_centroid.begin(),d_X_centroid.end(),blitz::TinyVector<double,NDIM>(0.0));
    for (int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            d_X_centroid[m] += d_X_perimeter[m](n);
        }
        d_X_centroid[m] /= double(d_num_perimeter_nodes[m]);
    }

    // Determine the maximum distance from perimeter nodes to centroids.
    std::vector<double> r_max(d_num_meters,0.0);
    for (int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            const blitz::TinyVector<double,NDIM> r(d_X_perimeter[m](n)-d_X_centroid[m]);
            r_max[m] = std::max(r_max[m],sqrt(blitz::dot(r,r)));
        }
    }

    // Determine the finest grid spacing in the Cartesian grid hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domainXLower = grid_geom->getXLower();
    const double* const domainXUpper = grid_geom->getXUpper();
    const double* const dx_coarsest = grid_geom->getDx();
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
    const SAMRAI::hier::Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];

    const SAMRAI::hier::IntVector<NDIM>& ratio_to_level_zero = hierarchy->getPatchLevel(finest_ln)->getRatio();
    std::vector<double> dx_finest(NDIM,0.0);
    for (int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d]/double(ratio_to_level_zero(d));
    }
    const double h_finest = *(std::min_element(dx_finest.begin(),dx_finest.end()));

    // Build the meter web patch centroids and area elements.
    //
    // Note that we set the number of web nodes in each meter to that the
    // spacing is approximately half a meshwidth on the finest level of the
    // Cartesian grid hierarchy.
    for (int m = 0; m < d_num_meters; ++m)
    {
        const int num_web_nodes = 2*int(ceil(r_max[m]/h_finest));
        d_X_web [m].resize(d_num_perimeter_nodes[m],num_web_nodes);
        d_dA_web[m].resize(d_num_perimeter_nodes[m],num_web_nodes);
        init_meter_elements(d_X_web[m],d_dA_web[m],d_X_perimeter[m],d_X_centroid[m]);
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
    d_web_patch_map.resize(finest_ln+1);
    d_web_centroid_map.clear();
    d_web_centroid_map.resize(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& ratio = level->getRatio();
        const SAMRAI::hier::Box<NDIM> domain_box_level = SAMRAI::hier::Box<NDIM>::refine(domain_box, ratio);
        const SAMRAI::hier::Index<NDIM>& domain_box_level_lower = domain_box_level.lower();
        const SAMRAI::hier::Index<NDIM>& domain_box_level_upper = domain_box_level.upper();
        double dx[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d]/double(ratio(d));
        }

        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > finer_level = (ln < finest_ln ? hierarchy->getPatchLevel(ln+1) : SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL));
        const SAMRAI::hier::IntVector<NDIM>& finer_ratio = (ln < finest_ln ? finer_level->getRatio() : SAMRAI::hier::IntVector<NDIM>(1));
        const SAMRAI::hier::Box<NDIM> finer_domain_box_level = SAMRAI::hier::Box<NDIM>::refine(domain_box, finer_ratio);
        const SAMRAI::hier::Index<NDIM>& finer_domain_box_level_lower = finer_domain_box_level.lower();
        const SAMRAI::hier::Index<NDIM>& finer_domain_box_level_upper = finer_domain_box_level.upper();
        double finer_dx[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            finer_dx[d] = dx_coarsest[d]/double(finer_ratio(d));
        }

        for (int l = 0; l < d_num_meters; ++l)
        {
            // Setup the web patch mapping.
            for (int m = 0; m < d_X_web[l].extent(0); ++m)
            {
                for (int n = 0; n < d_X_web[l].extent(1); ++n)
                {
                    const double* const X = d_X_web[l](m,n).data();
                    const SAMRAI::hier::Index<NDIM> i = IBTK::IndexUtilities::getCellIndex(X, domainXLower, domainXUpper, dx, domain_box_level_lower, domain_box_level_upper);
                    const SAMRAI::hier::Index<NDIM> finer_i = IBTK::IndexUtilities::getCellIndex(X, domainXLower, domainXUpper, finer_dx, finer_domain_box_level_lower, finer_domain_box_level_upper);
                    if (level->getBoxes().contains(i) && (ln == finest_ln || !finer_level->getBoxes().contains(finer_i)))
                    {
                        WebPatch p;
                        p.meter_num = l;
                        p.X  = &d_X_web [l](m,n);
                        p.dA = &d_dA_web[l](m,n);
                        d_web_patch_map[ln].insert(std::make_pair(i,p));
                    }
                }
            }

            // Setup the web centroid mapping.
            const double* const X = d_X_centroid[l].data();
            const SAMRAI::hier::Index<NDIM> i = IBTK::IndexUtilities::getCellIndex(X, domainXLower, domainXUpper, dx, domain_box_level_lower, domain_box_level_upper);
            const SAMRAI::hier::Index<NDIM> finer_i = IBTK::IndexUtilities::getCellIndex(X, domainXLower, domainXUpper, finer_dx, finer_domain_box_level_lower, finer_domain_box_level_upper);
            if (level->getBoxes().contains(i) && (ln == finest_ln || !finer_level->getBoxes().contains(finer_i)))
            {
                WebCentroid c;
                c.meter_num = l;
                c.X = &d_X_centroid[l];
                d_web_centroid_map[ln].insert(std::make_pair(i,c));
            }
        }
    }

    t_initialize_hierarchy_dependent_data->stop();
    return;
}// initializeHierarchyDependentData

void
IBInstrumentPanel::readInstrumentData(
    const int U_data_idx,
    const int P_data_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    IBTK::LDataManager* const lag_manager,
    const int timestep_num,
    const double data_time)
{
    if (d_num_meters == 0) return;

    t_read_instrument_data->start();

    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Ensure we are collecting data at the
    if (timestep_num != d_instrument_read_timestep_num)
    {
        TBOX_ERROR(d_object_name << "::readInstrumentData():\n"
                   << "  time step number: " << timestep_num
                   << " is != instrumentation time step number: " << d_instrument_read_timestep_num
                   << std::endl);
    }

    if (!SAMRAI::tbox::MathUtilities<double>::equalEps(data_time, d_instrument_read_time))
    {
        TBOX_ERROR(d_object_name << "::readInstrumentData():\n"
                   << "  data read time: " << data_time
                   << " is != instrumentation data read time: " << d_instrument_read_time
                   << std::endl);
    }

    // Reset the instrument values.
    std::fill(d_flow_values      .begin(), d_flow_values      .end(), 0.0);
    std::fill(d_mean_pres_values .begin(), d_mean_pres_values .end(), 0.0);
    std::fill(d_point_pres_values.begin(), d_point_pres_values.end(), 0.0);
    std::vector<double> A(d_num_meters, 0.0);

    // Compute the local contributions to the flux of U through the flow meter,
    // the average value of P in the flow meter, and the pointwise value of P at
    // the centroid of the meter.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
            const SAMRAI::hier::Index<NDIM>& patch_upper = patch_box.upper();

            const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const xLower = pgeom->getXLower();
            const double* const xUpper = pgeom->getXUpper();
            const double* const dx = pgeom->getDx();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_cc_data = patch->getPatchData(U_data_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > U_sc_data = patch->getPatchData(U_data_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > P_cc_data = patch->getPatchData(P_data_idx);

            for (SAMRAI::hier::Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const SAMRAI::hier::Index<NDIM>& i = b();
                std::pair<WebPatchMap::const_iterator,WebPatchMap::const_iterator> patch_range = d_web_patch_map[ln].equal_range(i);
                if (patch_range.first != patch_range.second)
                {
                    const blitz::TinyVector<double,NDIM> X_cell(xLower[0] + dx[0]*(double(i(0)-patch_lower(0))+0.5),
                                                                xLower[1] + dx[1]*(double(i(1)-patch_lower(1))+0.5)
#if (NDIM == 3)
                                                                ,
                                                                xLower[2] + dx[2]*(double(i(2)-patch_lower(2))+0.5)
#endif
                                                                );
                    if (!U_cc_data.isNull())
                    {
                        for (WebPatchMap::const_iterator it = patch_range.first; it != patch_range.second; ++it)
                        {
                            const int& meter_num = (*it).second.meter_num;
                            const blitz::TinyVector<double,NDIM>& X = *((*it).second.X);
                            const blitz::TinyVector<double,NDIM>& dA = *((*it).second.dA);
                            const blitz::TinyVector<double,NDIM> U = linear_interp<NDIM>(X, i, X_cell, *U_cc_data, patch_lower, patch_upper, xLower, xUpper, dx);
                            d_flow_values[meter_num] += blitz::dot(U,dA);
                        }
                    }
                    if (!U_sc_data.isNull())
                    {
                        for (WebPatchMap::const_iterator it = patch_range.first; it != patch_range.second; ++it)
                        {
                            const int& meter_num = (*it).second.meter_num;
                            const blitz::TinyVector<double,NDIM>& X = *((*it).second.X);
                            const blitz::TinyVector<double,NDIM>& dA = *((*it).second.dA);
                            const blitz::TinyVector<double,NDIM> U = linear_interp(X, i, X_cell, *U_sc_data, patch_lower, patch_upper, xLower, xUpper, dx);
                            d_flow_values[meter_num] += blitz::dot(U,dA);
                        }
                    }
                    if (!P_cc_data.isNull())
                    {
                        for (WebPatchMap::const_iterator it = patch_range.first; it != patch_range.second; ++it)
                        {
                            const int& meter_num = (*it).second.meter_num;
                            const blitz::TinyVector<double,NDIM>& X = *((*it).second.X);
                            const blitz::TinyVector<double,NDIM>& dA = *((*it).second.dA);
                            const blitz::TinyVector<double,1> P = linear_interp<1>(X, i, X_cell, *P_cc_data, patch_lower, patch_upper, xLower, xUpper, dx);
                            d_mean_pres_values[meter_num] += P(0)*norm(dA);
                            A                 [meter_num] += norm(dA);
                        }
                    }
                }

                std::pair<WebCentroidMap::const_iterator,WebCentroidMap::const_iterator> centroid_range = d_web_centroid_map[ln].equal_range(i);
                if (centroid_range.first != centroid_range.second)
                {
                    const blitz::TinyVector<double,NDIM> X_cell(xLower[0] + dx[0]*(double(i(0)-patch_lower(0))+0.5),
                                                                xLower[1] + dx[1]*(double(i(1)-patch_lower(1))+0.5)
#if (NDIM == 3)
                                                                ,
                                                                xLower[2] + dx[2]*(double(i(2)-patch_lower(2))+0.5)
#endif
                                                                );
                    if (!P_cc_data.isNull())
                    {
                        for (WebCentroidMap::const_iterator it = centroid_range.first; it != centroid_range.second; ++it)
                        {
                            const int& meter_num = (*it).second.meter_num;
                            const blitz::TinyVector<double,NDIM>& X = *((*it).second.X);
                            const blitz::TinyVector<double,1> P = linear_interp<1>(X, i, X_cell, *P_cc_data, patch_lower, patch_upper, xLower, xUpper, dx);
                            d_point_pres_values[meter_num] = P(0);
                        }
                    }
                }
            }
        }
    }

    // Synchronize the values across all processes.
    SAMRAI::tbox::SAMRAI_MPI::sumReduction(&d_flow_values      [0],d_num_meters);
    SAMRAI::tbox::SAMRAI_MPI::sumReduction(&d_mean_pres_values [0],d_num_meters);
    SAMRAI::tbox::SAMRAI_MPI::sumReduction(&d_point_pres_values[0],d_num_meters);
    SAMRAI::tbox::SAMRAI_MPI::sumReduction(&A                  [0],d_num_meters);

    // Normalize the mean pressure.
    for (int m = 0; m < d_num_meters; ++m)
    {
        d_mean_pres_values[m] /= A[m];
    }

    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager-> getLNodeIndexPatchDescriptorIndex();

    // Loop over all local nodes to determine the velocities of the local
    // perimeter nodes.
    std::vector<blitz::Array<blitz::TinyVector<double,NDIM>,1> > U_perimeter(d_num_meters);
    for (int m = 0; m < d_num_meters; ++m)
    {
        U_perimeter[m].resize(d_num_perimeter_nodes[m]);
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            U_perimeter[m](n) = 0.0;
        }
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (lag_manager->levelContainsLagrangianData(ln))
        {
            // Extract the local velocity array.
            SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> U_data = lag_manager->getLNodeLevelData(IBTK::LDataManager::VEL_DATA_NAME,ln);
            Vec U_vec = U_data->getGlobalVec();
            double* U_arr;
            int ierr = VecGetArray(U_vec, &U_arr);  IBTK_CHKERRQ(ierr);

            // Store the local velocities of the perimeter nodes.
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
                for (IBTK::LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
                     it != idx_data->lnode_index_end(); ++it)
                {
                    const IBTK::LNodeIndex& node_idx = *it;
                    SAMRAI::tbox::Pointer<IBInstrumentationSpec> spec = node_idx.getStashData<IBInstrumentationSpec>();
                    if (!spec.isNull())
                    {
                        const int& petsc_idx = node_idx.getLocalPETScIndex();
                        const double* const U = &U_arr[NDIM*petsc_idx];
                        const int m = spec->getMeterIndex();
                        const int n = spec->getNodeIndex();
                        std::copy(U,U+NDIM,U_perimeter[m](n).data());
                    }
                }
            }

            // Restore the local position array.
            ierr = VecGetArray(U_vec, &U_arr);  IBTK_CHKERRQ(ierr);
        }
    }

    // Set the velocities of all perimeter nodes on all processes.
    std::vector<double> U_perimeter_flattened;
    for (int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            U_perimeter_flattened.insert(U_perimeter_flattened.end(),U_perimeter[m](n).data(),U_perimeter[m](n).data()+NDIM);
        }
    }
    SAMRAI::tbox::SAMRAI_MPI::sumReduction(&U_perimeter_flattened[0],U_perimeter_flattened.size());
    for (int m = 0, k = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n, ++k)
        {
            std::copy(&U_perimeter_flattened[NDIM*k],(&U_perimeter_flattened[NDIM*k])+NDIM,U_perimeter[m](n).data());
        }
    }

    // Determine the velocity of the centroid of each perimeter.
    std::vector<blitz::TinyVector<double,NDIM> > U_centroid(d_num_meters,blitz::TinyVector<double,NDIM>(0.0));
    for (int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            U_centroid[m] += U_perimeter[m](n);
        }
        U_centroid[m] /= double(d_num_perimeter_nodes[m]);
    }

    // Correct for the relative motion of the flow meters.
    for (int m = 0; m < d_num_meters; ++m)
    {
        d_flow_values[m] -= compute_flow_correction(U_perimeter[m],U_centroid[m],d_X_perimeter[m],d_X_centroid[m]);
    }

    // Output meter data.
    SAMRAI::tbox::plog << std::string(d_max_instrument_name_len+94,'*') << "\n";
    SAMRAI::tbox::plog << d_object_name << "::readInstrumentData():\n";
    SAMRAI::tbox::plog << std::string(d_max_instrument_name_len,' ')
                       << "  time       "
                       << "  x_centroid "
                       << "  y_centroid "
                       << (NDIM == 3 ? "  z_centroid " : "")
                       << "  flow rate   "
                       << "  mean pres.  "
                       << "  point pres. "
                       << "\n";

    outputLogData(SAMRAI::tbox::plog);
    if (d_flow_units != "") SAMRAI::tbox::plog << "flow units: " << d_flow_units << "    ";
    if (d_pres_units != "") SAMRAI::tbox::plog << "pressure units: " << d_pres_units;
    if (d_flow_units !="" || d_pres_units != "") SAMRAI::tbox::plog << "\n";
    SAMRAI::tbox::plog << std::string(d_max_instrument_name_len+94,'*') << "\n";

    if (d_output_log_file && SAMRAI::tbox::SAMRAI_MPI::getRank() == 0)
    {
        outputLogData(d_log_file_stream);
        d_log_file_stream.flush();
    }

    t_read_instrument_data->stop();
    return;
}// readInstrumentData

void
IBInstrumentPanel::writePlotData(
    const int timestep_num,
    const double simulation_time)
{
    if (d_num_meters == 0) return;

    t_write_plot_data->start();
#if HAVE_LIBSILO
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(timestep_num >= 0);
    TBOX_ASSERT(!d_plot_directory_name.empty());
#endif

    if (timestep_num != d_instrument_read_timestep_num)
    {
        TBOX_ERROR(d_object_name << "::writePlotData():\n"
                   << "  time step number: " << timestep_num
                   << " is != last time step number: " << d_instrument_read_timestep_num
                   << std::endl);
    }

    if (d_plot_directory_name.empty())
    {
        TBOX_ERROR(d_object_name << "::writePlotData():\n"
                   << "  dump directory name is empty" << std::endl);
    }

    char temp_buf[SILO_NAME_BUFSIZE];
    std::string current_file_name;
    DBfile* dbfile;
    const int mpi_rank  = SAMRAI::tbox::SAMRAI_MPI::getRank();
    const int mpi_nodes = SAMRAI::tbox::SAMRAI_MPI::getNodes();

    // Create the working directory.
    sprintf(temp_buf, "%06d", d_instrument_read_timestep_num);
    std::string current_dump_directory_name = SILO_DUMP_DIR_PREFIX + temp_buf;
    std::string dump_dirname = d_plot_directory_name + "/" + current_dump_directory_name;

    SAMRAI::tbox::Utilities::recursiveMkdir(dump_dirname);

    // Create one local DBfile per MPI process.
    sprintf(temp_buf, "%04d", mpi_rank);
    current_file_name = dump_dirname + "/" + SILO_PROCESSOR_FILE_PREFIX;
    current_file_name += temp_buf;
    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

    if ((dbfile = DBCreate(current_file_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB)) == NULL)
    {
        TBOX_ERROR(d_object_name + "::writePlotData():\n"
                   << "  Could not create DBfile named " << current_file_name << std::endl);
    }

    // Output the web data on the available MPI processes.
    for (int meter = 0; meter < d_num_meters; ++meter)
    {
        if (meter%mpi_nodes == mpi_rank)
        {
            std::string dirname = d_instrument_names[meter];
            if (DBMkDir(dbfile, dirname.c_str()) == -1)
            {
                TBOX_ERROR(d_object_name + "::writePlotData():\n"
                           << "  Could not create directory named "
                           << dirname << std::endl);
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
        std::string summary_file_name = dump_dirname + "/" + SILO_SUMMARY_FILE_PREFIX + temp_buf + SILO_SUMMARY_FILE_POSTFIX;
        if ((dbfile = DBCreate(summary_file_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB)) == NULL)
        {
            TBOX_ERROR(d_object_name + "::writePlotData():\n"
                       << "  Could not create DBfile named " << summary_file_name << std::endl);
        }

        int    cycle = timestep_num;
        float  time  = simulation_time;
        double dtime = simulation_time;

        static const int MAX_OPTS = 3;
        DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
        DBAddOption(optlist, DBOPT_CYCLE, &cycle);
        DBAddOption(optlist, DBOPT_TIME , &time );
        DBAddOption(optlist, DBOPT_DTIME, &dtime);

        for (int meter = 0; meter < d_num_meters; ++meter)
        {
            const int proc = meter%mpi_nodes;
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
                           << meter_name << std::endl);
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
        std::string file = current_dump_directory_name + "/" + SILO_SUMMARY_FILE_PREFIX + temp_buf + SILO_SUMMARY_FILE_POSTFIX;
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

    SAMRAI::tbox::SAMRAI_MPI::barrier();
#else
    TBOX_WARNING("IBInstrumentPanel::writePlotData(): SILO is not installed; cannot write data." << std::endl);
#endif //if HAVE_LIBSILO
    t_write_plot_data->stop();
    return;
}// writePlotData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBInstrumentPanel::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    d_plot_directory_name = db->getStringWithDefault("plot_directory_name", d_plot_directory_name);
    d_output_log_file = db->getBoolWithDefault("output_log_file", d_output_log_file);
    if (d_output_log_file)
    {
        d_log_file_name = db->getStringWithDefault("log_file_name", d_log_file_name);
    }
    d_flow_conv = db->getDoubleWithDefault("flow_conv", d_flow_conv);
    d_pres_conv = db->getDoubleWithDefault("pres_conv", d_pres_conv);
    d_flow_units = db->getStringWithDefault("flow_units", d_flow_units);
    d_pres_units = db->getStringWithDefault("pres_units", d_pres_units);
    return;
}// getFromInput

void
IBInstrumentPanel::outputLogData(
    std::ostream& os)
{
    for (int m = 0; m < d_num_meters; ++m)
    {
        std::string meter_name(d_instrument_names[m]);
        meter_name.resize(d_max_instrument_name_len,' ');
        os << meter_name;

        os.setf(std::ios_base::scientific); os.precision(5);
        os  << "  " << d_instrument_read_time;

        os.setf(std::ios_base::scientific); os.precision(5);
        os << "  " << d_X_centroid[m](0);
        os.setf(std::ios_base::scientific); os.precision(5);
        os << "  " << d_X_centroid[m](1);
#if (NDIM == 3)
        os.setf(std::ios_base::scientific); os.precision(5);
        os << "  " << d_X_centroid[m](2);
#endif

        os.setf(std::ios_base::scientific); os.setf(std::ios_base::showpos); os.precision(5);
        os << "  " << d_flow_conv*d_flow_values[m];

        os.setf(std::ios_base::scientific); os.setf(std::ios_base::showpos); os.precision(5);
        os << "  " << d_pres_conv*d_mean_pres_values[m];

        os.setf(std::ios_base::scientific); os.setf(std::ios_base::showpos); os.precision(5);
        os << "  " << d_pres_conv*d_point_pres_values[m];

        os << "\n";

        os.unsetf(std::ios_base::scientific); os.unsetf(std::ios_base::showpos);
    }
    return;
}// outputLogData

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBInstrumentPanel>;

//////////////////////////////////////////////////////////////////////////////
