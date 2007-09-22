// Filename: LagM3DDataWriter.C
// Last modified: <17.Sep.2007 19:44:23 griffith@box221.cims.nyu.edu>
// Created on 26 Apr 2005 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)

#include "LagM3DDataWriter.h"

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
#include <ibamr/IBMarker.h>

// STOOLS INCLUDES
#include <stools/PETSC_SAMRAI_ERROR.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <IndexData.h>
#include <tbox/MPI.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// The rank of the root MPI process and the MPI tag number.
static const int M3D_MPI_ROOT = 0;

// The name of the myocardial3D dumps and database filenames.
static const int M3D_NFG_MAX = 999;
static const int M3D_BUFSIZE = 128;

void
build_local_marker_cloud(
    std::ostream& os,
    const int& nmarks,
    const double* const X,
    const int& node_offset,
    const int& cloud_number)
{
    for (int k = 0; k < nmarks; ++k)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            os << std::setw(7) << std::fixed << std::setprecision(3) << X[NDIM*k+d] << " ";
        }
        os << std::setw(4) << node_offset+k+1 << " " << std::setw(2) << cloud_number+1 << "\n";
    }
    return;
}//build_local_marker_cloud

void
build_local_cart_block(
    std::ostream& os,
    const SAMRAI::hier::IntVector<NDIM>& nelem,
    const SAMRAI::hier::IntVector<NDIM>& periodic,
    const double* const X,
    const int& fiber_offset,
    const int& group_offset,
    const int& layer_number)
{
    int group_counter = 0;
    if ((nelem[0] == 1 && nelem[1] == 1) || (nelem[0] == 1 && nelem[2] == 1) || (nelem[1] == 1 && nelem[2] == 1))
    {
        // Output a single fiber.
        const int fiber_number = fiber_offset+1;

        os << std::setw(7) << fiber_number << " " << std::setw(7) << nelem.getProduct() << " = FIBER POINTS\n";
        for (int k = 0; k < nelem.getProduct(); ++k)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                os << std::setw(7) << std::fixed << std::setprecision(3) << X[NDIM*k+d] << " ";
            }
            os << std::setw(4) << k+1 << " " << std::setw(4) << fiber_number << " " << std::setw(4) << group_offset+group_counter+1 << " " << std::setw(4) << layer_number << "\n";
        }
    }
    else if ((nelem[0] == 1) || (nelem[1] == 1) || (nelem[2] == 1))
    {
        // Output a 2D sheet of fibers.
        int fiber_counter = 0;
        for (int d0 = 0; d0 < NDIM; ++d0)
        {
            if (nelem[d0] > 1)
            {
                // Find the other nontrivial dimension.
                for (int d1 = 0; d1 < NDIM; ++d1)
                {
                    if (d1 != d0 && nelem[d1] > 1)
                    {
                        for (int j = 0; j < nelem[d0]; ++j)
                        {
                            const int fiber_number = fiber_offset+fiber_counter+1;
                            fiber_counter += 1;

                            os << std::setw(7) << fiber_number << " " << std::setw(7) << nelem[d1] + (periodic[d1] ? 1 : 0) << " = FIBER POINTS\n";
                            for (int k = 0; k < nelem[d1] + (periodic[d1] ? 1 : 0); ++k)
                            {
                                bool end_of_fiber = false;
                                if (periodic[d1] && k == nelem[d1])
                                {
                                    k = 0;
                                    end_of_fiber = true;
                                }

                                int idx[NDIM] = {0 , 0 , 0};
                                idx[d0] = j;
                                idx[d1] = k;
                                const int offset = idx[0] + idx[1]*nelem[0] + idx[2]*nelem[0]*nelem[1];

                                for (int d = 0; d < NDIM; ++d)
                                {
                                    os << std::setw(7) << std::fixed << std::setprecision(3) << X[NDIM*offset+d] << " ";
                                }
                                os << std::setw(4) << k+1 << " " << std::setw(4) << fiber_number << " " << std::setw(4) << group_offset+group_counter+1 << " " << std::setw(4) << layer_number << "\n";

                                if (end_of_fiber) break;
                            }
                        }
                        group_counter += 1;
                    }
                }
            }
        }
    }
    else
    {
        // Output a 3D volume of fibers.
        int fiber_counter = 0;

        // Loop over all pairs of dimensions.
        for (int d0 = 0; d0 < NDIM; ++d0)
        {
            const int d1 = (d0+1)%NDIM;
            const int d2 = (d1+1)%NDIM;

            // myocardial3D cannot cleanly handle headers specifying groups with
            // more than 999 fibers.
            int nfibers_per_group = nelem[d0]*nelem[d1];
            while(nfibers_per_group > M3D_NFG_MAX)
            {
                nfibers_per_group /= 2;
            }

            for (int i = 0; i < nelem[d0]; ++i)
            {
                for (int j = 0; j < nelem[d1]; ++j)
                {
                    const int fiber_number = fiber_offset+fiber_counter+1;
                    fiber_counter += 1;
                    os << std::setw(7) << fiber_number << " " << std::setw(7) << nelem[d2] + (periodic[d2] ? 1 : 0) << " = FIBER POINTS\n";
                    for (int k = 0; k < nelem[d2] + (periodic[d2] ? 1 : 0); ++k)
                    {
                        bool end_of_fiber = false;
                        if (periodic[d2] && k == nelem[d2])
                        {
                            k = 0;
                            end_of_fiber = true;
                        }

                        int idx[NDIM] = {0 , 0 , 0};
                        idx[d0] = i;
                        idx[d1] = j;
                        idx[d2] = k;
                        const int offset = idx[0] + idx[1]*nelem[0] + idx[2]*nelem[0]*nelem[1];

                        for (int d = 0; d < NDIM; ++d)
                        {
                            os << std::setw(7) << std::fixed << std::setprecision(3) << X[NDIM*offset+d] << " ";
                        }
                        os << std::setw(4) << k+1 << " " << std::setw(4) << fiber_number << " " << std::setw(4) << group_offset+group_counter+1 << " " << std::setw(4) << layer_number << "\n";

                        if (end_of_fiber) break;
                    }
                }
            }
            group_counter += 1;
        }
    }
    return;
}//build_local_cart_block
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
      d_coords_data(d_finest_ln+1,SAMRAI::tbox::Pointer<LNodeLevelData>(NULL)),
      d_nvars(d_finest_ln+1,0),
      d_var_names(d_finest_ln+1),
      d_var_depths(d_finest_ln+1),
      d_var_data(d_finest_ln+1),
      d_ao(d_finest_ln+1),
      d_build_vec_scatters(d_finest_ln+1),
      d_src_vec(d_finest_ln+1),
      d_dst_vec(d_finest_ln+1),
      d_vec_scatter(d_finest_ln+1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_experiment_name.size() == 3);
    assert((0 <= d_experiment_number) && (d_experiment_number <= 9999));
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
            Vec& v = (*it).second;
            if (v)
            {
                ierr = VecDestroy(v);
                PETSC_SAMRAI_ERROR(ierr);
            }
        }
        for (std::map<int,VecScatter>::iterator it = d_vec_scatter[ln].begin();
             it != d_vec_scatter[ln].end(); ++it)
        {
            VecScatter& vs = (*it).second;
            if (vs)
            {
                ierr = VecScatterDestroy(vs);
                PETSC_SAMRAI_ERROR(ierr);
            }
        }
    }
    return;
}// ~LagM3DDataWriter

void
LagM3DDataWriter::setPatchHierarchy(
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert(hierarchy->getFinestLevelNumber() >= d_finest_ln);
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
    assert((coarsest_ln >= 0) && (finest_ln >= coarsest_ln));
    if (!d_hierarchy.isNull())
    {
        assert(finest_ln <= d_hierarchy->getFinestLevelNumber());
    }
#endif
    // Destroy any un-needed PETSc objects.
    int ierr;
    for (int ln = SAMRAI::tbox::Utilities::imax(d_coarsest_ln,0);
         (ln <= d_finest_ln) && (ln < coarsest_ln); ++ln)
    {
        for (std::map<int,Vec>::iterator it = d_dst_vec[ln].begin();
             it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = (*it).second;
            if (v != static_cast<Vec>(NULL))
            {
                ierr = VecDestroy(v);  PETSC_SAMRAI_ERROR(ierr);
            }
        }
        for (std::map<int,VecScatter>::iterator it = d_vec_scatter[ln].begin();
             it != d_vec_scatter[ln].end(); ++it)
        {
            VecScatter& vs = (*it).second;
            if (vs != static_cast<VecScatter>(NULL))
            {
                ierr = VecScatterDestroy(vs);  PETSC_SAMRAI_ERROR(ierr);
            }
        }
    }

    for (int ln = finest_ln+1; ln <= d_finest_ln; ++ln)
    {
        for (std::map<int,Vec>::iterator it = d_dst_vec[ln].begin();
             it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = (*it).second;
            if (v != static_cast<Vec>(NULL))
            {
                ierr = VecDestroy(v);  PETSC_SAMRAI_ERROR(ierr);
            }
        }
        for (std::map<int,VecScatter>::iterator it = d_vec_scatter[ln].begin();
             it != d_vec_scatter[ln].end(); ++it)
        {
            VecScatter& vs = (*it).second;
            if (vs != static_cast<VecScatter>(NULL))
            {
                ierr = VecScatterDestroy(vs);  PETSC_SAMRAI_ERROR(ierr);
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
    d_nvars      .resize(d_finest_ln+1,0);
    d_var_names  .resize(d_finest_ln+1);
    d_var_depths .resize(d_finest_ln+1);
    d_var_data   .resize(d_finest_ln+1);

    d_ao                .resize(d_finest_ln+1);
    d_build_vec_scatters.resize(d_finest_ln+1);
    d_src_vec           .resize(d_finest_ln+1);
    d_dst_vec           .resize(d_finest_ln+1);
    d_vec_scatter       .resize(d_finest_ln+1);

    return;
}// resetLevels

void
LagM3DDataWriter::registerIBMarkerPatchDataIndex(
    const int mark_idx)
{
    d_mark_idx = mark_idx;
    return;
}// registerIBMarkerPatchDataIndex

void
LagM3DDataWriter::registerMarkerCloud(
    const std::string& name,
    const int nmarks,
    const int first_mark_idx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(nmarks > 0);
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
    const SAMRAI::hier::IntVector<NDIM>& nelem,
    const SAMRAI::hier::IntVector<NDIM>& periodic,
    const int first_lag_idx,
    const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(min(level_number,d_coarsest_ln),max(level_number,d_finest_ln));
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 0; d < NDIM; ++d)
    {
        assert(nelem[d] > 0);
        assert(periodic(d) == 0 || periodic(d) == 1);
    }
    assert(d_coarsest_ln <= level_number &&
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
    SAMRAI::tbox::Pointer<LNodeLevelData> coords_data,
    const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(min(level_number,d_coarsest_ln),max(level_number,d_finest_ln));
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!coords_data.isNull());
    assert(coords_data->getDepth() == NDIM);
    assert(d_coarsest_ln <= level_number &&
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
        resetLevels(min(level_number,d_coarsest_ln),max(level_number,d_finest_ln));
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_coarsest_ln <= level_number &&
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
    assert(coarsest_ln <= finest_ln);
#endif

    if (coarsest_ln < d_coarsest_ln || finest_ln > d_finest_ln)
    {
        resetLevels(min(coarsest_ln,d_coarsest_ln),max(finest_ln,d_finest_ln));
    }

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_coarsest_ln <= coarsest_ln && finest_ln <= d_finest_ln);
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
    assert(time_step_number >= 0);
    assert(!d_dump_directory_name.empty());
#endif

    if (time_step_number <= d_time_step_number)
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                   << "  data writer with name " << d_object_name << "\n"
                   << "  time step number: " << time_step_number
                   << " is <= last time step number: " << d_time_step_number
                   << endl);
    }
    d_time_step_number = time_step_number;

    if (d_dump_directory_name.empty())
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                   << "  data writer with name " << d_object_name << "\n"
                   << "  dump directory name is empty" << endl);
    }

    int ierr;
    const int mpi_rank = SAMRAI::tbox::MPI::getRank();
    const int mpi_size = SAMRAI::tbox::MPI::getNodes();

    // Determine the marker cloud and fiber offsets.
    int num_local_marker_nodes = std::accumulate(d_cloud_nmarks.begin(), d_cloud_nmarks.end(), 0);
    std::vector<int> num_marker_nodes_proc(mpi_size,0);
    SAMRAI::tbox::MPI::allGather(num_local_marker_nodes, &num_marker_nodes_proc[0]);
    const int marker_node_offset = std::accumulate(
        num_marker_nodes_proc.begin(), num_marker_nodes_proc.begin()+mpi_rank, 0);

    const int num_local_marker_clouds = d_nclouds;
    std::vector<int> num_marker_clouds_proc(mpi_size,0);
    SAMRAI::tbox::MPI::allGather(num_local_marker_clouds, &num_marker_clouds_proc[0]);
    const int marker_cloud_offset = std::accumulate(
        num_marker_clouds_proc.begin(), num_marker_clouds_proc.begin()+mpi_rank, 0);
    const int num_marker_clouds = std::accumulate(
        num_marker_clouds_proc.begin()+mpi_rank, num_marker_clouds_proc.end(), marker_cloud_offset);

    int local_marker_node_counter, local_marker_cloud_counter;

    int num_local_fibers = 0;
    for (std::vector<std::vector<int> >::const_iterator cit = d_block_nfibers.begin();
         cit != d_block_nfibers.end(); ++cit)
    {
        num_local_fibers += std::accumulate((*cit).begin(), (*cit).end(), 0);
    }
    std::vector<int> num_fibers_proc(mpi_size,0);
    SAMRAI::tbox::MPI::allGather(num_local_fibers, &num_fibers_proc[0]);
    const int fiber_offset = std::accumulate(
        num_fibers_proc.begin(), num_fibers_proc.begin()+mpi_rank, 0);
    const int num_fibers = std::accumulate(
        num_fibers_proc.begin()+mpi_rank, num_fibers_proc.end(), fiber_offset);

    int num_local_groups = 0;
    for (std::vector<std::vector<int> >::const_iterator cit = d_block_ngroups.begin();
         cit != d_block_ngroups.end(); ++cit)
    {
        num_local_groups += std::accumulate((*cit).begin(), (*cit).end(), 0);
    }
    std::vector<int> num_groups_proc(mpi_size,0);
    SAMRAI::tbox::MPI::allGather(num_local_groups, &num_groups_proc[0]);
    const int group_offset = std::accumulate(
        num_groups_proc.begin(), num_groups_proc.begin()+mpi_rank, 0);

    const int num_local_layers = std::accumulate(d_nblocks.begin(), d_nblocks.end(), 0);
    std::vector<int> num_layers_proc(mpi_size,0);
    SAMRAI::tbox::MPI::allGather(num_local_layers, &num_layers_proc[0]);
    const int layer_offset = std::accumulate(
        num_layers_proc.begin(), num_layers_proc.begin()+mpi_rank, 0);
    const int num_layers = std::accumulate(
        num_layers_proc.begin()+mpi_rank, num_layers_proc.end(), layer_offset);

    int local_fiber_counter, local_group_counter, local_layer_counter;

    // Construct the VecScatter objectss required to write the plot data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (d_build_vec_scatters[ln])
        {
            buildVecScatters(d_ao[ln], ln);
        }
        d_build_vec_scatters[ln] = false;
    }

    // Create the dump directory.
    SAMRAI::tbox::Utilities::recursiveMkdir(d_dump_directory_name);

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

    // Local data streams.
    std::ofstream local_marker_stream(std::string(d_dump_directory_name + "/" + local_marker_file_name).c_str(), std::ios::out);
    std::ofstream local_fiber_stream(std::string(d_dump_directory_name + "/" + local_fiber_file_name).c_str(), std::ios::out);

    // Gather marker data.
    int max_marker_idx = 0;
    for (int cloud = 0; cloud < d_nclouds; ++cloud)
    {
        max_marker_idx = std::max(max_marker_idx,d_cloud_first_mark_idx[cloud]+d_cloud_nmarks[cloud]-1);
    }
    SAMRAI::tbox::MPI::maxReduction(&max_marker_idx);

    const int total_num_marks = max_marker_idx+1;
    std::vector<double> X_marker(NDIM*total_num_marks,0.0);
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            if (d_mark_idx != -1)
            {
                SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<NDIM,IBMarker> > mark_data =
                    patch->getPatchData(d_mark_idx);
                for (SAMRAI::pdat::IndexData<NDIM,IBMarker>::Iterator it(*mark_data); it; it++)
                {
                    const IBMarker& mark = it();
                    const std::vector<double>& X = mark.getPositions();
                    const std::vector<int>& idx = mark.getIndices();
                    for (size_t k = 0; k < X.size()/NDIM; ++k)
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X_marker[NDIM*idx[k]+d] = X[NDIM*k+d];  // XXXX add index error checking!
                        }
                    }
                }
            }
        }
    }
    SAMRAI::tbox::MPI::sumReduction(&X_marker[0], NDIM*max_marker_idx);

    // Add the local clouds to the local marker file.
    local_marker_node_counter = 0;
    local_marker_cloud_counter = 0;
    for (int cloud = 0; cloud < d_nclouds; ++cloud)
    {
        const int nmarks = d_cloud_nmarks[cloud];
        const int first_mark_idx = d_cloud_first_mark_idx[cloud];
        const double* const X = &X_marker[NDIM*first_mark_idx];
        build_local_marker_cloud(local_marker_stream, nmarks, X,
                                 marker_node_offset+local_marker_node_counter,
                                 marker_cloud_offset+local_marker_cloud_counter);
        local_marker_node_counter += nmarks;
        local_marker_cloud_counter += 1;
    }

    // Add the local fibers to the local fiber file.
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
            PETSC_SAMRAI_ERROR(ierr);

            Vec& global_X_vec = d_coords_data[ln]->getGlobalVec();
            ierr = VecScatterBegin(d_vec_scatter[ln][NDIM], global_X_vec, local_X_vec,
                                   INSERT_VALUES, SCATTER_FORWARD);
            PETSC_SAMRAI_ERROR(ierr);
            ierr = VecScatterEnd(d_vec_scatter[ln][NDIM], global_X_vec, local_X_vec,
                                 INSERT_VALUES, SCATTER_FORWARD);
            PETSC_SAMRAI_ERROR(ierr);

            double* local_X_arr;
            ierr = VecGetArray(local_X_vec, &local_X_arr);
            PETSC_SAMRAI_ERROR(ierr);

            std::vector<Vec> local_v_vecs;
            std::vector<double*> local_v_arrs;

            for (int v = 0; v < d_nvars[ln]; ++v)
            {
                const int var_depth = d_var_depths[ln][v];
                Vec local_v_vec;
                ierr = VecDuplicate(d_dst_vec[ln][var_depth], &local_v_vec);
                PETSC_SAMRAI_ERROR(ierr);

                Vec& global_v_vec = d_var_data[ln][v]->getGlobalVec();
                ierr = VecScatterBegin(d_vec_scatter[ln][var_depth], global_v_vec, local_v_vec,
                                       INSERT_VALUES, SCATTER_FORWARD);
                PETSC_SAMRAI_ERROR(ierr);
                ierr = VecScatterEnd(d_vec_scatter[ln][var_depth], global_v_vec, local_v_vec,
                                     INSERT_VALUES, SCATTER_FORWARD);
                PETSC_SAMRAI_ERROR(ierr);

                double* local_v_arr;
                ierr = VecGetArray(local_v_vec, &local_v_arr);
                PETSC_SAMRAI_ERROR(ierr);

                local_v_vecs.push_back(local_v_vec);
                local_v_arrs.push_back(local_v_arr);
            }

            // Keep track of the current offset in the local Vec data.
            int offset = 0;

            // Add the local blocks to the local fiber file.
            for (int block = 0; block < d_nblocks[ln]; ++block)
            {
                const SAMRAI::hier::IntVector<NDIM>& nelem    = d_block_nelems  [ln][block];
                const SAMRAI::hier::IntVector<NDIM>& periodic = d_block_periodic[ln][block];
                const int ntot = nelem.getProduct();
                const double* const X = local_X_arr + NDIM*offset;
                const int layer_number = layer_offset+local_layer_counter+1;
                build_local_cart_block(local_fiber_stream, nelem, periodic, X,
                                       fiber_offset+local_fiber_counter,
                                       group_offset+local_group_counter, layer_number);
                offset += ntot;
                local_fiber_counter += d_block_nfibers[ln][block];
                local_group_counter += d_block_ngroups[ln][block];
                local_layer_counter += 1;
            }

            // Clean up allocated data.
            ierr = VecRestoreArray(local_X_vec, &local_X_arr);
            PETSC_SAMRAI_ERROR(ierr);
            ierr = VecDestroy(local_X_vec);
            PETSC_SAMRAI_ERROR(ierr);
            for (int v = 0; v < d_nvars[ln]; ++v)
            {
                ierr = VecRestoreArray(local_v_vecs[v], &local_v_arrs[v]);
                PETSC_SAMRAI_ERROR(ierr);
                ierr = VecDestroy(local_v_vecs[v]);
                PETSC_SAMRAI_ERROR(ierr);
            }
        }
    }

    // Determine the length of the computational domain.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const grid_xLower = grid_geom->getXLower();
    const double* const grid_xUpper = grid_geom->getXUpper();
    const double L = std::max(std::max(grid_xUpper[0]-grid_xLower[0], grid_xUpper[1]-grid_xLower[1]), grid_xUpper[2]-grid_xLower[2]);

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
                    const SAMRAI::hier::IntVector<NDIM>& nelem = d_block_nelems[ln][block];
                    const SAMRAI::hier::IntVector<NDIM>& periodic = d_block_periodic[ln][block];
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
        SAMRAI::tbox::MPI::barrier();
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
            SAMRAI::tbox::MPI::barrier();
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
            SAMRAI::tbox::MPI::barrier();
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
        SAMRAI::tbox::MPI::barrier();
    }

    // Create list file and cat script.
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

        // Setup cat commands to create a unified marker file.
        cat_file_stream << "cat " << marker_header_file_name;
        for (int rank = 0; rank < mpi_size; ++rank)
        {
            std::ostringstream stream;
            stream << std::setfill('0') << std::setw(4) << mpi_rank;
            const std::string mpi_rank_string = stream.str();
            cat_file_stream << " " << marker_file_name + "." + mpi_rank_string;
        }
        cat_file_stream << " > " << marker_file_name << "\n";

        // Setup cat commands to create a unified fiber file.
        cat_file_stream << "cat " << fiber_header_file_name;
        for (int rank = 0; rank < mpi_size; ++rank)
        {
            std::ostringstream stream;
            stream << std::setfill('0') << std::setw(4) << mpi_rank;
            const std::string mpi_rank_string = stream.str();
            cat_file_stream << " " << fiber_file_name + "." + mpi_rank_string;
        }
        cat_file_stream << " > " << fiber_file_name << "\n";
    }
    SAMRAI::tbox::MPI::barrier();
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
        const SAMRAI::hier::IntVector<NDIM>& nelem = d_block_nelems[level_number][block];
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
        (!ref_is_idxs.empty() ? static_cast<int>(ref_is_idxs.size()) : static_cast<int>(ao_dummy.size())),
        (!ref_is_idxs.empty() ? &ref_is_idxs[0]                      : &ao_dummy[0]));
    PETSC_SAMRAI_ERROR(ierr);

    // Setup IS indices for all necessary data depths.
    std::map<int,std::vector<int> > src_is_idxs;

    src_is_idxs[NDIM].resize(ref_is_idxs.size());
    std::transform(ref_is_idxs.begin(), ref_is_idxs.end(),
                   src_is_idxs[NDIM].begin(),
                   std::bind2nd(std::multiplies<int>(),NDIM));
    d_src_vec[level_number][NDIM] = d_coords_data[level_number]->getGlobalVec();

    for (int v = 0; v < d_nvars[level_number]; ++v)
    {
        const int var_depth = d_var_depths[level_number][v];
        if (src_is_idxs.find(var_depth) == src_is_idxs.end())
        {
            src_is_idxs[var_depth].resize(ref_is_idxs.size());
            std::transform(ref_is_idxs.begin(), ref_is_idxs.end(),
                           src_is_idxs[var_depth].begin(),
                           std::bind2nd(std::multiplies<int>(),var_depth));
            d_src_vec[level_number][var_depth] = d_var_data[level_number][v]->getGlobalVec();
        }
    }

    // Create the VecScatters to scatter data from the global PETSc Vec to
    // contiguous local subgrids.  VecScatter objects are individually created
    // for data depths as necessary.
    for (std::map<int,std::vector<int> >::iterator it = src_is_idxs.begin();
         it != src_is_idxs.end(); ++it)
    {
        const int depth = (*it).first;
        const std::vector<int>& idxs = (*it).second;

        IS src_is;
        ierr = ISCreateBlock(PETSC_COMM_WORLD, depth, idxs.size(),
                             &idxs[0], &src_is);
        PETSC_SAMRAI_ERROR(ierr);

        Vec& src_vec = d_src_vec[level_number][depth];
        Vec& dst_vec = d_dst_vec[level_number][depth];
        if (dst_vec)
        {
            ierr = VecDestroy(dst_vec);
            PETSC_SAMRAI_ERROR(ierr);
        }
        ierr = VecCreateMPI(PETSC_COMM_WORLD, depth*idxs.size(),
                            PETSC_DETERMINE, &dst_vec);
        PETSC_SAMRAI_ERROR(ierr);

        ierr = VecSetBlockSize(dst_vec, depth);
        PETSC_SAMRAI_ERROR(ierr);

        VecScatter& vec_scatter = d_vec_scatter[level_number][depth];
        if (vec_scatter)
        {
            ierr = VecScatterDestroy(vec_scatter);
            PETSC_SAMRAI_ERROR(ierr);
        }
        ierr = VecScatterCreate(src_vec, src_is, dst_vec, PETSC_NULL,
                                &vec_scatter);
        PETSC_SAMRAI_ERROR(ierr);

        ierr = ISDestroy(src_is);  PETSC_SAMRAI_ERROR(ierr);
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

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LagM3DDataWriter>;

//////////////////////////////////////////////////////////////////////////////
