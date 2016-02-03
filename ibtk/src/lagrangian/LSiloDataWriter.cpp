// Filename: LSiloDataWriter.cpp
// Created on 26 Apr 2005 by Boyce Griffith
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "IBTK_config.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "mpi.h"
#include "petscao.h"
#include "petscis.h"
#include "petsclog.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"
// IWYU pragma: no_include "petsc-private/vecimpl.h"

#if defined(IBTK_HAVE_SILO)
#include "silo.h"
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// The rank of the root MPI process and the MPI tag number.
static const int SILO_MPI_ROOT = 0;
static const int SILO_MPI_TAG = 0;

// The name of the Silo dumps and database filenames.
static const int SILO_NAME_BUFSIZE = 128;
static const std::string VISIT_DUMPS_FILENAME = "lag_data.visit";
static const std::string SILO_DUMP_DIR_PREFIX = "lag_data.cycle_";
static const std::string SILO_SUMMARY_FILE_PREFIX = "lag_data.cycle_";
static const std::string SILO_SUMMARY_FILE_POSTFIX = ".summary.silo";
static const std::string SILO_PROCESSOR_FILE_PREFIX = "lag_data.proc_";
static const std::string SILO_PROCESSOR_FILE_POSTFIX = ".silo";

// Version of LSiloDataWriter restart file data.
static const int LAG_SILO_DATA_WRITER_VERSION = 1;

#if defined(IBTK_HAVE_SILO)
/*!
 * \brief Build a local mesh database entry corresponding to a cloud of marker
 * points.
 */
void
build_local_marker_cloud(DBfile* dbfile,
                         std::string& dirname,
                         const int nmarks,
                         const double* const X,
                         const int nvars,
                         const std::vector<std::string>& varnames,
                         const std::vector<int>& varstartdepths,
                         const std::vector<int>& varplotdepths,
                         const std::vector<int>& vardepths,
                         const std::vector<const double*> varvals,
                         const int time_step,
                         const double simulation_time)
{
    std::vector<float> block_X(NDIM * nmarks);
    std::vector<std::vector<float> > block_varvals(nvars);
    for (int v = 0; v < nvars; ++v)
    {
        const int varplotdepth = varplotdepths[v];
        block_varvals[v].resize(varplotdepth * nmarks);
    }

    for (int i = 0; i < nmarks; ++i)
    {
        // Get the coordinate data.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            block_X[d * nmarks + i] = float(X[NDIM * i + d]);
        }

        // Get the variable data.
        for (int v = 0; v < nvars; ++v)
        {
            const int varstartdepth = varstartdepths[v];
            const int varplotdepth = varplotdepths[v];
            const int vardepth = vardepths[v];
            for (int d = 0; d < varplotdepth; ++d)
            {
                block_varvals[v][d * nmarks + i] = float(varvals[v][vardepth * i + varstartdepth + d]);
            }
        }
    }

    // Set the working directory in the Silo database.
    if (DBSetDir(dbfile, dirname.c_str()) == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_marker_cloud()\n"
                   << "  Could not set directory "
                   << dirname
                   << std::endl);
    }

    // Write out the variables.
    int cycle = time_step;
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
        coords[d] = nmarks > 0 ? &block_X[d * nmarks] : NULL;
    }

    int ndims = NDIM;

    DBPutPointmesh(dbfile, meshname, ndims, &coords[0], nmarks, DB_FLOAT, optlist);

    for (int v = 0; v < nvars; ++v)
    {
        const char* varname = varnames[v].c_str();
        const int varplotdepth = varplotdepths[v];

        std::vector<float*> vars(varplotdepth);
        for (int d = 0; d < varplotdepth; ++d)
        {
            vars[d] = nmarks > 0 ? &block_varvals[v][d * nmarks] : NULL;
        }

        if (varplotdepth == 1)
        {
            DBPutPointvar1(dbfile, varname, meshname, vars[0], nmarks, DB_FLOAT, optlist);
        }
        else
        {
            DBPutPointvar(dbfile, varname, meshname, varplotdepth, &vars[0], nmarks, DB_FLOAT, optlist);
        }
    }

    DBFreeOptlist(optlist);

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_marker_cloud()\n"
                   << "  Could not return to the base directory from subdirectory "
                   << dirname
                   << std::endl);
    }
    return;
} // build_local_marker_cloud

/*!
 * \brief Build a local mesh database entry corresponding to a quadrilateral
 * curvilinear block.
 */
void
build_local_curv_block(DBfile* dbfile,
                       std::string& dirname,
                       const IntVector<NDIM>& nelem_in,
                       const IntVector<NDIM>& periodic,
                       const double* const X,
                       const int nvars,
                       const std::vector<std::string>& varnames,
                       const std::vector<int>& varstartdepths,
                       const std::vector<int>& varplotdepths,
                       const std::vector<int>& vardepths,
                       const std::vector<const double*> varvals,
                       const int time_step,
                       const double simulation_time)
{
    // Check for co-dimension 1 or 2 data.
    IntVector<NDIM> nelem, degenerate;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (nelem_in(d) == 1)
        {
            nelem(d) = 2;
            degenerate(d) = 1;
        }
        else
        {
            nelem(d) = nelem_in(d);
            degenerate(d) = 0;
        }
    }

    // Rearrange the data into the format required by Silo.
    const int ntot = 1 * (periodic(0) ? nelem(0) + 1 : nelem(0))
#if (NDIM > 1)
                     * (periodic(1) ? nelem(1) + 1 : nelem(1))
#if (NDIM > 2)
                     * (periodic(2) ? nelem(2) + 1 : nelem(2))
#endif
#endif
        ;

    std::vector<float> block_X(NDIM * ntot);
    std::vector<std::vector<float> > block_varvals(nvars);
    for (int v = 0; v < nvars; ++v)
    {
        const int varplotdepth = varplotdepths[v];
        block_varvals[v].resize(varplotdepth * ntot);
    }

    int offset = 0;
#if (NDIM > 2)
    for (int k = 0; k < nelem(2) + (periodic(2) ? 1 : 0); ++k)
    {
#endif
#if (NDIM > 1)
        for (int j = 0; j < nelem(1) + (periodic(1) ? 1 : 0); ++j)
        {
#endif
            for (int i = 0; i < nelem(0) + (periodic(0) ? 1 : 0); ++i)
            {
                const int idx = +(degenerate(0) ? 0 : (i % nelem(0)))
#if (NDIM > 1)
                                + (degenerate(1) ? 0 : (j % nelem(1)) * nelem(0))
#if (NDIM > 2)
                                + (degenerate(2) ? 0 : (k % nelem(2)) * nelem(1) * nelem(0))
#endif
#endif
                    ;

                // Get the coordinate data.
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    block_X[d * ntot + offset] = float(X[NDIM * idx + d]);
                }

                // Get the variable data.
                for (int v = 0; v < nvars; ++v)
                {
                    const int varstartdepth = varstartdepths[v];
                    const int varplotdepth = varplotdepths[v];
                    const int vardepth = vardepths[v];
                    for (int d = 0; d < varplotdepth; ++d)
                    {
                        block_varvals[v][d * ntot + offset] = float(varvals[v][vardepth * idx + varstartdepth + d]);
                    }
                }

                // Increment the counter.
                ++offset;
            }
#if (NDIM > 1)
        }
#endif
#if (NDIM > 2)
    }
#endif

    // Set the working directory in the Silo database.
    if (DBSetDir(dbfile, dirname.c_str()) == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_curv_block()\n"
                   << "  Could not set directory "
                   << dirname
                   << std::endl);
    }

    // Write out the variables.
    int cycle = time_step;
    float time = float(simulation_time);
    double dtime = simulation_time;

    static const int MAX_OPTS = 3;
    DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    DBAddOption(optlist, DBOPT_TIME, &time);
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

    const char* meshname = "mesh";
    const char* coordnames[3] = { "xcoords", "ycoords", "zcoords" };
    std::vector<float*> coords(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        coords[d] = ntot > 0 ? &block_X[d * ntot] : 0;
    }

    int ndims = NDIM;
    std::vector<int> dims(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        dims[d] = nelem(d) + (periodic(d) ? 1 : 0);
    }

    DBPutQuadmesh(dbfile,
                  meshname,
                  const_cast<char**>(coordnames),
                  &coords[0],
                  &dims[0],
                  ndims,
                  DB_FLOAT,
                  DB_NONCOLLINEAR,
                  optlist);

    for (int v = 0; v < nvars; ++v)
    {
        const char* varname = varnames[v].c_str();
        const int varplotdepth = varplotdepths[v];
        std::vector<char*> compnames(varplotdepth);
        for (int d = 0; d < varplotdepth; ++d)
        {
            std::ostringstream stream;
            stream << "_" << d;
            const std::string compname = varnames[v] + stream.str();
            compnames[d] = strdup(compname.c_str());
        }

        std::vector<float*> vars(varplotdepth);
        for (int d = 0; d < varplotdepth; ++d)
        {
            vars[d] = ntot > 0 ? &block_varvals[v][d * ntot] : NULL;
        }

        if (varplotdepth == 1)
        {
            DBPutQuadvar1(dbfile, varname, meshname, vars[0], &dims[0], ndims, NULL, 0, DB_FLOAT, DB_NODECENT, optlist);
        }
        else
        {
            DBPutQuadvar(dbfile,
                         varname,
                         meshname,
                         varplotdepth,
                         &compnames[0],
                         &vars[0],
                         &dims[0],
                         ndims,
                         NULL,
                         0,
                         DB_FLOAT,
                         DB_NODECENT,
                         optlist);
        }

        for (int d = 0; d < varplotdepth; ++d)
        {
            free(compnames[d]);
        }
    }

    DBFreeOptlist(optlist);

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_curv_block()\n"
                   << "  Could not return to the base directory from subdirectory "
                   << dirname
                   << std::endl);
    }
    return;
} // build_local_curv_block

/*!
 * \brief Build a local mesh database entry corresponding to an unstructured
 * mesh.
 */
void
build_local_ucd_mesh(DBfile* dbfile,
                     std::string& dirname,
                     const std::set<int>& vertices,
                     const std::multimap<int, std::pair<int, int> >& edge_map,
                     const double* const X,
                     const int nvars,
                     const std::vector<std::string>& varnames,
                     const std::vector<int>& varstartdepths,
                     const std::vector<int>& varplotdepths,
                     const std::vector<int>& vardepths,
                     const std::vector<const double*> varvals,
                     const int time_step,
                     const double simulation_time)
{
    // Rearrange the data into the format required by Silo.
    const int ntot = static_cast<int>(vertices.size());

    std::vector<float> block_X(NDIM * ntot);
    std::vector<std::vector<float> > block_varvals(nvars);
    for (int v = 0; v < nvars; ++v)
    {
        const int varplotdepth = varplotdepths[v];
        block_varvals[v].resize(varplotdepth * ntot);
    }

    int offset = 0;
    std::map<int, int> local_vertex_map;
    for (std::set<int>::const_iterator it = vertices.begin(); it != vertices.end(); ++it)
    {
        const int idx = (*it);
        local_vertex_map[idx] = offset;

        // Get the coordinate data.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            block_X[d * ntot + offset] = float(X[NDIM * offset + d]);
        }

        // Get the variable data.
        for (int v = 0; v < nvars; ++v)
        {
            const int varstartdepth = varstartdepths[v];
            const int varplotdepth = varplotdepths[v];
            const int vardepth = vardepths[v];
            for (int d = 0; d < varplotdepth; ++d)
            {
                block_varvals[v][d * ntot + offset] = float(varvals[v][vardepth * offset + varstartdepth + d]);
            }
        }

        // Increment the counter.
        ++offset;
    }

    // Prune duplicate edges.
    std::set<std::pair<int, int> > local_edge_set;
    for (std::multimap<int, std::pair<int, int> >::const_iterator it = edge_map.begin(); it != edge_map.end(); ++it)
    {
        std::pair<int, int> e = it->second;
#if !defined(NDEBUG)
        TBOX_ASSERT(vertices.count(e.first) == 1);
        TBOX_ASSERT(vertices.count(e.second) == 1);
#endif
        if (e.first > e.second)
        {
            std::swap<int>(e.first, e.second);
        }
        local_edge_set.insert(e);
    }

    // Create an edge map corresponding to the pruned edge list.
    std::multimap<int, int> local_edge_map;
    for (std::set<std::pair<int, int> >::const_iterator it = local_edge_set.begin(); it != local_edge_set.end(); ++it)
    {
        const int e1 = it->first;
        const int e2 = it->second;
        local_edge_map.insert(std::make_pair(local_vertex_map[e1], local_vertex_map[e2]));
    }

    // Set the working directory in the Silo database.
    if (DBSetDir(dbfile, dirname.c_str()) == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_ucd_mesh()\n"
                   << "  Could not set directory "
                   << dirname
                   << std::endl);
    }

    // Node coordinates.
    int ndims = NDIM;

    int cycle = time_step;
    float time = float(simulation_time);
    double dtime = simulation_time;

    static const int MAX_OPTS = 3;
    DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    DBAddOption(optlist, DBOPT_TIME, &time);
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

    const char* meshname = "mesh";
    const char* coordnames[3] = { "xcoords", "ycoords", "zcoords" };
    std::vector<float*> coords(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        coords[d] = ntot > 0 ? &block_X[d * ntot] : NULL;
    }
    const int nnodes = ntot;

    // Connectivity.
    std::vector<int> nodelist;
    nodelist.reserve(2 * local_edge_map.size());

    for (std::multimap<int, int>::const_iterator it = local_edge_map.begin(); it != local_edge_map.end(); ++it)
    {
        nodelist.push_back(it->first);
        nodelist.push_back(it->second);
    }
    int lnodelist = static_cast<int>(nodelist.size());
    int nshapetypes = 1;
    int shapecnt[] = { static_cast<int>(local_edge_map.size()) };
    int shapesize[] = { 2 };
    int shapetype[] = { DB_ZONETYPE_BEAM };
    int nzones = static_cast<int>(local_edge_map.size());

    // Write out connectivity information.
    const int origin = 0;
    const int lo_offset = 0;
    const int hi_offset = 0;

    // Write out connectivity information.
    DBPutZonelist2(dbfile,
                   "zonelist",
                   nzones,
                   ndims,
                   lnodelist > 0 ? &nodelist[0] : NULL,
                   lnodelist,
                   origin,
                   lo_offset,
                   hi_offset,
                   shapetype,
                   shapesize,
                   shapecnt,
                   nshapetypes,
                   optlist);

    // Write an unstructured mesh.
    DBPutUcdmesh(dbfile,
                 meshname,
                 ndims,
                 const_cast<char**>(coordnames),
                 &coords[0],
                 nnodes,
                 nzones,
                 "zonelist",
                 NULL,
                 DB_FLOAT,
                 NULL);

    // Write the variables defined on the unstructured mesh.
    for (int v = 0; v < nvars; ++v)
    {
        const char* varname = varnames[v].c_str();
        const int varplotdepth = varplotdepths[v];
        std::vector<char*> compnames(varplotdepth);
        for (int d = 0; d < varplotdepth; ++d)
        {
            std::ostringstream stream;
            stream << "_" << d;
            const std::string compname = varnames[v] + stream.str();
            compnames[d] = strdup(compname.c_str());
        }

        std::vector<float*> vars(varplotdepth);
        for (int d = 0; d < varplotdepth; ++d)
        {
            vars[d] = ntot > 0 ? &block_varvals[v][d * ntot] : NULL;
        }

        if (varplotdepth == 1)
        {
            DBPutUcdvar1(dbfile, varname, meshname, vars[0], nnodes, NULL, 0, DB_FLOAT, DB_NODECENT, optlist);
        }
        else
        {
            DBPutUcdvar(dbfile,
                        varname,
                        meshname,
                        varplotdepth,
                        &compnames[0],
                        &vars[0],
                        nnodes,
                        NULL,
                        0,
                        DB_FLOAT,
                        DB_NODECENT,
                        optlist);
        }

        for (int d = 0; d < varplotdepth; ++d)
        {
            free(compnames[d]);
        }
    }

    DBFreeOptlist(optlist);

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        TBOX_ERROR("LSiloDataWriter::build_local_ucd_mesh()\n"
                   << "  Could not return to the base directory from subdirectory "
                   << dirname
                   << std::endl);
    }
    return;
} // build_local_ucd_mesh
#endif // if defined(IBTK_HAVE_SILO)
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

LSiloDataWriter::LSiloDataWriter(const std::string& object_name,
                                 const std::string& dump_directory_name,
                                 bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_dump_directory_name(dump_directory_name),
      d_time_step_number(-1),
      d_hierarchy(),
      d_coarsest_ln(0),
      d_finest_ln(0),
      d_nclouds(d_finest_ln + 1, 0),
      d_cloud_names(d_finest_ln + 1),
      d_cloud_nmarks(d_finest_ln + 1),
      d_cloud_first_lag_idx(d_finest_ln + 1),
      d_nblocks(d_finest_ln + 1, 0),
      d_block_names(d_finest_ln + 1),
      d_block_nelems(d_finest_ln + 1),
      d_block_periodic(d_finest_ln + 1),
      d_block_first_lag_idx(d_finest_ln + 1),
      d_nmbs(d_finest_ln + 1, 0),
      d_mb_names(d_finest_ln + 1),
      d_mb_nblocks(d_finest_ln + 1),
      d_mb_nelems(d_finest_ln + 1),
      d_mb_periodic(d_finest_ln + 1),
      d_mb_first_lag_idx(d_finest_ln + 1),
      d_nucd_meshes(d_finest_ln + 1, 0),
      d_ucd_mesh_names(d_finest_ln + 1),
      d_ucd_mesh_vertices(d_finest_ln + 1),
      d_ucd_mesh_edge_maps(d_finest_ln + 1),
      d_coords_data(d_finest_ln + 1, Pointer<LData>(NULL)),
      d_nvars(d_finest_ln + 1, 0),
      d_var_names(d_finest_ln + 1),
      d_var_start_depths(d_finest_ln + 1),
      d_var_plot_depths(d_finest_ln + 1),
      d_var_depths(d_finest_ln + 1),
      d_var_data(d_finest_ln + 1),
      d_ao(d_finest_ln + 1),
      d_build_vec_scatters(d_finest_ln + 1),
      d_src_vec(d_finest_ln + 1),
      d_dst_vec(d_finest_ln + 1),
      d_vec_scatter(d_finest_ln + 1)
{
#if defined(IBTK_HAVE_SILO)
// intentionally blank
#else
    TBOX_WARNING("LSiloDataWriter::LSiloDataWriter(): SILO is not installed; cannot write data." << std::endl);
#endif
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the restart database.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    return;
} // LSiloDataWriter

LSiloDataWriter::~LSiloDataWriter()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    // Destroy any remaining PETSc objects.
    int ierr;
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        for (std::map<int, Vec>::iterator it = d_dst_vec[ln].begin(); it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = it->second;
            if (v)
            {
                ierr = VecDestroy(&v);
                IBTK_CHKERRQ(ierr);
            }
        }
        for (std::map<int, VecScatter>::iterator it = d_vec_scatter[ln].begin(); it != d_vec_scatter[ln].end(); ++it)
        {
            VecScatter& vs = it->second;
            if (vs)
            {
                ierr = VecScatterDestroy(&vs);
                IBTK_CHKERRQ(ierr);
            }
        }
    }
    return;
} // ~LSiloDataWriter

void
LSiloDataWriter::setPatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(hierarchy->getFinestLevelNumber() >= d_finest_ln);
#endif
    // Reset the hierarchy.
    d_hierarchy = hierarchy;
    return;
} // setPatchHierarchy

void
LSiloDataWriter::resetLevels(const int coarsest_ln, const int finest_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT((coarsest_ln >= 0) && (finest_ln >= coarsest_ln));
    if (d_hierarchy)
    {
        TBOX_ASSERT(finest_ln <= d_hierarchy->getFinestLevelNumber());
    }
#endif
    // Destroy any unneeded PETSc objects.
    int ierr;
    for (int ln = std::max(d_coarsest_ln, 0); (ln <= d_finest_ln) && (ln < coarsest_ln); ++ln)
    {
        for (std::map<int, Vec>::iterator it = d_dst_vec[ln].begin(); it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = it->second;
            if (v)
            {
                ierr = VecDestroy(&v);
                IBTK_CHKERRQ(ierr);
            }
        }
        for (std::map<int, VecScatter>::iterator it = d_vec_scatter[ln].begin(); it != d_vec_scatter[ln].end(); ++it)
        {
            VecScatter& vs = it->second;
            if (vs)
            {
                ierr = VecScatterDestroy(&vs);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    for (int ln = finest_ln + 1; ln <= d_finest_ln; ++ln)
    {
        for (std::map<int, Vec>::iterator it = d_dst_vec[ln].begin(); it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = it->second;
            if (v)
            {
                ierr = VecDestroy(&v);
                IBTK_CHKERRQ(ierr);
            }
        }
        for (std::map<int, VecScatter>::iterator it = d_vec_scatter[ln].begin(); it != d_vec_scatter[ln].end(); ++it)
        {
            VecScatter& vs = it->second;
            if (vs)
            {
                ierr = VecScatterDestroy(&vs);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    // Reset the level numbers.
    d_coarsest_ln = coarsest_ln;
    d_finest_ln = finest_ln;

    // Resize some arrays.
    d_nclouds.resize(d_finest_ln + 1, 0);
    d_cloud_names.resize(d_finest_ln + 1);
    d_cloud_nmarks.resize(d_finest_ln + 1);
    d_cloud_first_lag_idx.resize(d_finest_ln + 1);

    d_nblocks.resize(d_finest_ln + 1, 0);
    d_block_names.resize(d_finest_ln + 1);
    d_block_nelems.resize(d_finest_ln + 1);
    d_block_periodic.resize(d_finest_ln + 1);
    d_block_first_lag_idx.resize(d_finest_ln + 1);

    d_nmbs.resize(d_finest_ln + 1, 0);
    d_mb_nblocks.resize(d_finest_ln + 1);
    d_mb_names.resize(d_finest_ln + 1);
    d_mb_nelems.resize(d_finest_ln + 1);
    d_mb_periodic.resize(d_finest_ln + 1);
    d_mb_first_lag_idx.resize(d_finest_ln + 1);

    d_nucd_meshes.resize(d_finest_ln + 1, 0);
    d_ucd_mesh_names.resize(d_finest_ln + 1);
    d_ucd_mesh_vertices.resize(d_finest_ln + 1);
    d_ucd_mesh_edge_maps.resize(d_finest_ln + 1);

    d_coords_data.resize(d_finest_ln + 1, NULL);
    d_nvars.resize(d_finest_ln + 1, 0);
    d_var_names.resize(d_finest_ln + 1);
    d_var_start_depths.resize(d_finest_ln + 1);
    d_var_plot_depths.resize(d_finest_ln + 1);
    d_var_depths.resize(d_finest_ln + 1);
    d_var_data.resize(d_finest_ln + 1);

    d_ao.resize(d_finest_ln + 1);
    d_build_vec_scatters.resize(d_finest_ln + 1);
    d_src_vec.resize(d_finest_ln + 1);
    d_dst_vec.resize(d_finest_ln + 1);
    d_vec_scatter.resize(d_finest_ln + 1);
    return;
} // resetLevels

void
LSiloDataWriter::registerMarkerCloud(const std::string& name,
                                     const int nmarks,
                                     const int first_lag_idx,
                                     const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(std::min(level_number, d_coarsest_ln), std::max(level_number, d_finest_ln));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(nmarks > 0);
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
#endif

    // Check to see if the cloud name has already been registered.
    if (find(d_block_names[level_number].begin(), d_block_names[level_number].end(), name) !=
        d_block_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerMarkerCloud()\n"
                                 << "  marker clouds must have unique names.\n"
                                 << "  a Cartesian block named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    if (find(d_mb_names[level_number].begin(), d_mb_names[level_number].end(), name) != d_mb_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerMarkerCloud()\n"
                                 << "  marker clouds must have unique names.\n"
                                 << "  a Cartesian multiblock named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    if (find(d_ucd_mesh_names[level_number].begin(), d_ucd_mesh_names[level_number].end(), name) !=
        d_ucd_mesh_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerMarkerCloud()\n"
                                 << "  marker clouds must have unique names.\n"
                                 << "  an unstructured mesh named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    // Check to see if we are updating a previously registered cloud.
    for (int k = 0; k < d_nclouds[level_number]; ++k)
    {
        if (d_cloud_names[level_number][k] == name)
        {
            d_cloud_nmarks[level_number][k] = nmarks;
            d_cloud_first_lag_idx[level_number][k] = first_lag_idx;
            return;
        }
    }

    // Record the layout of the marker cloud.
    ++d_nclouds[level_number];
    d_cloud_names[level_number].push_back(name);
    d_cloud_nmarks[level_number].push_back(nmarks);
    d_cloud_first_lag_idx[level_number].push_back(first_lag_idx);
    return;
} // registerMarkerCloud

void
LSiloDataWriter::registerLogicallyCartesianBlock(const std::string& name,
                                                 const IntVector<NDIM>& nelem,
                                                 const IntVector<NDIM>& periodic,
                                                 const int first_lag_idx,
                                                 const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(std::min(level_number, d_coarsest_ln), std::max(level_number, d_finest_ln));
    }

#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(nelem(d) > 0);
        TBOX_ASSERT(periodic(d) == 0 || periodic(d) == 1);
    }
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
#endif

    // Check to see if the block name has already been registered.
    if (find(d_cloud_names[level_number].begin(), d_cloud_names[level_number].end(), name) !=
        d_cloud_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                                 << "  Cartesian blocks must have unique names.\n"
                                 << "  a marker cloud named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    if (find(d_mb_names[level_number].begin(), d_mb_names[level_number].end(), name) != d_mb_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                                 << "  Cartesian blocks must have unique names.\n"
                                 << "  a Cartesian multiblock named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    if (find(d_ucd_mesh_names[level_number].begin(), d_ucd_mesh_names[level_number].end(), name) !=
        d_ucd_mesh_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                                 << "  Cartesian blocks must have unique names.\n"
                                 << "  an unstructured mesh named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    // Check to see if we are updating a previously registered block.
    for (int k = 0; k < d_nblocks[level_number]; ++k)
    {
        if (d_block_names[level_number][k] == name)
        {
            d_block_nelems[level_number][k] = nelem;
            d_block_periodic[level_number][k] = periodic;
            d_block_first_lag_idx[level_number][k] = first_lag_idx;
            return;
        }
    }

    // Record the layout of the logically Cartesian block.
    ++d_nblocks[level_number];
    d_block_names[level_number].push_back(name);
    d_block_nelems[level_number].push_back(nelem);
    d_block_periodic[level_number].push_back(periodic);
    d_block_first_lag_idx[level_number].push_back(first_lag_idx);
    return;
} // registerLogicallyCartesianBlock

void
LSiloDataWriter::registerLogicallyCartesianMultiblock(const std::string& name,
                                                      const std::vector<IntVector<NDIM> >& nelem,
                                                      const std::vector<IntVector<NDIM> >& periodic,
                                                      const std::vector<int>& first_lag_idx,
                                                      const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(std::min(level_number, d_coarsest_ln), std::max(level_number, d_finest_ln));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(periodic.size() == nelem.size());
    TBOX_ASSERT(first_lag_idx.size() == nelem.size());
    size_t sz = nelem.size();
    for (size_t i = 0; i < sz; ++i)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(nelem[i](d) > 0);
            TBOX_ASSERT(periodic[i](d) == 0 || periodic[i](d) == 1);
        }
    }
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
#endif

    // Check to see if the multiblock name has already been registered.
    if (find(d_cloud_names[level_number].begin(), d_cloud_names[level_number].end(), name) !=
        d_cloud_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianMultiblock()\n"
                                 << "  Cartesian multiblocks must have unique names.\n"
                                 << "  a marker cloud named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    if (find(d_block_names[level_number].begin(), d_block_names[level_number].end(), name) !=
        d_block_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianMultiblock()\n"
                                 << "  Cartesian multiblocks must have unique names.\n"
                                 << "  a Cartesian block named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    if (find(d_ucd_mesh_names[level_number].begin(), d_ucd_mesh_names[level_number].end(), name) !=
        d_ucd_mesh_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianMultiblock()\n"
                                 << "  Cartesian multiblocks must have unique names.\n"
                                 << "  an unstructured mesh named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    // Check to see if we are updating a previously registered multiblock.
    for (int k = 0; k < d_nmbs[level_number]; ++k)
    {
        if (d_mb_names[level_number][k] == name)
        {
            d_mb_nblocks[level_number][k] = static_cast<int>(nelem.size());
            d_mb_nelems[level_number][k] = nelem;
            d_mb_periodic[level_number][k] = periodic;
            d_mb_first_lag_idx[level_number][k] = first_lag_idx;
            return;
        }
    }

    // Record the layout of the logically Cartesian multiblock.
    ++d_nmbs[level_number];
    d_mb_names[level_number].push_back(name);
    d_mb_nblocks[level_number].push_back(static_cast<int>(nelem.size()));
    d_mb_nelems[level_number].push_back(nelem);
    d_mb_periodic[level_number].push_back(periodic);
    d_mb_first_lag_idx[level_number].push_back(first_lag_idx);
    return;
} // registerLogicallyCartesianMultiblock

void
LSiloDataWriter::registerUnstructuredMesh(const std::string& name,
                                          const std::multimap<int, std::pair<int, int> >& edge_map,
                                          const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(std::min(level_number, d_coarsest_ln), std::max(level_number, d_finest_ln));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
#endif

    // Check to see if the unstructured mesh name has already been registered.
    if (find(d_cloud_names[level_number].begin(), d_cloud_names[level_number].end(), name) !=
        d_cloud_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerUnstructuredMesh()\n"
                                 << "  unstructured meshes must have unique names.\n"
                                 << "  a marker cloud named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    if (find(d_block_names[level_number].begin(), d_block_names[level_number].end(), name) !=
        d_block_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerUnstructuredMesh()\n"
                                 << "  unstructured meshes must have unique names.\n"
                                 << "  a Cartesian block named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    if (find(d_mb_names[level_number].begin(), d_mb_names[level_number].end(), name) != d_mb_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerUnstructuredMesh()\n"
                                 << "  unstructured meshes must have unique names.\n"
                                 << "  a Cartesian multiblock named ``"
                                 << name
                                 << "'' has already been registered.\n");
    }

    // Extract the list of vertices from the list of edges.
    std::set<int> vertices;
    for (std::multimap<int, std::pair<int, int> >::const_iterator it = edge_map.begin(); it != edge_map.end(); ++it)
    {
        const std::pair<int, int>& e = it->second;
        vertices.insert(e.first);
        vertices.insert(e.second);
    }

    // Check to see if we are updating a previously registered unstructured
    // mesh.
    for (int k = 0; k < d_nucd_meshes[level_number]; ++k)
    {
        if (d_ucd_mesh_names[level_number][k] == name)
        {
            d_ucd_mesh_vertices[level_number][k] = vertices;
            d_ucd_mesh_edge_maps[level_number][k] = edge_map;
            return;
        }
    }

    // Record the layout of the unstructured mesh.
    ++d_nucd_meshes[level_number];
    d_ucd_mesh_names[level_number].push_back(name);
    d_ucd_mesh_vertices[level_number].push_back(vertices);
    d_ucd_mesh_edge_maps[level_number].push_back(edge_map);
    return;
} // registerUnstructuredMesh

void
LSiloDataWriter::registerCoordsData(Pointer<LData> coords_data, const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(std::min(level_number, d_coarsest_ln), std::max(level_number, d_finest_ln));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(coords_data);
    TBOX_ASSERT(coords_data->getDepth() == NDIM);
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
#endif
    d_coords_data[level_number] = coords_data;
    return;
} // registerCoordsData

void
LSiloDataWriter::registerVariableData(const std::string& var_name, Pointer<LData> var_data, const int level_number)
{
    const int start_depth = 0;
    const int var_depth = var_data->getDepth();
    registerVariableData(var_name, var_data, start_depth, var_depth, level_number);
    return;
} // registerVariableData

void
LSiloDataWriter::registerVariableData(const std::string& var_name,
                                      Pointer<LData> var_data,
                                      const int start_depth,
                                      const int var_depth,
                                      const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(std::min(level_number, d_coarsest_ln), std::max(level_number, d_finest_ln));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(!var_name.empty());
    TBOX_ASSERT(var_data);
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
#endif
    if (find(d_var_names[level_number].begin(), d_var_names[level_number].end(), var_name) !=
        d_var_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerVariableData()\n"
                                 << "  variable with name "
                                 << var_name
                                 << " already registered for plotting\n"
                                 << "  on patch level "
                                 << level_number
                                 << std::endl);
    }
    ++d_nvars[level_number];
    d_var_names[level_number].push_back(var_name);
    d_var_start_depths[level_number].push_back(start_depth);
    d_var_plot_depths[level_number].push_back(var_depth);
    d_var_depths[level_number].push_back(var_data->getDepth());
    d_var_data[level_number].push_back(var_data);
    return;
} // registerVariableData

void
LSiloDataWriter::registerLagrangianAO(AO& ao, const int level_number)
{
    if (level_number < d_coarsest_ln || level_number > d_finest_ln)
    {
        resetLevels(std::min(level_number, d_coarsest_ln), std::max(level_number, d_finest_ln));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(d_coarsest_ln <= level_number && d_finest_ln >= level_number);
#endif
    d_ao[level_number] = ao;
    d_build_vec_scatters[level_number] = true;
    return;
} // registerLagrangianAO

void
LSiloDataWriter::registerLagrangianAO(std::vector<AO>& ao, const int coarsest_ln, const int finest_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(coarsest_ln <= finest_ln);
#endif

    if (coarsest_ln < d_coarsest_ln || finest_ln > d_finest_ln)
    {
        resetLevels(std::min(coarsest_ln, d_coarsest_ln), std::max(finest_ln, d_finest_ln));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(d_coarsest_ln <= coarsest_ln && finest_ln <= d_finest_ln);
#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        registerLagrangianAO(ao[ln], ln);
    }
    return;
} // registerLagrangianAO

void
LSiloDataWriter::writePlotData(const int time_step_number, const double simulation_time)
{
#if defined(IBTK_HAVE_SILO)
#if !defined(NDEBUG)
    TBOX_ASSERT(time_step_number >= 0);
    TBOX_ASSERT(!d_dump_directory_name.empty());
#endif

    if (time_step_number <= d_time_step_number)
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                 << "  data writer with name "
                                 << d_object_name
                                 << "\n"
                                 << "  time step number: "
                                 << time_step_number
                                 << " is <= last time step number: "
                                 << d_time_step_number
                                 << std::endl);
    }
    d_time_step_number = time_step_number;

    if (d_dump_directory_name.empty())
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                 << "  data writer with name "
                                 << d_object_name
                                 << "\n"
                                 << "  dump directory name is empty"
                                 << std::endl);
    }

    int ierr;
    char temp_buf[SILO_NAME_BUFSIZE];
    std::string current_file_name;
    DBfile* dbfile;
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int mpi_nodes = SAMRAI_MPI::getNodes();

    // Construct the VecScatter objects required to write the plot data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (d_build_vec_scatters[ln])
        {
            buildVecScatters(d_ao[ln], ln);
        }
        d_build_vec_scatters[ln] = false;
    }

    // Create the working directory.
    sprintf(temp_buf, "%06d", d_time_step_number);
    std::string current_dump_directory_name = SILO_DUMP_DIR_PREFIX + temp_buf;
    std::string dump_dirname = d_dump_directory_name + "/" + current_dump_directory_name;

    Utilities::recursiveMkdir(dump_dirname);

    // Create one local DBfile per MPI process.
    sprintf(temp_buf, "%04d", mpi_rank);
    current_file_name = dump_dirname + "/" + SILO_PROCESSOR_FILE_PREFIX;
    current_file_name += temp_buf;
    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

    if (!(dbfile = DBCreate(current_file_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB)))
    {
        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                 << "  Could not create DBfile named "
                                 << current_file_name
                                 << std::endl);
    }

    std::vector<std::vector<int> > meshtype(d_finest_ln + 1), vartype(d_finest_ln + 1);
    std::vector<std::vector<std::vector<int> > > multimeshtype(d_finest_ln + 1), multivartype(d_finest_ln + 1);

    // Set the local data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (d_coords_data[ln])
        {
            // Scatter the data from "global" to "local" form.
            Vec local_X_vec;
            ierr = VecDuplicate(d_dst_vec[ln][NDIM], &local_X_vec);
            IBTK_CHKERRQ(ierr);

            Vec global_X_vec = d_coords_data[ln]->getVec();
            ierr = VecScatterBegin(d_vec_scatter[ln][NDIM], global_X_vec, local_X_vec, INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterEnd(d_vec_scatter[ln][NDIM], global_X_vec, local_X_vec, INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);

            double* local_X_arr;
            ierr = VecGetArray(local_X_vec, &local_X_arr);
            IBTK_CHKERRQ(ierr);

            std::vector<Vec> local_v_vecs;
            std::vector<double*> local_v_arrs;

            for (int v = 0; v < d_nvars[ln]; ++v)
            {
                const int var_depth = d_var_depths[ln][v];
                Vec local_v_vec;
                ierr = VecDuplicate(d_dst_vec[ln][var_depth], &local_v_vec);
                IBTK_CHKERRQ(ierr);

                Vec global_v_vec = d_var_data[ln][v]->getVec();
                ierr = VecScatterBegin(
                    d_vec_scatter[ln][var_depth], global_v_vec, local_v_vec, INSERT_VALUES, SCATTER_FORWARD);
                IBTK_CHKERRQ(ierr);
                ierr = VecScatterEnd(
                    d_vec_scatter[ln][var_depth], global_v_vec, local_v_vec, INSERT_VALUES, SCATTER_FORWARD);
                IBTK_CHKERRQ(ierr);

                double* local_v_arr;
                ierr = VecGetArray(local_v_vec, &local_v_arr);
                IBTK_CHKERRQ(ierr);

                local_v_vecs.push_back(local_v_vec);
                local_v_arrs.push_back(local_v_arr);
            }

            // Keep track of the current offset in the local Vec data.
            int offset = 0;

            // Add the local clouds to the local DBfile.
            for (int cloud = 0; cloud < d_nclouds[ln]; ++cloud)
            {
                const int nmarks = d_cloud_nmarks[ln][cloud];

                std::ostringstream stream;
                stream << "level_" << ln << "_cloud_" << cloud;
                std::string dirname = stream.str();

                if (DBMkDir(dbfile, dirname.c_str()) == -1)
                {
                    TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                             << "  Could not create directory named "
                                             << dirname
                                             << std::endl);
                }

                const double* const X = local_X_arr + NDIM * offset;
                std::vector<const double*> var_vals(d_nvars[ln]);
                for (int v = 0; v < d_nvars[ln]; ++v)
                {
                    var_vals[v] = local_v_arrs[v] + d_var_depths[ln][v] * offset;
                }

                build_local_marker_cloud(dbfile,
                                         dirname,
                                         nmarks,
                                         X,
                                         d_nvars[ln],
                                         d_var_names[ln],
                                         d_var_start_depths[ln],
                                         d_var_plot_depths[ln],
                                         d_var_depths[ln],
                                         var_vals,
                                         time_step_number,
                                         simulation_time);

                offset += nmarks;
            }

            // Add the local blocks to the local DBfile.
            for (int block = 0; block < d_nblocks[ln]; ++block)
            {
                const IntVector<NDIM>& nelem = d_block_nelems[ln][block];
                const IntVector<NDIM>& periodic = d_block_periodic[ln][block];
                const int ntot = nelem.getProduct();

                std::ostringstream stream;
                stream << "level_" << ln << "_block_" << block;
                std::string dirname = stream.str();

                if (DBMkDir(dbfile, dirname.c_str()) == -1)
                {
                    TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                             << "  Could not create directory named "
                                             << dirname
                                             << std::endl);
                }

                const double* const X = local_X_arr + NDIM * offset;
                std::vector<const double*> var_vals(d_nvars[ln]);
                for (int v = 0; v < d_nvars[ln]; ++v)
                {
                    var_vals[v] = local_v_arrs[v] + d_var_depths[ln][v] * offset;
                }

                build_local_curv_block(dbfile,
                                       dirname,
                                       nelem,
                                       periodic,
                                       X,
                                       d_nvars[ln],
                                       d_var_names[ln],
                                       d_var_start_depths[ln],
                                       d_var_plot_depths[ln],
                                       d_var_depths[ln],
                                       var_vals,
                                       time_step_number,
                                       simulation_time);
                meshtype[ln].push_back(DB_QUAD_CURV);
                vartype[ln].push_back(DB_QUADVAR);

                offset += ntot;
            }

            // Add the local multiblocks to the local DBfile.
            multimeshtype[ln].resize(d_nmbs[ln]);
            multivartype[ln].resize(d_nmbs[ln]);
            for (int mb = 0; mb < d_nmbs[ln]; ++mb)
            {
                for (int block = 0; block < d_mb_nblocks[ln][mb]; ++block)
                {
                    const IntVector<NDIM>& nelem = d_mb_nelems[ln][mb][block];
                    const IntVector<NDIM>& periodic = d_mb_periodic[ln][mb][block];
                    const int ntot = nelem.getProduct();

                    std::ostringstream stream;
                    stream << "level_" << ln << "_mb_" << mb << "_block_" << block;
                    std::string dirname = stream.str();

                    if (DBMkDir(dbfile, dirname.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                                 << "  Could not create directory named "
                                                 << dirname
                                                 << std::endl);
                    }

                    const double* const X = local_X_arr + NDIM * offset;
                    std::vector<const double*> var_vals(d_nvars[ln]);
                    for (int v = 0; v < d_nvars[ln]; ++v)
                    {
                        var_vals[v] = local_v_arrs[v] + d_var_depths[ln][v] * offset;
                    }

                    build_local_curv_block(dbfile,
                                           dirname,
                                           nelem,
                                           periodic,
                                           X,
                                           d_nvars[ln],
                                           d_var_names[ln],
                                           d_var_start_depths[ln],
                                           d_var_plot_depths[ln],
                                           d_var_depths[ln],
                                           var_vals,
                                           time_step_number,
                                           simulation_time);
                    multimeshtype[ln][mb].push_back(DB_QUAD_CURV);
                    multivartype[ln][mb].push_back(DB_QUADVAR);

                    offset += ntot;
                }
            }

            // Add the local UCD meshes to the local DBfile.
            for (int mesh = 0; mesh < d_nucd_meshes[ln]; ++mesh)
            {
                const std::set<int>& vertices = d_ucd_mesh_vertices[ln][mesh];
                const std::multimap<int, std::pair<int, int> >& edge_map = d_ucd_mesh_edge_maps[ln][mesh];
                const size_t ntot = vertices.size();

                std::ostringstream stream;
                stream << "level_" << ln << "_mesh_" << mesh;
                std::string dirname = stream.str();

                if (DBMkDir(dbfile, dirname.c_str()) == -1)
                {
                    TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                             << "  Could not create directory named "
                                             << dirname
                                             << std::endl);
                }

                const double* const X = local_X_arr + NDIM * offset;
                std::vector<const double*> var_vals(d_nvars[ln]);
                for (int v = 0; v < d_nvars[ln]; ++v)
                {
                    var_vals[v] = local_v_arrs[v] + d_var_depths[ln][v] * offset;
                }

                build_local_ucd_mesh(dbfile,
                                     dirname,
                                     vertices,
                                     edge_map,
                                     X,
                                     d_nvars[ln],
                                     d_var_names[ln],
                                     d_var_start_depths[ln],
                                     d_var_plot_depths[ln],
                                     d_var_depths[ln],
                                     var_vals,
                                     time_step_number,
                                     simulation_time);

                offset += ntot;
            }

            // Clean up allocated data.
            ierr = VecRestoreArray(local_X_vec, &local_X_arr);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&local_X_vec);
            IBTK_CHKERRQ(ierr);
            for (int v = 0; v < d_nvars[ln]; ++v)
            {
                ierr = VecRestoreArray(local_v_vecs[v], &local_v_arrs[v]);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&local_v_vecs[v]);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    DBClose(dbfile);

    // Send data to the root MPI process required to create the multimesh and
    // multivar objects.
    std::vector<std::vector<int> > nclouds_per_proc, nblocks_per_proc, nmbs_per_proc, nucd_meshes_per_proc;
    std::vector<std::vector<std::vector<int> > > meshtypes_per_proc, vartypes_per_proc, mb_nblocks_per_proc;
    std::vector<std::vector<std::vector<std::vector<int> > > > multimeshtypes_per_proc, multivartypes_per_proc;
    std::vector<std::vector<std::vector<std::string> > > cloud_names_per_proc, block_names_per_proc, mb_names_per_proc,
        ucd_mesh_names_per_proc;

    if (mpi_rank == SILO_MPI_ROOT)
    {
        nclouds_per_proc.resize(d_finest_ln + 1);
        nblocks_per_proc.resize(d_finest_ln + 1);
        nmbs_per_proc.resize(d_finest_ln + 1);
        nucd_meshes_per_proc.resize(d_finest_ln + 1);
        meshtypes_per_proc.resize(d_finest_ln + 1);
        vartypes_per_proc.resize(d_finest_ln + 1);
        mb_nblocks_per_proc.resize(d_finest_ln + 1);
        multimeshtypes_per_proc.resize(d_finest_ln + 1);
        multivartypes_per_proc.resize(d_finest_ln + 1);
        cloud_names_per_proc.resize(d_finest_ln + 1);
        block_names_per_proc.resize(d_finest_ln + 1);
        mb_names_per_proc.resize(d_finest_ln + 1);
        ucd_mesh_names_per_proc.resize(d_finest_ln + 1);
    }

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (mpi_rank == SILO_MPI_ROOT)
        {
            nclouds_per_proc[ln].resize(mpi_nodes);
            nblocks_per_proc[ln].resize(mpi_nodes);
            nmbs_per_proc[ln].resize(mpi_nodes);
            nucd_meshes_per_proc[ln].resize(mpi_nodes);
            meshtypes_per_proc[ln].resize(mpi_nodes);
            vartypes_per_proc[ln].resize(mpi_nodes);
            mb_nblocks_per_proc[ln].resize(mpi_nodes);
            multimeshtypes_per_proc[ln].resize(mpi_nodes);
            multivartypes_per_proc[ln].resize(mpi_nodes);
            cloud_names_per_proc[ln].resize(mpi_nodes);
            block_names_per_proc[ln].resize(mpi_nodes);
            mb_names_per_proc[ln].resize(mpi_nodes);
            ucd_mesh_names_per_proc[ln].resize(mpi_nodes);
        }

        // Set the values for the root process.
        if (mpi_rank == SILO_MPI_ROOT)
        {
            nclouds_per_proc[ln][mpi_rank] = d_nclouds[ln];
            nblocks_per_proc[ln][mpi_rank] = d_nblocks[ln];
            nmbs_per_proc[ln][mpi_rank] = d_nmbs[ln];
            nucd_meshes_per_proc[ln][mpi_rank] = d_nucd_meshes[ln];
            meshtypes_per_proc[ln][mpi_rank] = meshtype[ln];
            vartypes_per_proc[ln][mpi_rank] = vartype[ln];
            mb_nblocks_per_proc[ln][mpi_rank] = d_mb_nblocks[ln];
            multimeshtypes_per_proc[ln][mpi_rank] = multimeshtype[ln];
            multivartypes_per_proc[ln][mpi_rank] = multivartype[ln];
            cloud_names_per_proc[ln][mpi_rank] = d_cloud_names[ln];
            block_names_per_proc[ln][mpi_rank] = d_block_names[ln];
            mb_names_per_proc[ln][mpi_rank] = d_mb_names[ln];
            ucd_mesh_names_per_proc[ln][mpi_rank] = d_ucd_mesh_names[ln];
        }

        // Get the values for the non-root processes.
        int one = 1;
        for (int proc = 0; proc < mpi_nodes; ++proc)
        {
            // Skip the root process; we already have those values.
            if (proc == SILO_MPI_ROOT)
            {
                proc += 1;
                if (proc >= mpi_nodes) break;
            }

            if (mpi_rank == proc)
            {
                SAMRAI_MPI::send(&d_nclouds[ln], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                SAMRAI_MPI::recv(&nclouds_per_proc[ln][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nclouds[ln] > 0)
            {
                int num_bytes;
                for (int cloud = 0; cloud < d_nclouds[ln]; ++cloud)
                {
                    num_bytes = static_cast<int>((d_cloud_names[ln][cloud].size() + 1) * sizeof(char));
                    SAMRAI_MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);
                    SAMRAI_MPI::sendBytes(
                        static_cast<const void*>(d_cloud_names[ln][cloud].c_str()), num_bytes, SILO_MPI_ROOT);
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nclouds_per_proc[ln][proc] > 0)
            {
                cloud_names_per_proc[ln][proc].resize(nclouds_per_proc[ln][proc]);

                int num_bytes;
                char* name;
                for (int cloud = 0; cloud < nclouds_per_proc[ln][proc]; ++cloud)
                {
                    SAMRAI_MPI::recv(&num_bytes, one, proc, false);
                    name = new char[num_bytes / sizeof(char)];
                    SAMRAI_MPI::recvBytes(static_cast<void*>(name), num_bytes);
                    cloud_names_per_proc[ln][proc][cloud].assign(name);
                    delete[] name;
                }
            }

            if (mpi_rank == proc)
            {
                SAMRAI_MPI::send(&d_nblocks[ln], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                SAMRAI_MPI::recv(&nblocks_per_proc[ln][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nblocks[ln] > 0)
            {
                SAMRAI_MPI::send(&meshtype[ln][0], d_nblocks[ln], SILO_MPI_ROOT, false);
                SAMRAI_MPI::send(&vartype[ln][0], d_nblocks[ln], SILO_MPI_ROOT, false);

                int num_bytes;
                for (int block = 0; block < d_nblocks[ln]; ++block)
                {
                    num_bytes = static_cast<int>((d_block_names[ln][block].size() + 1) * sizeof(char));
                    SAMRAI_MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);
                    SAMRAI_MPI::sendBytes(
                        static_cast<const void*>(d_block_names[ln][block].c_str()), num_bytes, SILO_MPI_ROOT);
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nblocks_per_proc[ln][proc] > 0)
            {
                meshtypes_per_proc[ln][proc].resize(nblocks_per_proc[ln][proc]);
                vartypes_per_proc[ln][proc].resize(nblocks_per_proc[ln][proc]);
                block_names_per_proc[ln][proc].resize(nblocks_per_proc[ln][proc]);

                SAMRAI_MPI::recv(&meshtypes_per_proc[ln][proc][0], nblocks_per_proc[ln][proc], proc, false);
                SAMRAI_MPI::recv(&vartypes_per_proc[ln][proc][0], nblocks_per_proc[ln][proc], proc, false);

                int num_bytes;
                char* name;
                for (int block = 0; block < nblocks_per_proc[ln][proc]; ++block)
                {
                    SAMRAI_MPI::recv(&num_bytes, one, proc, false);
                    name = new char[num_bytes / sizeof(char)];
                    SAMRAI_MPI::recvBytes(static_cast<void*>(name), num_bytes);
                    block_names_per_proc[ln][proc][block].assign(name);
                    delete[] name;
                }
            }

            if (mpi_rank == proc)
            {
                SAMRAI_MPI::send(&d_nmbs[ln], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                SAMRAI_MPI::recv(&nmbs_per_proc[ln][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nmbs[ln] > 0)
            {
                SAMRAI_MPI::send(&d_mb_nblocks[ln][0], d_nmbs[ln], SILO_MPI_ROOT, false);

                int num_bytes;
                for (int mb = 0; mb < d_nmbs[ln]; ++mb)
                {
                    SAMRAI_MPI::send(&multimeshtype[ln][mb][0], d_mb_nblocks[ln][mb], SILO_MPI_ROOT, false);
                    SAMRAI_MPI::send(&multivartype[ln][mb][0], d_mb_nblocks[ln][mb], SILO_MPI_ROOT, false);

                    num_bytes = static_cast<int>(d_mb_names[ln][mb].size()) + 1;
                    SAMRAI_MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);

                    const void* mb_name_ptr = d_mb_names[ln][mb].c_str();
                    MPI_Send(const_cast<void*>(mb_name_ptr),
                             num_bytes,
                             MPI_CHAR,
                             SILO_MPI_ROOT,
                             SILO_MPI_TAG,
                             SAMRAI_MPI::commWorld);
                    const int tree = SAMRAI_MPI::getTreeDepth();
                    SAMRAI_MPI::updateOutgoingStatistics(tree, num_bytes * sizeof(char));
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nmbs_per_proc[ln][proc] > 0)
            {
                mb_nblocks_per_proc[ln][proc].resize(nmbs_per_proc[ln][proc]);
                multimeshtypes_per_proc[ln][proc].resize(nmbs_per_proc[ln][proc]);
                multivartypes_per_proc[ln][proc].resize(nmbs_per_proc[ln][proc]);
                mb_names_per_proc[ln][proc].resize(nmbs_per_proc[ln][proc]);

                SAMRAI_MPI::recv(&mb_nblocks_per_proc[ln][proc][0], nmbs_per_proc[ln][proc], proc, false);

                int num_bytes;
                char* name;
                for (int mb = 0; mb < nmbs_per_proc[ln][proc]; ++mb)
                {
                    multimeshtypes_per_proc[ln][proc][mb].resize(mb_nblocks_per_proc[ln][proc][mb]);
                    multivartypes_per_proc[ln][proc][mb].resize(mb_nblocks_per_proc[ln][proc][mb]);

                    SAMRAI_MPI::recv(
                        &multimeshtypes_per_proc[ln][proc][mb][0], mb_nblocks_per_proc[ln][proc][mb], proc, false);
                    SAMRAI_MPI::recv(
                        &multivartypes_per_proc[ln][proc][mb][0], mb_nblocks_per_proc[ln][proc][mb], proc, false);

                    SAMRAI_MPI::recv(&num_bytes, one, proc, false);
                    name = new char[num_bytes];

                    MPI_Status status;
                    MPI_Recv(name, num_bytes, MPI_CHAR, proc, SILO_MPI_TAG, SAMRAI_MPI::commWorld, &status);
                    const int tree = SAMRAI_MPI::getTreeDepth();
                    SAMRAI_MPI::updateIncomingStatistics(tree, num_bytes * sizeof(char));

                    mb_names_per_proc[ln][proc][mb].assign(name);
                    delete[] name;
                }
            }

            if (mpi_rank == proc)
            {
                SAMRAI_MPI::send(&d_nucd_meshes[ln], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                SAMRAI_MPI::recv(&nucd_meshes_per_proc[ln][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nucd_meshes[ln] > 0)
            {
                int num_bytes;
                for (int mesh = 0; mesh < d_nucd_meshes[ln]; ++mesh)
                {
                    num_bytes = static_cast<int>((d_ucd_mesh_names[ln][mesh].size() + 1) * sizeof(char));
                    SAMRAI_MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);
                    SAMRAI_MPI::sendBytes(
                        static_cast<const void*>(d_ucd_mesh_names[ln][mesh].c_str()), num_bytes, SILO_MPI_ROOT);
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nucd_meshes_per_proc[ln][proc] > 0)
            {
                ucd_mesh_names_per_proc[ln][proc].resize(nucd_meshes_per_proc[ln][proc]);
                int num_bytes;
                char* name;
                for (int mesh = 0; mesh < nucd_meshes_per_proc[ln][proc]; ++mesh)
                {
                    SAMRAI_MPI::recv(&num_bytes, one, proc, false);
                    name = new char[num_bytes / sizeof(char)];
                    SAMRAI_MPI::recvBytes(static_cast<void*>(name), num_bytes);
                    ucd_mesh_names_per_proc[ln][proc][mesh].assign(name);
                    delete[] name;
                }
            }

            SAMRAI_MPI::barrier();
        }
    }

    if (mpi_rank == SILO_MPI_ROOT)
    {
        // Create and initialize the multimesh Silo database on the root MPI
        // process.
        sprintf(temp_buf, "%06d", d_time_step_number);
        std::string summary_file_name =
            dump_dirname + "/" + SILO_SUMMARY_FILE_PREFIX + temp_buf + SILO_SUMMARY_FILE_POSTFIX;
        if (!(dbfile = DBCreate(summary_file_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB)))
        {
            TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                     << "  Could not create DBfile named "
                                     << summary_file_name
                                     << std::endl);
        }

        int cycle = time_step_number;
        float time = float(simulation_time);
        double dtime = simulation_time;

        static const int MAX_OPTS = 3;
        DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
        DBAddOption(optlist, DBOPT_CYCLE, &cycle);
        DBAddOption(optlist, DBOPT_TIME, &time);
        DBAddOption(optlist, DBOPT_DTIME, &dtime);

        for (int proc = 0; proc < mpi_nodes; ++proc)
        {
            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                for (int cloud = 0; cloud < nclouds_per_proc[ln][proc]; ++cloud)
                {
                    sprintf(temp_buf, "%04d", proc);
                    current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                    current_file_name += temp_buf;
                    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                    std::ostringstream stream;
                    stream << current_file_name << ":level_" << ln << "_cloud_" << cloud << "/mesh";
                    std::string meshname = stream.str();
                    char* meshname_ptr = const_cast<char*>(meshname.c_str());
                    int meshtype = DB_POINTMESH;

                    std::string& cloud_name = cloud_names_per_proc[ln][proc][cloud];

                    DBPutMultimesh(dbfile, cloud_name.c_str(), 1, &meshname_ptr, &meshtype, optlist);

                    if (DBMkDir(dbfile, cloud_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                                 << "  Could not create directory named "
                                                 << cloud_name
                                                 << std::endl);
                    }
                }

                for (int block = 0; block < nblocks_per_proc[ln][proc]; ++block)
                {
                    sprintf(temp_buf, "%04d", proc);
                    current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                    current_file_name += temp_buf;
                    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                    std::ostringstream stream;
                    stream << current_file_name << ":level_" << ln << "_block_" << block << "/mesh";
                    std::string meshname = stream.str();
                    char* meshname_ptr = const_cast<char*>(meshname.c_str());
                    int meshtype = meshtypes_per_proc[ln][proc][block];

                    std::string& block_name = block_names_per_proc[ln][proc][block];

                    DBPutMultimesh(dbfile, block_name.c_str(), 1, &meshname_ptr, &meshtype, optlist);

                    if (DBMkDir(dbfile, block_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                                 << "  Could not create directory named "
                                                 << block_name
                                                 << std::endl);
                    }
                }

                for (int mb = 0; mb < nmbs_per_proc[ln][proc]; ++mb)
                {
                    sprintf(temp_buf, "%04d", proc);
                    current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                    current_file_name += temp_buf;
                    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                    const int nblocks = mb_nblocks_per_proc[ln][proc][mb];
                    char** meshnames = new char*[nblocks];

                    for (int block = 0; block < nblocks; ++block)
                    {
                        std::ostringstream stream;
                        stream << current_file_name << ":level_" << ln << "_mb_" << mb << "_block_" << block << "/mesh";
                        meshnames[block] = strdup(stream.str().c_str());
                    }

                    std::string& mb_name = mb_names_per_proc[ln][proc][mb];

                    DBPutMultimesh(dbfile,
                                   mb_name.c_str(),
                                   nblocks,
                                   meshnames,
                                   &multimeshtypes_per_proc[ln][proc][mb][0],
                                   optlist);

                    if (DBMkDir(dbfile, mb_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                                 << "  Could not create directory named "
                                                 << mb_name
                                                 << std::endl);
                    }

                    for (int block = 0; block < nblocks; ++block)
                    {
                        free(meshnames[block]);
                    }
                    delete[] meshnames;
                }

                for (int mesh = 0; mesh < nucd_meshes_per_proc[ln][proc]; ++mesh)
                {
                    sprintf(temp_buf, "%04d", proc);
                    current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                    current_file_name += temp_buf;
                    current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                    std::ostringstream stream;
                    stream << current_file_name << ":level_" << ln << "_mesh_" << mesh << "/mesh";
                    std::string meshname = stream.str();
                    char* meshname_ptr = const_cast<char*>(meshname.c_str());
                    int meshtype = DB_UCDMESH;

                    std::string& mesh_name = ucd_mesh_names_per_proc[ln][proc][mesh];

                    DBPutMultimesh(dbfile, mesh_name.c_str(), 1, &meshname_ptr, &meshtype, optlist);

                    if (DBMkDir(dbfile, mesh_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name << "::writePlotData()\n"
                                                 << "  Could not create directory named "
                                                 << mesh_name
                                                 << std::endl);
                    }
                }

                for (int v = 0; v < d_nvars[ln]; ++v)
                {
                    for (int cloud = 0; cloud < nclouds_per_proc[ln][proc]; ++cloud)
                    {
                        sprintf(temp_buf, "%04d", proc);
                        current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                        current_file_name += temp_buf;
                        current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                        std::ostringstream varname_stream;
                        varname_stream << current_file_name << ":level_" << ln << "_cloud_" << cloud << "/"
                                       << d_var_names[ln][v];
                        std::string varname = varname_stream.str();
                        char* varname_ptr = const_cast<char*>(varname.c_str());
                        int vartype = DB_POINTVAR;

                        std::string& cloud_name = cloud_names_per_proc[ln][proc][cloud];

                        std::ostringstream stream;
                        stream << cloud_name << "/" << d_var_names[ln][v];
                        std::string var_name = stream.str();

                        DBPutMultivar(dbfile, var_name.c_str(), 1, &varname_ptr, &vartype, optlist);
                    }

                    for (int block = 0; block < nblocks_per_proc[ln][proc]; ++block)
                    {
                        sprintf(temp_buf, "%04d", proc);
                        current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                        current_file_name += temp_buf;
                        current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                        std::ostringstream varname_stream;
                        varname_stream << current_file_name << ":level_" << ln << "_block_" << block << "/"
                                       << d_var_names[ln][v];
                        std::string varname = varname_stream.str();
                        char* varname_ptr = const_cast<char*>(varname.c_str());
                        int vartype = vartypes_per_proc[ln][proc][block];

                        std::string& block_name = block_names_per_proc[ln][proc][block];

                        std::ostringstream stream;
                        stream << block_name << "/" << d_var_names[ln][v];
                        std::string var_name = stream.str();

                        DBPutMultivar(dbfile, var_name.c_str(), 1, &varname_ptr, &vartype, optlist);
                    }

                    for (int mb = 0; mb < nmbs_per_proc[ln][proc]; ++mb)
                    {
                        sprintf(temp_buf, "%04d", proc);
                        current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                        current_file_name += temp_buf;
                        current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                        const int nblocks = mb_nblocks_per_proc[ln][proc][mb];
                        char** varnames = new char*[nblocks];

                        for (int block = 0; block < nblocks; ++block)
                        {
                            std::ostringstream varname_stream;
                            varname_stream << current_file_name << ":level_" << ln << "_mb_" << mb << "_block_" << block
                                           << d_var_names[ln][v];
                            varnames[block] = strdup(varname_stream.str().c_str());
                        }

                        std::string& mb_name = mb_names_per_proc[ln][proc][mb];

                        std::ostringstream stream;
                        stream << mb_name << "/" << d_var_names[ln][v];
                        std::string var_name = stream.str();

                        DBPutMultivar(dbfile,
                                      var_name.c_str(),
                                      nblocks,
                                      varnames,
                                      &multivartypes_per_proc[ln][proc][mb][0],
                                      optlist);

                        for (int block = 0; block < nblocks; ++block)
                        {
                            free(varnames[block]);
                        }
                        delete[] varnames;
                    }

                    for (int mesh = 0; mesh < nucd_meshes_per_proc[ln][proc]; ++mesh)
                    {
                        sprintf(temp_buf, "%04d", proc);
                        current_file_name = SILO_PROCESSOR_FILE_PREFIX;
                        current_file_name += temp_buf;
                        current_file_name += SILO_PROCESSOR_FILE_POSTFIX;

                        std::ostringstream varname_stream;
                        varname_stream << current_file_name << ":level_" << ln << "_mesh_" << mesh << "/"
                                       << d_var_names[ln][v];
                        std::string varname = varname_stream.str();
                        char* varname_ptr = const_cast<char*>(varname.c_str());
                        int vartype = DB_UCDVAR;

                        std::string& mesh_name = ucd_mesh_names_per_proc[ln][proc][mesh];

                        std::ostringstream stream;
                        stream << mesh_name << "/" << d_var_names[ln][v];
                        std::string var_name = stream.str();

                        DBPutMultivar(dbfile, var_name.c_str(), 1, &varname_ptr, &vartype, optlist);
                    }
                }
            }
        }

        DBClose(dbfile);

        // Create or update the dumps file on the root MPI process.
        static bool summary_file_opened = false;
        std::string path = d_dump_directory_name + "/" + VISIT_DUMPS_FILENAME;
        sprintf(temp_buf, "%06d", d_time_step_number);
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
    SAMRAI_MPI::barrier();
#else
    TBOX_WARNING("LSiloDataWriter::writePlotData(): SILO is not installed; cannot write data." << std::endl);
#endif // if defined(IBTK_HAVE_SILO)
    return;
} // writePlotData

void
LSiloDataWriter::putToDatabase(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    db->putInteger("LAG_SILO_DATA_WRITER_VERSION", LAG_SILO_DATA_WRITER_VERSION);

    db->putInteger("d_coarsest_ln", d_coarsest_ln);
    db->putInteger("d_finest_ln", d_finest_ln);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        std::ostringstream ln_stream;
        ln_stream << "_" << ln;
        const std::string ln_string = ln_stream.str();

        db->putInteger("d_nclouds" + ln_string, d_nclouds[ln]);
        if (d_nclouds[ln] > 0)
        {
            db->putStringArray(
                "d_cloud_names" + ln_string, &d_cloud_names[ln][0], static_cast<int>(d_cloud_names[ln].size()));
            db->putIntegerArray(
                "d_cloud_nmarks" + ln_string, &d_cloud_nmarks[ln][0], static_cast<int>(d_cloud_nmarks[ln].size()));
            db->putIntegerArray("d_cloud_first_lag_idx" + ln_string,
                                &d_cloud_first_lag_idx[ln][0],
                                static_cast<int>(d_cloud_first_lag_idx[ln].size()));
        }

        db->putInteger("d_nblocks" + ln_string, d_nblocks[ln]);
        if (d_nblocks[ln] > 0)
        {
            db->putStringArray(
                "d_block_names" + ln_string, &d_block_names[ln][0], static_cast<int>(d_block_names[ln].size()));

            std::vector<int> flattened_block_nelems;
            flattened_block_nelems.reserve(NDIM * d_block_nelems.size());
            for (std::vector<IntVector<NDIM> >::const_iterator cit = d_block_nelems[ln].begin();
                 cit != d_block_nelems[ln].end();
                 ++cit)
            {
                flattened_block_nelems.insert(flattened_block_nelems.end(), &(*cit)[0], &(*cit)[0] + NDIM);
            }
            db->putIntegerArray("flattened_block_nelems" + ln_string,
                                &flattened_block_nelems[0],
                                static_cast<int>(flattened_block_nelems.size()));

            std::vector<int> flattened_block_periodic;
            flattened_block_periodic.reserve(NDIM * d_block_periodic.size());
            for (std::vector<IntVector<NDIM> >::const_iterator cit = d_block_periodic[ln].begin();
                 cit != d_block_periodic[ln].end();
                 ++cit)
            {
                flattened_block_periodic.insert(flattened_block_periodic.end(), &(*cit)[0], &(*cit)[0] + NDIM);
            }
            db->putIntegerArray("flattened_block_periodic" + ln_string,
                                &flattened_block_periodic[0],
                                static_cast<int>(flattened_block_periodic.size()));

            db->putIntegerArray("d_block_first_lag_idx" + ln_string,
                                &d_block_first_lag_idx[ln][0],
                                static_cast<int>(d_block_first_lag_idx[ln].size()));
        }

        db->putInteger("d_nmbs" + ln_string, d_nmbs[ln]);
        if (d_nmbs[ln] > 0)
        {
            db->putStringArray("d_mb_names" + ln_string, &d_mb_names[ln][0], static_cast<int>(d_mb_names[ln].size()));

            for (int mb = 0; mb < d_nmbs[ln]; ++mb)
            {
                std::ostringstream mb_stream;
                mb_stream << "_" << mb;
                const std::string mb_string = mb_stream.str();

                db->putInteger("d_mb_nblocks" + ln_string + mb_string, d_mb_nblocks[ln][mb]);
                if (d_mb_nblocks[ln][mb] > 0)
                {
                    std::vector<int> flattened_mb_nelems;
                    flattened_mb_nelems.reserve(NDIM * d_mb_nelems.size());
                    for (std::vector<IntVector<NDIM> >::const_iterator cit = d_mb_nelems[ln][mb].begin();
                         cit != d_mb_nelems[ln][mb].end();
                         ++cit)
                    {
                        flattened_mb_nelems.insert(flattened_mb_nelems.end(), &(*cit)[0], &(*cit)[0] + NDIM);
                    }
                    db->putIntegerArray("flattened_mb_nelems" + ln_string + mb_string,
                                        &flattened_mb_nelems[0],
                                        static_cast<int>(flattened_mb_nelems.size()));

                    std::vector<int> flattened_mb_periodic;
                    flattened_mb_periodic.reserve(NDIM * d_mb_periodic.size());
                    for (std::vector<IntVector<NDIM> >::const_iterator cit = d_mb_periodic[ln][mb].begin();
                         cit != d_mb_periodic[ln][mb].end();
                         ++cit)
                    {
                        flattened_mb_periodic.insert(flattened_mb_periodic.end(), &(*cit)[0], &(*cit)[0] + NDIM);
                    }
                    db->putIntegerArray("flattened_mb_periodic" + ln_string + mb_string,
                                        &flattened_mb_periodic[0],
                                        static_cast<int>(flattened_mb_periodic.size()));

                    db->putIntegerArray("d_mb_first_lag_idx" + ln_string + mb_string,
                                        &d_mb_first_lag_idx[ln][mb][0],
                                        static_cast<int>(d_mb_first_lag_idx[ln][mb].size()));
                }
            }
        }

        db->putInteger("d_nucd_meshes" + ln_string, d_nucd_meshes[ln]);
        if (d_nucd_meshes[ln] > 0)
        {
            db->putStringArray("d_ucd_mesh_names" + ln_string,
                               &d_ucd_mesh_names[ln][0],
                               static_cast<int>(d_ucd_mesh_names[ln].size()));

            for (int mesh = 0; mesh < d_nucd_meshes[ln]; ++mesh)
            {
                std::ostringstream mesh_stream;
                mesh_stream << "_" << mesh;
                const std::string mesh_string = mesh_stream.str();

                std::vector<int> ucd_mesh_vertices_vector;
                ucd_mesh_vertices_vector.reserve(d_ucd_mesh_vertices[ln][mesh].size());
                for (std::set<int>::const_iterator cit = d_ucd_mesh_vertices[ln][mesh].begin();
                     cit != d_ucd_mesh_vertices[ln][mesh].end();
                     ++cit)
                {
                    ucd_mesh_vertices_vector.push_back(*cit);
                }
                db->putInteger("ucd_mesh_vertices_vector.size()" + ln_string + mesh_string,
                               static_cast<int>(ucd_mesh_vertices_vector.size()));
                db->putIntegerArray("ucd_mesh_vertices_vector" + ln_string + mesh_string,
                                    &ucd_mesh_vertices_vector[0],
                                    static_cast<int>(ucd_mesh_vertices_vector.size()));

                std::vector<int> ucd_mesh_edge_maps_vector;
                ucd_mesh_edge_maps_vector.reserve(3 * d_ucd_mesh_edge_maps[ln][mesh].size());
                for (std::multimap<int, std::pair<int, int> >::const_iterator cit =
                         d_ucd_mesh_edge_maps[ln][mesh].begin();
                     cit != d_ucd_mesh_edge_maps[ln][mesh].end();
                     ++cit)
                {
                    const int i = cit->first;
                    std::pair<int, int> e = cit->second;
                    ucd_mesh_edge_maps_vector.push_back(i);
                    ucd_mesh_edge_maps_vector.push_back(e.first);
                    ucd_mesh_edge_maps_vector.push_back(e.second);
                }
                db->putInteger("ucd_mesh_edge_maps_vector.size()" + ln_string + mesh_string,
                               static_cast<int>(ucd_mesh_edge_maps_vector.size()));
                db->putIntegerArray("ucd_mesh_edge_maps_vector" + ln_string + mesh_string,
                                    &ucd_mesh_edge_maps_vector[0],
                                    static_cast<int>(ucd_mesh_edge_maps_vector.size()));
            }
        }
    }
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LSiloDataWriter::buildVecScatters(AO& ao, const int level_number)
{
    if (!d_coords_data[level_number]) return;

    int ierr;

    // Setup the IS data used to generate the VecScatters that redistribute the
    // distributed data into local marker clouds, local logically Cartesian
    // blocks, and local UCD meshes.
    std::vector<int> ref_is_idxs;
    for (int cloud = 0; cloud < d_nclouds[level_number]; ++cloud)
    {
        const int nmarks = d_cloud_nmarks[level_number][cloud];
        const int first_lag_idx = d_cloud_first_lag_idx[level_number][cloud];
        ref_is_idxs.reserve(ref_is_idxs.size() + nmarks);

        for (int idx = first_lag_idx; idx < first_lag_idx + nmarks; ++idx)
        {
            ref_is_idxs.push_back(idx);
        }
    }

    for (int block = 0; block < d_nblocks[level_number]; ++block)
    {
        const IntVector<NDIM>& nelem = d_block_nelems[level_number][block];
        const int ntot = nelem.getProduct();
        const int first_lag_idx = d_block_first_lag_idx[level_number][block];
        ref_is_idxs.reserve(ref_is_idxs.size() + ntot);

        for (int idx = first_lag_idx; idx < first_lag_idx + ntot; ++idx)
        {
            ref_is_idxs.push_back(idx);
        }
    }

    for (int mb = 0; mb < d_nmbs[level_number]; ++mb)
    {
        for (int block = 0; block < d_mb_nblocks[level_number][mb]; ++block)
        {
            const IntVector<NDIM>& nelem = d_mb_nelems[level_number][mb][block];
            const int ntot = nelem.getProduct();
            const int first_lag_idx = d_mb_first_lag_idx[level_number][mb][block];
            ref_is_idxs.reserve(ref_is_idxs.size() + ntot);

            for (int idx = first_lag_idx; idx < first_lag_idx + ntot; ++idx)
            {
                ref_is_idxs.push_back(idx);
            }
        }
    }

    for (int mesh = 0; mesh < d_nucd_meshes[level_number]; ++mesh)
    {
        ref_is_idxs.insert(ref_is_idxs.end(),
                           d_ucd_mesh_vertices[level_number][mesh].begin(),
                           d_ucd_mesh_vertices[level_number][mesh].end());
    }

    // Map Lagrangian indices to PETSc indices.
    std::vector<int> ao_dummy(1, -1);
    ierr = AOApplicationToPetsc(
        ao,
        (!ref_is_idxs.empty() ? static_cast<int>(ref_is_idxs.size()) : static_cast<int>(ao_dummy.size())),
        (!ref_is_idxs.empty() ? &ref_is_idxs[0] : &ao_dummy[0]));
    IBTK_CHKERRQ(ierr);

    // Setup IS indices for all necessary data depths.
    std::map<int, std::vector<int> > src_is_idxs;

    src_is_idxs[NDIM] = ref_is_idxs;
    d_src_vec[level_number][NDIM] = d_coords_data[level_number]->getVec();

    for (int v = 0; v < d_nvars[level_number]; ++v)
    {
        const int var_depth = d_var_depths[level_number][v];
        if (src_is_idxs.find(var_depth) == src_is_idxs.end())
        {
            src_is_idxs[var_depth] = ref_is_idxs;
            d_src_vec[level_number][var_depth] = d_var_data[level_number][v]->getVec();
        }
    }

    // Create the VecScatters to scatter data from the global PETSc Vec to
    // contiguous local subgrids.  VecScatter objects are individually created
    // for data depths as necessary.
    for (std::map<int, std::vector<int> >::iterator it = src_is_idxs.begin(); it != src_is_idxs.end(); ++it)
    {
        const int depth = it->first;
        const std::vector<int>& idxs = it->second;
        const int idxs_sz = static_cast<int>(idxs.size());

        IS src_is;
        ierr = ISCreateBlock(
            PETSC_COMM_WORLD, depth, idxs_sz, (idxs.empty() ? NULL : &idxs[0]), PETSC_COPY_VALUES, &src_is);
        IBTK_CHKERRQ(ierr);

        Vec& src_vec = d_src_vec[level_number][depth];
        Vec& dst_vec = d_dst_vec[level_number][depth];
        if (dst_vec)
        {
            ierr = VecDestroy(&dst_vec);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecCreateMPI(PETSC_COMM_WORLD, depth * idxs_sz, PETSC_DETERMINE, &dst_vec);
        IBTK_CHKERRQ(ierr);

        VecScatter& vec_scatter = d_vec_scatter[level_number][depth];
        if (vec_scatter)
        {
            ierr = VecScatterDestroy(&vec_scatter);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecScatterCreate(src_vec, src_is, dst_vec, NULL, &vec_scatter);
        IBTK_CHKERRQ(ierr);

        ierr = ISDestroy(&src_is);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // buildVecScatters

void
LSiloDataWriter::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("LAG_SILO_DATA_WRITER_VERSION");
    if (ver != LAG_SILO_DATA_WRITER_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                                 << "Restart file version different than class version.");
    }

    const int coarsest_ln = db->getInteger("d_coarsest_ln");
    const int finest_ln = db->getInteger("d_finest_ln");
    resetLevels(coarsest_ln, finest_ln);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        std::ostringstream ln_stream;
        ln_stream << "_" << ln;
        const std::string ln_string = ln_stream.str();

        d_nclouds[ln] = db->getInteger("d_nclouds" + ln_string);
        if (d_nclouds[ln] > 0)
        {
            d_cloud_names[ln].resize(d_nclouds[ln]);
            db->getStringArray(
                "d_cloud_names" + ln_string, &d_cloud_names[ln][0], static_cast<int>(d_cloud_names[ln].size()));

            d_cloud_nmarks[ln].resize(d_nclouds[ln]);
            db->getIntegerArray(
                "d_cloud_nmarks" + ln_string, &d_cloud_nmarks[ln][0], static_cast<int>(d_cloud_nmarks[ln].size()));

            d_cloud_first_lag_idx[ln].resize(d_nclouds[ln]);
            db->getIntegerArray("d_cloud_first_lag_idx" + ln_string,
                                &d_cloud_first_lag_idx[ln][0],
                                static_cast<int>(d_cloud_first_lag_idx[ln].size()));
        }

        d_nblocks[ln] = db->getInteger("d_nblocks" + ln_string);
        if (d_nblocks[ln] > 0)
        {
            d_block_names[ln].resize(d_nblocks[ln]);
            db->getStringArray(
                "d_block_names" + ln_string, &d_block_names[ln][0], static_cast<int>(d_block_names[ln].size()));

            d_block_nelems[ln].resize(d_nblocks[ln]);
            std::vector<int> flattened_block_nelems;
            flattened_block_nelems.resize(NDIM * d_block_nelems[ln].size());
            db->getIntegerArray("flattened_block_nelems" + ln_string,
                                &flattened_block_nelems[0],
                                static_cast<int>(flattened_block_nelems.size()));
            for (unsigned int l = 0; l < d_block_nelems[ln].size(); ++l)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    d_block_nelems[ln][l](d) = flattened_block_nelems[NDIM * l + d];
                }
            }

            d_block_periodic[ln].resize(d_nblocks[ln]);
            std::vector<int> flattened_block_periodic;
            flattened_block_periodic.resize(NDIM * d_block_periodic[ln].size());
            db->getIntegerArray("flattened_block_periodic" + ln_string,
                                &flattened_block_periodic[0],
                                static_cast<int>(flattened_block_periodic.size()));
            for (unsigned int l = 0; l < d_block_periodic[ln].size(); ++l)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    d_block_periodic[ln][l](d) = flattened_block_periodic[NDIM * l + d];
                }
            }

            d_block_first_lag_idx[ln].resize(d_nblocks[ln]);
            db->getIntegerArray("d_block_first_lag_idx" + ln_string,
                                &d_block_first_lag_idx[ln][0],
                                static_cast<int>(d_block_first_lag_idx[ln].size()));
        }

        d_nmbs[ln] = db->getInteger("d_nmbs" + ln_string);
        if (d_nmbs[ln] > 0)
        {
            d_mb_names[ln].resize(d_nmbs[ln]);
            db->getStringArray("d_mb_names" + ln_string, &d_mb_names[ln][0], static_cast<int>(d_mb_names[ln].size()));

            d_mb_nblocks.resize(d_nmbs[ln]);
            d_mb_nelems.resize(d_nmbs[ln]);
            d_mb_periodic.resize(d_nmbs[ln]);
            d_mb_first_lag_idx.resize(d_nmbs[ln]);
            for (int mb = 0; mb < d_nmbs[ln]; ++mb)
            {
                std::ostringstream mb_stream;
                mb_stream << "_" << mb;
                const std::string mb_string = mb_stream.str();

                d_mb_nblocks[ln][mb] = db->getInteger("d_mb_nblocks" + ln_string + mb_string);
                if (d_mb_nblocks[ln][mb] > 0)
                {
                    d_mb_nelems[ln][mb].resize(d_mb_nblocks[ln][mb]);
                    std::vector<int> flattened_mb_nelems;
                    flattened_mb_nelems.resize(NDIM * d_mb_nelems[ln][mb].size());
                    db->getIntegerArray("flattened_mb_nelems" + ln_string + mb_string,
                                        &flattened_mb_nelems[0],
                                        static_cast<int>(flattened_mb_nelems.size()));
                    for (unsigned int l = 0; l < d_mb_nelems[ln][mb].size(); ++l)
                    {
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            d_mb_nelems[ln][mb][l](d) = flattened_mb_nelems[NDIM * l + d];
                        }
                    }

                    d_mb_periodic[ln][mb].resize(d_mb_nblocks[ln][mb]);
                    std::vector<int> flattened_mb_periodic;
                    flattened_mb_periodic.resize(NDIM * d_mb_periodic[ln][mb].size());
                    db->getIntegerArray("flattened_mb_periodic" + ln_string + mb_string,
                                        &flattened_mb_periodic[0],
                                        static_cast<int>(flattened_mb_periodic.size()));
                    for (unsigned int l = 0; l < d_mb_periodic[ln][mb].size(); ++l)
                    {
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            d_mb_periodic[ln][mb][l](d) = flattened_mb_periodic[NDIM * l + d];
                        }
                    }

                    d_mb_first_lag_idx[ln][mb].resize(d_mb_nblocks[ln][mb]);
                    db->getIntegerArray("d_mb_first_lag_idx" + ln_string + mb_string,
                                        &d_mb_first_lag_idx[ln][mb][0],
                                        static_cast<int>(d_mb_first_lag_idx[ln][mb].size()));
                }
            }
        }

        d_nucd_meshes[ln] = db->getInteger("d_nucd_meshes" + ln_string);
        if (d_nucd_meshes[ln] > 0)
        {
            d_ucd_mesh_names[ln].resize(d_nucd_meshes[ln]);
            db->getStringArray("d_ucd_mesh_names" + ln_string,
                               &d_ucd_mesh_names[ln][0],
                               static_cast<int>(d_ucd_mesh_names[ln].size()));

            d_ucd_mesh_vertices[ln].resize(d_nucd_meshes[ln]);
            d_ucd_mesh_edge_maps[ln].resize(d_nucd_meshes[ln]);
            for (int mesh = 0; mesh < d_nucd_meshes[ln]; ++mesh)
            {
                std::ostringstream mesh_stream;
                mesh_stream << "_" << mesh;
                const std::string mesh_string = mesh_stream.str();

                const int ucd_mesh_vertices_vector_size =
                    db->getInteger("ucd_mesh_vertices_vector.size()" + ln_string + mesh_string);
                std::vector<int> ucd_mesh_vertices_vector(ucd_mesh_vertices_vector_size);
                db->getIntegerArray("ucd_mesh_vertices_vector" + ln_string + mesh_string,
                                    &ucd_mesh_vertices_vector[0],
                                    static_cast<int>(ucd_mesh_vertices_vector.size()));
                d_ucd_mesh_vertices[ln][mesh].insert(ucd_mesh_vertices_vector.begin(), ucd_mesh_vertices_vector.end());

                const int ucd_mesh_edge_maps_vector_size =
                    db->getInteger("ucd_mesh_edge_maps_vector.size()" + ln_string + mesh_string);
                std::vector<int> ucd_mesh_edge_maps_vector(ucd_mesh_edge_maps_vector_size);
                db->getIntegerArray("ucd_mesh_edge_maps_vector" + ln_string + mesh_string,
                                    &ucd_mesh_edge_maps_vector[0],
                                    static_cast<int>(ucd_mesh_edge_maps_vector.size()));
                for (int l = 0; l < ucd_mesh_edge_maps_vector_size / 3; ++l)
                {
                    const int idx1 = ucd_mesh_edge_maps_vector[3 * l];
                    const std::pair<int, int> e(ucd_mesh_edge_maps_vector[3 * l + 1],
                                                ucd_mesh_edge_maps_vector[3 * l + 2]);
#if !defined(NDEBUG)
                    TBOX_ASSERT(idx1 == e.first);
#endif
                    d_ucd_mesh_edge_maps[ln][mesh].insert(std::make_pair(idx1, e));
                }
            }
        }
    }
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
