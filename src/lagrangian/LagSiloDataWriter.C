// Filename: LagSiloDataWriter.C
// Created on 26 Apr 2005 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)
// Last modified: <02.Oct.2006 14:44:45 boyce@boyce-griffiths-powerbook-g4-15.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "LagSiloDataWriter.h"

// IBAMR INCLUDES
#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#endif

// STOOLS INCLUDES
#include <stools/PETSC_SAMRAI_ERROR.h>

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

#include <tbox/MPI.h>
#include <tbox/Utilities.h>

// SILO INCLUDES
#if HAVE_LIBSILO
extern "C" {
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
// The rank of the root MPI process and the MPI tag number.
static const int SILO_MPI_ROOT = 0;
static const int SILO_MPI_TAG  = 0;

// The name of the Silo dumps and database filenames.
static const int SILO_NAME_BUFSIZE = 128;
static const string SILO_DUMPS_FILENAME = "lag_dumps.visit";
static const string SILO_DB_FILENAME = "lag_data.summary.pdb";

#if HAVE_LIBSILO
/*!
 * @brief Build a local mesh database entry corresponding to a cloud
 * of marker points.
 */
void buildLocalMarkerCloud(
    DBfile* dbfile,
    string& dirname,
    const int nmarks,
    const double* const X,
    const int time_step,
    const double simulation_time)
{
    vector<float> block_X(NDIM*nmarks);

    for (int i = 0; i < nmarks; ++i)
    {
        // Get the coordinate data.
        for (int d = 0; d < NDIM; ++d)
        {
            block_X[d*nmarks+i] = static_cast<float>(X[NDIM*i + d]);
        }
    }

    // Set the working directory in the Silo database.
    if (DBSetDir(dbfile, dirname.c_str()) == -1)
    {
        TBOX_ERROR("LagSiloDataWriter::buildLocalMarkerCloud()\n"
                   << "  Could not set directory " << dirname << endl);
    }

    // Write out the variables.
    int    cycle = time_step;
    float  time  = static_cast<float>(simulation_time);
    double dtime = simulation_time;

    static const int MAX_OPTS = 3;
    DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    DBAddOption(optlist, DBOPT_TIME , &time);
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

    const char* meshname = "mesh";
    vector<float*> coords(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        coords[d] = &block_X[d*nmarks];
    }

    int ndims = NDIM;

    DBPutPointmesh(dbfile, meshname, ndims, &coords[0], nmarks,
                   DB_FLOAT, optlist);

    DBFreeOptlist(optlist);

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        TBOX_ERROR("LagSiloDataWriter::buildLocalMarkerCloud()\n"
                   << "  Could not return to the base directory from subdirectory " << dirname << endl);
    }
    return;
}// buildLocalMarkerCloud

/*!
 * @brief Build a local mesh database entry corresponding to a
 * quadrilateral curvilinear block.
 */
void buildLocalCurvBlock(
    DBfile* dbfile,
    string& dirname,
    const SAMRAI::hier::IntVector<NDIM>& nelem_in,
    const SAMRAI::hier::IntVector<NDIM>& periodic,
    const double* const X,
    const int nvars,
    const vector<string>& varnames,
    const vector<int>& vardepths,
    const vector<const double*> varvals,
    const int time_step,
    const double simulation_time)
{
    // Check for codimension 1 or 2 data.
    SAMRAI::hier::IntVector<NDIM> nelem, degenerate;
    for (int d = 0; d < NDIM; ++d)
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
    const int ntot = 1
        *(periodic(0) ? nelem(0)+1 : nelem(0))
#if (NDIM > 1)
        *(periodic(1) ? nelem(1)+1 : nelem(1))
#if (NDIM > 2)
        *(periodic(2) ? nelem(2)+1 : nelem(2))
#endif
#endif
        ;

    vector<float> block_X(NDIM*ntot);
    vector<vector<float> > block_varvals(nvars);
    for (int v = 0; v < nvars; ++v)
    {
        const int vardepth = vardepths[v];
        block_varvals[v].resize(vardepth*ntot);
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
                const int idx =
                    + (degenerate(0) ? 0 : (i%nelem(0)))
#if (NDIM > 1)
                    + (degenerate(1) ? 0 : (j%nelem(1))*nelem(0))
#if (NDIM > 2)
                    + (degenerate(2) ? 0 : (k%nelem(2))*nelem(1)*nelem(0))
#endif
#endif
                    ;

                // Get the coordinate data.
                for (int d = 0; d < NDIM; ++d)
                {
                    block_X[d*ntot+offset] = static_cast<float>(X[NDIM*idx + d]);
                }

                // Get the variable data.
                for (int v = 0; v < nvars; ++v)
                {
                    const int vardepth = vardepths[v];
                    for (int d = 0; d < vardepth; ++d)
                    {
                        block_varvals[v][d*ntot+offset] = static_cast<float>(varvals[v][vardepth*idx + d]);
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
        TBOX_ERROR("LagSiloDataWriter::buildLocalCurvBlock()\n"
                   << "  Could not set directory " << dirname << endl);
    }

    // Write out the variables.
    int    cycle = time_step;
    float  time  = static_cast<float>(simulation_time);
    double dtime = simulation_time;

    static const int MAX_OPTS = 3;
    DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    DBAddOption(optlist, DBOPT_TIME , &time);
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

    static const int MAX_NDIM = 3;  // We may need to increase this value if
                                    // string theory really takes off.
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(NDIM <= MAX_NDIM);
#endif
    const char* meshname = "mesh";
    char* coordnames[MAX_NDIM] = {
        const_cast<char*>("xcoords"),
        const_cast<char*>("ycoords"),
        const_cast<char*>("zcoords") };
    vector<float*> coords(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        coords[d] = &block_X[d*ntot];
    }

    int ndims = NDIM;
    vector<int> dims(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        dims[d] = nelem(d) + (periodic(d) ? 1 : 0);
    }

    DBPutQuadmesh(dbfile, meshname, coordnames, &coords[0], &dims[0], ndims,
                  DB_FLOAT, DB_NONCOLLINEAR, optlist);

    for (int v = 0; v < nvars; ++v)
    {
        const char* varname = varnames[v].c_str();
        const int vardepth = vardepths[v];
        vector<char*> compnames(vardepth);
        for (int d = 0; d < vardepth; ++d)
        {
            ostringstream stream;
            stream << "_" << d;
            const string compname = varnames[v] + stream.str();
            compnames[d] = strdup(compname.c_str());
        }

        vector<float*> vars(vardepth);
        for (int d = 0; d < vardepth; ++d)
        {
            vars[d] = &block_varvals[v][d*ntot];
        }

        if (vardepth == 1)
        {
            DBPutQuadvar1(dbfile, varname, meshname, vars[0], &dims[0], ndims,
                          NULL, 0, DB_FLOAT, DB_NODECENT, optlist);
        }
        else
        {
            DBPutQuadvar(dbfile, varname, meshname, vardepth, &compnames[0], &vars[0], &dims[0], ndims,
                         NULL, 0, DB_FLOAT, DB_NODECENT, optlist);
        }

        for (int d = 0; d < vardepth; ++d)
        {
            free(compnames[d]);
        }
    }

    DBFreeOptlist(optlist);

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        TBOX_ERROR("LagSiloDataWriter::buildLocalCurvBlock()\n"
                   << "  Could not return to the base directory from subdirectory " << dirname << endl);
    }
    return;
}// buildLocalCurvBlock
#endif //if HAVE_LIBSILO

}

/////////////////////////////// PUBLIC ///////////////////////////////////////

LagSiloDataWriter::LagSiloDataWriter(
    const string& object_name,
    const string& dump_directory_name)
    : d_object_name(object_name),
      d_dump_directory_name(dump_directory_name)
{
    // Initialize some default values.
    d_time_step_number = -1;
    d_coarsest_ln = -1;
    d_finest_ln = -1;

    return;
}// LagSiloDataWriter

LagSiloDataWriter::~LagSiloDataWriter()
{
    // Destroy any remaining PETSc objects.
    int ierr;
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        for (map<int,Vec>::iterator it = d_dst_vec[ln].begin();
             it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = (*it).second;
            if (v)
            {
                ierr = VecDestroy(v);
                PETSC_SAMRAI_ERROR(ierr);
            }
        }
        for (map<int,VecScatter>::iterator it = d_vec_scatter[ln].begin();
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
}// ~LagSiloDataWriter

void
LagSiloDataWriter::resetLevels(
    const int coarsest_ln,
    const int finest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(coarsest_ln >= 0 && finest_ln >= coarsest_ln);
#endif
    // Destroy any un-needed PETSc objects.
    int ierr;
    for (int ln = SAMRAI::tbox::Utilities::imax(d_coarsest_ln,0);
         (ln <= d_finest_ln) && (ln < coarsest_ln); ++ln)
    {
        for (map<int,Vec>::iterator it = d_dst_vec[ln].begin();
             it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = (*it).second;
            if (v)
            {
                ierr = VecDestroy(v);
                PETSC_SAMRAI_ERROR(ierr);
            }
        }
        for (map<int,VecScatter>::iterator it = d_vec_scatter[ln].begin();
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

    for (int ln = finest_ln+1; ln <= d_finest_ln; ++ln)
    {
        for (map<int,Vec>::iterator it = d_dst_vec[ln].begin();
             it != d_dst_vec[ln].end(); ++it)
        {
            Vec& v = (*it).second;
            if (v)
            {
                ierr = VecDestroy(v);
                PETSC_SAMRAI_ERROR(ierr);
            }
        }
        for (map<int,VecScatter>::iterator it = d_vec_scatter[ln].begin();
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

    // Reset the level numbers.
    d_coarsest_ln = coarsest_ln;
    d_finest_ln   = finest_ln;

    // Resize some arrays.
    d_nclouds            .resize(d_finest_ln+1,0);
    d_cloud_names        .resize(d_finest_ln+1);
    d_cloud_nmarks       .resize(d_finest_ln+1);
    d_cloud_first_lag_idx.resize(d_finest_ln+1);

    d_nblocks            .resize(d_finest_ln+1,0);
    d_block_names        .resize(d_finest_ln+1);
    d_block_nelems       .resize(d_finest_ln+1);
    d_block_periodic     .resize(d_finest_ln+1);
    d_block_first_lag_idx.resize(d_finest_ln+1);

    d_nmbs            .resize(d_finest_ln+1,0);
    d_mb_nblocks      .resize(d_finest_ln+1);
    d_mb_names        .resize(d_finest_ln+1);
    d_mb_nelems       .resize(d_finest_ln+1);
    d_mb_periodic     .resize(d_finest_ln+1);
    d_mb_first_lag_idx.resize(d_finest_ln+1);

    d_coords_data.resize(d_finest_ln+1);
    d_nvars      .resize(d_finest_ln+1);
    d_var_names  .resize(d_finest_ln+1);
    d_var_depths .resize(d_finest_ln+1);
    d_var_data   .resize(d_finest_ln+1);

    d_src_vec    .resize(d_finest_ln+1);
    d_dst_vec    .resize(d_finest_ln+1);
    d_vec_scatter.resize(d_finest_ln+1);

    return;
}// resetLevels

void
LagSiloDataWriter::registerMarkerCloud(
    const string& name,
    const int nmarks,
    const int first_lag_idx,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(nmarks > 0);
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
#endif

    // Check to see if the cloud name has already been registered.
    if (find(d_cloud_names[level_number].begin(), d_cloud_names[level_number].end(),
             name) != d_cloud_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerMarkerCloud()\n"
                   << "  marker clouds must have unique names.\n"
                   << "  a marker cloud named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_block_names[level_number].begin(), d_block_names[level_number].end(),
             name) != d_block_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerMarkerCloud()\n"
                   << "  marker clouds must have unique names.\n"
                   << "  a Cartesian block named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_mb_names[level_number].begin(), d_mb_names[level_number].end(),
             name) != d_mb_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerMarkerCloud()\n"
                   << "  marker clouds must have unique names.\n"
                   << "  a Cartesian multiblock named ``" << name << "'' has already been registered.\n");
    }

    // Record the layout of the marker cloud.
    ++d_nclouds[level_number];
    d_cloud_names        [level_number].push_back(name);
    d_cloud_nmarks       [level_number].push_back(nmarks);
    d_cloud_first_lag_idx[level_number].push_back(first_lag_idx);

    return;
}// registerMarkerCloud

void
LagSiloDataWriter::registerLogicallyCartesianBlock(
    const string& name,
    const SAMRAI::hier::IntVector<NDIM>& nelem,
    const SAMRAI::hier::IntVector<NDIM>& periodic,
    const int first_lag_idx,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int d = 0; d < NDIM; ++d)
    {
        assert(nelem(d) > 0);
        assert(periodic(d) == 0 || periodic(d) == 1);
    }
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
#endif

    // Check to see if the block name has already been registered.
    if (find(d_cloud_names[level_number].begin(), d_cloud_names[level_number].end(),
             name) != d_cloud_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                   << "  Cartesian blocks must have unique names.\n"
                   << "  a marker cloud named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_block_names[level_number].begin(), d_block_names[level_number].end(),
             name) != d_block_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                   << "  Cartesian blocks must have unique names.\n"
                   << "  a Cartesian block named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_mb_names[level_number].begin(), d_mb_names[level_number].end(),
             name) != d_mb_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianBlock()\n"
                   << "  Cartesian blocks must have unique names.\n"
                   << "  a Cartesian multiblock named ``" << name << "'' has already been registered.\n");
    }

    // Record the layout of the logically Cartesian block.
    ++d_nblocks[level_number];
    d_block_names        [level_number].push_back(name);
    d_block_nelems       [level_number].push_back(nelem);
    d_block_periodic     [level_number].push_back(periodic);
    d_block_first_lag_idx[level_number].push_back(first_lag_idx);

    return;
}// registerLogicallyCartesianBlock

void
LagSiloDataWriter::registerLogicallyCartesianMultiblock(
    const string& name,
    const vector<SAMRAI::hier::IntVector<NDIM> >& nelem,
    const vector<SAMRAI::hier::IntVector<NDIM> >& periodic,
    const vector<int>& first_lag_idx,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(periodic     .size() == nelem.size());
    assert(first_lag_idx.size() == nelem.size());
    int sz = nelem.size();
    for (int i = 0; i < sz; ++i)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            assert(nelem[i](d) > 0);
            assert(periodic[i](d) == 0 || periodic[i](d) == 1);
        }
    }
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
#endif

    // Check to see if the multiblock name has already been registered.
    if (find(d_cloud_names[level_number].begin(), d_cloud_names[level_number].end(),
             name) != d_cloud_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianMultiblock()\n"
                   << "  Cartesian multiblocks must have unique names.\n"
                   << "  a marker cloud named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_block_names[level_number].begin(), d_block_names[level_number].end(),
             name) != d_block_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianMultiblock()\n"
                   << "  Cartesian multiblocks must have unique names.\n"
                   << "  a Cartesian block named ``" << name << "'' has already been registered.\n");
    }

    if (find(d_mb_names[level_number].begin(), d_mb_names[level_number].end(),
             name) != d_mb_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerLogicallyCartesianMultiblock()\n"
                   << "  Cartesian multiblocks must have unique names.\n"
                   << "  a Cartesian multiblock named ``" << name << "'' has already been registered.\n");
    }

    // Record the layout of the logically Cartesian multiblock.
    ++d_nmbs[level_number];
    d_mb_names        [level_number].push_back(name);
    d_mb_nblocks      [level_number].push_back(nelem.size());
    d_mb_nelems       [level_number].push_back(nelem);
    d_mb_periodic     [level_number].push_back(periodic);
    d_mb_first_lag_idx[level_number].push_back(first_lag_idx);

    return;
}// registerLogicallyCartesianMultiblock

void
LagSiloDataWriter::registerCoordsData(
    SAMRAI::tbox::Pointer<LNodeLevelData> coords_data,
    const int level_number)
{
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
LagSiloDataWriter::registerVariableData(
    const string& var_name,
    SAMRAI::tbox::Pointer<LNodeLevelData> var_data,
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!var_name.empty());
    assert(!var_data.isNull());
    assert(d_coarsest_ln <= level_number &&
           d_finest_ln   >= level_number);
#endif
    if (find(d_var_names[level_number].begin(),
             d_var_names[level_number].end(),
             var_name) != d_var_names[level_number].end())
    {
        TBOX_ERROR(d_object_name << "::registerVariableData()\n"
                   << "  variable with name " << var_name << " already registered for plotting\n"
                   << "  on patch level " << level_number << endl);
    }
    ++d_nvars[level_number];
    d_var_names [level_number].push_back(var_name);
    d_var_depths[level_number].push_back(var_data->getDepth());
    d_var_data  [level_number].push_back(var_data);
    return;
}// registerVariableData

void
LagSiloDataWriter::setLagrangianAO(
    vector<AO>& ao,
    const int coarsest_ln,
    const int finest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(coarsest_ln <= finest_ln);
    assert(d_coarsest_ln <= coarsest_ln && finest_ln <= d_finest_ln);
#endif

    int ierr;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_coords_data[ln].isNull())
        {
            // Setup the IS data used to generate the VecScatters that
            // redistribute the distributed data into local marker
            // clouds and local logically Cartesian blocks.
            vector<int> ref_is_idxs;
            for (int cloud = 0; cloud < d_nclouds[ln]; ++cloud)
            {
                const int nmarks = d_cloud_nmarks[ln][cloud];
                const int first_lag_idx = d_cloud_first_lag_idx[ln][cloud];
                ref_is_idxs.reserve(ref_is_idxs.size()+nmarks);

                for (int idx = first_lag_idx; idx < first_lag_idx + nmarks;
                     ++idx)
                {
                    ref_is_idxs.push_back(idx);
                }
            }

            for (int block = 0; block < d_nblocks[ln]; ++block)
            {
                const SAMRAI::hier::IntVector<NDIM>& nelem = d_block_nelems[ln][block];
                const int ntot = nelem.getProduct();
                const int first_lag_idx = d_block_first_lag_idx[ln][block];
                ref_is_idxs.reserve(ref_is_idxs.size()+ntot);

                for (int idx = first_lag_idx; idx < first_lag_idx + ntot; ++idx)
                {
                    ref_is_idxs.push_back(idx);
                }
            }

            for (int mb = 0; mb < d_nmbs[ln]; ++mb)
            {
                for (int block = 0; block < d_mb_nblocks[ln][mb]; ++block)
                {
                    const SAMRAI::hier::IntVector<NDIM>& nelem = d_mb_nelems[ln][mb][block];
                    const int ntot = nelem.getProduct();
                    const int first_lag_idx = d_mb_first_lag_idx[ln][mb][block];
                    ref_is_idxs.reserve(ref_is_idxs.size()+ntot);

                    for (int idx = first_lag_idx; idx < first_lag_idx + ntot; ++idx)
                    {
                        ref_is_idxs.push_back(idx);
                    }
                }
            }

            // Map Lagrangian indices to PETSc indices.
            vector<int> ao_dummy(1,-1);
            ierr = AOApplicationToPetsc(
                ao[ln],
                (!ref_is_idxs.empty() ? static_cast<int>(ref_is_idxs.size()) : static_cast<int>(ao_dummy.size())),
                (!ref_is_idxs.empty() ? &ref_is_idxs[0]                      : &ao_dummy[0]));
            PETSC_SAMRAI_ERROR(ierr);

            // Setup IS indices for all necessary data depths.
            map<int,vector<int> > src_is_idxs;

            src_is_idxs[NDIM].resize(ref_is_idxs.size());
            transform(ref_is_idxs.begin(), ref_is_idxs.end(),
                      src_is_idxs[NDIM].begin(),
                      bind2nd(multiplies<int>(),NDIM));
            d_src_vec[ln][NDIM] = d_coords_data[ln]->getGlobalVec();

            for (int v = 0; v < d_nvars[ln]; ++v)
            {
                const int var_depth = d_var_depths[ln][v];
                if (src_is_idxs.find(var_depth) == src_is_idxs.end())
                {
                    src_is_idxs[var_depth].resize(ref_is_idxs.size());
                    transform(ref_is_idxs.begin(), ref_is_idxs.end(),
                              src_is_idxs[var_depth].begin(),
                              bind2nd(multiplies<int>(),var_depth));
                    d_src_vec[ln][var_depth] = d_var_data[ln][v]->getGlobalVec();
                }
            }

            // Create the VecScatters to scatter data from the global
            // PETSc Vec to contiguous local subgrids.  VecScatter
            // objects are individually created for data depths as
            // necessary.
            for (map<int,vector<int> >::iterator it = src_is_idxs.begin();
                 it != src_is_idxs.end(); ++it)
            {
                const int depth = (*it).first;
                const vector<int>& idxs = (*it).second;

                IS src_is;
                ierr = ISCreateBlock(PETSC_COMM_WORLD, depth, idxs.size(),
                                     &idxs[0], &src_is);
                PETSC_SAMRAI_ERROR(ierr);

                Vec& src_vec = d_src_vec[ln][depth];
                Vec& dst_vec = d_dst_vec[ln][depth];
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

                VecScatter& vec_scatter = d_vec_scatter[ln][depth];
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
        }
    }

    return;
}// setLagrangianAO

void
LagSiloDataWriter::writePlotData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int time_step_number,
    const double simulation_time)
{
#if HAVE_LIBSILO
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(time_step_number >= 0);
    assert(!d_dump_directory_name.empty());
#endif

    if (time_step_number <= d_time_step_number)
    {
        TBOX_ERROR("LagSiloDataWriter::writePlotData()\n"
                   << "  data writer with name " << d_object_name << "\n"
                   << "  time step number: " << time_step_number
                   << " is <= last time step number: " << d_time_step_number
                   << endl);
    }
    d_time_step_number = time_step_number;

    if (d_dump_directory_name.empty())
    {
        TBOX_ERROR("LagSiloDataWriter::writePlotData()\n"
                   << "  data writer with name " << d_object_name << "\n"
                   << "  dump directory name is empty" << endl);
    }

    int ierr;
    char temp_buf[SILO_NAME_BUFSIZE];
    string current_file_name;
    DBfile* dbfile;
    const int mpi_rank  = SAMRAI::tbox::MPI::getRank();
    const int mpi_nodes = SAMRAI::tbox::MPI::getNodes();

    // Create the working directory.
    sprintf(temp_buf, "%05d", d_time_step_number);
    string current_dump_directory_name = "silo_dump.";
    current_dump_directory_name += temp_buf;
    string dump_dirname = d_dump_directory_name + "/";
    dump_dirname += current_dump_directory_name;

    SAMRAI::tbox::Utilities::recursiveMkdir(dump_dirname);

    // Create one local DBfile per MPI process.
    sprintf(temp_buf, "%05d", mpi_rank);
    current_file_name = dump_dirname + "/" + "lag_data.";
    current_file_name += temp_buf;
    current_file_name += ".pdb";

    if ((dbfile = DBCreate(current_file_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB))
        == NULL)
    {
        TBOX_ERROR(d_object_name + "::writePlotData()\n"
                   << "  Could not create DBfile named " << current_file_name << endl);
    }

    vector<vector<int> > meshtype(d_finest_ln+1), vartype(d_finest_ln+1);
    vector<vector<vector<int> > > multimeshtype(d_finest_ln+1), multivartype(d_finest_ln+1);

    // Set the local data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (!d_coords_data[ln].isNull())
        {
            // Scatter the data from "global" to "local" form.
            Vec local_X_vec;
            ierr = VecDuplicate(d_dst_vec[ln][NDIM], &local_X_vec);
            PETSC_SAMRAI_ERROR(ierr);

            Vec& global_X_vec = d_coords_data[ln]->getGlobalVec();
            ierr = VecScatterBegin(global_X_vec, local_X_vec, INSERT_VALUES,
                                   SCATTER_FORWARD, d_vec_scatter[ln][NDIM]);
            PETSC_SAMRAI_ERROR(ierr);
            ierr = VecScatterEnd(global_X_vec, local_X_vec, INSERT_VALUES,
                                 SCATTER_FORWARD, d_vec_scatter[ln][NDIM]);
            PETSC_SAMRAI_ERROR(ierr);

            double* local_X_arr;
            ierr = VecGetArray(local_X_vec, &local_X_arr);
            PETSC_SAMRAI_ERROR(ierr);

            vector<Vec> local_v_vecs;
            vector<double*> local_v_arrs;

            for (int v = 0; v < d_nvars[ln]; ++v)
            {
                const int var_depth = d_var_depths[ln][v];
                Vec local_v_vec;
                ierr = VecDuplicate(d_dst_vec[ln][var_depth], &local_v_vec);
                PETSC_SAMRAI_ERROR(ierr);

                Vec& global_v_vec = d_var_data[ln][v]->getGlobalVec();
                ierr = VecScatterBegin(global_v_vec, local_v_vec, INSERT_VALUES,
                                       SCATTER_FORWARD, d_vec_scatter[ln][var_depth]);
                PETSC_SAMRAI_ERROR(ierr);
                ierr = VecScatterEnd(global_v_vec, local_v_vec, INSERT_VALUES,
                                     SCATTER_FORWARD, d_vec_scatter[ln][var_depth]);
                PETSC_SAMRAI_ERROR(ierr);

                double* local_v_arr;
                ierr = VecGetArray(local_v_vec, &local_v_arr);
                PETSC_SAMRAI_ERROR(ierr);

                local_v_vecs.push_back(local_v_vec);
                local_v_arrs.push_back(local_v_arr);
            }

            // Keep track of the current offset in the local Vec data.
            int offset = 0;

            // Add the local clouds to the local DBfile.
            for (int cloud = 0; cloud < d_nclouds[ln]; ++cloud)
            {
                const int nmarks = d_cloud_nmarks[ln][cloud];

                ostringstream stream;
                stream << "level_" << ln << "_cloud_" << cloud;
                string dirname = stream.str();

                if (DBMkDir(dbfile, dirname.c_str()) == -1)
                {
                    TBOX_ERROR(d_object_name + "::writePlotData()\n"
                               << "  Could not create directory named "
                               << dirname << endl);
                }

                const double* const X = local_X_arr + NDIM*offset;

                buildLocalMarkerCloud(dbfile, dirname, nmarks, X,
                                      time_step_number, simulation_time);

                offset += nmarks;
            }

            // Add the local blocks to the local DBfile.
            for (int block = 0; block < d_nblocks[ln]; ++block)
            {
                const SAMRAI::hier::IntVector<NDIM>& nelem    = d_block_nelems  [ln][block];
                const SAMRAI::hier::IntVector<NDIM>& periodic = d_block_periodic[ln][block];
                const int ntot = nelem.getProduct();

                ostringstream stream;
                stream << "level_" << ln << "_block_" << block;
                string dirname = stream.str();

                if (DBMkDir(dbfile, dirname.c_str()) == -1)
                {
                    TBOX_ERROR(d_object_name + "::writePlotData()\n"
                               << "  Could not create directory named "
                               << dirname << endl);
                }

                const double* const X = local_X_arr + NDIM*offset;

                vector<const double*> var_vals(d_nvars[ln]);
                for (int v = 0; v < d_nvars[ln]; ++v)
                {
                    var_vals[v] = local_v_arrs[v] + d_var_depths[ln][v]*offset;
                }

                buildLocalCurvBlock(dbfile, dirname, nelem, periodic, X,
                                    d_nvars[ln], d_var_names[ln], d_var_depths[ln], var_vals,
                                    time_step_number, simulation_time);
                meshtype[ln].push_back(DB_QUAD_CURV);
                vartype [ln].push_back(DB_QUADVAR);

                offset += ntot;
            }

            // Add the local multiblocks to the local DBfile.
            multimeshtype[ln].resize(d_nmbs[ln]);
            multivartype [ln].resize(d_nmbs[ln]);
            for (int mb = 0; mb < d_nmbs[ln]; ++mb)
            {
                for (int block = 0; block < d_mb_nblocks[ln][mb]; ++block)
                {
                    const SAMRAI::hier::IntVector<NDIM>& nelem    = d_mb_nelems  [ln][mb][block];
                    const SAMRAI::hier::IntVector<NDIM>& periodic = d_mb_periodic[ln][mb][block];
                    const int ntot = nelem.getProduct();

                    ostringstream stream;
                    stream << "level_" << ln << "_mb_" << mb << "_block_" << block;
                    string dirname = stream.str();

                    if (DBMkDir(dbfile, dirname.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name + "::writePlotData()\n"
                                   << "  Could not create directory named "
                                   << dirname << endl);
                    }

                    const double* const X = local_X_arr + NDIM*offset;

                    vector<const double*> var_vals(d_nvars[ln]);
                    for (int v = 0; v < d_nvars[ln]; ++v)
                    {
                        var_vals[v] = local_v_arrs[v] + d_var_depths[ln][v]*offset;
                    }

                    buildLocalCurvBlock(dbfile, dirname, nelem, periodic, X,
                                        d_nvars[ln], d_var_names[ln], d_var_depths[ln], var_vals,
                                        time_step_number, simulation_time);
                    multimeshtype[ln][mb].push_back(DB_QUAD_CURV);
                    multivartype [ln][mb].push_back(DB_QUADVAR);

                    offset += ntot;
                }
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

    DBClose(dbfile);

    // Send data to the root MPI process required to create the
    // multimesh and multivar objects.
    vector<vector<int> > nclouds_per_proc, nblocks_per_proc, nmbs_per_proc;
    vector<vector<vector<int> > > meshtypes_per_proc, vartypes_per_proc, mb_nblocks_per_proc;
    vector<vector<vector<vector<int> > > > multimeshtypes_per_proc, multivartypes_per_proc;
    vector<vector<vector<string> > > cloud_names_per_proc, block_names_per_proc, mb_names_per_proc;

    if (mpi_rank == SILO_MPI_ROOT)
    {
        nclouds_per_proc       .resize(d_finest_ln+1);
        nblocks_per_proc       .resize(d_finest_ln+1);
        nmbs_per_proc          .resize(d_finest_ln+1);
        meshtypes_per_proc     .resize(d_finest_ln+1);
        vartypes_per_proc      .resize(d_finest_ln+1);
        mb_nblocks_per_proc    .resize(d_finest_ln+1);
        multimeshtypes_per_proc.resize(d_finest_ln+1);
        multivartypes_per_proc .resize(d_finest_ln+1);
        cloud_names_per_proc   .resize(d_finest_ln+1);
        block_names_per_proc   .resize(d_finest_ln+1);
        mb_names_per_proc      .resize(d_finest_ln+1);
    }

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (mpi_rank == SILO_MPI_ROOT)
        {
            nclouds_per_proc       [ln].resize(mpi_nodes);
            nblocks_per_proc       [ln].resize(mpi_nodes);
            nmbs_per_proc          [ln].resize(mpi_nodes);
            meshtypes_per_proc     [ln].resize(mpi_nodes);
            vartypes_per_proc      [ln].resize(mpi_nodes);
            mb_nblocks_per_proc    [ln].resize(mpi_nodes);
            multimeshtypes_per_proc[ln].resize(mpi_nodes);
            multivartypes_per_proc [ln].resize(mpi_nodes);
            cloud_names_per_proc   [ln].resize(mpi_nodes);
            block_names_per_proc   [ln].resize(mpi_nodes);
            mb_names_per_proc      [ln].resize(mpi_nodes);
        }

        int one = 1;
        for (int proc = 0; proc < mpi_nodes; ++proc)
        {
            if (mpi_rank == proc)
            {
                SAMRAI::tbox::MPI::send(&d_nclouds[ln], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                SAMRAI::tbox::MPI::recv(&nclouds_per_proc[ln][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nclouds[ln] > 0)
            {
                int num_bytes;
                for (int cloud = 0; cloud < d_nclouds[ln]; ++cloud)
                {
                    num_bytes = (d_cloud_names[ln][cloud].size()+1)*sizeof(char);
                    SAMRAI::tbox::MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);
                    SAMRAI::tbox::MPI::sendBytes((void*)d_cloud_names[ln][cloud].c_str(), num_bytes, SILO_MPI_ROOT);
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nclouds_per_proc[ln][proc] > 0)
            {
                cloud_names_per_proc[ln][proc].resize(nclouds_per_proc[ln][proc]);

                int num_bytes;
                char* name;
                for (int cloud = 0; cloud < nclouds_per_proc[ln][proc]; ++cloud)
                {
                    SAMRAI::tbox::MPI::recv(&num_bytes, one, proc, false);
                    name = new char[num_bytes/sizeof(char)];
                    SAMRAI::tbox::MPI::recvBytes((void*)name, num_bytes);
                    cloud_names_per_proc[ln][proc][cloud].assign(name);
                    delete[] name;
                }
            }

            if (mpi_rank == proc)
            {
                SAMRAI::tbox::MPI::send(&d_nblocks[ln], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                SAMRAI::tbox::MPI::recv(&nblocks_per_proc[ln][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nblocks[ln] > 0)
            {
                SAMRAI::tbox::MPI::send(&meshtype[ln][0], d_nblocks[ln], SILO_MPI_ROOT, false);
                SAMRAI::tbox::MPI::send(&vartype [ln][0], d_nblocks[ln], SILO_MPI_ROOT, false);

                int num_bytes;
                for (int block = 0; block < d_nblocks[ln]; ++block)
                {
                    num_bytes = (d_block_names[ln][block].size()+1)*sizeof(char);
                    SAMRAI::tbox::MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);
                    SAMRAI::tbox::MPI::sendBytes((void*)d_block_names[ln][block].c_str(), num_bytes, SILO_MPI_ROOT);
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nblocks_per_proc[ln][proc] > 0)
            {
                meshtypes_per_proc  [ln][proc].resize(nblocks_per_proc[ln][proc]);
                vartypes_per_proc   [ln][proc].resize(nblocks_per_proc[ln][proc]);
                block_names_per_proc[ln][proc].resize(nblocks_per_proc[ln][proc]);

                SAMRAI::tbox::MPI::recv(&meshtypes_per_proc[ln][proc][0],
                                nblocks_per_proc   [ln][proc], proc, false);
                SAMRAI::tbox::MPI::recv(&vartypes_per_proc [ln][proc][0],
                                nblocks_per_proc   [ln][proc], proc, false);

                int num_bytes;
                char* name;
                for (int block = 0; block < nblocks_per_proc[ln][proc]; ++block)
                {
                    SAMRAI::tbox::MPI::recv(&num_bytes, one, proc, false);
                    name = new char[num_bytes/sizeof(char)];
                    SAMRAI::tbox::MPI::recvBytes((void*)name, num_bytes);
                    block_names_per_proc[ln][proc][block].assign(name);
                    delete[] name;
                }
            }

            if (mpi_rank == proc)
            {
                SAMRAI::tbox::MPI::send(&d_nmbs[ln], one, SILO_MPI_ROOT, false);
            }
            if (mpi_rank == SILO_MPI_ROOT)
            {
                SAMRAI::tbox::MPI::recv(&nmbs_per_proc[ln][proc], one, proc, false);
            }

            if (mpi_rank == proc && d_nmbs[ln] > 0)
            {
                SAMRAI::tbox::MPI::send(&d_mb_nblocks[ln][0], d_nmbs[ln], SILO_MPI_ROOT, false);

                int num_bytes;
                for (int mb = 0; mb < d_nmbs[ln]; ++mb)
                {
                    SAMRAI::tbox::MPI::send(&multimeshtype[ln][mb][0], d_mb_nblocks[ln][mb], SILO_MPI_ROOT, false);
                    SAMRAI::tbox::MPI::send(&multivartype [ln][mb][0], d_mb_nblocks[ln][mb], SILO_MPI_ROOT, false);

                    num_bytes = d_mb_names[ln][mb].size()+1;
                    SAMRAI::tbox::MPI::send(&num_bytes, one, SILO_MPI_ROOT, false);

                    (void)MPI_Send((void*)d_mb_names[ln][mb].c_str(), num_bytes, MPI_CHAR,
                                   SILO_MPI_ROOT, SILO_MPI_TAG, SAMRAI::tbox::MPI::commWorld);
                    const int tree = SAMRAI::tbox::MPI::getTreeDepth();
                    SAMRAI::tbox::MPI::updateOutgoingStatistics(tree, num_bytes*sizeof(char));
                }
            }
            if (mpi_rank == SILO_MPI_ROOT && nmbs_per_proc[ln][proc] > 0)
            {
                mb_nblocks_per_proc    [ln][proc].resize(nmbs_per_proc[ln][proc]);
                multimeshtypes_per_proc[ln][proc].resize(nmbs_per_proc[ln][proc]);
                multivartypes_per_proc [ln][proc].resize(nmbs_per_proc[ln][proc]);
                mb_names_per_proc      [ln][proc].resize(nmbs_per_proc[ln][proc]);

                SAMRAI::tbox::MPI::recv(&mb_nblocks_per_proc[ln][proc][0],
                                nmbs_per_proc       [ln][proc], proc, false);

                int num_bytes;
                char* name;
                for (int mb = 0; mb < nmbs_per_proc[ln][proc]; ++mb)
                {
                    multimeshtypes_per_proc[ln][proc][mb].resize(mb_nblocks_per_proc[ln][proc][mb]);
                    multivartypes_per_proc [ln][proc][mb].resize(mb_nblocks_per_proc[ln][proc][mb]);

                    SAMRAI::tbox::MPI::recv(&multimeshtypes_per_proc[ln][proc][mb][0],
                                    mb_nblocks_per_proc     [ln][proc][mb], proc, false);
                    SAMRAI::tbox::MPI::recv(&multivartypes_per_proc [ln][proc][mb][0],
                                    mb_nblocks_per_proc     [ln][proc][mb], proc, false);

                    SAMRAI::tbox::MPI::recv(&num_bytes, one, proc, false);
                    name = new char[num_bytes];

                    MPI_Status status;
                    (void)MPI_Recv(name, num_bytes, MPI_CHAR,
                                   proc, SILO_MPI_TAG, SAMRAI::tbox::MPI::commWorld, &status);
                    const int tree = SAMRAI::tbox::MPI::getTreeDepth();
                    SAMRAI::tbox::MPI::updateIncomingStatistics(tree, num_bytes*sizeof(char));

                    mb_names_per_proc[ln][proc][mb].assign(name);
                    delete[] name;
                }
            }
            SAMRAI::tbox::MPI::barrier();
        }
    }

    if (mpi_rank == SILO_MPI_ROOT)
    {
        // Create and initialize the multimesh Silo database on the
        // root MPI process.
        string summary_file_name = dump_dirname + "/" + SILO_DB_FILENAME;
        if ((dbfile = DBCreate(summary_file_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB))
            == NULL)
        {
            TBOX_ERROR(d_object_name + "::writePlotData()\n"
                       << "  Could not create DBfile named " << summary_file_name << endl);
        }

        int    cycle = time_step_number;
        float  time  = static_cast<float>(simulation_time);
        double dtime = simulation_time;

        static const int MAX_OPTS = 3;
        DBoptlist* optlist = DBMakeOptlist(MAX_OPTS);
        DBAddOption(optlist, DBOPT_CYCLE, &cycle);
        DBAddOption(optlist, DBOPT_TIME , &time );
        DBAddOption(optlist, DBOPT_DTIME, &dtime);

        for (int proc = 0; proc < mpi_nodes; ++proc)
        {
            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                for (int cloud = 0; cloud < nclouds_per_proc[ln][proc]; ++cloud)
                {
                    sprintf(temp_buf, "%05d", proc);
                    current_file_name = "lag_data.";
                    current_file_name += temp_buf;
                    current_file_name += ".pdb";

                    ostringstream stream;
                    stream << current_file_name << ":level_" << ln << "_cloud_" << cloud << "/mesh";
                    string meshname = stream.str();
                    char* meshname_ptr = const_cast<char*>(meshname.c_str());
                    int meshtype = DB_POINTMESH;

                    string& cloud_name = cloud_names_per_proc[ln][proc][cloud];

                    DBPutMultimesh(dbfile, cloud_name.c_str(), 1, &meshname_ptr, &meshtype, optlist);

                    if (DBMkDir(dbfile, cloud_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name + "::writePlotData()\n"
                                   << "  Could not create directory named "
                                   << cloud_name << endl);
                    }
                }

                for (int block = 0; block < nblocks_per_proc[ln][proc]; ++block)
                {
                    sprintf(temp_buf, "%05d", proc);
                    current_file_name = "lag_data.";
                    current_file_name += temp_buf;
                    current_file_name += ".pdb";

                    ostringstream stream;
                    stream << current_file_name << ":level_" << ln << "_block_" << block << "/mesh";
                    string meshname = stream.str();
                    char* meshname_ptr = const_cast<char*>(meshname.c_str());
                    int meshtype = meshtypes_per_proc[ln][proc][block];

                    string& block_name = block_names_per_proc[ln][proc][block];

                    DBPutMultimesh(dbfile, block_name.c_str(), 1, &meshname_ptr, &meshtype, optlist);

                    if (DBMkDir(dbfile, block_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name + "::writePlotData()\n"
                                   << "  Could not create directory named "
                                   << block_name << endl);
                    }
                }

                for (int mb = 0; mb < nmbs_per_proc[ln][proc]; ++mb)
                {
                    sprintf(temp_buf, "%05d", proc);
                    current_file_name = "lag_data.";
                    current_file_name += temp_buf;
                    current_file_name += ".pdb";

                    const int nblocks = mb_nblocks_per_proc[ln][proc][mb];
                    char** meshnames = new char*[nblocks];

                    for (int block = 0; block < nblocks; ++block)
                    {
                        ostringstream stream;
                        stream << current_file_name << ":level_" << ln << "_mb_" << mb << "_block_" << block << "/mesh";
                        meshnames[block] = strdup(stream.str().c_str());
                    }

                    string& mb_name = mb_names_per_proc[ln][proc][mb];

                    DBPutMultimesh(dbfile, mb_name.c_str(), nblocks, meshnames, &multimeshtypes_per_proc[ln][proc][mb][0], optlist);

                    if (DBMkDir(dbfile, mb_name.c_str()) == -1)
                    {
                        TBOX_ERROR(d_object_name + "::writePlotData()\n"
                                   << "  Could not create directory named "
                                   << mb_name << endl);
                    }

                    for (int block = 0; block < nblocks; ++block)
                    {
                        free(meshnames[block]);
                    }
                    delete[] meshnames;
                }

                for (int v = 0; v < d_nvars[ln]; ++v)
                {
                    for (int block = 0; block < nblocks_per_proc[ln][proc]; ++block)
                    {
                        sprintf(temp_buf, "%05d", proc);
                        current_file_name = "lag_data.";
                        current_file_name += temp_buf;
                        current_file_name += ".pdb";

                        ostringstream varname_stream;
                        varname_stream << current_file_name << ":level_" << ln << "_block_" << block << "/" << d_var_names[ln][v];
                        string varname = varname_stream.str();
                        char* varname_ptr = const_cast<char*>(varname.c_str());
                        int vartype = vartypes_per_proc[ln][proc][block];

                        string& block_name = block_names_per_proc[ln][proc][block];

                        ostringstream stream;
                        stream << block_name << "/" << d_var_names[ln][v];
                        string var_name = stream.str();

                        DBPutMultivar(dbfile, var_name.c_str(), 1, &varname_ptr, &vartype, optlist);
                    }

                    for (int mb = 0; mb < nmbs_per_proc[ln][proc]; ++mb)
                    {
                        sprintf(temp_buf, "%05d", proc);
                        current_file_name = "lag_data.";
                        current_file_name += temp_buf;
                        current_file_name += ".pdb";

                        const int nblocks = mb_nblocks_per_proc[ln][proc][mb];
                        char** varnames = new char*[nblocks];

                        for (int block = 0; block < nblocks; ++block)
                        {
                            ostringstream varname_stream;
                            varname_stream << current_file_name << ":level_" << ln << "_mb_" << mb << "_block_" << block << d_var_names[ln][v];
                            varnames[block] = strdup(varname_stream.str().c_str());
                        }

                        string& mb_name = mb_names_per_proc[ln][proc][mb];

                        ostringstream stream;
                        stream << mb_name << "/" << d_var_names[ln][v];
                        string var_name = stream.str();

                        DBPutMultivar(dbfile, var_name.c_str(), nblocks, varnames, &multivartypes_per_proc[ln][proc][mb][0], optlist);

                        for (int block = 0; block < nblocks; ++block)
                        {
                            free(varnames[block]);
                        }
                        delete[] varnames;
                    }
                }
            }
        }

        DBClose(dbfile);

        // Create or update the dumps file on the root MPI process.
        static bool summary_file_opened = false;
        string path = d_dump_directory_name + "/" + SILO_DUMPS_FILENAME;
        string file = current_dump_directory_name + "/" + SILO_DB_FILENAME;

        if (!summary_file_opened)
        {
            summary_file_opened = true;
            ofstream sfile(path.c_str(), ios::out);
            sfile << file << "\n";
            sfile.close();
        }
        else
        {
            ofstream sfile(path.c_str(), ios::app);
            sfile << file << "\n";
            sfile.close();
        }
    }
    SAMRAI::tbox::MPI::barrier();
#else
    TBOX_WARNING("LagSiloDataWriter::writePlotData(): SILO is not installed; cannot write data." << endl);
#endif //if HAVE_LIBSILO
    return;
}// writePlotData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LagSiloDataWriter>;

//////////////////////////////////////////////////////////////////////////////
