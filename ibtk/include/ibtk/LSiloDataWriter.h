// Filename: LSiloDataWriter.h
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

#ifndef included_LSiloDataWriter
#define included_LSiloDataWriter

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "petscao.h"
#include "petscvec.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

namespace IBTK
{
class LData;
} // namespace IBTK
namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LSiloDataWriter provides functionality to output Lagrangian
 * data for visualization via the <A HREF="http://www.llnl.gov/visit">VisIt
 * visualization tool</A> in the Silo data format.
 *
 * For more information about Silo, see the Silo manual <A
 * HREF="http://www.llnl.gov/bdiv/meshtv/manuals/silo.pdf">here</A>.
 */
class LSiloDataWriter : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name           String used for error reporting.
     * \param dump_directory_name   String indicating the directory where visualization data is
     *to
     *be written.
     * \param register_for_restart  Boolean indicating whether to register this object with the
     *restart manager.
     */
    LSiloDataWriter(const std::string& object_name,
                    const std::string& dump_directory_name,
                    bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    ~LSiloDataWriter();

    /*!
     * \name Methods to set the hierarchy and range of levels.
     */
    //\{

    /*!
     * \brief Reset the patch hierarchy over which operations occur.
     */
    void setPatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Reset range of patch levels over which operations occur.
     */
    void resetLevels(int coarsest_ln, int finest_ln);

    //\}

    /*!
     * \brief Register or update a range of Lagrangian indices that are to be
     * visualized as a cloud of marker particles.
     *
     * \note This method is not collective over all MPI processes.  A particular
     * cloud of markers must be registered on only \em one MPI process.
     */
    void registerMarkerCloud(const std::string& name, int nmarks, int first_lag_idx, int level_number);

    /*!
     * \brief Register or update a range of Lagrangian indices that are to be
     * treated as a logically Cartesian block.
     *
     * \note This method is not collective over all MPI processes.  A particular
     * block of indices must be registered on only \em one MPI process.
     */
    void registerLogicallyCartesianBlock(const std::string& name,
                                         const SAMRAI::hier::IntVector<NDIM>& nelem,
                                         const SAMRAI::hier::IntVector<NDIM>& periodic,
                                         int first_lag_idx,
                                         int level_number);

    /*!
     * \brief Register or update several ranges of Lagrangian indices that are
     * to be treated as logically Cartesian blocks.
     *
     * \note This method is not collective over all MPI processes.  A particular
     * block of indices must be registered on only \em one MPI process.
     */
    void registerLogicallyCartesianMultiblock(const std::string& name,
                                              const std::vector<SAMRAI::hier::IntVector<NDIM> >& nelem,
                                              const std::vector<SAMRAI::hier::IntVector<NDIM> >& periodic,
                                              const std::vector<int>& first_lag_idx,
                                              int level_number);

    /*!
     * \brief Register or update an unstructured mesh.
     *
     * \note This method is not collective over all MPI processes.  A particular
     * collection of indices must be registered on only \em one MPI process.
     */
    void registerUnstructuredMesh(const std::string& name,
                                  const std::multimap<int, std::pair<int, int> >& edge_map,
                                  int level_number);

    /*!
     * \brief Register the coordinates of the curvilinear mesh with the Silo
     * data writer.
     */
    void registerCoordsData(SAMRAI::tbox::Pointer<LData> coords_data, int level_number);

    /*!
     * \brief Register a variable for plotting with the Silo data writer.
     */
    void registerVariableData(const std::string& var_name, SAMRAI::tbox::Pointer<LData> var_data, int level_number);

    /*!
     * \brief Register a variable for plotting with the Silo data writer with a
     * specified starting depth and data depth.
     */
    void registerVariableData(const std::string& var_name,
                              SAMRAI::tbox::Pointer<LData> var_data,
                              int start_depth,
                              int var_depth,
                              int level_number);

    /*!
     * \brief Register or update a single Lagrangian AO (application ordering)
     * objects with the Silo data writer.
     *
     * These AO objects are used to map between (fixed) Lagrangian indices and
     * (time-dependent) PETSc indices.  Each time that the AO objects are reset
     * (e.g., during adaptive regridding), the new AO objects must be supplied
     * to the Silo data writer.
     */
    void registerLagrangianAO(AO& ao, int level_number);

    /*!
     * \brief Register or update a collection of Lagrangian AO (application
     * ordering) objects with the Silo data writer.
     *
     * These AO objects are used to map between (fixed) Lagrangian indices and
     * (time-dependent) PETSc indices.  Each time that the AO objects are reset
     * (e.g., during adaptive regridding), the new AO objects must be supplied
     * to the Silo data writer.
     */
    void registerLagrangianAO(std::vector<AO>& ao, int coarsest_ln, int finest_ln);

    /*!
     * \brief Write the plot data to disk.
     */
    void writePlotData(int time_step_number, double simulation_time);

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database pointer must be non-null.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LSiloDataWriter();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LSiloDataWriter(const LSiloDataWriter& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LSiloDataWriter& operator=(const LSiloDataWriter& that);

    /*!
     * \brief Build the VecScatter objects required to communicate data for
     * plotting.
     */
    void buildVecScatters(AO& ao, int level_number);

    /*!
     * Read object state from the restart file and initialize class data
     * members.  The database from which the restart data is read is determined
     * by the object_name specified in the constructor.
     *
     * Unrecoverable Errors:
     *
     *    -   The database corresponding to object_name is not found in the
     *        restart file.
     *
     *    -   The class version number and restart version number do not match.
     *
     */
    void getFromRestart();

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.  The boolean is used to control restart
     * file writing operations.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*
     * The directory where data is to be dumped.
     */
    std::string d_dump_directory_name;

    /*
     * Time step number (passed in by user).
     */
    int d_time_step_number;

    /*
     * Grid hierarchy information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    /*
     * Information about the indices in the local marker clouds.
     */
    std::vector<int> d_nclouds;
    std::vector<std::vector<std::string> > d_cloud_names;
    std::vector<std::vector<int> > d_cloud_nmarks, d_cloud_first_lag_idx;

    /*
     * Information about the indices in the logically Cartesian subgrids.
     */
    std::vector<int> d_nblocks;
    std::vector<std::vector<std::string> > d_block_names;
    std::vector<std::vector<SAMRAI::hier::IntVector<NDIM> > > d_block_nelems;
    std::vector<std::vector<SAMRAI::hier::IntVector<NDIM> > > d_block_periodic;
    std::vector<std::vector<int> > d_block_first_lag_idx;

    /*
     * Information about the indices in the logically Cartesian multiblock
     * subgrids.
     */
    std::vector<int> d_nmbs;
    std::vector<std::vector<std::string> > d_mb_names;
    std::vector<std::vector<int> > d_mb_nblocks;
    std::vector<std::vector<std::vector<SAMRAI::hier::IntVector<NDIM> > > > d_mb_nelems;
    std::vector<std::vector<std::vector<SAMRAI::hier::IntVector<NDIM> > > > d_mb_periodic;
    std::vector<std::vector<std::vector<int> > > d_mb_first_lag_idx;

    /*
     * Information about the indices in the unstructured meshes.
     */
    std::vector<int> d_nucd_meshes;
    std::vector<std::vector<std::string> > d_ucd_mesh_names;
    std::vector<std::vector<std::set<int> > > d_ucd_mesh_vertices;
    std::vector<std::vector<std::multimap<int, std::pair<int, int> > > > d_ucd_mesh_edge_maps;

    /*
     * Coordinates and variable data for plotting.
     */
    std::vector<SAMRAI::tbox::Pointer<LData> > d_coords_data;

    std::vector<int> d_nvars;
    std::vector<std::vector<std::string> > d_var_names;
    std::vector<std::vector<int> > d_var_start_depths, d_var_plot_depths, d_var_depths;
    std::vector<std::vector<SAMRAI::tbox::Pointer<LData> > > d_var_data;

    /*
     * Data for obtaining local data.
     */
    std::vector<AO> d_ao;
    std::vector<bool> d_build_vec_scatters;
    std::vector<std::map<int, Vec> > d_src_vec, d_dst_vec;
    std::vector<std::map<int, VecScatter> > d_vec_scatter;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LSiloDataWriter
