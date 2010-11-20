// Filename: IBHDF5Initializer.h
// Created on 26 Sep 2006 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_IBHDF5Initializer
#define included_IBHDF5Initializer

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/IBBeamForceSpec.h>
#include <ibamr/IBInstrumentationSpec.h>
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/IBTargetPointForceSpec.h>

// IBTK INCLUDES
#include <ibtk/LNodeInitStrategy.h>

// IBTK INCLUDES
#include <ibtk/Streamable.h>

// C++ STDLIB INCLUDES
#include <map>
#include <set>
#include <vector>

// HDF5 INCLUDES
#include <hdf5.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBHDF5Initializer is a concrete LNodeInitStrategy that
 * initializes the configuration of one or more Lagrangian structures from HDF5
 * input files.
 */
class IBHDF5Initializer
    : public IBTK::LNodeInitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBHDF5Initializer(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual
    ~IBHDF5Initializer();

    /*!
     * \brief Determine whether there are any Lagrangian nodes on the specified
     * patch level.
     *
     * \return A boolean value indicating whether Lagrangian data is associated
     * with the given level in the patch hierarchy.
     */
    virtual bool
    getLevelHasLagrangianData(
        const int level_number,
        const bool can_be_refined) const;

    /*!
     * \brief Determine the number of local nodes on the specified patch level.
     *
     * \return The number of local nodes on the specified level.
     */
    virtual int
    getLocalNodeCountOnPatchLevel(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time);

    /*!
     * \brief Initialize the LNodeIndex and LNodeLevel data needed to specify
     * the configuration of the curvilinear mesh on the patch level.
     *
     * \return The number of local nodes initialized on the patch level.
     */
    virtual int
    initializeDataOnPatchLevel(
        const int lag_node_index_idx,
        const int global_index_offset,
        const int local_index_offset,
        SAMRAI::tbox::Pointer<IBTK::LNodeLevelData>& X_data,
        SAMRAI::tbox::Pointer<IBTK::LNodeLevelData>& U_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        IBTK::LDataManager* const lag_manager);

    /*!
     * \brief Initialize the LNodeLevel data needed to specify the mass and
     * spring constant data required by the penalty IB method.
     *
     * \return The number of local nodes initialized on the patch level.
     */
    virtual int
    initializeMassDataOnPatchLevel(
        const int global_index_offset,
        const int local_index_offset,
        SAMRAI::tbox::Pointer<IBTK::LNodeLevelData>& M_data,
        SAMRAI::tbox::Pointer<IBTK::LNodeLevelData>& K_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        IBTK::LDataManager* const lag_manager);

    /*!
     * \brief Tag cells for initial refinement.
     *
     * When the patch hierarchy is being constructed at the initial simulation
     * time, it is necessary to instruct the gridding algorithm where to place
     * local refinement in order to accommodate portions of the curvilinear mesh
     * that will reside in any yet-to-be-constructed level(s) of the patch
     * hierarchy.
     */
    virtual void
    tagCellsForInitialRefinement(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index);

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBHDF5Initializer();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBHDF5Initializer(
        const IBHDF5Initializer& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBHDF5Initializer&
    operator=(
        const IBHDF5Initializer& that);

    /*!
     * \brief Compute the cell indices of all local vertices assigned to the
     * specified level of the patch hierarchy.
     */
    void
    findLocalPatchIndices(
        std::vector<SAMRAI::hier::Index<NDIM> >& cell_idxs,
        std::vector<int>& patch_nums,
        const int level_number,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level) const;

    /*!
     * \brief Compute the cell indices of all local vertices provided by the
     * specified HDF5 file.
     */
    void
    findLocalPatchIndicesFromHDF5(
        std::vector<SAMRAI::hier::Index<NDIM> >& cell_idxs,
        std::vector<int>& patch_nums,
        const hid_t file_id,
        const std::string& base_group_name,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
        const std::string& filename) const;

    /*!
     * \brief Cache all local data that is associated with the specified level
     * of the patch hierarchy.
     */
    void
    buildLevelDataCache(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time);

    /*!
     * \brief Cache all local data provided by the specified HDF5 file that is
     * associated with the specified level of the patch hierarchy.
     */
    void
    buildLevelVertexDataCacheFromHDF5(
        int& num_vertex,
        int& num_local_vertex,
        std::vector<std::vector<double> >& posns,
        std::vector<std::pair<int,int> >& vertex_idxs,
        std::vector<SAMRAI::hier::Index<NDIM> >& cell_idxs,
        std::vector<int>& patch_nums,
        std::set<int>& local_vertex_idx_set,
        const hid_t file_id,
        const std::string& base_group_name,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
        const std::string& filename,
        const int file_number,
        const int num_files) const;

    void
    buildLevelSpringDataCacheFromHDF5(
        int& num_spring,
        int& num_local_spring,
        std::map<int,SAMRAI::tbox::Pointer<IBSpringForceSpec> >& spring_data_map,
        const std::set<int>& local_vertex_idx_set,
        const hid_t file_id,
        const std::string& base_group_name,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
        const std::string& filename,
        const int file_number,
        const int num_files) const;

    void
    buildLevelBeamDataCacheFromHDF5(
        int& num_beam,
        int& num_local_beam,
        std::map<int,SAMRAI::tbox::Pointer<IBBeamForceSpec> >& beam_data_map,
        const std::set<int>& local_vertex_idx_set,
        const hid_t file_id,
        const std::string& base_group_name,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
        const std::string& filename,
        const int file_number,
        const int num_files) const;

    void
    buildLevelTargetPointDataCacheFromHDF5(
        int& num_target_point,
        int& num_local_target_point,
        std::map<int,SAMRAI::tbox::Pointer<IBTargetPointForceSpec> >& target_point_data_map,
        const std::set<int>& local_vertex_idx_set,
        const hid_t file_id,
        const std::string& base_group_name,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
        const std::string& filename,
        const int file_number,
        const int num_files) const;

    void
    buildLevelInstrumentationDataCacheFromHDF5(
        std::vector<std::string>& instrument_names,
        int& num_inst_point,
        int& num_local_inst_point,
        std::map<int,SAMRAI::tbox::Pointer<IBInstrumentationSpec> >& inst_point_data_map,
        const std::set<int>& local_vertex_idx_set,
        const hid_t file_id,
        const std::string& base_group_name,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
        const std::string& filename,
        const int file_number,
        const int num_files) const;

    /*!
     * \brief Clear all cached level data.
     */
    void
    clearLevelDataCache();

    /*!
     * \return The canonical Lagrangian index of the specified vertex.
     */
    int
    getCanonicalLagrangianIndex(
        const std::pair<int,int>& global_vertex_idx,
        const int global_index_offset) const;

    /*!
     * \return The specification objects associated with the specified vertex.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::Streamable> >
    initializeSpecs(
        const std::pair<int,int>& local_vertex_idx,
        const std::pair<int,int>& global_vertex_idx,
        const int global_index_offset);

    /*!
     * Read input values, indicated above, from given database.
     *
     * When assertion checking is active, the database pointer must be non-null.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    std::string d_object_name;

    /*
     * The boolean value determines whether file read batons are employed to
     * prevent multiple MPI processes from accessing the same input files
     * simultaneously.
     */
    bool d_use_file_batons;

    /*
     * The maximum number of levels in the Cartesian grid patch hierarchy and a
     * vector of boolean values indicating whether a particular level has been
     * initialized yet.
     */
    int d_max_levels;
    std::vector<bool> d_level_is_initialized;

    /*
     * The filenames of the structures.
     */
    std::vector<std::vector<std::string> > d_filenames;

    /*
     * Spring information.
     */
    std::vector<std::vector<bool> >   d_enable_springs;
    std::vector<std::vector<bool> >   d_using_uniform_spring_stiffness;
    std::vector<std::vector<double> > d_uniform_spring_stiffness;
    std::vector<std::vector<bool> >   d_using_uniform_spring_rest_length;
    std::vector<std::vector<double> > d_uniform_spring_rest_length;
    std::vector<std::vector<bool> >   d_using_uniform_spring_force_fcn_idx;
    std::vector<std::vector<int> >    d_uniform_spring_force_fcn_idx;

    /*
     * Beam information.
     */
    std::vector<std::vector<bool> >   d_enable_beams;
    std::vector<std::vector<bool> >   d_using_uniform_beam_bend_rigidity;
    std::vector<std::vector<double> > d_uniform_beam_bend_rigidity;

    /*
     * Target point information.
     */
    std::vector<std::vector<bool> >   d_enable_target_points;
    std::vector<std::vector<bool> >   d_using_uniform_target_stiffness;
    std::vector<std::vector<double> > d_uniform_target_stiffness;
    std::vector<std::vector<bool> >   d_using_uniform_target_damping;
    std::vector<std::vector<double> > d_uniform_target_damping;

    /*
     * Instrumentation information.
     */
    std::vector<std::vector<bool> > d_enable_instrumentation;
    std::vector<std::vector<std::vector<std::string> > > d_instrument_names;

    /*
     * Cached level data.
     */
    int d_cache_level_number;

    std::vector<int> d_level_num_vertex, d_level_num_local_vertex, d_level_vertex_offset;
    std::vector<std::vector<std::vector<double> > >       d_level_posns;
    std::vector<std::vector<std::pair<int,int> > >        d_level_vertex_idxs;
    std::vector<std::vector<SAMRAI::hier::Index<NDIM> > > d_level_cell_idxs;
    std::vector<std::vector<int> >                        d_level_patch_nums;

    std::vector<std::set<int> > d_level_reset_specs_set;

    std::vector<int> d_level_num_spring, d_level_num_local_spring;
    std::vector<std::map<int,SAMRAI::tbox::Pointer<IBSpringForceSpec> > > d_level_spring_data_map;

    std::vector<int> d_level_num_beam, d_level_num_local_beam;
    std::vector<std::map<int,SAMRAI::tbox::Pointer<IBBeamForceSpec> > > d_level_beam_data_map;

    std::vector<int> d_level_num_target_point, d_level_num_local_target_point;
    std::vector<std::map<int,SAMRAI::tbox::Pointer<IBTargetPointForceSpec> > > d_level_target_point_data_map;

    std::vector<int> d_level_num_inst_point, d_level_num_local_inst_point;
    std::vector<std::map<int,SAMRAI::tbox::Pointer<IBInstrumentationSpec> > > d_level_inst_point_data_map;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "IBHDF5Initializer.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBHDF5Initializer
