// Filename: LM3DDataWriter.h
// Created on 11 Sep 2007 by Boyce Griffith
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

#ifndef included_LM3DDataWriter
#define included_LM3DDataWriter

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LData.h>

// SAMRAI INCLUDES
#include <IntVector.h>
#include <PatchHierarchy.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

// PETSc INCLUDES
#include <petscvec.h>
#include <petscao.h>

// C++ STDLIB INCLUDES
#include <map>
#include <set>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LM3DDataWriter provides functionality to output Lagrangian
 * data for visualization via myocardial3D (referred to as m3D throughout the
 * class documentation).
 */
class LM3DDataWriter
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name          String used for error reporting.
     * \param dump_directory_name  String indicating the directory where visualization data is to be written.
     * \param experiment_name      String indicating the three letter "experiment name" used by myocardial3D.  (default="amr")
     * \param experiment_number    Integer indicating the "experiment number" used by myocardial3D.  (default=000)
     *
     * \note The value for \p experiment_name must be a three-letter string, and
     * the value for \p experiment_num must be an integer between 0 and 9999.
     */
    LM3DDataWriter(
        const std::string& object_name,
        const std::string& dump_directory_name,
        const std::string& experiment_name="amr",
        const int& experiment_number=0);

    /*!
     * \brief Destructor.
     */
    ~LM3DDataWriter();

    /*!
     * \name Methods to set the hierarchy and range of levels.
     */
    //\{

    /*!
     * \brief Reset the patch hierarchy over which operations occur.
     */
    void
    setPatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Reset range of patch levels over which operations occur.
     */
    void
    resetLevels(
        const int coarsest_ln,
        const int finest_ln);

    //\}

    /*!
     * \brief Register the patch data descriptor index corresponding to the
     * SAMRAI::pdat::IndexData containing the marker data to be visualized.
     */
    void
    registerLMarkerPatchDataIndex(
        const int mark_idx);

    /*!
     * \brief Register a range of marker indices that are to be visualized as a
     * cloud of marker particles.
     *
     * \note This method is not collective over all MPI processes.  A particular
     * cloud of markers must be registered on only \em one MPI process.
     */
    void
    registerMarkerCloud(
        const std::string& name,
        const int nmarks,
        const int first_marker_idx);

    /*!
     * \brief Register a range of Lagrangian indices that are to be treated as a
     * logically Cartesian block.
     *
     * \note This method is not collective over all MPI processes.  A particular
     * block of indices must be registered on only \em one MPI process.
     */
    void
    registerLogicallyCartesianBlock(
        const std::string& name,
        const SAMRAI::hier::IntVector<NDIM>& nelem,
        const SAMRAI::hier::IntVector<NDIM>& periodic,
        const int first_lag_idx,
        const int level_number);

    /*!
     * \brief Register the coordinates of the curvilinear mesh with the m3D data
     * writer.
     */
    void
    registerCoordsData(
        SAMRAI::tbox::Pointer<LData> coords_data,
        const int level_number);

    /*!
     * \brief Register a single Lagrangian AO (application ordering) objects
     * with the m3D data writer.
     *
     * These AO objects are used to map between (fixed) Lagrangian indices and
     * (time-dependent) PETSc indices.  Each time that the AO objects are reset
     * (e.g., during adaptive regridding), the new AO objects must be supplied
     * to the m3D data writer.
     */
    void
    registerLagrangianAO(
        AO& ao,
        const int level_number);

    /*!
     * \brief Register a collection of Lagrangian AO (application ordering)
     * objects with the m3D data writer.
     *
     * These AO objects are used to map between (fixed) Lagrangian indices and
     * (time-dependent) PETSc indices.  Each time that the AO objects are reset
     * (e.g., during adaptive regridding), the new AO objects must be supplied
     * to the m3D data writer.
     */
    void
    registerLagrangianAO(
        std::vector<AO>& ao,
        const int coarsest_ln,
        const int finest_ln);

    /*!
     * \brief Write the plot data to disk.
     */
    void
    writePlotData(
        const int time_step_number,
        const double simulation_time);

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LM3DDataWriter();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LM3DDataWriter(
        const LM3DDataWriter& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LM3DDataWriter&
    operator=(
        const LM3DDataWriter& that);

    /*!
     * \brief Build the VecScatter objects required to communicate data for
     * plotting.
     */
    void
    buildVecScatters(
        AO& ao,
        const int level_number);

    /*!
     * \brief Return the base-name of the marker output file corresponding to
     * the given timestep number.
     */
    std::string
    getMarkerFileName(
        const int& time_step_number) const;

    /*!
     * \brief Return the base-name of the fiber output file corresponding to the
     * given timestep number.
     */
    std::string
    getFiberFileName(
        const int& time_step_number) const;

    /*!
     * \brief Return the name of the menu file.
     */
    std::string
    getMenuFileName() const;

    /*!
     * \brief Return the name of the list file.
     */
    std::string
    getListFileName() const;

    /*!
     * \brief Return the name of the cat script file.
     */
    std::string
    getCatScriptFileName() const;

    /*
     * The object name is used for error reporting purposes.
     */
    std::string d_object_name;

    /*
     * The directory where data is to be dumped.
     */
    std::string d_dump_directory_name;

    /*
     * The experiment name, number, and corresponding filename prefix.
     */
    const std::string d_experiment_name;
    const int d_experiment_number;
    std::string d_file_prefix;

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
    int d_mark_idx, d_nclouds;
    std::vector<std::string> d_cloud_names;
    std::vector<int> d_cloud_nmarks, d_cloud_first_mark_idx;

    /*
     * Information about the indices in the logically Cartesian subgrids.
     */
    std::vector<int> d_nblocks;
    std::vector<std::vector<std::string> > d_block_names;
    std::vector<std::vector<SAMRAI::hier::IntVector<NDIM> > > d_block_nelems;
    std::vector<std::vector<SAMRAI::hier::IntVector<NDIM> > > d_block_periodic;
    std::vector<std::vector<int> > d_block_nfibers, d_block_ngroups, d_block_first_lag_idx;

    /*
     * Coordinate data for plotting.
     */
    std::vector<SAMRAI::tbox::Pointer<LData> > d_coords_data;

    /*
     * Data for obtaining local data.
     */
    std::vector<AO> d_ao;
    std::vector<bool> d_build_vec_scatters;
    std::vector<std::map<int,Vec> > d_src_vec, d_dst_vec;
    std::vector<std::map<int,VecScatter> > d_vec_scatter;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/LM3DDataWriter.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LM3DDataWriter
