#ifndef included_LagSiloDataWriter
#define included_LagSiloDataWriter

// Filename: LagSiloDataWriter.h
// Created on 26 Apr 2005 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)
// Last modified: <27.Nov.2006 01:46:00 boyce@bigboy.nyconnect.com>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodeLevelData.h>

// SAMRAI INCLUDES
#include <IntVector.h>
#include <PatchHierarchy.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

// PETSC INCLUDES
#include <petscvec.h>
#include <petscao.h>

// C++ STDLIB INCLUDES
#include <map>
#include <set>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LagSiloDataWriter provides functionality to output
 * Lagrangian data for visualization in the Silo data format.
 *
 * For more information about Silo, see the Silo manual <A
 * HREF="http://www.llnl.gov/bdiv/meshtv/manuals/silo.pdf">here</A>.
 * Silo data may be visualized using the <A
 * HREF="http://www.llnl.gov/visit">VisIt visualization tool</A>.
 */
class LagSiloDataWriter
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name std::string used for error reporting
     * purposes
     * \param dump_directory_name std::string used to specify the
     * directory where visualization data is to be written
     */
    LagSiloDataWriter(
        const std::string& object_name,
        const std::string& dump_directory_name);

    /*!
     * \brief Destructor.
     */
    ~LagSiloDataWriter();

    /*!
     * \name Methods to set the hierarchy and range of levels.
     */
    //\{

    /*!
     * \brief Reset the patch hierarchy over which operations occur.
     */
    void setPatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Reset range of patch levels over which operations occur.
     *
     * The specified levels must exist in the hierarchy or an
     * assertion will result.
     */
    void resetLevels(
        const int coarsest_ln,
        const int finest_ln);

    //\}

    /*!
     * \brief Register a range of Lagrangian indices that are to be
     * visualized as a cloud of marker particles.
     *
     * \note This method is not collective over all MPI processes.  A
     * particular cloud of markers must be registered on only \em one
     * MPI process.
     */
    void registerMarkerCloud(
        const std::string& name,
        const int nmarks,
        const int first_lag_idx,
        const int level_number);

    /*!
     * \brief Register a range of Lagrangian indices that are to be
     * treated as a logically Cartesian block.
     *
     * \note This method is not collective over all MPI processes.  A
     * particular block of indices must be registered on only \em one
     * MPI process.
     */
    void registerLogicallyCartesianBlock(
        const std::string& name,
        const SAMRAI::hier::IntVector<NDIM>& nelem,
        const SAMRAI::hier::IntVector<NDIM>& periodic,
        const int first_lag_idx,
        const int level_number);

    /*!
     * \brief Register several ranges of Lagrangian indices that are
     * to be treated as logically Cartesian blocks.
     *
     * \note This method is not collective over all MPI processes.  A
     * particular block of indices must be registered on only \em one
     * MPI process.
     */
    void registerLogicallyCartesianMultiblock(
        const std::string& name,
        const std::vector<SAMRAI::hier::IntVector<NDIM> >& nelem,
        const std::vector<SAMRAI::hier::IntVector<NDIM> >& periodic,
        const std::vector<int>& first_lag_idx,
        const int level_number);

    /*!
     * \brief Register an unstructured mesh.
     *
     * \note This method is not collective over all MPI processes.  A
     * particular collection of indices must be registered on only \em
     * one MPI process.
     */
    void registerUnstructuredMesh(
        const std::string& name,
        const std::multimap<int,std::pair<int,int> > edge_map,
        const int level_number);

    /*!
     * \brief Register the coordinates of the curvilinear mesh with
     * the Silo data writer.
     */
    void registerCoordsData(
        SAMRAI::tbox::Pointer<LNodeLevelData> coords_data,
        const int level_number);

    /*!
     * \brief Register a variable for plotting with the Silo data
     * writer.
     */
    void registerVariableData(
        const std::string& var_name,
        SAMRAI::tbox::Pointer<LNodeLevelData> var_data,
        const int level_number);

    /*!
     * \brief Register a collection of Lagrangian AO (application
     * ordering) objects with the Silo data writer.
     *
     * These AO objects are used to map between (fixed) Lagrangian
     * indices and (time-dependent) PETSc indices.  Each time that the
     * AO objects are reset (e.g., during adaptive regridding), the
     * new AO objects must be supplied to the Silo data writer.
     */
    void setLagrangianAO(
        std::vector<AO>& ao,
        const int coarsest_ln,
        const int finest_ln);

    /*!
     * \brief Write the plot data to disk.
     */
    void writePlotData(
        const int time_step_number,
        const double simulation_time);

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    LagSiloDataWriter();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    LagSiloDataWriter(
        const LagSiloDataWriter& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LagSiloDataWriter& operator=(
        const LagSiloDataWriter& that);

    /*
     * The object name is used for error reporting purposes.
     */
    std::string d_object_name;

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
    std::vector<std::vector<string> > d_cloud_names;
    std::vector<std::vector<int> > d_cloud_nmarks, d_cloud_first_lag_idx;

    /*
     * Information about the indices in the logically Cartesian
     * subgrids.
     */
    std::vector<int> d_nblocks;
    std::vector<std::vector<string> > d_block_names;
    std::vector<std::vector<SAMRAI::hier::IntVector<NDIM> > > d_block_nelems;
    std::vector<std::vector<SAMRAI::hier::IntVector<NDIM> > > d_block_periodic;
    std::vector<std::vector<int> > d_block_first_lag_idx;

    /*
     * Information about the indices in the logically Cartesian
     * multiblock subgrids.
     */
    std::vector<int> d_nmbs;
    std::vector<std::vector<string> > d_mb_names;
    std::vector<std::vector<int> > d_mb_nblocks;
    std::vector<std::vector<std::vector<SAMRAI::hier::IntVector<NDIM> > > > d_mb_nelems;
    std::vector<std::vector<std::vector<SAMRAI::hier::IntVector<NDIM> > > > d_mb_periodic;
    std::vector<std::vector<std::vector<int> > > d_mb_first_lag_idx;

    /*
     * Information about the indices in the unstructured meshes.
     */
    std::vector<int> d_nucd_meshes;
    std::vector<std::vector<string> > d_ucd_mesh_names;
    std::vector<std::vector<std::set<int> > > d_ucd_mesh_vertices;
    std::vector<std::vector<std::multimap<int,std::pair<int,int> > > > d_ucd_mesh_edge_maps;

    /*
     * Coordinates and variable data for plotting.
     */
    std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > d_coords_data;

    std::vector<int> d_nvars;
    std::vector<std::vector<string > > d_var_names;
    std::vector<std::vector<int> > d_var_depths;
    std::vector<std::vector<SAMRAI::tbox::Pointer<LNodeLevelData> > > d_var_data;

    /*
     * Data for obtaining local data.
     */
    std::vector<std::map<int,Vec> > d_src_vec, d_dst_vec;
    std::vector<std::map<int,VecScatter> > d_vec_scatter;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/LagSiloDataWriter.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LagSiloDataWriter
