#ifndef included_LagSiloDataWriter
#define included_LagSiloDataWriter

// Filename: LagSiloDataWriter.h
// Created on 26 Apr 2005 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)
// Last modified: <04.Oct.2006 19:51:55 boyce@boyce-griffiths-powerbook-g4-15.local>

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
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * @brief Description of class.
 */
class LagSiloDataWriter
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * @brief Constructor.
     */
    LagSiloDataWriter(
        const std::string& object_name,
        const std::string& dump_directory_name);

    /*!
     * @brief Destructor.
     */
    ~LagSiloDataWriter();

    //@{ @name Methods to set the range of levels.

    /*!
     * @brief Reset range of patch levels over which operations occur.
     *
     * The levels must exist in the hierarchy or an assertion failure
     * will result.
     */
    void resetLevels(
        const int coarsest_ln,
        const int finest_ln);

    //@}

    /*!
     * @brief XXXX
     */
    void registerMarkerCloud(
        const std::string& name,
        const int nmarks,
        const int first_lag_idx,
        const int level_number);

    void registerLogicallyCartesianBlock(
        const std::string& name,
        const SAMRAI::hier::IntVector<NDIM>& nelem,
        const SAMRAI::hier::IntVector<NDIM>& periodic,
        const int first_lag_idx,
        const int level_number);

    void registerLogicallyCartesianMultiblock(
        const std::string& name,
        const std::vector<SAMRAI::hier::IntVector<NDIM> >& nelem,
        const std::vector<SAMRAI::hier::IntVector<NDIM> >& periodic,
        const std::vector<int>& first_lag_idx,
        const int level_number);

    void registerCoordsData(
        SAMRAI::tbox::Pointer<LNodeLevelData> coords_data,
        const int level_number);

    void registerVariableData(
        const std::string& var_name,
        SAMRAI::tbox::Pointer<LNodeLevelData> var_data,
        const int level_number);

    void setLagrangianAO(
        std::vector<AO>& ao,
        const int coarsest_ln,
        const int finest_ln);

    void writePlotData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int time_step_number,
        const double simulation_time);

protected:

private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    LagSiloDataWriter();

    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    LagSiloDataWriter(
        const LagSiloDataWriter& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
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
