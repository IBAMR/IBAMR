//
// LagSiloDataWriter.h
//
// Created on 26 Apr 2005
//         by Boyce Griffith (boyce@mstu1.cims.nyu.edu).
//
// Last modified: <25.Jun.2005 00:26:35 boyce@mstu1.cims.nyu.edu>
//

#ifndef included_LagSiloDataWriter
#define included_LagSiloDataWriter

// STL INCLUDES
//
#include <map>
#include <vector>

// SILO INCLUDES
//
extern "C" {
#include <silo.h>
}

// SAMRAI-tools INCLUDES
//
#include "LNodeLevelData.h"

// SAMRAI INCLUDES
//
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

using namespace SAMRAI;
using namespace std;

//
//  Include "petscvec.h" so that we can use PETSc Vecs.  Note that
//  this file automatically includes:
//  
//     petsc.h       - base PETSc routines   petscis.h     - index sets
//     petscsys.h    - system routines       petscviewer.h - viewers
//
// Include "petscao.h" allows use of the AO (application ordering)
// object.
//
#include "petscvec.h"
#include "petscao.h"

// CLASS DEFINITION
//

/*!
 * @brief Description of class.
 */
class LagSiloDataWriter
    : public tbox::DescribedClass
{
public:
    /*!
     * @brief Constructor.
     */
    LagSiloDataWriter(
        const string& object_name,
        const string& dump_directory_name);
    
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
        const string& name,
        const int nmarks,
        const int first_lag_idx,
        const int level_number);

    void registerLogicallyCartesianBlock(
        const string& name,
        const hier::IntVector<NDIM>& nelem,
        const hier::IntVector<NDIM>& periodic,
        const int first_lag_idx,
        const int level_number);

    void registerLogicallyCartesianMultiblock(
        const string& name,
        const vector<hier::IntVector<NDIM> >& nelem,
        const vector<hier::IntVector<NDIM> >& periodic,
        const vector<int>& first_lag_idx,
        const int level_number);

    void registerCoordsData(
        tbox::Pointer<LNodeLevelData> coords_data,
        const int level_number);

    void registerVariableData(
        const string& var_name,
        tbox::Pointer<LNodeLevelData> var_data,
        const int level_number);

    void setLagrangianAO(
        vector<AO>& ao,
        const int coarsest_ln,
        const int finest_ln);

    void writePlotData(
        const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
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

    /*!
     * @brief Build a local mesh database entry corresponding to a
     * cloud of marker points.
     */
    static void buildLocalMarkerCloud(
        DBfile* dbfile,
        string& dirname,
        const int nmarks,
        const double* const X,
        const int time_step,
        const double simulation_time);

    /*!
     * @brief Build a local mesh database entry corresponding to a
     * quadrilateral curvilinear block.
     */
    static void buildLocalCurvBlock(
        DBfile* dbfile,
        string& dirname,
        const hier::IntVector<NDIM>& nelem,
        const hier::IntVector<NDIM>& periodic,
        const double* const X,
        const int nvars,
        const vector<string>& varnames,
        const vector<int>& vardepths,
        const vector<const double*> varvals,
        const int time_step,
        const double simulation_time);

    /*
     * The object name is used for error reporting purposes.
     */
    string d_object_name;

    /*
     * The directory where data is to be dumped.
     */
    string d_dump_directory_name;
    
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
    vector<int> d_nclouds;
    vector<vector<string> > d_cloud_names;
    vector<vector<int> > d_cloud_nmarks, d_cloud_first_lag_idx;

    /*
     * Information about the indices in the logically Cartesian
     * subgrids.
     */
    vector<int> d_nblocks;
    vector<vector<string> > d_block_names;
    vector<vector<hier::IntVector<NDIM> > > d_block_nelems;
    vector<vector<hier::IntVector<NDIM> > > d_block_periodic;
    vector<vector<int> > d_block_first_lag_idx;

    /*
     * Information about the indices in the logically Cartesian
     * multiblock subgrids.
     */
    vector<int> d_nmbs;
    vector<vector<string> > d_mb_names;
    vector<vector<int> > d_mb_nblocks;
    vector<vector<vector<hier::IntVector<NDIM> > > > d_mb_nelems;
    vector<vector<vector<hier::IntVector<NDIM> > > > d_mb_periodic;
    vector<vector<vector<int> > > d_mb_first_lag_idx;

    /*
     * Coordinates and variable data for plotting.
     */
    vector<tbox::Pointer<LNodeLevelData> > d_coords_data;

    vector<int> d_nvars;
    vector<vector<string > > d_var_names;
    vector<vector<int> > d_var_depths;
    vector<vector<tbox::Pointer<LNodeLevelData> > > d_var_data;

    /*
     * Data for obtaining local data.
     */
    vector<map<int,Vec> > d_src_vec, d_dst_vec;
    vector<map<int,VecScatter> > d_vec_scatter;
};

// INLINED FUNCTION DEFINITIONS
//
//#ifndef DEBUG_NO_INLINE
//#include "LagSiloDataWriter.I"
//#endif

#endif //#ifndef included_LagSiloDataWriter

//////////////////////////////////////////////////////////////////////////////
