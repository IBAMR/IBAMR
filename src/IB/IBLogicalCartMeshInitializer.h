#ifndef included_IBLogicalCartMeshInitializer
#define included_IBLogicalCartMeshInitializer

// Filename: IBLogicalCartMeshInitializer.h
// Last modified: <03.Oct.2006 13:16:37 boyce@boyce-griffiths-powerbook-g4-15.local>
// Created on 06 Dec 2005 by Boyce Griffith (boyce@boyce.cims.nyu.edu).

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodePosnInitStrategy.h>

// SAMRAI INCLUDES
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <map>
#include <set>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * @brief A concrete LNodePosnInitStrategy for initializing Lagrangian
 * meshes that are comprised of one or more logically Cartesian
 * components.
 */
class IBLogicalCartMeshInitializer
    : public LNodePosnInitStrategy
{
public:
    /*!
     * @brief Constructor.
     */
    IBLogicalCartMeshInitializer(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * @brief Destructor.
     */
    ~IBLogicalCartMeshInitializer();

    /*!
     * @return A boolean value indicating whether Lagrangian data is
     * associated with the given level in the patch hierarchy.
     */
    bool getLevelHasLagrangianData(
        const int level_number,
        const bool can_be_refined) const;

    /*!
     * @return The number of local nodes on the patch level.
     */
    int getLocalNodeCountOnPatchLevel(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time);

    /*!
     * @brief Initialize the LNodeIndex and LNodeLevel data on the
     * patch level.
     */
    void initializeDataOnPatchLevel(
        const int lag_node_index_idx,
        SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time);

    /*!
     * @brief Provide cell tagging for the initial configuration of
     * the Lagrangian mesh.
     *
     * When the patch hierarchy is being constructed at the initial
     * simulation time, it is necessary that the gridding algorithm be
     * instructed where to place local refinement in order to
     * accomodate portions of the curvilinear mesh that will reside in
     * the yet-to-be-constructed level(s) of the patch hierarchy.
     */
    void tagCellsForInitialRefinement(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index);

protected:

private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    IBLogicalCartMeshInitializer();

    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    IBLogicalCartMeshInitializer(
        const IBLogicalCartMeshInitializer& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    IBLogicalCartMeshInitializer& operator=(
        const IBLogicalCartMeshInitializer& that);

    /*!
     * @brief Count the number of nodes in the specified component
     * that are initially within the local portion of the specified
     * PatchLevel.
     */
    int getMeshLocalNodeCount(
        const std::string& mesh_name,
        const int comp_num,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& patch_level);

    /*!
     * @brief Initialize the data corresponding to the portions of the
     * specified component that are initially within the local portion
     * of the specified PatchLevel.
     */
    void initializeLocalData(
        const std::string& mesh_name,
        const int comp_num,
        const int lag_node_index_idx,
        SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& patch_level);

    /*!
     * @brief Tag cells for refinement whenever they contain local
     * mesh nodes that are assigned to finer levels of the patch
     * hierarchy than the present level.
     */
    void tagLocalCellsForRefinement(
        const std::string& mesh_name,
        const int comp_num,
        const int tag_index,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& patch_level);

    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*
     * Data to define the logically Cartesian components.
     */
    std::set<std::string> d_mesh_name;
    std::map<std::string,std::string> d_mesh_file;
    std::map<std::string,int> d_ncomp, d_mesh_offset;
    std::map<std::string,std::vector<int> > d_mesh_rank, d_mesh_level;
    std::map<std::string,std::vector<std::vector<int> > > d_mesh_dim;
    std::map<std::string,std::vector<std::vector<bool> > > d_mesh_periodic;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "IBLogicalCartMeshInitializer.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBLogicalCartMeshInitializer
