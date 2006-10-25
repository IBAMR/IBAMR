#ifndef included_XInit
#define included_XInit

// Filename: XInit.h
// Created on 12 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <25.Oct.2006 12:32:56 boyce@bigboy.nyconnect.com>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodePosnInitStrategy.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <tbox/Database.h>
#include <tbox/DescribedClass.h>
#include <GridGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <tbox/Pointer.h>

// NAMESPACE
using namespace IBAMR;
using namespace SAMRAI;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * @brief Method to initialize the position of the Lagrangian fiber
 * mesh.
 */
class XInit
    : public LNodePosnInitStrategy
{
public:
    /*!
     * @brief Default constructor.
     */
    XInit(
        const string& object_name,
        tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
        tbox::Pointer<tbox::Database> input_db);

    /*!
     * @brief Destructor.
     */
    ~XInit();

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
        const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
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
        tbox::Pointer<LNodeLevelData>& X_data,
        const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
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
     *
     * NOTE: A default empty implementation is provided when support
     * for local mesh refinement is not required.
     */
    void tagCellsForInitialRefinement(
        const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
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
    XInit();

    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    XInit(
        const XInit& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    XInit& operator=(
        const XInit& that);

    /*!
     * Get the location of node l.
     */
    void getNodePosn(
        std::vector<double>& X,
        const int l) const;

    /*!
     * Read input values, indicated above, from given database.
     */
    void getFromInput(
        tbox::Pointer<tbox::Database> db);

    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.
     */
    string d_object_name;

    /*
     * The grid geometry.
     */
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geom;

    /*
     * The number of nodes and some other parameters.
     */
    int d_num_nodes;
    tbox::Array<double> d_center;
    double d_radius, d_stiffness;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "XInit.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_XInit
