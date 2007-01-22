#ifndef included_TargetPointInitializer
#define included_TargetPointInitializer

// Filename: TargetPointInitializer.h
// Last modified: <09.Nov.2006 14:36:08 griffith@box221.cims.nyu.edu>
// Created on 26 Oct 2006 by Boyce Griffith (boyce@bigboy.nyconnect.com)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodePosnInitStrategy.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Concrete LNodePosnInitStrategy to specify the initial
 * positions of a collection of target points.
 *
 * \todo document input database entries
 */
class TargetPointInitializer
    : public LNodePosnInitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    TargetPointInitializer(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual ~TargetPointInitializer();

    /*!
     * \brief Determine whether there are any Lagrangian nodes on the
     * specified patch level.
     *
     * \return A boolean value indicating whether Lagrangian data is
     * associated with the given level in the patch hierarchy.
     */
    virtual bool getLevelHasLagrangianData(
        const int level_number,
        const bool can_be_refined) const;

    /*!
     * \brief Determine the number of local nodes on the specified
     * patch level.
     *
     * \return The number of local nodes on the specified level.
     */
    virtual int getLocalNodeCountOnPatchLevel(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time);

    /*!
     * \brief Initialize the LNodeIndex and LNodeLevel data needed to
     * specify the configuration of the curvilinear mesh on the patch
     * level.
     *
     * \return The number of local nodes initialized on the patch
     * level.
     */
    virtual int initializeDataOnPatchLevel(
        const int lag_node_index_idx,
        const int global_index_offset,
        const int local_index_offset,
        SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time);

    /*!
     * \brief Tag cells for initial refinement.
     *
     * When the patch hierarchy is being constructed at the initial
     * simulation time, it is necessary to instruct the gridding
     * algorithm where to place local refinement in order to
     * accomodate portions of the curvilinear mesh that will reside in
     * any yet-to-be-constructed level(s) of the patch hierarchy.
     */
    virtual void tagCellsForInitialRefinement(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index);

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    TargetPointInitializer();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    TargetPointInitializer(
        const TargetPointInitializer& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    TargetPointInitializer& operator=(
        const TargetPointInitializer& that);

    /*!
     * \brief Determine the indices of the target points initially
     * located within the specified patch.
     */
    void getPatchTargetPointIndices(
        std::vector<int>& point_indices,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
        const int level_number,
        const bool can_be_refined) const;

    /*!
     * \return The initial position of the specified target point.
     */
    std::vector<double> getTargetPointPosn(
        const int point_index) const;

    /*!
     * \return The spring stiffness associated with the specified
     * target point.
     */
    double getTargetPointStiffness(
        const int point_index) const;

    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.
     */
    std::string d_object_name;

    /*
     * The number of target points, their initial locations, and the
     * corresponding spring stiffnesses.
     */
    int d_num_points;
    std::vector<double> d_initial_posns;
    std::vector<double> d_stiffnesses;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "TargetPointInitializer.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_TargetPointInitializer
