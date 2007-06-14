#ifndef included_LNodeInitStrategy
#define included_LNodeInitStrategy

// Filename: LNodeInitStrategy.h
// Last modified: <13.Jun.2007 17:48:02 griffith@box221.cims.nyu.edu>
// Created on 11 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodeLevelData.h>

// SAMRAI INCLUDES
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LNodeInitStrategy provides a mechanism for specifying the
 * initial configuration of the curvilinear mesh.
 */
class LNodeInitStrategy
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    LNodeInitStrategy();

    /*!
     * \brief Destructor.
     */
    virtual
    ~LNodeInitStrategy();

    /*!
     * \return A boolean value indicating whether Lagrangian data is associated
     * with the given level in the patch hierarchy.
     */
    virtual bool
    getLevelHasLagrangianData(
        const int level_number,
        const bool can_be_refined) const = 0;

    /*!
     * \return The number of local nodes on the patch level.
     */
    virtual int
    getLocalNodeCountOnPatchLevel(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time) = 0;

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
        SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
        SAMRAI::tbox::Pointer<LNodeLevelData>& U_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        LDataManager* const lag_manager) = 0;

    /*!
     * \brief Initialize the LNodeLevel data needed to specify the mass and
     * spring constant data required by the penalty IB method.
     *
     * \return The number of local nodes initialized on the patch level.
     *
     * \note A default empty implementation is provided when support for massive
     * boundaries is not required.
     */
    virtual int
    initializeMassDataOnPatchLevel(
        const int global_index_offset,
        const int local_index_offset,
        SAMRAI::tbox::Pointer<LNodeLevelData>& M_data,
        SAMRAI::tbox::Pointer<LNodeLevelData>& K_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        LDataManager* const lag_manager);

    /*!
     * \brief Provide cell tagging for the initial configuration of the
     * Lagrangian mesh.
     *
     * When the patch hierarchy is being constructed at the initial simulation
     * time, it is necessary that the gridding algorithm be instructed where to
     * place local refinement in order to accomodate portions of the curvilinear
     * mesh that will reside in the yet-to-be-constructed level(s) of the patch
     * hierarchy.
     *
     * \note A default empty implementation is provided when support for local
     * mesh refinement is not required.
     */
    virtual void
    tagCellsForInitialRefinement(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LNodeInitStrategy(
        const LNodeInitStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeInitStrategy&
    operator=(
        const LNodeInitStrategy& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/LNodeInitStrategy.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeInitStrategy
