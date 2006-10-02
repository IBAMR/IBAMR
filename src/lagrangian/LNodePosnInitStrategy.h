//
// LNodePosnInitStrategy.h
//
// Created on 11 Jul 2004
//         by Boyce Griffith (boyce@trasnaform.speakeasy.net).
//
// Last modified: <17.Jun.2005 16:05:41 boyce@bigboy.verizon.net>
//

#ifndef included_LNodePosnInitStrategy
#define included_LNodePosnInitStrategy

// SAMRAI-tools INCLUDES
//
#include "LNodeLevelData.h"

// SAMRAI INCLUDES
//
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

using namespace SAMRAI;
using namespace std;

// CLASS DEFINITION
//

/*!
 * @brief Class LNodePosnInitStrategy provides a mechanism for
 * specifying the initial configuration of the curvilinear mesh.
 */
class LNodePosnInitStrategy
    : public virtual tbox::DescribedClass
{
public:
    /*!
     * @brief Default constructor.
     */
    LNodePosnInitStrategy();
    
    /*!
     * @brief Destructor.
     */
    virtual ~LNodePosnInitStrategy();

    /*!
     * @return A boolean value indicating whether Lagrangian data is
     * associated with the given level in the patch hierarchy.
     */
    virtual bool getLevelHasLagrangianData(
        const int level_number,
        const bool can_be_refined) const = 0;

    /*!
     * @return The number of local nodes on the patch level.
     */
    virtual int getLocalNodeCountOnPatchLevel(
        const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time) = 0;

    /*!
     * @brief Initialize the LNodeIndex and LNodeLevel data needed to
     * specify the configuration of the curvilinear mesh on the patch
     * level.
     */
    virtual void initializeDataOnPatchLevel(
        const int lag_node_index_idx,
        tbox::Pointer<LNodeLevelData>& X_data,
        const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time) = 0;

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
    virtual void tagCellsForInitialRefinement(
        const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
        const int level_number,
        const double error_data_time,
        const int tag_index);
    
private:
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    LNodePosnInitStrategy(
        const LNodePosnInitStrategy& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     * 
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    LNodePosnInitStrategy& operator=(
        const LNodePosnInitStrategy& that);
};

// INLINED FUNCTION DEFINITIONS
//
//#ifndef DEBUG_NO_INLINE
//#include "LNodePosnInitStrategy.I"
//#endif

#endif //#ifndef included_LNodePosnInitStrategy

//////////////////////////////////////////////////////////////////////////////
