//
// LNodeJacobianStrategy.h
//
// Created on 29 May 2006
//         by Boyce Griffith (boyce@boyce.cims.nyu.edu).
//
// Last modified: <09.Jun.2006 15:44:15 boyce@boyce.cims.nyu.edu>
//

#ifndef included_LNodeJacobianStrategy
#define included_LNodeJacobianStrategy

// SAMRAI-tools INCLUDES
//
#include "LNodeLevelData.h"

// SAMRAI INCLUDES
//
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "PatchHierarchy.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

using namespace SAMRAI;
using namespace std;

// CLASS DEFINITION
//

/*!
 * @brief Class LNodeJacobianStrategy provides a mechanism for
 * specifying Jacobian determinants on the curvilinear mesh.
 */
class LNodeJacobianStrategy
    : public virtual tbox::DescribedClass
{
public:
    /*!
     * @brief Default constructor.
     */
    LNodeJacobianStrategy();
    
    /*!
     * @brief Destructor.
     */
    virtual ~LNodeJacobianStrategy();

    /*!
     * @brief Initialize the Jacobian determinant for the given
     * configuration of the curvilinear mesh.
     */
    virtual void initializeJacobianDet(
        const int lag_node_index_idx,
        tbox::Pointer<LNodeLevelData> J_data,
        tbox::Pointer<LNodeLevelData> X_data,
        const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time) = 0;

private:
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    LNodeJacobianStrategy(
        const LNodeJacobianStrategy& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     * 
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    LNodeJacobianStrategy& operator=(
        const LNodeJacobianStrategy& that);
};

// INLINED FUNCTION DEFINITIONS
//
//#ifndef DEBUG_NO_INLINE
//#include "LNodeJacobianStrategy.I"
//#endif

#endif //#ifndef included_LNodeJacobianStrategy

//////////////////////////////////////////////////////////////////////////////
