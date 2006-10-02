#ifndef included_LNodeJacobianInitStrategy
#define included_LNodeJacobianInitStrategy

// Filename: LNodeJacobianInitStrategy.h
// Created on 29 May 2006 by Boyce Griffith (boyce@boyce.cims.nyu.edu)
// Last modified: <02.Oct.2006 14:27:11 boyce@boyce-griffiths-powerbook-g4-15.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodeLevelData.h>

// SAMRAI INCLUDES
#include <PatchHierarchy.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * @brief Class LNodeJacobianInitStrategy provides a mechanism for
 * specifying Jacobian determinants on the curvilinear mesh.
 */
class LNodeJacobianInitStrategy
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * @brief Default constructor.
     */
    LNodeJacobianInitStrategy();

    /*!
     * @brief Destructor.
     */
    virtual ~LNodeJacobianInitStrategy();

    /*!
     * @brief Initialize the Jacobian determinant for the given
     * configuration of the curvilinear mesh.
     */
    virtual void initializeJacobianDet(
        const int lag_node_index_idx,
        SAMRAI::tbox::Pointer<LNodeLevelData> J_data,
        SAMRAI::tbox::Pointer<LNodeLevelData> X_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
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
    LNodeJacobianInitStrategy(
        const LNodeJacobianInitStrategy& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    LNodeJacobianInitStrategy& operator=(
        const LNodeJacobianInitStrategy& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "LNodeJacobianInitStrategy.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeJacobianInitStrategy
