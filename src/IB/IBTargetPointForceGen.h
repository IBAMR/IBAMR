#ifndef included_IBTargetPointForceGen
#define included_IBTargetPointForceGen

// Filename: IBTargetPointForceGen.h
// Last modified: <21.Mar.2007 23:47:44 griffith@box221.cims.nyu.edu>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/IBLagrangianForceStrategy.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{

/*!
 * \brief Class IBTargetPointForceGen computes the penalty forces
 * associated with a collection of target points.
 */
class IBTargetPointForceGen
    : public IBLagrangianForceStrategy
{
public:
    /*!
     * \brief Default constructor.
     */
    IBTargetPointForceGen(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db=NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBTargetPointForceGen();

    /*!
     * \brief Compute the penalty forces determined with the present
     * configuration of the Lagrangian mesh.
     *
     * \note Nodal forces computed by this method are \em added to the
     * force vector.
     */
    virtual void computeLagrangianForce(
        SAMRAI::tbox::Pointer<LNodeLevelData> F_data,
        SAMRAI::tbox::Pointer<LNodeLevelData> X_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double data_time,
        const LDataManager* const lag_manager);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    IBTargetPointForceGen(
        const IBTargetPointForceGen& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBTargetPointForceGen& operator=(
        const IBTargetPointForceGen& that);

    /*!
     * \brief Read input values, indicated above, from given database.
     *
     * The database pointer may be null.
     */
    void getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBTargetPointForceGen.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTargetPointForceGen
