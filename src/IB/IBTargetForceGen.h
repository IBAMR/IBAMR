#ifndef included_IBTargetForceGen
#define included_IBTargetForceGen

// Filename: IBTargetForceGen.h
// Last modified: <21.Mar.2007 20:13:26 griffith@box221.cims.nyu.edu>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/IBLagrangianForceStrategy.h>
#include <ibamr/LDataManager.h>
#include <ibamr/LNodeLevelData.h>

// SAMRAI INCLUDES
#include <PatchHierarchy.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

// PETSc INCLUDES
#include <petscmat.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{

/*!
 * \brief Class IBTargetForceGen computes the penalty forces
 * associated with a collection of target points.
 */
class IBTargetForceGen
    : public IBLagrangianForceStrategy
{
public:
    /*!
     * \brief Default constructor.
     */
    IBTargetForceGen(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db=NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBTargetForceGen();

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
    IBTargetForceGen(
        const IBTargetForceGen& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBTargetForceGen& operator=(
        const IBTargetForceGen& that);

    /*!
     * \brief Read input values, indicated above, from given database.
     *
     * When assertion checking is active, the database pointer must be
     * non-null.
     */
    void getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBTargetForceGen.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTargetForceGen
