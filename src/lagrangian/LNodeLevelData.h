#ifndef included_LNodeLevelData
#define included_LNodeLevelData

// Filename: LNodeLevelData.h
// Last modified: <17.Apr.2007 19:35:44 griffith@box221.cims.nyu.edu>
// Created on 08 Mar 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Serializable.h>

// PETSc INCLUDES
#include <petscvec.h>
#include <petscao.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBAMR
{
class LDataManager;
}// namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LNodeLevelData provides storage for a single variable defined at
 * the nodes of the Lagrangian mesh.
 */
class LNodeLevelData
    : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Destructor.
     */
    ~LNodeLevelData();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeLevelData&
    operator=(
        const LNodeLevelData& that);

    /*!
     * \brief This operator() method provides a non-const refrence to the
     * Lagrangian nodal quantity at local index i at the specified depth.
     *
     * \return A non-const refrence to the Lagrangian nodal quantity at local
     * index i at the specified depth.
     */
    double&
    operator()(
        const int i,
        const int depth=0);

    /*!
     * \return The total number of nodes of the Lagrangian data.
     */
    int
    getGlobalNodeCount();

    /*!
     * \return The number of local (i.e., on processor) nodes of the Lagrangian
     * data.
     *
     * \note This count does not include ghost nodes.
     */
    int
    getLocalNodeCount();

    /*!
     * \return The depth (i.e., number of components per node) of the Lagrangian
     * data.
     */
    int
    getDepth() const;

    /*!
     * \brief Begin updating the ghost nodes.
     */
    void
    beginGhostUpdate();

    /*!
     * \brief End updating the ghost nodes.
     */
    void
    endGhostUpdate();

    /*!
     * \return The global PETSc Vec object which contains the level data.
     */
    Vec&
    getGlobalVec();

    /*!
     * \return The ghosted local form of the global PETSc Vec object which
     * contains the level data.
     */
    Vec&
    getLocalFormVec();

    /*!
     * \brief Restore the ghosted local form of the global PETSc Vec object.
     */
    void
    restoreLocalFormVec();

    /*!
     * \return The array corresponding to the PETSc Vec object in ghosted local
     * form.
     */
    PetscScalar*
    getLocalFormArray();

    /*!
     * \brief Restore the array to the ghosted local form of the global PETSc
     * Vec object.
     */
    void
    restoreLocalFormArray();

    ///
    ///  The following routines:
    ///
    ///      putToDatabase()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::tbox::Serializable abstract base class.
    ///

    /*!
     * Store the contents of the class in a database.
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    friend class LDataManager;

    /*!
     * \brief Constructor.
     */
    LNodeLevelData(
        const std::string& name,
        const int num_local_nodes,
        const int depth,
        const std::vector<int>& nonlocal_petsc_indices=std::vector<int>(0));

    /*!
     * \brief Constructor.
     */
    LNodeLevelData(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    LNodeLevelData(
        const LNodeLevelData& from);

    /*!
     * \brief Reset the PETSc Vec and related data.
     */
    void
    resetData(
        Vec& new_global_vec,
        const std::vector<int>& new_nonlocal_petsc_indices);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LNodeLevelData();

    /*!
     * \brief Assign that to this.
     */
    void
    assignThatToThis(
        const LNodeLevelData& that);

    /*
     * The name of the LNodeLevelData.
     */
    std::string d_name;

    /*
     * The depth (i.e. number of components per node) of the LNodeLevelData.
     */
    int d_depth;

    /*
     * The nonlocal PETSc indices available in the ghosted local form.
     */
    std::vector<int> d_nonlocal_petsc_indices;

    /*
     * The global PETSc Vec object which contains the level data and its ghosted
     * local form.
     */
    Vec d_global_vec, d_local_vec;
    bool d_in_local_form;

    /*
     * The array corresponding to the PETSc Vec object in local form.
     */
    PetscScalar* d_local_vec_array;
    bool d_extracted_local_array;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/LNodeLevelData.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeLevelData
