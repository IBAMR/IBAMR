#ifndef included_PhysicalBCDataStrategy
#define included_PhysicalBCDataStrategy

// Filename: PhysicalBCDataStrategy.h
// Last modified: <04.Oct.2006 19:49:35 boyce@boyce-griffiths-powerbook-g4-15.local>
// Created on 15 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <Variable.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <ostream>
#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Interface to allow the specification of physical boundary
 * conditions for a single quantity.
 *
 * The interface specified by this class is very similar to that of
 * the SetDataStrategy class.  The different interfaces are provided
 * for the sake of code clarity, i.e., to separate classes which set
 * internal patch data from those which set boundary data.
 * Additionally, note that the interfaces required of concrete
 * implementations of both strategy class are different.
 *
 * \see SetDataStrategy
 */
class PhysicalBCDataStrategy
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief The constructor sets the name of the boundary conditions
     * strategy.
     */
    PhysicalBCDataStrategy(
        const std::string& object_name);

    /*!
     * \brief Empty virtual destructor.
     */
    virtual ~PhysicalBCDataStrategy();

    /*!
     *  \name Methods to set patch boundary data.
     */
    //\{

    /*!
     * \brief Set data in ghost cells corresponding to physical
     * boundary conditions on the specified levels of the patch
     * hierarchy using the virtual function
     * setPhysicalBoundaryConditionsOnPatch().
     *
     * \see setPhysicalBoundaryConditionsOnPatch
     */
    void setPhysicalBoundaryConditionsOnPatchHierarchy(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const double fill_time,
        const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill,
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * \brief Set data in ghost cells corresponding to physical
     * boundary conditions on the specified level of the patch
     * hierarchy using the virtual function
     * setPhysicalBoundaryConditionsOnPatch().
     *
     * \see setPhysicalBoundaryConditionsOnPatch
     */
    void setPhysicalBoundaryConditionsOnPatchLevel(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
        const double fill_time,
        const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill);

    /*!
     * \brief Pure virtual function to set data in ghost cells
     * corresponding to physical boundary conditions.
     */
    virtual void setPhysicalBoundaryConditionsOnPatch(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::hier::Patch<NDIM>& patch,
        const double fill_time,
        const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) = 0;

    //\}

    /*!
     * \brief Print the name of the boundary conditions strategy.
     */
    virtual void printClassData(
        std::ostream& os) const;

private:
    /*!
     * \brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    PhysicalBCDataStrategy();

    /*!
     * \brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    PhysicalBCDataStrategy(
        const PhysicalBCDataStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PhysicalBCDataStrategy& operator=(
        const PhysicalBCDataStrategy& that);

    /*
     * The object name is used for error/warning reporting.
     */
    std::string d_object_name;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/PhysicalBCDataStrategy.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PhysicalBCDataStrategy
