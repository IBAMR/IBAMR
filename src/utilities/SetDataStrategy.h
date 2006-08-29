#ifndef included_SetDataStrategy
#define included_SetDataStrategy

// Filename: SetDataStrategy.h
// Last modified: <28.Aug.2006 21:39:20 boyce@bigboy.nyconnect.com>
// Created on 15 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
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
 * \brief Interface to allow the specification of data for a single
 * quantity at a specified time.
 *
 * The interface specified by this class is very similar to that of
 * the PhysicalBCDataStrategy class.  The different interfaces are
 * provided for the sake of code clarity, i.e., to separate classes
 * which set internal patch data from those which set boundary data.
 * Additionally, note that the interfaces required of concrete
 * implementations of both strategy class are different.
 *
 * \see PhysicalBCDataStrategy
 */
class SetDataStrategy
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief The constructor sets the name of the strategy object.
     */
    SetDataStrategy(
        const std::string& object_name);
    
    /*!
     * \brief Empty virtual destructor.
     */
    virtual ~SetDataStrategy();
    
    /*!
     * \name Methods to set patch interior data.
     */
    //\{
    
    /*!
     * \brief Indicates whether the concrete SetDataStrategy object is
     * time-dependent.
     */
    virtual bool isTimeDependent() const = 0;
    
    /*!
     * \brief Set data on the patch interiors on the specified levels
     * of the patch hierarchy using the virtual function
     * setDataOnPatch().
     * 
     * \see setDataOnPatch
     */
    void setDataOnPatchHierarchy(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const double data_time,
        const bool initial_time=false,
        const int coarsest_ln=-1,
        const int finest_ln=-1);
    
    /*!
     * \brief Set data on the patch interiors on the specified level
     * of the patch hierarchy using the virtual function
     * setDataOnPatch().
     *
     * \see setDataOnPatch
     */
    void setDataOnPatchLevel(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
        const double data_time,
        const bool initial_time=false);
    
    /*!
     * \brief Pure virtual function to set data on the patch interior.
     */
    virtual void setDataOnPatch(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false) = 0;
    
    //\}

    /*!
     * \brief Print the name of the strategy object.
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
    SetDataStrategy();
    
    /*!
     * \brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     * 
     * \param from The value to copy to this object.
     */
    SetDataStrategy(
        const SetDataStrategy& from);
    
    /*!
     * \brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     * 
     * \param that The value to assign to this object.
     * 
     * \return A reference to this object.
     */
    SetDataStrategy& operator=(
        const SetDataStrategy& that);

    /*
     * The object name is used for error/warning reporting.
     */
    std::string d_object_name;    
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "SetDataStrategy.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_SetDataStrategy
