#ifndef included_SetDataStrategy
#define included_SetDataStrategy

// Filename: SetDataStrategy.h
// Last modified: <17.Aug.2006 15:55:24 boyce@bigboy.nyconnect.com>
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
#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * @brief Interface to allow the specification of data for a single
 * quantity at a specified time.
 *
 * The interface specified by this class is very similar to that of
 * the PhysicalBCDataStrategy class.  The different interfaces are
 * provided for the sake of code clarity---i.e. in order to easily
 * separate classes which set internal patch data from those which set
 * boundary data.  Additionally, note that the interfaces required of
 * concrete implementations of both strategy class are different.
 *
 * @see PhysicalBCDataStrategy
 */
namespace IBAMR
{
class SetDataStrategy
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * @brief The constructor sets the name of the strategy object.
     */
    SetDataStrategy(
        const std::string& object_name);
    
    /*!
     * @brief Empty virtual destructor.
     */
    virtual ~SetDataStrategy();
    
    //@{ @name Methods to set the data.

    /*!
     * Indicates whether the concrete SetDataStrategy object is time
     * dependent.
     */
    virtual bool isTimeDependent() const = 0;
    
    /*!
     * Set data on the patch interiors on the specified levels of the
     * patch hierarchy using the virtual function setDataOnPatch().
     * 
     * @see setDataOnPatch
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
     * Set data on the patch interiors on the specified level of the
     * patch hierarchy using the virtual function setDataOnPatch().
     *
     * @see setDataOnPatch
     */
    void setDataOnPatchLevel(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
        const double data_time,
        const bool initial_time=false);
    
    /*!
     * Set the data on the patch interior.
     */
    virtual void setDataOnPatch(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false) = 0;
    
    //@}    

    /*!
     * Prints the name of the strategy object.
     */
    virtual void printClassData(
        ostream& os) const;

private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    SetDataStrategy();
    
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     * 
     * @param from The value to copy to this object.
     */
    SetDataStrategy(
        const SetDataStrategy& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     * 
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
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

//#ifndef DEBUG_NO_INLINE
//#include "SetDataStrategy.I"
//#endif

#endif //#ifndef included_SetDataStrategy
