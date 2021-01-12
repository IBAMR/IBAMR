// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_CartGridFunction
#define included_IBTK_CartGridFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "PatchLevel.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <string>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CartGridFunction provides an abstract interface for objects for
 * evaluating functions to set values in SAMRAI::hier::PatchData objects.
 */
class CartGridFunction : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief The default constructor sets the name of the strategy object.
     */
    CartGridFunction(std::string object_name = "");

    /*!
     * \brief Empty virtual destructor.
     */
    virtual ~CartGridFunction() = default;

    /*!
     * \name Methods to set patch interior data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CartGridFunction object is
     * time-dependent.
     */
    virtual bool isTimeDependent() const = 0;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy using the virtual function
     * setDataOnPatch().
     *
     * \see setDataOnPatch
     */
    virtual void setDataOnPatchHierarchy(int data_idx,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                         double data_time,
                                         bool initial_time = false,
                                         int coarsest_ln = -1,
                                         int finest_ln = -1);

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * level of the patch hierarchy using the virtual function setDataOnPatch().
     *
     * \see setDataOnPatch
     */
    virtual void setDataOnPatchLevel(int data_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                                     double data_time,
                                     bool initial_time = false);

    /*!
     * \brief Pure virtual function to evaluate the function on the patch
     * interior.
     */
    virtual void setDataOnPatch(int data_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                double data_time,
                                bool initial_time = false,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL)) = 0;

    //\}

protected:
    /*
     * The object name is used for error/warning reporting.
     */
    std::string d_object_name;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CartGridFunction(const CartGridFunction& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartGridFunction& operator=(const CartGridFunction& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CartGridFunction
