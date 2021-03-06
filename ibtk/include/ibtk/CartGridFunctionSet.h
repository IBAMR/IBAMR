// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_CartGridFunctionSet
#define included_IBTK_CartGridFunctionSet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/CartGridFunction.h"

#include "IntVector.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

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
 * \brief Class CartGridFunctionSet is a concrete CartGridFunction that is used
 * to allow multiple CartGridFunction objects to act as a single function.
 */
class CartGridFunctionSet : public CartGridFunction
{
public:
    /*!
     * \brief The default constructor sets the name of the strategy object and
     * sets the collection of functions to be empty.
     */
    CartGridFunctionSet(std::string object_name = "");

    /*!
     * \brief Empty virtual destructor.
     */
    virtual ~CartGridFunctionSet() = default;

    /*!
     * \brief Add a CartGridFunction to the set of functions grouped together by
     * this object.
     */
    void addFunction(SAMRAI::tbox::Pointer<CartGridFunction> fcn);

    /*!
     * \name Methods to set patch interior data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CartGridFunctionSet object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy using the implementations of
     * setDataOnPatchHierarchy() provided by the component function objects.
     */
    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1) override;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * level of the patch hierarchy using the implementations of
     * setDataOnPatchLevel() provided by the component function objects.
     */
    void setDataOnPatchLevel(int data_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                             double data_time,
                             bool initial_time = false) override;

    /*!
     * \brief Evaluate the function on the patch interior using the
     * implementations of setDataOnPatch() provided by the component function
     * objects.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL)) override;

    //\}

protected:
    /*
     * The collection of function objects.
     */
    std::vector<SAMRAI::tbox::Pointer<CartGridFunction> > d_fcns;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CartGridFunctionSet(const CartGridFunctionSet& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartGridFunctionSet& operator=(const CartGridFunctionSet& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CartGridFunctionSet
