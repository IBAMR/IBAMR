// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// INCLUDE GAURDS
#ifndef included_CartGridBodyForce
#define included_CartGridBodyForce

//////////////////////////////////INCLUDES///////////////////////////////////////////////

#include "ibtk/CartGridFunction.h"

#include "Patch.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "tbox/Pointer.h"

namespace IBTK
{
/*!
 * \brief CartGridBodyForce class.
 *
 * This class is a concrete implementation of IBTK::CartGridFunction class and is used to set the
 * body force on patches for the fluid solver.
 */

class CartGridBodyForce : public IBTK::CartGridFunction
{
private:
    /*!
     * Patch index for body force var.
     */
    int d_body_force_idx;

    /*!
     * \note The default ctor is not implemented and should not be used.
     */
    CartGridBodyForce();

    /*!
     * \note The default copy ctor is not implemented and should not be used.
     */
    CartGridBodyForce(const CartGridBodyForce& from);

    /*!
     * \note The default assignment operator is not implemented and should not be used.
     */
    CartGridBodyForce& operator=(const CartGridBodyForce& that);

public:
    /*!
     * \brief Constructor
     *
     * \param body_force_idx is the array index of the patch where an external body force for INS is calculated.
     */
    CartGridBodyForce(const int body_force_idx);

    /*!
     * Returns if the body force cartesian grid function is time dependent or not.
     */
    virtual bool isTimeDependent() const;

    /*!
     * This method copies the body force calculated a \textit{piori} on the  patch interior passed to this method.
     */
    virtual void setDataOnPatch(const int data_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                const double data_time,
                                const bool initial_time = false,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

}; // CartGridBodyForce

} // namespace IBTK
#endif // #ifndef included_CartGridBodyForce
