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

#ifndef included_IBTK_LIndexSetVariable
#define included_IBTK_LIndexSetVariable

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "Variable.h"

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LIndexSetVariable provides a SAMRAI::hier::Variable class
 * corresponding to patch data of type LIndexSetData.
 */
template <class T>
class LIndexSetVariable : public SAMRAI::hier::Variable<NDIM>
{
public:
    /*!
     * Create an LIndexSetVariable object with the specified name.
     */
    LIndexSetVariable(std::string name);

    /*!
     * Virtual destructor for LIndexSetVariable objects.
     */
    virtual ~LIndexSetVariable() = default;

    /*!
     * Return false since the LIndexSet data index space matches the
     * cell-centered index space for AMR patches.  Thus, LIndexSet data does not
     * live on patch borders.
     */
    bool dataLivesOnPatchBorder() const override;

    /*!
     * Return true so that the LIndexSet data quantities will always be treated
     * as though fine values represent them on coarse-fine interfaces.  Note
     * that this is really artificial since the LIndexSet data index space
     * matches the cell-centered index space for AMR patches.  Thus, LIndexSet
     * data does not live on patch borders and so there is no ambiguity
     * regarding coarse-fine interface values.
     */
    bool fineBoundaryRepresentsVariable() const override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LIndexSetVariable() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LIndexSetVariable(const LIndexSetVariable<T>& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LIndexSetVariable<T>& operator=(const LIndexSetVariable<T>& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LIndexSetVariable
