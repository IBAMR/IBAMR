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

#ifndef included_IBTK_LSetVariable
#define included_IBTK_LSetVariable

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "Variable.h"

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LSetVariable provides a SAMRAI::hier::Variable class
 * corresponding to patch data of type LSetData.
 */
template <class T>
class LSetVariable : public SAMRAI::hier::Variable<NDIM>
{
public:
    /*!
     * Create an LSetVariable object with the specified name.
     */
    LSetVariable(std::string name);

    /*!
     * Virtual destructor for LSetVariable objects.
     */
    virtual ~LSetVariable() = default;

    /*!
     * Return false since the LSet data index space matches the
     * cell-centered index space for AMR patches.  Thus, LSet data does
     * not live on patch borders.
     */
    bool dataLivesOnPatchBorder() const override;

    /*!
     * Return true so that the LSet data quantities will always be treated
     * as though fine values represent them on coarse-fine interfaces.  Note
     * that this is really artificial since the LSet data index space
     * matches the cell-centered index space for AMR patches.  Thus, LSet
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
    LSetVariable() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LSetVariable(const LSetVariable<T>& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LSetVariable<T>& operator=(const LSetVariable<T>& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LSetVariable
