// Filename: LIndexSetVariable.h
// Created on 13 May 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_LIndexSetVariable
#define included_LIndexSetVariable

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "Variable.h"

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
    LIndexSetVariable(const std::string& name);

    /*!
     * Virtual destructor for LIndexSetVariable objects.
     */
    virtual ~LIndexSetVariable();

    /*!
     * Return false since the LIndexSet data index space matches the
     * cell-centered index space for AMR patches.  Thus, LIndexSet data does not
     * live on patch borders.
     */
    bool dataLivesOnPatchBorder() const;

    /*!
     * Return true so that the LIndexSet data quantities will always be treated
     * as though fine values represent them on coarse-fine interfaces.  Note
     * that this is really artificial since the LIndexSet data index space
     * matches the cell-centered index space for AMR patches.  Thus, LIndexSet
     * data does not live on patch borders and so there is no ambiguity
     * regarding coarse-fine interface values.
     */
    bool fineBoundaryRepresentsVariable() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LIndexSetVariable();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LIndexSetVariable(const LIndexSetVariable<T>& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LIndexSetVariable<T>& operator=(const LIndexSetVariable<T>& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LIndexSetVariable
