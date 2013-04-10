// Filename: VecCellVariable.h
// Created on 09 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_VecCellVariable
#define included_VecCellVariable

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <Variable.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * Class VecCellVariable<NDIM> is a templated variable class used to define
 * vector-valued cell-centered quantities on an AMR mesh.  It is a subclass of
 * hier::Variable and is templated on the type of the underlying data (e.g.,
 * double, int, bool, etc.).
 *
 * See header file for VecCellData<NDIM> class for a more detailed description of
 * the data layout.
 *
 * \see VecCellData
 * \see VecCellDataFactory
 * \see SAMRAI::hier::Variable
 */
template<class TYPE>
class VecCellVariable
    : public SAMRAI::hier::Variable<NDIM>
{
public:
    /*!
     * \brief Create a cell-centered variable object with the given name and
     * depth (i.e., number of data values at each cell index location).
     */
    VecCellVariable(
        const std::string& name,
        int depth);

    /*!
     * \brief Virtual destructor for cell variable objects.
     */
    virtual
    ~VecCellVariable<TYPE>();

    /*!
     * \brief Return true indicating that cell data quantities will always be
     * treated as though fine values take precedence on coarse-fine interfaces.
     * Note that this is really artificial since the cell data index space
     * matches the cell-centered index space for AMR patches.  However, some
     * value must be supplied for communication operations.
     */
    bool
    fineBoundaryRepresentsVariable() const;

    /*!
     * \brief Return false indicating that cell data on a patch interior does not
     * exist on the patch boundary.
     */
    bool
    dataLivesOnPatchBorder() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VecCellVariable();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VecCellVariable(
        const VecCellVariable& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VecCellVariable&
    operator=(
        const VecCellVariable& that);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/VecCellVariable.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VecCellVariable
