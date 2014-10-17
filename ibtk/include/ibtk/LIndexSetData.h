// Filename: LIndexSetData.h
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

#ifndef included_LIndexSetData
#define included_LIndexSetData

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "Box.h"
#include "IntVector.h"
#include "ibtk/LSetData.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LIndexSetData is a specialization of the templated class
 * LSetData that is intended to be used with Lagrangian data objects that
 * provide Lagrangian and PETSc indexing information.
 *
 * \see LSetData
 * \see SAMRAI::pdat::IndexData
 */
template <class T>
class LIndexSetData : public LSetData<T>
{
public:
    /*!
     * The constructor for an LIndexSetData object.  The box describes the
     * interior of the index space and the ghosts vector describes the ghost
     * nodes in each coordinate direction.
     */
    LIndexSetData(const SAMRAI::hier::Box<NDIM>& box, const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * The virtual destructor for an LIndexSetData object.
     */
    virtual ~LIndexSetData();

    /*!
     * \brief Update the cached indexing data.
     */
    void cacheLocalIndices(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           const SAMRAI::hier::IntVector<NDIM>& periodic_shift);

    /*!
     * \return A constant reference to the set of Lagrangian data indices that
     * lie in the patch (including the ghost cell region).
     */
    const std::vector<int>& getLagrangianIndices() const;

    /*!
     * \return A constant reference to the set of Lagrangian data indices that
     * lie in the patch interior.
     */
    const std::vector<int>& getInteriorLagrangianIndices() const;

    /*!
     * \return A constant reference to the set of Lagrangian data indices that
     * lie in the ghost cell region of the patch data object.
     */
    const std::vector<int>& getGhostLagrangianIndices() const;

    /*!
     * \return A constant reference to the set of global PETSc data indices that
     * lie in the patch (including the ghost cell region).
     */
    const std::vector<int>& getGlobalPETScIndices() const;

    /*!
     * \return A constant reference to the set of global PETSc data indices that
     * lie in the patch interior.
     */
    const std::vector<int>& getInteriorGlobalPETScIndices() const;

    /*!
     * \return A constant reference to the set of global PETSc data indices that
     * lie in the ghost cell region of the patch data object.
     */
    const std::vector<int>& getGhostGlobalPETScIndices() const;

    /*!
     * \return A constant reference to the set of local PETSc data indices that
     * lie in the patch (including the ghost cell region).
     */
    const std::vector<int>& getLocalPETScIndices() const;

    /*!
     * \return A constant reference to the set of local PETSc data indices that
     * lie in the patch interior.
     */
    const std::vector<int>& getInteriorLocalPETScIndices() const;

    /*!
     * \return A constant reference to the set of local PETSc data indices that
     * lie in the ghost cell region of the patch data object.
     */
    const std::vector<int>& getGhostLocalPETScIndices() const;

    /*!
     * \return A constant reference to the periodic shifts for the indices that
     * lie in the patch (including the ghost cell region).
     */
    const std::vector<double>& getPeriodicShifts() const;

    /*!
     * \return A constant reference to the periodic shifts for the indices that
     * lie in the patch interior.
     */
    const std::vector<double>& getInteriorPeriodicShifts() const;

    /*!
     * \return A constant reference to the periodic shifts for the indices that
     * lie in the ghost cell region of the patch data object.
     */
    const std::vector<double>& getGhostPeriodicShifts() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LIndexSetData();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LIndexSetData(const LIndexSetData<T>& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LIndexSetData& operator=(const LIndexSetData<T>& that);

    std::vector<int> d_lag_indices, d_interior_lag_indices, d_ghost_lag_indices;
    std::vector<int> d_global_petsc_indices, d_interior_global_petsc_indices, d_ghost_global_petsc_indices;
    std::vector<int> d_local_petsc_indices, d_interior_local_petsc_indices, d_ghost_local_petsc_indices;
    std::vector<double> d_periodic_shifts, d_interior_periodic_shifts, d_ghost_periodic_shifts;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LIndexSetData-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LIndexSetData
