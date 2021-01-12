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

#ifndef included_IBTK_LIndexSetData
#define included_IBTK_LIndexSetData

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LSetData.h"

#include "Box.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

#include <vector>

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
    LIndexSetData(SAMRAI::hier::Box<NDIM> box, SAMRAI::hier::IntVector<NDIM> ghosts);

    /*!
     * The virtual destructor for an LIndexSetData object.
     */
    virtual ~LIndexSetData() = default;

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
    LIndexSetData() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LIndexSetData(const LIndexSetData<T>& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LIndexSetData& operator=(const LIndexSetData<T>& that) = delete;

    std::vector<int> d_lag_indices, d_interior_lag_indices, d_ghost_lag_indices;
    std::vector<int> d_global_petsc_indices, d_interior_global_petsc_indices, d_ghost_global_petsc_indices;
    std::vector<int> d_local_petsc_indices, d_interior_local_petsc_indices, d_ghost_local_petsc_indices;
    std::vector<double> d_periodic_shifts, d_interior_periodic_shifts, d_ghost_periodic_shifts;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LIndexSetData-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LIndexSetData
