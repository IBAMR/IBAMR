// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2022 by the IBAMR developers
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

#ifndef included_IBTK_LNodeIndex
#define included_IBTK_LNodeIndex

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/ibtk_utilities.h"

#include "IntVector.h"
#include "tbox/DescribedClass.h"
#include "tbox/MathUtilities.h"
#include "tbox/Utilities.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <functional>
#include <ostream>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Index;
} // namespace hier
namespace tbox
{
class AbstractStream;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LNodeIndex provides Lagrangian and <A
 * HREF="http://www.mcs.anl.gov/petsc">PETSc</A> indexing information for a
 * single node of a Lagrangian mesh.
 */
class LNodeIndex : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    LNodeIndex(int lagrangian_nidx = -1,
               int global_petsc_nidx = -1,
               int local_petsc_nidx = -1,
               const SAMRAI::hier::IntVector<NDIM>& initial_periodic_offset = SAMRAI::hier::IntVector<NDIM>(0),
               const SAMRAI::hier::IntVector<NDIM>& current_periodic_offset = SAMRAI::hier::IntVector<NDIM>(0),
               const Vector& initial_periodic_displacement = Vector::Zero(),
               const Vector& current_periodic_displacement = Vector::Zero());

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    LNodeIndex(const LNodeIndex& from);

    /*!
     * \brief Constructor that unpacks data from an input stream.
     */
    LNodeIndex(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~LNodeIndex();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeIndex& operator=(const LNodeIndex& that);

    /*!
     * \return The Lagrangian index referenced by this LNodeIndex object.
     */
    int getLagrangianIndex() const;

    /*!
     * \brief Reset the Lagrangian index referenced by this LNodeIndex object.
     */
    void setLagrangianIndex(int lagrangian_nidx);

    /*!
     * \return The global PETSc index referenced by this LNodeIndex object.
     */
    int getGlobalPETScIndex() const;

    /*!
     * \brief Reset the global PETSc index referenced by this LNodeIndex object.
     */
    void setGlobalPETScIndex(int global_petsc_nidx);

    /*!
     * \return The local PETSc index referenced by this LNodeIndex object.
     */
    int getLocalPETScIndex() const;

    /*!
     * \brief Reset the local PETSc index referenced by this LNodeIndex object.
     */
    void setLocalPETScIndex(int local_petsc_nidx);

    /*!
     * \brief Indicate that the LNodeIndex object has been shifted across a
     * periodic boundary.
     */
    virtual void registerPeriodicShift(const SAMRAI::hier::IntVector<NDIM>& offset, const Vector& displacement);

    /*!
     * \brief Get the initial (t = 0) periodic offset.
     */
    virtual const SAMRAI::hier::IntVector<NDIM>& getInitialPeriodicOffset() const;

    /*!
     * \brief Get the periodic offset.
     */
    virtual const SAMRAI::hier::IntVector<NDIM>& getPeriodicOffset() const;

    /*!
     * \brief Get the initial (t = 0) periodic displacement.
     */
    virtual const Vector& getInitialPeriodicDisplacement() const;

    /*!
     * \brief Get the periodic displacement.
     */
    virtual const Vector& getPeriodicDisplacement() const;

    /*!
     * \brief Copy data from the source.
     *
     * \note The cell index of the destination object is src_index + src_offset.
     */
    virtual void copySourceItem(const SAMRAI::hier::Index<NDIM>& src_index,
                                const SAMRAI::hier::IntVector<NDIM>& src_offset,
                                const LNodeIndex& src_item);

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    virtual size_t getDataStreamSize() const;

    /*!
     * \brief Pack data into the output stream.
     */
    virtual void packStream(SAMRAI::tbox::AbstractStream& stream);

    /*!
     * \brief Unpack data from the input stream.
     */
    virtual void unpackStream(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset);

private:
    /*!
     * Assign that to this.
     */
    void assignThatToThis(const LNodeIndex& that);

    int d_lagrangian_nidx;   // the fixed global Lagrangian index
    int d_global_petsc_nidx; // the global PETSc index
    int d_local_petsc_nidx;  // the local PETSc index

    // the periodic offset and displacement
    SAMRAI::hier::IntVector<NDIM> d_offset_0, d_offset;
    Vector d_displacement_0, d_displacement;
};

/*!
 * \brief Comparison functor to order on the physical location of the Lagrangian
 * node.
 */
class LNodeIndexPosnComp
{
public:
    LNodeIndexPosnComp(const boost::multi_array_ref<double, 2>& X_ghosted_local_form_array)
        : d_X_ghosted_local_form_array(&X_ghosted_local_form_array)
    {
        // intentionally blank
        return;
    }

    ~LNodeIndexPosnComp()
    {
        // intentionally blank
        return;
    }

    inline bool operator()(const LNodeIndex& lhs, const LNodeIndex& rhs) const
    {
#if !defined(NDEBUG)
#if ((NDIM > 3) || (NDIM < 1))
        TBOX_ERROR("operator<(const LNodeIndex&,const LNodeIndex&): not implemented for NDIM==" << NDIM << endl);
#endif
#endif
        const double* const X_lhs = &(*d_X_ghosted_local_form_array)[lhs.getLocalPETScIndex()][0];
        const double* const X_rhs = &(*d_X_ghosted_local_form_array)[rhs.getLocalPETScIndex()][0];
        return ((X_lhs[0] < X_rhs[0])) ||
#if (NDIM > 1)
               ((X_lhs[0] == X_rhs[0]) && (X_lhs[1] < X_rhs[1])) ||
#if (NDIM > 2)
               ((X_lhs[0] == X_rhs[0]) && (X_lhs[1] == X_rhs[1]) && (X_lhs[2] < X_rhs[2])) ||
#endif
#endif
               ((X_lhs[0] == X_rhs[0]) &&
#if (NDIM > 1)
                (X_lhs[1] == X_rhs[1]) &&
#if (NDIM > 2)
                (X_lhs[2] == X_rhs[2]) &&
#endif
#endif
                (lhs.getLagrangianIndex() < rhs.getLagrangianIndex()));
    } // operator()

    inline bool operator()(const LNodeIndex* lhs, const LNodeIndex* rhs) const
    {
        return (*this)(*lhs, *rhs);
    } // operator()

private:
    const boost::multi_array_ref<double, 2>* const d_X_ghosted_local_form_array;
};

/*!
 * \brief Comparison functor to order on the Lagrangian index of the Lagrangian
 * node.
 */
struct LNodeIndexLagrangianIndexComp
{
    inline bool operator()(const LNodeIndex& lhs, const LNodeIndex& rhs) const
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(lhs.getLagrangianIndex() >= 0);
        TBOX_ASSERT(rhs.getLagrangianIndex() >= 0);
#endif
        return lhs.getLagrangianIndex() < rhs.getLagrangianIndex();
    } // operator()

    inline bool operator()(const LNodeIndex* lhs, const LNodeIndex* rhs) const
    {
        return (*this)(*lhs, *rhs);
    } // operator()
};

/*!
 * \brief Comparison functor to order on the global PETSc index of the
 * Lagrangian node.
 */
struct LNodeIndexGlobalPETScIndexComp
{
    inline bool operator()(const LNodeIndex& lhs, const LNodeIndex& rhs) const
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(lhs.getGlobalPETScIndex() >= 0);
        TBOX_ASSERT(rhs.getGlobalPETScIndex() >= 0);
#endif
        return lhs.getGlobalPETScIndex() < rhs.getGlobalPETScIndex();
    } // operator()

    inline bool operator()(const LNodeIndex* lhs, const LNodeIndex* rhs) const
    {
        return (*this)(*lhs, *rhs);
    } // operator()
};

/*!
 * \brief Comparison functor to order on the local PETSc index of the
 * Lagrangian node.
 */
struct LNodeIndexLocalPETScIndexComp
{
    inline bool operator()(const LNodeIndex& lhs, const LNodeIndex& rhs) const
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(lhs.getLocalPETScIndex() >= 0);
        TBOX_ASSERT(rhs.getLocalPETScIndex() >= 0);
#endif
        return lhs.getLocalPETScIndex() < rhs.getLocalPETScIndex();
    } // operator()

    inline bool operator()(const LNodeIndex* lhs, const LNodeIndex* rhs) const
    {
        return (*this)(*lhs, *rhs);
    } // operator()
};

/*!
 * \brief Comparison functor to check for equality between LNodeIndex objects
 * based on their positions.
 */
class LNodeIndexPosnEqual
{
public:
    LNodeIndexPosnEqual(const boost::multi_array_ref<double, 2>& X_ghosted_local_form_array)
        : d_X_ghosted_local_form_array(&X_ghosted_local_form_array)
    {
        // intentionally blank
        return;
    }

    ~LNodeIndexPosnEqual()
    {
        // intentionally blank
        return;
    }

    inline bool operator()(const LNodeIndex& lhs, const LNodeIndex& rhs) const
    {
        const double* const X_lhs = &(*d_X_ghosted_local_form_array)[lhs.getLocalPETScIndex()][0];
        const double* const X_rhs = &(*d_X_ghosted_local_form_array)[rhs.getLocalPETScIndex()][0];
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (!IBTK::rel_equal_eps(X_lhs[d], X_rhs[d])) return false;
        }
        return true;
    } // operator()

    inline bool operator()(const LNodeIndex* lhs, const LNodeIndex* rhs) const
    {
        return (*this)(*lhs, *rhs);
    } // operator()

private:
    const boost::multi_array_ref<double, 2>* const d_X_ghosted_local_form_array;
};

/*!
 * \brief Comparison functor to check for equality between LNodeIndex objects
 * based on their Lagrangian indices.
 */
struct LNodeIndexLagrangianIndexEqual
{
    inline bool operator()(const LNodeIndex& lhs, const LNodeIndex& rhs) const
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(lhs.getLagrangianIndex() >= 0);
        TBOX_ASSERT(rhs.getLagrangianIndex() >= 0);
#endif
        return lhs.getLagrangianIndex() == rhs.getLagrangianIndex();
    } // operator()

    inline bool operator()(const LNodeIndex* lhs, const LNodeIndex* rhs) const
    {
        return (*this)(*lhs, *rhs);
    } // operator()
};

/*!
 * \brief Comparison functor to check for equality between between LNodeIndex
 * objects based on their global PETSc indices.
 */
struct LNodeIndexGlobalPETScIndexEqual
{
    inline bool operator()(const LNodeIndex& lhs, const LNodeIndex& rhs) const
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(lhs.getGlobalPETScIndex() >= 0);
        TBOX_ASSERT(rhs.getGlobalPETScIndex() >= 0);
#endif
        return lhs.getGlobalPETScIndex() == rhs.getGlobalPETScIndex();
    } // operator()

    inline bool operator()(const LNodeIndex* lhs, const LNodeIndex* rhs) const
    {
        return (*this)(*lhs, *rhs);
    } // operator()
};

/*!
 * \brief Comparison functor to check for equality between LNodeIndex objects
 * based on their local PETSc indices.
 */
struct LNodeIndexLocalPETScIndexEqual
{
    inline bool operator()(const LNodeIndex& lhs, const LNodeIndex& rhs) const
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(lhs.getLocalPETScIndex() >= 0);
        TBOX_ASSERT(rhs.getLocalPETScIndex() >= 0);
#endif
        return lhs.getLocalPETScIndex() == rhs.getLocalPETScIndex();
    } // operator()

    inline bool operator()(const LNodeIndex* lhs, const LNodeIndex* rhs) const
    {
        return (*this)(*lhs, *rhs);
    } // operator()
};

} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LNodeIndex-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LNodeIndex
