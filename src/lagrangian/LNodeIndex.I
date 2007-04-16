// Filename: LNodeIndex.I
// Created on 28 Feb 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <16.Apr.2007 02:03:27 boyce@trasnaform2.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/StashableManager.h>

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
LNodeIndex::LNodeIndex(
    const int lagrangian_nidx,
    const int local_petsc_nidx,
    double* const X_ptr,
    const std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data)
    : d_lagrangian_nidx(lagrangian_nidx),
      d_local_petsc_nidx(local_petsc_nidx),
      d_X_ptr(X_ptr),
      d_stash_data(stash_data)
{
    // intentionally blank
    return;
}// LNodeIndex

inline
LNodeIndex::LNodeIndex(
    const LNodeIndex& from)
    : d_lagrangian_nidx(from.d_lagrangian_nidx),
      d_local_petsc_nidx(from.d_local_petsc_nidx),
      d_X_ptr(from.d_X_ptr),
      d_stash_data(from.d_stash_data)
{
    // intentionally blank
    return;
}// LNodeIndex

inline
LNodeIndex::~LNodeIndex()
{
    // intentionally blank
    return;
}// ~LNodeIndex

inline LNodeIndex&
LNodeIndex::operator=(
    const LNodeIndex& that)
{
    if (this != &that)
    {
        assignThatToThis(that);
    }
    return *this;
}// operator=

inline int
LNodeIndex::getLagrangianIndex() const
{
    return d_lagrangian_nidx;
}// getLagrangianIndex

inline void
LNodeIndex::setLagrangianIndex(
    const int lagrangian_nidx)
{
    d_lagrangian_nidx = lagrangian_nidx;
    return;
}// setLagrangianIndex

inline int
LNodeIndex::getLocalPETScIndex() const
{
    return d_local_petsc_nidx;
}// getLocalPETScIndex

inline void
LNodeIndex::setLocalPETScIndex(
    const int local_petsc_nidx)
{
    d_local_petsc_nidx = local_petsc_nidx;
    return;
}// setLocalPETScIndex

inline double*
LNodeIndex::getNodeLocation() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_X_ptr != NULL);
#endif
    return d_X_ptr;
}// getNodeLocation

inline void
LNodeIndex::setNodeLocation(
    double* const X_ptr)
{
    d_X_ptr = X_ptr;
    return;
}// setNodeLocation

inline const std::vector<SAMRAI::tbox::Pointer<Stashable> >&
LNodeIndex::getStashData() const
{
    return d_stash_data;
}// getStashData

inline vector<SAMRAI::tbox::Pointer<Stashable> >&
LNodeIndex::getStashData()
{
    return d_stash_data;
}// getStashData

inline void
LNodeIndex::copySourceItem(
    const SAMRAI::hier::Index<NDIM>& src_index,
    const SAMRAI::hier::IntVector<NDIM>& src_offset,
    const LNodeIndex& src_item)
{
    assignThatToThis(src_item);
    return;
}// copySourceItem

inline size_t
LNodeIndex::getDataStreamSize() const
{
    return 2*SAMRAI::tbox::AbstractStream::sizeofInt()
        + StashableManager::getManager()->getDataStreamSize(d_stash_data);
}// getDataStreamSize

inline void
LNodeIndex::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_lagrangian_nidx,1);
    stream.pack(&d_local_petsc_nidx, 1);

    StashableManager::getManager()->packStream(stream, d_stash_data);
    return;
}// packStream

inline void
LNodeIndex::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    stream.unpack(&d_lagrangian_nidx,1);
    stream.unpack(&d_local_petsc_nidx, 1);
    d_X_ptr = NULL;

    d_stash_data.clear();
    StashableManager::getManager()->unpackStream(stream, offset, d_stash_data);
    return;
}// unpackStream

/////////////////////////////// PRIVATE //////////////////////////////////////

inline void
LNodeIndex::assignThatToThis(
    const LNodeIndex& that)
{
    d_lagrangian_nidx = that.d_lagrangian_nidx;
    d_local_petsc_nidx = that.d_local_petsc_nidx;
    d_X_ptr = that.d_X_ptr;
    d_stash_data = that.d_stash_data;
    return;
}// assignThatToThis

//////////////////////////////////////////////////////////////////////////////

inline bool
operator<(
    const LNodeIndex& lhs,
    const LNodeIndex& rhs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(lhs.d_X_ptr != NULL);
    assert(rhs.d_X_ptr != NULL);
#if ((NDIM>3) || (NDIM<1))
    TBOX_ERROR("operator<(const LNodeIndex&,const LNodeIndex&): not implemented for NDIM=="
               << NDIM << endl);
#endif
#endif

    return (lhs.d_X_ptr[0]<rhs.d_X_ptr[0]) ||
#if (NDIM>1)
        ((lhs.d_X_ptr[0]==rhs.d_X_ptr[0])&&
         (lhs.d_X_ptr[1]< rhs.d_X_ptr[1])) ||
#if (NDIM>2)
        ((lhs.d_X_ptr[0]==rhs.d_X_ptr[0])&&
         (lhs.d_X_ptr[1]==rhs.d_X_ptr[1])&&
         (lhs.d_X_ptr[2]< rhs.d_X_ptr[2])) ||
#endif
#endif
        ((lhs.d_X_ptr[0]==rhs.d_X_ptr[0])&&
#if (NDIM>1)
         (lhs.d_X_ptr[1]==rhs.d_X_ptr[1])&&
#if (NDIM>2)
         (lhs.d_X_ptr[2]==rhs.d_X_ptr[2])&&
#endif
#endif
         (lhs.d_lagrangian_nidx<rhs.d_lagrangian_nidx));
}// operator<

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
