// Filename: IBAnchorPointSpec.C
// Last modified: <15.Dec.2009 19:13:06 griffith@boyce-griffiths-mac-pro.local>
// Created on 18 Aug 2008 by Boyce Griffith (boyce@dm-linux.maths.gla.ac.uk)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBAnchorPointSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline void*
IBAnchorPointSpec::operator new(
    std::size_t size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(size == sizeof(IBAnchorPointSpec));
#endif
    return s_arena.alloc(size);
}// operator new

inline void
IBAnchorPointSpec::operator delete(
    void* ptr,
    std::size_t size)
{
    if (ptr == NULL) return;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(size == sizeof(IBAnchorPointSpec));
#endif
    s_arena.free(ptr);
    return;
}// operator delete

inline
IBAnchorPointSpec::IBAnchorPointSpec(
    const int node_idx)
    : d_node_idx(node_idx)
{
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBAnchorPointSpec::IBAnchorPointSpec():\n"
                   << "  must call IBAnchorPointSpec::registerWithStashableManager() before\n"
                   << "  creating any IBAnchorPointSpec objects.\n");
    }
    return;
}// IBAnchorPointSpec

inline
IBAnchorPointSpec::~IBAnchorPointSpec()
{
    // intentionally blank
    return;
}// ~IBAnchorPointSpec

inline const int&
IBAnchorPointSpec::getNodeIndex() const
{
    return d_node_idx;
}// getNodeIndex

inline int&
IBAnchorPointSpec::getNodeIndex()
{
    return d_node_idx;
}// getNodeIndex

inline int
IBAnchorPointSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBAnchorPointSpec::getDataStreamSize() const
{
    return SAMRAI::tbox::AbstractStream::sizeofInt();
}// getDataStreamSize

inline void
IBAnchorPointSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_node_idx,1);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
