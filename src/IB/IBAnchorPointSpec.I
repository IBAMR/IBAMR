// Filename: IBAnchorPointSpec.C
// Last modified: <18.Aug.2008 13:35:17 boyce@dm-linux.maths.gla.ac.uk>
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
