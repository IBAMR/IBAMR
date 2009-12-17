// Filename: IBInstrumentationSpec.I
// Last modified: <15.Dec.2009 19:17:45 griffith@boyce-griffiths-mac-pro.local>
// Created on 11 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBInstrumentationSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

inline const std::vector<std::string>&
IBInstrumentationSpec::getInstrumentNames()
{
    return s_instrument_names;
}// getInstrumentNames

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline void*
IBInstrumentationSpec::operator new(
    std::size_t size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(size == sizeof(IBInstrumentationSpec));
#endif
    return s_arena.alloc(size);
}// operator new

inline void
IBInstrumentationSpec::operator delete(
    void* ptr,
    std::size_t size)
{
    if (ptr == NULL) return;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(size == sizeof(IBInstrumentationSpec));
#endif
    s_arena.free(ptr);
    return;
}// operator delete

inline
IBInstrumentationSpec::IBInstrumentationSpec(
    const int master_idx,
    const int meter_idx,
    const int node_idx)
    : d_master_idx(master_idx),
      d_meter_idx(meter_idx),
      d_node_idx(node_idx)
{
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBInstrumentationSpec::IBInstrumentationSpec():\n"
                   << "  must call IBInstrumentationSpec::registerWithStashableManager() before\n"
                   << "  creating any IBInstrumentationSpec objects.\n");
    }
    return;
}// IBInstrumentationSpec

inline
IBInstrumentationSpec::~IBInstrumentationSpec()
{
    // intentionally blank
    return;
}// ~IBInstrumentationSpec

inline const int&
IBInstrumentationSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

inline int&
IBInstrumentationSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

inline const int&
IBInstrumentationSpec::getMeterIndex() const
{
    return d_meter_idx;
}// getMeterIndex

inline int&
IBInstrumentationSpec::getMeterIndex()
{
    return d_meter_idx;
}// getMeterIndex

inline const int&
IBInstrumentationSpec::getNodeIndex() const
{
    return d_node_idx;
}// getNodeIndex

inline int&
IBInstrumentationSpec::getNodeIndex()
{
    return d_node_idx;
}// getNodeIndex

inline int
IBInstrumentationSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBInstrumentationSpec::getDataStreamSize() const
{
    return 5*SAMRAI::tbox::AbstractStream::sizeofInt();
}// getDataStreamSize

inline void
IBInstrumentationSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_master_idx,1);
    stream.pack(&d_meter_idx,1);
    stream.pack(&d_node_idx,1);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// NAMESPACE ////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
