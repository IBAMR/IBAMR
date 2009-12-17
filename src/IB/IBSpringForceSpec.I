// Filename: IBSpringForceSpec.I
// Last modified: <15.Dec.2009 19:21:46 griffith@boyce-griffiths-mac-pro.local>
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
IBSpringForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline void*
IBSpringForceSpec::operator new(
    std::size_t size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(size == sizeof(IBSpringForceSpec));
#endif
    return s_arena.alloc(size);
}// operator new

inline void
IBSpringForceSpec::operator delete(
    void* ptr,
    std::size_t size)
{
    if (ptr == NULL) return;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(size == sizeof(IBSpringForceSpec));
#endif
    s_arena.free(ptr);
    return;
}// operator delete

inline
IBSpringForceSpec::IBSpringForceSpec(
    const int master_idx,
    const std::vector<int>& slave_idxs,
    const std::vector<int>& force_fcn_idxs,
    const std::vector<double>& stiffnesses,
    const std::vector<double>& rest_lengths)
    : d_master_idx(master_idx),
      d_slave_idxs(slave_idxs),
      d_force_fcn_idxs(force_fcn_idxs),
      d_stiffnesses(stiffnesses),
      d_rest_lengths(rest_lengths)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    const int num_springs = d_slave_idxs.size();
    TBOX_ASSERT(num_springs == int(d_force_fcn_idxs.size()));
    TBOX_ASSERT(num_springs == int(d_stiffnesses.size()));
    TBOX_ASSERT(num_springs == int(d_rest_lengths.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBSpringForceSpec::IBSpringForceSpec():\n"
                   << "  must call IBSpringForceSpec::registerWithStashableManager() before\n"
                   << "  creating any IBSpringForceSpec objects.\n");
    }
    return;
}// IBSpringForceSpec

inline
IBSpringForceSpec::~IBSpringForceSpec()
{
    // intentionally blank
    return;
}// ~IBSpringForceSpec

inline unsigned
IBSpringForceSpec::getNumberOfSprings() const
{
    const int num_springs = d_slave_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_springs == int(d_force_fcn_idxs.size()));
    TBOX_ASSERT(num_springs == int(d_stiffnesses.size()));
    TBOX_ASSERT(num_springs == int(d_rest_lengths.size()));
#endif
    return num_springs;
}// getNumberOfSprings

inline const int&
IBSpringForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

inline int&
IBSpringForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

inline const std::vector<int>&
IBSpringForceSpec::getSlaveNodeIndices() const
{
    return d_slave_idxs;
}// getSlaveNodeIndices

inline std::vector<int>&
IBSpringForceSpec::getSlaveNodeIndices()
{
    return d_slave_idxs;
}// getSlaveNodeIndices

inline const std::vector<int>&
IBSpringForceSpec::getForceFunctionIndices() const
{
    return d_force_fcn_idxs;
}// getForceFunctionIndices

inline std::vector<int>&
IBSpringForceSpec::getForceFunctionIndices()
{
    return d_force_fcn_idxs;
}// getForceFunctionIndices

inline const std::vector<double>&
IBSpringForceSpec::getStiffnesses() const
{
    return d_stiffnesses;
}// getStiffnesses

inline std::vector<double>&
IBSpringForceSpec::getStiffnesses()
{
    return d_stiffnesses;
}// getStiffnesses

inline const std::vector<double>&
IBSpringForceSpec::getRestingLengths() const
{
    return d_rest_lengths;
}// getRestingLengths

inline std::vector<double>&
IBSpringForceSpec::getRestingLengths()
{
    return d_rest_lengths;
}// getRestingLengths

inline int
IBSpringForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBSpringForceSpec::getDataStreamSize() const
{
    const int num_springs = d_slave_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_springs == int(d_force_fcn_idxs.size()));
    TBOX_ASSERT(num_springs == int(d_stiffnesses.size()));
    TBOX_ASSERT(num_springs == int(d_rest_lengths.size()));
#endif
    return ((2+2*num_springs)*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (  2*num_springs)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

inline void
IBSpringForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    const int num_springs = d_slave_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_springs == int(d_force_fcn_idxs.size()));
    TBOX_ASSERT(num_springs == int(d_stiffnesses.size()));
    TBOX_ASSERT(num_springs == int(d_rest_lengths.size()));
#endif
    stream.pack(&num_springs,1);
    stream.pack(&d_master_idx,1);
    stream.pack(&d_slave_idxs[0],num_springs);
    stream.pack(&d_force_fcn_idxs[0],num_springs);
    stream.pack(&d_stiffnesses[0],num_springs);
    stream.pack(&d_rest_lengths[0],num_springs);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
