// Filename: IBSpringForceSpec.I
// Last modified: <12.Feb.2008 21:17:37 griffith@box221.cims.nyu.edu>
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

inline
IBSpringForceSpec::IBSpringForceSpec(
    const int master_idx,
    const std::vector<int>& slave_idxs,
    const std::vector<int>& force_fcn_idxs,
    const std::vector<double>& stiffnesses,
    const std::vector<double>& rest_lengths)
    : d_num_springs(static_cast<int>(slave_idxs.size())),
      d_master_idx(master_idx),
      d_slave_idxs(slave_idxs),
      d_force_fcn_idxs(force_fcn_idxs),
      d_stiffnesses(stiffnesses),
      d_rest_lengths(rest_lengths)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_slave_idxs.size()));
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_force_fcn_idxs.size()));
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_stiffnesses.size()));
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_rest_lengths.size()));
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
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_slave_idxs.size()));
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_force_fcn_idxs.size()));
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_stiffnesses.size()));
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_rest_lengths.size()));
#endif
    return d_num_springs;
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
    return ((2+2*d_num_springs)*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (  2*d_num_springs)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

inline void
IBSpringForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_slave_idxs.size()));
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_force_fcn_idxs.size()));
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_stiffnesses.size()));
    TBOX_ASSERT(d_num_springs == static_cast<int>(d_rest_lengths.size()));
#endif
    stream.pack(&d_num_springs,1);
    stream.pack(&d_master_idx,1);
    stream.pack(&d_slave_idxs[0],d_num_springs);
    stream.pack(&d_force_fcn_idxs[0],d_num_springs);
    stream.pack(&d_stiffnesses[0],d_num_springs);
    stream.pack(&d_rest_lengths[0],d_num_springs);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
