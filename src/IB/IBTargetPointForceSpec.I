// Filename: IBTargetPointForceSpec.C
// Last modified: <12.Mar.2008 22:44:33 griffith@box221.cims.nyu.edu>
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
IBTargetPointForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IBTargetPointForceSpec::IBTargetPointForceSpec(
    const int master_idx,
    const double& kappa_target,
    const std::vector<double>& X_target)
    : d_master_idx(master_idx),
      d_kappa_target(kappa_target),
      d_X_target(X_target)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(NDIM == static_cast<int>(d_X_target.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBTargetPointForceSpec::IBTargetPointForceSpec():\n"
                   << "  must call IBTargetPointForceSpec::registerWithIBTK::StashableManager() before\n"
                   << "  creating any IBTargetPointForceSpec objects.\n");
    }
    return;
}// IBTargetPointForceSpec

inline
IBTargetPointForceSpec::~IBTargetPointForceSpec()
{
    // intentionally blank
    return;
}// ~IBTargetPointForceSpec

inline const int&
IBTargetPointForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

inline int&
IBTargetPointForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

inline const double&
IBTargetPointForceSpec::getStiffness() const
{
    return d_kappa_target;
}// getStiffness

inline double&
IBTargetPointForceSpec::getStiffness()
{
    return d_kappa_target;
}// getStiffness

inline const std::vector<double>&
IBTargetPointForceSpec::getTargetPointPosition() const
{
    return d_X_target;
}// getTargetPointPosition

inline std::vector<double>&
IBTargetPointForceSpec::getTargetPointPosition()
{
    return d_X_target;
}// getTargetPointPosition

inline int
IBTargetPointForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBTargetPointForceSpec::getDataStreamSize() const
{
    return ((1     )*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (1+NDIM)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

inline void
IBTargetPointForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(NDIM == static_cast<int>(d_X_target.size()));
#endif
    stream.pack(&d_master_idx,1);
    stream.pack(&d_kappa_target,1);
    stream.pack(&d_X_target[0],NDIM);
    return;
}// packStream

inline void
IBTargetPointForceSpec::registerPeriodicShift(
    const SAMRAI::hier::IntVector<NDIM>& offset,
    const std::vector<double>& displacement)
{
    for (int d = 0; d < NDIM; ++d)
    {
        d_X_target[d] += displacement[d];
    }
    return;
}// registerPeriodicShift

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
