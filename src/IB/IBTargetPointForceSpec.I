// Filename: IBTargetPointForceSpec.C
// Last modified: <15.Dec.2009 19:23:25 griffith@boyce-griffiths-mac-pro.local>
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

inline void*
IBTargetPointForceSpec::operator new(
    std::size_t size)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(size == sizeof(IBTargetPointForceSpec));
#endif
    return s_arena.alloc(size);
}// operator new

inline void
IBTargetPointForceSpec::operator delete(
    void* ptr,
    std::size_t size)
{
    if (ptr == NULL) return;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(size == sizeof(IBTargetPointForceSpec));
#endif
    s_arena.free(ptr);
    return;
}// operator delete

inline
IBTargetPointForceSpec::IBTargetPointForceSpec(
    const int master_idx,
    const double& kappa_target,
    const double& eta_target,
    const std::vector<double>& X_target)
    : d_master_idx(master_idx),
      d_kappa_target(kappa_target),
      d_eta_target(eta_target),
      d_X_target(X_target)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(NDIM == int(d_X_target.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBTargetPointForceSpec::IBTargetPointForceSpec():\n"
                   << "  must call IBTargetPointForceSpec::registerWithStashableManager() before\n"
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

inline const double&
IBTargetPointForceSpec::getDamping() const
{
    return d_eta_target;
}// getDamping

inline double&
IBTargetPointForceSpec::getDamping()
{
    return d_eta_target;
}// getDamping

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
            (2+NDIM)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

inline void
IBTargetPointForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(NDIM == int(d_X_target.size()));
#endif
    stream.pack(&d_master_idx,1);
    stream.pack(&d_kappa_target,1);
    stream.pack(&d_eta_target,1);
    stream.pack(&d_X_target[0],NDIM);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
