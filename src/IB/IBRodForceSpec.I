// Filename: IBRodForceSpec.I
// Last modified: <23.Jun.2010 16:42:25 griffith@boyce-griffiths-mac-pro.local>
// Created on 23 Jun 2010 by Boyce Griffith (griffith@boyce-griffiths-mac-pro.local)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBRodForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IBRodForceSpec::IBRodForceSpec(
    const int master_idx,
    const std::vector<int>& next_idxs,
    const std::vector<std::vector<double> >& material_params)
    : d_master_idx(master_idx),
      d_next_idxs(next_idxs),
      d_material_params(material_params)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    const int num_rods = d_next_idxs.size();
    TBOX_ASSERT(num_rods == int(d_material_params.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBRodForceSpec::IBRodForceSpec():\n"
                   << "  must call IBRodForceSpec::registerWithStashableManager() before\n"
                   << "  creating any IBRodForceSpec objects.\n");
    }
    return;
}// IBRodForceSpec

inline
IBRodForceSpec::~IBRodForceSpec()
{
    // intentionally blank
    return;
}// ~IBRodForceSpec

inline unsigned
IBRodForceSpec::getNumberOfRods() const
{
    const int num_rods = d_next_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_rods == int(d_material_params.size()));
#endif
    return num_rods;
}// getNumberOfRods

inline const int&
IBRodForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

inline int&
IBRodForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

inline const std::vector<int>&
IBRodForceSpec::getNextNodeIndices() const
{
    return d_next_idxs;
}// getNextNodeIndices

inline std::vector<int>&
IBRodForceSpec::getNextNodeIndices()
{
    return d_next_idxs;
}// getNextNodeIndices

inline const std::vector<std::vector<double> >&
IBRodForceSpec::getMaterialParams() const
{
    return d_material_params;
}// getMaterialParams

inline std::vector<std::vector<double> >&
IBRodForceSpec::getMaterialParams()
{
    return d_material_params;
}// getMaterialParams

inline int
IBRodForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBRodForceSpec::getDataStreamSize() const
{
    const int num_rods = d_next_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_rods == int(d_material_params.size()));
#endif
    return (( 2+num_rods)*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (10*num_rods)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

inline void
IBRodForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    const int num_rods = d_next_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_rods == int(d_material_params.size()));
#endif
    stream.pack(&num_rods,1);
    stream.pack(&d_master_idx,1);
    stream.pack(&d_next_idxs[0],num_rods);
    for (int n = 0; n < num_rods; ++n)
    {
        stream.pack(&d_material_params[n][0],10);
    }
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
