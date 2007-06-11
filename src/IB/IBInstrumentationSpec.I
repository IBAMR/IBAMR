// Filename: IBInstrumentationSpec.I
// Last modified: <11.Jun.2007 17:37:23 griffith@box221.cims.nyu.edu>
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

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IBInstrumentationSpec::IBInstrumentationSpec(
    const int master_idx,
    const int flow_meter_idx,
    const int pressure_gauge_idx)
    : d_master_idx(master_idx),
      d_flow_meter_idx(flow_meter_idx),
      d_pressure_gauge_idx(pressure_gauge_idx)
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
IBInstrumentationSpec::getFlowMeterIndex() const
{
    return d_flow_meter_idx;
}// getFlowMeterIndex

inline int&
IBInstrumentationSpec::getFlowMeterIndex()
{
    return d_flow_meter_idx;
}// getFlowMeterIndex

inline const int&
IBInstrumentationSpec::getPressureGaugeIndex() const
{
    return d_pressure_gauge_idx;
}// getPressureGaugeIndex

inline int&
IBInstrumentationSpec::getPressureGaugeIndex()
{
    return d_pressure_gauge_idx;
}// getPressureGaugeIndex

inline int
IBInstrumentationSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBInstrumentationSpec::getDataStreamSize() const
{
    return 3*SAMRAI::tbox::AbstractStream::sizeofInt();
}// getDataStreamSize

inline void
IBInstrumentationSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_master_idx,1);
    stream.pack(&d_flow_meter_idx,1);
    stream.pack(&d_pressure_gauge_idx,1);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// NAMESPACE ////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
