// Filename: IBDataPostProcessor.C
// Last modified: <24.Sep.2008 17:20:07 griffith@box230.cims.nyu.edu>
// Created on 24 Sep 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "IBDataPostProcessor.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBDataPostProcessor::IBDataPostProcessor()
{
    // intentionally blank
    return;
}// IBDataPostProcessor

IBDataPostProcessor::~IBDataPostProcessor()
{
    // intentionally blank
    return;
}// ~IBDataPostProcessor

void
IBDataPostProcessor::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    IBTK::LDataManager* const lag_manager)
{
    // intentionally blank
    return;
}// initializeLevelData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBDataPostProcessor>;

//////////////////////////////////////////////////////////////////////////////
