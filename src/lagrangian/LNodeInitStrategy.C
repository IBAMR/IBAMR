// Filename: LNodeInitStrategy.C
// Created on 11 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <23.Jan.2007 02:55:54 boyce@bigboy.nyconnect.com>

#include "LNodeInitStrategy.h"

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

LNodeInitStrategy::LNodeInitStrategy()
{
    // intentionally blank
    return;
}// LNodeInitStrategy

LNodeInitStrategy::~LNodeInitStrategy()
{
    // intentionally blank
    return;
}// ~LNodeInitStrategy

int
LNodeInitStrategy::initializeMassDataOnPatchLevel(
    const int global_index_offset,
    const int local_index_offset,
    SAMRAI::tbox::Pointer<LNodeLevelData>& M_data,
    SAMRAI::tbox::Pointer<LNodeLevelData>& K_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    TBOX_WARNING("LNodeInitStrategy::initializeMassDataOnPatchLevel()\n"
                 << "  default implementation employed, no mass data initialized.\n");
    return 0;
}// initializeMassDataOnPatchLevel

void
LNodeInitStrategy::tagCellsForInitialRefinement(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    TBOX_WARNING("LNodeInitStrategy::tagCellsForInitialRefinement()\n"
                 << "  default implementation employed, no cells tagged for refinement.\n");
    return;
}// tagCellsForInitialRefinement

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LNodeInitStrategy>;

//////////////////////////////////////////////////////////////////////////////
