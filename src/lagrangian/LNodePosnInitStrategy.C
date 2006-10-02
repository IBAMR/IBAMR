// Filename: LNodePosnInitStrategy.C
// Created on 11 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <02.Oct.2006 14:21:38 boyce@boyce-griffiths-powerbook-g4-15.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "LNodePosnInitStrategy.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodePosnInitStrategy::LNodePosnInitStrategy()
{
    // intentionally blank
    return;
}// LNodePosnInitStrategy

LNodePosnInitStrategy::~LNodePosnInitStrategy()
{
    // intentionally blank
    return;
}// ~LNodePosnInitStrategy

void
LNodePosnInitStrategy::tagCellsForInitialRefinement(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    TBOX_WARNING("LNodePosnInitStrategy::tagCellsForInitialRefinement()\n"
                 << "  default implementation employed, no cells tagged for refinement.\n");
    return;
}// tagCellsForInitialRefinement

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LNodePosnInitStrategy>;

//////////////////////////////////////////////////////////////////////////////
