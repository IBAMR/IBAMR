//
// LNodePosnInitStrategy.C
//
// Created on 11 Jul 2004
//         by Boyce Griffith (boyce@trasnaform.speakeasy.net).
//
// Last modified: <02.May.2005 00:28:40 boyce@mstu1.cims.nyu.edu>
//

#include "LNodePosnInitStrategy.h"

/////////////////////////////// INLINE ///////////////////////////////////////

//#ifdef DEBUG_NO_INLINE
//#include "LNodePosnInitStrategy.I"
//#endif

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

void LNodePosnInitStrategy::tagCellsForInitialRefinement(
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    TBOX_WARNING("LNodePosnInitStrategy::tagCellsForInitialRefinement()\n"
                 << "  default implementation employed, no cells tagged for refinement.\n");
    return;
}// tagCellsForInitialRefinement
    
/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION

#include "tbox/Pointer.C"

//////////////////////////////////////////////////////////////////////
///
/// These declarations are required to use the LNodePosnInitStrategy
/// class.
///
//////////////////////////////////////////////////////////////////////

template class tbox::Pointer<LNodePosnInitStrategy>;

#endif

//////////////////////////////////////////////////////////////////////////////
