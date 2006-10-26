// Filename: SetDataStrategy.C
// Last modified: <25.Oct.2006 18:26:09 boyce@bigboy.nyconnect.com>
// Created on 15 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "SetDataStrategy.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

SetDataStrategy::SetDataStrategy(
    const std::string& object_name)
    : d_object_name(object_name)
{
    // intentionally blank
    return;
}// SetDataStrategy

SetDataStrategy::~SetDataStrategy()
{
    // intentionally blank
    return;
}// ~SetDataStrategy

void
SetDataStrategy::setDataOnPatchHierarchy(
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const double data_time,
    const bool initial_time,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
#endif
    const int coarsest_ln =
        (coarsest_ln_in == -1
         ? 0
         : coarsest_ln_in);
    const int finest_ln =
        (finest_ln_in == -1
         ? hierarchy->getFinestLevelNumber()
         : finest_ln_in);
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(
            data_idx, var, hierarchy->getPatchLevel(level_num),
            data_time, initial_time);
    }
    return;
}// setDataOnPatchHierarchy

void
SetDataStrategy::setDataOnPatchLevel(
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
    const double data_time,
    const bool initial_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!level.isNull());
#endif
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        setDataOnPatch(
            data_idx, var, *(level->getPatch(p())), data_time, initial_time);
    }
    return;
}// setDataOnPatchLevel

void
SetDataStrategy::printClassData(
    std::ostream& os) const
{
    os << "++++++++++++++++++++++++++++++++++++++++++++++++\n";
    os << "\nSetDataStrategy::printClassData...\n";
    os << "SetDataStrategy: this = " << const_cast<SetDataStrategy*>(this) << "\n";
    os << "d_object_name = " << d_object_name << "\n";
    os << "++++++++++++++++++++++++++++++++++++++++++++++++\n";
    return;
}// printClassData

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::SetDataStrategy>;

//////////////////////////////////////////////////////////////////////////////
