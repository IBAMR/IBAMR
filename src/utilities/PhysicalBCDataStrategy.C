// Filename: PhysicalBCDataStrategy.C
// Last modified: <03.Sep.2006 19:53:13 boyce@bigboy.nyconnect.com>
// Created on 15 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "PhysicalBCDataStrategy.h"

// IBAMR INCLUDES
#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#endif

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PhysicalBCDataStrategy::PhysicalBCDataStrategy(
    const std::string& object_name)
    : d_object_name(object_name)
{
    // intentionally blank
    return;
}// PhysicalBCDataStrategy

PhysicalBCDataStrategy::~PhysicalBCDataStrategy()
{
    // intentionally blank
    return;
}// ~PhysicalBCDataStrategy

void
PhysicalBCDataStrategy::setPhysicalBoundaryConditionsOnPatchHierarchy(
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const double fill_time,
    const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill,
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
        setPhysicalBoundaryConditionsOnPatchLevel(
            data_idx, var, hierarchy->getPatchLevel(level_num),
            fill_time, ghost_width_to_fill);
    }
    return;
}// setPhysicalBoundaryConditionsOnPatchHierarchy

void
PhysicalBCDataStrategy::setPhysicalBoundaryConditionsOnPatchLevel(
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
    const double fill_time,
    const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!level.isNull());
#endif
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        setPhysicalBoundaryConditionsOnPatch(
            data_idx, var, *(level->getPatch(p())),
            fill_time, ghost_width_to_fill);
    }
    return;
}// setPhysicalBoundaryConditionsOnPatchLevel

void
PhysicalBCDataStrategy::printClassData(
    std::ostream& os) const
{
    os << "++++++++++++++++++++++++++++++++++++++++++++++++\n";
    os << "\nPhysicalBCDataStrategy::printClassData...\n";
    os << "PhysicalBCDataStrategy: this = " << const_cast<PhysicalBCDataStrategy*>(this) << "\n";
    os << "d_object_name = " << d_object_name << "\n";
    os << "++++++++++++++++++++++++++++++++++++++++++++++++\n";
    return;
}// printClassData

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::PhysicalBCDataStrategy>;

//////////////////////////////////////////////////////////////////////////////
