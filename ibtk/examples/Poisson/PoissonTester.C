// Filename: PoissonTester.C
// Created on 12 Feb 2005 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "PoissonTester.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <Index.h>
#include <Patch.h>
#include <VariableDatabase.h>

/////////////////////////////// STATIC ///////////////////////////////////////

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

/////////////////////////////// PUBLIC ///////////////////////////////////////

PoissonTester::PoissonTester()
    : d_context(hier::VariableDatabase<NDIM>::getDatabase()->getContext("Poisson_TESTER::CONTEXT")),
      d_U_idx(hier::VariableDatabase<NDIM>::getDatabase()->
              registerVariableAndContext(new pdat::CellVariable<NDIM,double>("U"), d_context, CELLG)),
      d_V_idx(hier::VariableDatabase<NDIM>::getDatabase()->
              registerVariableAndContext(new pdat::CellVariable<NDIM,double>("V"), d_context, CELLG)),
      d_F_idx(hier::VariableDatabase<NDIM>::getDatabase()->
              registerVariableAndContext(new pdat::CellVariable<NDIM,double>("F"), d_context, CELLG))
{
    // intentionally blank
    return;
}// PoissonTester

PoissonTester::~PoissonTester()
{
    // intentionally blank
    return;
}// ~PoissonTester

///
///  The following routines:
///
///      initializeLevelData(),
///      resetHierarchyConfiguration()
///
///  are concrete implementations of functions declared in the
///  mesh::StandardTagAndInitStrategy abstract base class.
///

void
PoissonTester::initializeLevelData(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    // Allocate storage needed to initialize the level and fill data from
    // coarser levels in AMR hierarchy, if any.
    //
    // Since time gets set when we allocate data, re-stamp it to current time if
    // we don't need to allocate.
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    if (allocate_data)
    {
        level->allocatePatchData(d_U_idx, init_data_time);
        level->allocatePatchData(d_V_idx, init_data_time);
        level->allocatePatchData(d_F_idx, init_data_time);
    }
    else
    {
        level->setTime(init_data_time, d_U_idx);
        level->setTime(init_data_time, d_V_idx);
        level->setTime(init_data_time, d_F_idx);
    }

    // Initialize level data.
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
        const hier::Box<NDIM>& patch_box = patch->getBox();
        const hier::Index<NDIM>& patch_lower = patch_box.lower();
        tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const xLower = pgeom->getXLower();
        const double* const dx = pgeom->getDx();

        tbox::Pointer< pdat::CellData<NDIM,double> > U_data = patch->getPatchData(d_U_idx);
        tbox::Pointer< pdat::CellData<NDIM,double> > V_data = patch->getPatchData(d_V_idx);
        tbox::Pointer< pdat::CellData<NDIM,double> > F_data = patch->getPatchData(d_F_idx);

        U_data->fillAll(0.0);
        for (hier::Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const hier::Index<NDIM>& i = b();
            (*F_data)(i) = 4.0*M_PI*M_PI*double(NDIM);
            (*V_data)(i) = 1.0;

            (*F_data)(i) *= sin(2.0*M_PI*(xLower[0] + (i[0] - patch_lower[0] + 0.5)*dx[0]));
            (*V_data)(i) *= sin(2.0*M_PI*(xLower[0] + (i[0] - patch_lower[0] + 0.5)*dx[0]));

            for (int d = 1; d < NDIM; ++d)
            {
                (*F_data)(i) *= cos(2.0*M_PI*(xLower[d] + (i[d] - patch_lower[d] + 0.5)*dx[d]));
                (*V_data)(i) *= cos(2.0*M_PI*(xLower[d] + (i[d] - patch_lower[d] + 0.5)*dx[d]));
            }
            (*V_data)(i) += (i(0)+0.5)*dx[0];
        }
    }
    return;
}// initializeLevelData

void
PoissonTester::resetHierarchyConfiguration(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level <= finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    return;
}// resetHierarchyConfiguration

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class tbox::Pointer<PoissonTester>;

//////////////////////////////////////////////////////////////////////////////
