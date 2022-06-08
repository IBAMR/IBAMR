// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_TagLSRefinementCells
#define included_TagLSRefinementCells

///////////////////////////// INCLUDES ///////////////////////////////////

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <tbox/Pointer.h>

#include <BasePatchHierarchy.h>

/*!
 * Pre processing call back function to be hooked into IBAMR::HierarchyIntegrator class.
 *
 * This static member should be registered with an appropriate hierarchy integrator
 * via registerApplyGradientDetectorCallback
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to TagLSRefinementCells class object.
 *
 */

void
callTagSolidLSRefinementCellsCallbackFunction(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                              const int level_number,
                                              const double error_data_time,
                                              const int tag_index,
                                              const bool initial_time,
                                              const bool uses_richardson_extrapolation_too,
                                              void* ctx);

struct TagLSRefinementCells
{
    TagLSRefinementCells(SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                         double tag_value,
                         double tag_abs_thresh)
        : d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_tag_value(tag_value), d_tag_abs_thresh(tag_abs_thresh)
    {
        return;
    }

    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_var;
    double d_tag_value;
    double d_tag_abs_thresh;
};

#endif // #ifndef included_TagLSRefinementCells
