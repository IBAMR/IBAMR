// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
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

#include "BasePatchHierarchy.h"
#include "PatchHierarchy.h"
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

#include <string>

/*!
 * Pre processing call back function to be hooked into IBAMR::HierarchyIntegrator class.
 *
 * This static member should be registered with an appropriate hierarchy integrator
 * via registerApplyGradientDetectorCallback
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to TagLSRefinementCells class object.
 *
 * \TODO: Let's move this out of the global namespace and use "snake case" for static function names.
 */
void
callTagGasLSRefinementCellsCallbackFunction(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                            const int level_number,
                                            const double error_data_time,
                                            const int tag_index,
                                            const bool initial_time,
                                            const bool uses_richardson_extrapolation_too,
                                            void* ctx);
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
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_var;
    double d_tag_value;
    double d_tag_abs_thresh;
};

#endif // #ifndef included_TagLSRefinementCells
