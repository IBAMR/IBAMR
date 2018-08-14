// Filename: TagLSRefinementCells.h
// Created on Dec 28, 2017 by Nishant Nangia

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_TagLSRefinementCells
#define included_TagLSRefinementCells

///////////////////////////// INCLUDES ///////////////////////////////////
#include <string>

#include "BasePatchHierarchy.h"
#include "PatchHierarchy.h"
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

/*!
 * Pre processing call back function to be hooked into IBAMR::HierarchyIntegrator class.
 *
 * This static member should be registered with an appropriate hierarchy integrator
 * via registerApplyGradientDetectorCallback
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to TagLSRefinementCells class object.
 */

void callTagLSRefinementCellsCallbackFunction(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                              const int level_number,
                                              const double error_data_time,
                                              const int tag_index,
                                              const bool initial_time,
                                              const bool uses_richardson_extrapolation_too,
                                              void* ctx);

struct TagLSRefinementCells
{
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_gas_var;
    double d_tag_value;
    double d_tag_abs_thresh;
};

#endif // #ifndef included_TagLSRefinementCells
