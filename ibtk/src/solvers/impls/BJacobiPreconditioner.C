// Filename: BJacobiPreconditioner.C
// Created on 11 Apr 2005 by Boyce Griffith
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

#include "BJacobiPreconditioner.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/namespaces.h>

// C++ STDLIB INCLUDES
#include <ostream>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

BJacobiPreconditioner::BJacobiPreconditioner(
    Pointer<Database> input_db)
    : d_is_initialized(false),
      d_pc_map(),
      d_initial_guess_nonzero(false),
      d_rel_residual_tol(1.0e-6),
      d_abs_residual_tol(1.0e-30),
      d_max_iterations(1)
{
    // Get configuration data from the input database.
    if (!input_db.isNull())
    {
        // LinearSolver options.
        d_initial_guess_nonzero = input_db->getBoolWithDefault(
            "initial_guess_nonzero", d_initial_guess_nonzero);
        d_rel_residual_tol = input_db->getDoubleWithDefault(
            "rel_residual_tol", d_rel_residual_tol);
        d_abs_residual_tol = input_db->getDoubleWithDefault(
            "abs_residual_tol", d_abs_residual_tol);
        d_max_iterations = input_db->getIntegerWithDefault(
            "max_iterations", d_max_iterations);
    }
    return;
}// BJacobiPreconditioner()

BJacobiPreconditioner::~BJacobiPreconditioner()
{
    if (d_is_initialized) deallocateSolverState();
    return;
}// ~BJacobiPreconditioner()

void
BJacobiPreconditioner::setComponentPreconditioner(
    Pointer<LinearSolver> preconditioner,
    const int component)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!preconditioner.isNull());
    TBOX_ASSERT(component >= 0);
#endif
    d_pc_map[component] = preconditioner;
    return;
}// setComponentPreconditioner

bool
BJacobiPreconditioner::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    /*!
     * \todo Setup the solver so that it works correctly for multiple
     * iterations.
     */
    TBOX_ASSERT(d_max_iterations == 1);

    // Initialize the preconditioner, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x,b);

    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int   finest_ln = x.getFinestLevelNumber()  ;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(x.getNumberOfComponents() == b.getNumberOfComponents());
    TBOX_ASSERT(hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(  finest_ln == b.getFinestLevelNumber()  );
#endif
    const std::string& x_name = x.getName();
    const std::string& b_name = b.getName();
    bool ret_val = true;
    for (int comp = 0; comp < x.getNumberOfComponents(); ++comp)
    {
        // Setup a SAMRAIVectorReal to correspond to the individual vector
        // component.
        std::ostringstream str;
        str << comp;

        SAMRAIVectorReal<NDIM,double> x_comp(
            x_name+"_component_"+str.str(),
            hierarchy, coarsest_ln, finest_ln);
        x_comp.addComponent(x.getComponentVariable(comp),
                            x.getComponentDescriptorIndex(comp),
                            x.getControlVolumeIndex(comp));

        SAMRAIVectorReal<NDIM,double> b_comp(
            b_name+"_component_"+str.str(),
            hierarchy, coarsest_ln, finest_ln);
        b_comp.addComponent(b.getComponentVariable(comp),
                            b.getComponentDescriptorIndex(comp),
                            b.getControlVolumeIndex(comp));

        // Configure the component preconditioner.
        Pointer<LinearSolver> pc_comp = d_pc_map[comp];
        pc_comp->setInitialGuessNonzero(d_initial_guess_nonzero);
        pc_comp->setMaxIterations(d_max_iterations);
        pc_comp->setAbsoluteTolerance(d_abs_residual_tol);
        pc_comp->setRelativeTolerance(d_rel_residual_tol);

        // Apply the component preconditioner.
        const bool ret_val_comp = pc_comp->solveSystem(x_comp, b_comp);
        ret_val = ret_val && ret_val_comp;
    }

    // Deallocate the preconditioner, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    return ret_val;
}// solveSystem

void
BJacobiPreconditioner::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& b)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int finest_ln = x.getFinestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(finest_ln == b.getFinestLevelNumber());
    TBOX_ASSERT(x.getNumberOfComponents() == b.getNumberOfComponents());
#endif
    // Initialize the component preconditioners.
    const std::string& x_name = x.getName();
    const std::string& b_name = b.getName();
    for (std::map<int,Pointer<LinearSolver> >::iterator it = d_pc_map.begin();
         it != d_pc_map.end(); ++it)
    {
        const int comp = (*it).first;
        SAMRAIVectorReal<NDIM,double> x_comp(
            x_name+"_component", hierarchy, coarsest_ln, finest_ln);
        SAMRAIVectorReal<NDIM,double> b_comp(
            b_name+"_component", hierarchy, coarsest_ln, finest_ln);
        x_comp.addComponent(x.getComponentVariable(comp),
                            x.getComponentDescriptorIndex(comp),
                            x.getControlVolumeIndex(comp));
        b_comp.addComponent(b.getComponentVariable(comp),
                            b.getComponentDescriptorIndex(comp),
                            b.getControlVolumeIndex(comp));
        d_pc_map[comp]->initializeSolverState(x_comp, b_comp);
    }

    // Indicate that the preconditioner is initialized.
    d_is_initialized = true;

    return;
}// initializeSolverState

void
BJacobiPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    // Deallocate the component preconditioners.
    for (std::map<int,Pointer<LinearSolver> >::iterator it = d_pc_map.begin();
         it != d_pc_map.end(); ++it)
    {
        const int comp = (*it).first;
        d_pc_map[comp]->deallocateSolverState();
    }

    // Indicate that the preconditioner is NOT initialized.
    d_is_initialized = false;

    return;
}// deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::BJacobiPreconditioner>;

//////////////////////////////////////////////////////////////////////////////
