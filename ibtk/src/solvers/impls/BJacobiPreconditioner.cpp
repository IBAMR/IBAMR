// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/BJacobiPreconditioner.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <map>
#include <ostream>
#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

BJacobiPreconditioner::BJacobiPreconditioner(std::string object_name,
                                             Pointer<Database> input_db,
                                             const std::string& /*default_options_prefix*/)
{
    // Setup default options.
    GeneralSolver::init(std::move(object_name), /*homogeneous_bc*/ true);
    d_initial_guess_nonzero = false;
    d_max_iterations = 1;

    // Get configuration data from the input database.
    if (input_db)
    {
        // LinearSolver options.
        if (input_db->keyExists("initial_guess_nonzero"))
            setInitialGuessNonzero(input_db->getBool("initial_guess_nonzero"));
        if (input_db->keyExists("rel_residual_tol")) setRelativeTolerance(input_db->getDouble("rel_residual_tol"));
        if (input_db->keyExists("abs_residual_tol")) setAbsoluteTolerance(input_db->getDouble("abs_residual_tol"));
        if (input_db->keyExists("max_iterations")) setMaxIterations(input_db->getInteger("max_iterations"));
    }
    return;
} // BJacobiPreconditioner()

BJacobiPreconditioner::~BJacobiPreconditioner()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~BJacobiPreconditioner()

void
BJacobiPreconditioner::setComponentPreconditioner(Pointer<LinearSolver> preconditioner, const unsigned int component)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(preconditioner);
#endif
    d_pc_map[component] = preconditioner;
    return;
} // setComponentPreconditioner

bool
BJacobiPreconditioner::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    // Initialize the preconditioner, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int finest_ln = x.getFinestLevelNumber();
#if !defined(NDEBUG)
    TBOX_ASSERT(x.getNumberOfComponents() == b.getNumberOfComponents());
    TBOX_ASSERT(hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(finest_ln == b.getFinestLevelNumber());
#endif
    const std::string& x_name = x.getName();
    const std::string& b_name = b.getName();
    bool ret_val = true;

// Zero out the initial guess.
#if !defined(NDEBUG)
    TBOX_ASSERT(d_initial_guess_nonzero == false);
#endif
    x.setToScalar(0.0, /*interior_only*/ false);

    for (int comp = 0; comp < x.getNumberOfComponents(); ++comp)
    {
        // Setup a SAMRAIVectorReal to correspond to the individual vector
        // component.
        SAMRAIVectorReal<NDIM, double> x_comp(
            x_name + "_component_" + std::to_string(comp), hierarchy, coarsest_ln, finest_ln);
        x_comp.addComponent(
            x.getComponentVariable(comp), x.getComponentDescriptorIndex(comp), x.getControlVolumeIndex(comp));

        SAMRAIVectorReal<NDIM, double> b_comp(
            b_name + "_component_" + std::to_string(comp), hierarchy, coarsest_ln, finest_ln);
        b_comp.addComponent(
            b.getComponentVariable(comp), b.getComponentDescriptorIndex(comp), b.getControlVolumeIndex(comp));

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
} // solveSystem

void
BJacobiPreconditioner::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                             const SAMRAIVectorReal<NDIM, double>& b)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int finest_ln = x.getFinestLevelNumber();
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(finest_ln == b.getFinestLevelNumber());
    TBOX_ASSERT(x.getNumberOfComponents() == b.getNumberOfComponents());
#endif
    // Initialize the component preconditioners.
    const std::string& x_name = x.getName();
    const std::string& b_name = b.getName();
    for (const auto& linearSolver_pair : d_pc_map)
    {
        const int comp = linearSolver_pair.first;
        SAMRAIVectorReal<NDIM, double> x_comp(x_name + "_component", hierarchy, coarsest_ln, finest_ln);
        x_comp.addComponent(
            x.getComponentVariable(comp), x.getComponentDescriptorIndex(comp), x.getControlVolumeIndex(comp));
        SAMRAIVectorReal<NDIM, double> b_comp(b_name + "_component", hierarchy, coarsest_ln, finest_ln);
        b_comp.addComponent(
            b.getComponentVariable(comp), b.getComponentDescriptorIndex(comp), b.getControlVolumeIndex(comp));
        d_pc_map[comp]->initializeSolverState(x_comp, b_comp);
    }

    // Indicate that the preconditioner is initialized.
    d_is_initialized = true;
    return;
} // initializeSolverState

void
BJacobiPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    // Deallocate the component preconditioners.
    for (const auto& linearSolver_pair : d_pc_map)
    {
        const int comp = linearSolver_pair.first;
        d_pc_map[comp]->deallocateSolverState();
    }

    // Indicate that the preconditioner is NOT initialized.
    d_is_initialized = false;
    return;
} // deallocateSolverState

void
BJacobiPreconditioner::setInitialGuessNonzero(bool initial_guess_nonzero)
{
    if (initial_guess_nonzero)
    {
        TBOX_ERROR("BJacobiPreconditioner::setInitialGuessNonzero()\n"
                   << "  class IBTK::BJacobiPreconditioner requires a zero initial guess" << std::endl);
    }
    d_initial_guess_nonzero = initial_guess_nonzero;
    return;
} // setInitialGuessNonzero

void
BJacobiPreconditioner::setMaxIterations(int max_iterations)
{
    if (max_iterations > 1)
    {
        TBOX_ERROR("BJacobiPreconditioner::setMaxIterations()\n"
                   << "  class IBTK::BJacobiPreconditioner requires max_iterations == 1" << std::endl);
    }
    d_max_iterations = max_iterations;
    return;
} // setMaxIterations

int
BJacobiPreconditioner::getNumIterations() const
{
    IBTK_DO_ONCE(TBOX_WARNING("BJacobiPreconditioner::getNumIterations() not supported" << std::endl););
    return 0;
} // getNumIterations

double
BJacobiPreconditioner::getResidualNorm() const
{
    IBTK_DO_ONCE(TBOX_WARNING("BJacobiPreconditioner::getResidualNorm() not supported" << std::endl););
    return 0.0;
} // getResidualNorm

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
