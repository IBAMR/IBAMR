// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#include "ibtk/BGaussSeidelPreconditioner.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/SAMRAIScopedVectorCopy.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "tbox/ConstPointer.h"
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

BGaussSeidelPreconditioner::BGaussSeidelPreconditioner(std::string object_name,
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
        // Block Gauss-Seidel options.
        if (input_db->keyExists("symmetric_preconditioner"))
            d_symmetric_preconditioner = input_db->getBool("symmetric_preconditioner");
        if (input_db->keyExists("reverse_order")) d_reverse_order = input_db->getBool("reverse_order");

        // LinearSolver options.
        if (input_db->keyExists("initial_guess_nonzero"))
            setInitialGuessNonzero(input_db->getBool("initial_guess_nonzero"));
        if (input_db->keyExists("rel_residual_tol")) setRelativeTolerance(input_db->getDouble("rel_residual_tol"));
        if (input_db->keyExists("abs_residual_tol")) setAbsoluteTolerance(input_db->getDouble("abs_residual_tol"));
        if (input_db->keyExists("max_iterations")) setMaxIterations(input_db->getInteger("max_iterations"));
    }
    return;
} // BGaussSeidelPreconditioner()

BGaussSeidelPreconditioner::~BGaussSeidelPreconditioner()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~BGaussSeidelPreconditioner()

void
BGaussSeidelPreconditioner::setComponentPreconditioner(Pointer<LinearSolver> preconditioner,
                                                       const unsigned int component)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(preconditioner);
#endif
    d_pc_map[component] = preconditioner;
    return;
} // setComponentPreconditioner

void
BGaussSeidelPreconditioner::setComponentOperators(const std::vector<Pointer<LinearOperator> >& linear_ops,
                                                  const unsigned int component)
{
#if !defined(NDEBUG)
    for (unsigned int k = 0; k < linear_ops.size(); ++k)
    {
        if (k != component) TBOX_ASSERT(linear_ops[k]);
    }
#endif
    d_linear_ops_map[component] = linear_ops;
    return;
} // setComponentOperators

void
BGaussSeidelPreconditioner::setSymmetricPreconditioner(const bool symmetric_preconditioner)
{
    d_symmetric_preconditioner = symmetric_preconditioner;
    return;
} // setSymmetricPreconditioner

void
BGaussSeidelPreconditioner::setReversedOrder(const bool reverse_order)
{
    d_reverse_order = reverse_order;
    return;
} // setReversedOrder

bool
BGaussSeidelPreconditioner::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    // Initialize the preconditioner, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

#if !defined(NDEBUG)
    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int finest_ln = x.getFinestLevelNumber();
    TBOX_ASSERT(x.getNumberOfComponents() == b.getNumberOfComponents());
    TBOX_ASSERT(hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(finest_ln == b.getFinestLevelNumber());
#endif
    bool ret_val = true;

// Zero out the initial guess.
#if !defined(NDEBUG)
    TBOX_ASSERT(d_initial_guess_nonzero == false);
#endif
    x.setToScalar(0.0, /*interior_only*/ false);

    // Setup SAMRAIVectorReal objects to correspond to the individual vector
    // components.
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > > x_comps =
        getComponentVectors(Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > > b_comps =
        getComponentVectors(Pointer<SAMRAIVectorReal<NDIM, double> >(&b, false));

    // Setup the order in which the component preconditioner are to be applied.
    const int ncomps = x.getNumberOfComponents();
    std::vector<int> comps;
    comps.reserve(2 * ncomps - 1);
    if (!d_reverse_order)
    {
        // Standard order: Run from comp = 0 to comp = ncomp-1.
        for (int comp = 0; comp < ncomps; ++comp)
        {
            comps.push_back(comp);
        }
        if (d_symmetric_preconditioner)
        {
            for (int comp = ncomps - 2; comp >= 0; --comp)
            {
                comps.push_back(comp);
            }
        }
    }
    else
    {
        // Reversed order: Run from comp = ncomp-1 to comp = 0.
        for (int comp = ncomps - 1; comp >= 0; --comp)
        {
            comps.push_back(comp);
        }
        if (d_symmetric_preconditioner)
        {
            for (int comp = 1; comp < ncomps; ++comp)
            {
                comps.push_back(comp);
            }
        }
    }

    // Clone the right-hand-side vector to avoid modifying it during the
    // preconditioning operation.
    SAMRAIScopedVectorCopy<double> f(b);
    auto f_comps = f.getComponentVectors();

    // Apply the component preconditioners.
    int count = 0;
    for (auto it = comps.begin(); it != comps.end(); ++it, ++count)
    {
        const int comp = (*it);

        Pointer<SAMRAIVectorReal<NDIM, double> > x_comp = x_comps[comp];
        Pointer<SAMRAIVectorReal<NDIM, double> > b_comp = b_comps[comp];
        Pointer<SAMRAIVectorReal<NDIM, double> > f_comp = f_comps[comp];

        // Update the right-hand-side vector.
        f_comp->setToScalar(0.0);
        for (int c = 0; c < ncomps; ++c)
        {
            if (c == comp) continue;
            d_linear_ops_map[comp][c]->applyAdd(*x_comps[c], *f_comp, *f_comp);
        }
        f_comp->subtract(b_comp, f_comp);

        // Configure the component preconditioner.
        Pointer<LinearSolver> pc_comp = d_pc_map[comp];
        pc_comp->setInitialGuessNonzero(count >= ncomps);
        pc_comp->setMaxIterations(d_max_iterations);
        pc_comp->setAbsoluteTolerance(d_abs_residual_tol);
        pc_comp->setRelativeTolerance(d_rel_residual_tol);

        // Apply the component preconditioner.
        const bool ret_val_comp = pc_comp->solveSystem(*x_comp, *f_comp);
        ret_val = ret_val && ret_val_comp;
    }

    // Deallocate the preconditioner, when necessary.
    if (deallocate_after_solve) deallocateSolverState();
    return ret_val;
} // solveSystem

void
BGaussSeidelPreconditioner::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                                  const SAMRAIVectorReal<NDIM, double>& b)
{
#if !defined(NDEBUG)
    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int finest_ln = x.getFinestLevelNumber();
    TBOX_ASSERT(hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(finest_ln == b.getFinestLevelNumber());
    TBOX_ASSERT(x.getNumberOfComponents() == b.getNumberOfComponents());
#endif
    // Setup SAMRAIVectorReal objects to correspond to the individual vector
    // components.
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > > x_comps =
        getComponentVectors(ConstPointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > > b_comps =
        getComponentVectors(ConstPointer<SAMRAIVectorReal<NDIM, double> >(&b, false));

    // Initialize the component operators and preconditioners.
    const int ncomps = x.getNumberOfComponents();
    for (int comp = 0; comp < ncomps; ++comp)
    {
        for (int c = 0; c < ncomps; ++c)
        {
            // Skip the diagonal operators.
            if (c == comp) continue;
            d_linear_ops_map[comp][c]->initializeOperatorState(*x_comps[comp], *b_comps[comp]);
        }
        d_pc_map[comp]->initializeSolverState(*x_comps[comp], *b_comps[comp]);
    }

    // Indicate that the preconditioner is initialized.
    d_is_initialized = true;
    return;
} // initializeSolverState

void
BGaussSeidelPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    // Deallocate the component preconditioners.
    for (const auto& linearSolver_pair : d_pc_map)
    {
        linearSolver_pair.second->deallocateSolverState();
    }

    // Deallocate the component operators.
    for (auto& linearSolverVec_pair : d_linear_ops_map)
    {
        std::vector<Pointer<LinearOperator> >& comp_linear_ops = linearSolverVec_pair.second;
        for (auto& comp_linear_op : comp_linear_ops)
        {
            if (comp_linear_op) comp_linear_op->deallocateOperatorState();
        }
    }

    // Indicate that the preconditioner is NOT initialized.
    d_is_initialized = false;
    return;
} // deallocateSolverState

void
BGaussSeidelPreconditioner::setInitialGuessNonzero(bool initial_guess_nonzero)
{
    if (initial_guess_nonzero)
    {
        TBOX_ERROR(d_object_name << "::setInitialGuessNonzero()\n"
                                 << "  class IBTK::BGaussSeidelPreconditioner requires a zero initial guess"
                                 << std::endl);
    }
    d_initial_guess_nonzero = initial_guess_nonzero;
    return;
} // setInitialGuessNonzero

void
BGaussSeidelPreconditioner::setMaxIterations(int max_iterations)
{
    if (max_iterations > 1)
    {
        TBOX_ERROR(d_object_name << "::setMaxIterations()\n"
                                 << "  class IBTK::BGaussSeidelPreconditioner requires max_iterations == 1"
                                 << std::endl);
    }
    d_max_iterations = max_iterations;
    return;
} // setMaxIterations

int
BGaussSeidelPreconditioner::getNumIterations() const
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::getNumIterations() not supported" << std::endl););
    return 0;
} // getNumIterations

double
BGaussSeidelPreconditioner::getResidualNorm() const
{
    IBTK_DO_ONCE(TBOX_WARNING(d_object_name << "::getResidualNorm() not supported" << std::endl););
    return 0.0;
} // getResidualNorm

/////////////////////////////// PRIVATE //////////////////////////////////////

std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > >
BGaussSeidelPreconditioner::getComponentVectors(const ConstPointer<SAMRAIVectorReal<NDIM, double> > x)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = x->getPatchHierarchy();
    const int coarsest_ln = x->getCoarsestLevelNumber();
    const int finest_ln = x->getFinestLevelNumber();
    const std::string& x_name = x->getName();
    const int ncomps = x->getNumberOfComponents();

    // Setup SAMRAIVectorReal objects to correspond to the individual vector
    // components.
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > > x_comps(ncomps);
    for (int comp = 0; comp < ncomps; ++comp)
    {
        x_comps[comp] = new SAMRAIVectorReal<NDIM, double>(
            x_name + "_component_" + std::to_string(comp), hierarchy, coarsest_ln, finest_ln);
        x_comps[comp]->addComponent(
            x->getComponentVariable(comp), x->getComponentDescriptorIndex(comp), x->getControlVolumeIndex(comp));
    }
    return x_comps;
} // getComponentVectors

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
