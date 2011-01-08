// Filename: CCPoissonHypreSolver.I
// Last modified: <06.Sep.2006 01:02:49 boyce@bigboy.nyconnect.com>
// Created on 20 Jan 2006 by Boyce Griffith (boyce@boyce.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAI INCLUDES
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
void
CCPoissonHypreSolver::setSolDataDepth(
    const int depth)
{
    d_sol_depth = depth;
    return;
}// setSolDataDepth

inline
void
CCPoissonHypreSolver::setRhsDataDepth(
    const int depth)
{
    d_rhs_depth = depth;
    return;
}// setRhsDataDepth

inline
void
CCPoissonHypreSolver::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    d_initial_guess_nonzero = initial_guess_nonzero;
    return;
}// setInitialGuessNonzero

inline
bool
CCPoissonHypreSolver::getInitialGuessNonzero() const
{
    return d_initial_guess_nonzero;
}// getInitialGuessNonzero

inline
void
CCPoissonHypreSolver::setMaxIterations(
    int max_iterations)
{
    d_max_iterations = max_iterations;
    return;
}// setMaxIterations

inline
int
CCPoissonHypreSolver::getMaxIterations() const
{
    return d_max_iterations;
}// getMaxIterations

inline
void
CCPoissonHypreSolver::setAbsoluteTolerance(
    double abs_residual_tol)
{
    TBOX_WARNING("CCPoissonHypreSolver: " <<
                 "setAbsoluteTolerance() not implemented.\n");
    return;
}//setAbsoluteTolerance

inline
double
CCPoissonHypreSolver::getAbsoluteTolerance() const
{
    TBOX_WARNING("CCPoissonHypreSolver: " <<
                 "getAbsoluteTolerance() not implemented.\n");
    return 0.0;
}// getAbsoluteTolerance

inline
void
CCPoissonHypreSolver::setRelativeTolerance(
    double rel_residual_tol)
{
    d_rel_residual_tol = rel_residual_tol;
    return;
}//setRelativeTolerance

inline
double
CCPoissonHypreSolver::getRelativeTolerance() const
{
    return d_rel_residual_tol;
}// getRelativeTolerance

inline
int
CCPoissonHypreSolver::getNumIterations() const
{
    return d_current_its;
}// getNumIterations

inline
double
CCPoissonHypreSolver::getResidualNorm() const
{
    return d_current_residual_norm;
}// getResidualNorm

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
