// Filename: DiagonalScalingPreconditioner.I
// Created on 18 Apr 2009 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline void
DiagonalScalingPreconditioner::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::setInitialGuessNonzero() not supported" << std::endl);
    return;
}// setInitialGuessNonzero

inline bool
DiagonalScalingPreconditioner::getInitialGuessNonzero() const
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::getInitialGuessNonzero() not supported" << std::endl);
    return true;
}// getInitialGuessNonzero

inline void
DiagonalScalingPreconditioner::setMaxIterations(
    int max_iterations)
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::setMaxIterations() not supported" << std::endl);
    return;
}// setMaxIterations

inline int
DiagonalScalingPreconditioner::getMaxIterations() const
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::getMaxIterations() not supported" << std::endl);
    return 0;
}// getMaxIterations

inline void
DiagonalScalingPreconditioner::setAbsoluteTolerance(
    double abs_residual_tol)
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::setAbsoluteTolerance() not supported" << std::endl);
    return;
}//setAbsoluteTolerance

inline double
DiagonalScalingPreconditioner::getAbsoluteTolerance() const
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::getAbsoluteTolerance() not supported" << std::endl);
    return 0.0;
}// getAbsoluteTolerance

inline void
DiagonalScalingPreconditioner::setRelativeTolerance(
    double rel_residual_tol)
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::setRelativeTolerance() not supported" << std::endl);
    return;
}//setRelativeTolerance

inline double
DiagonalScalingPreconditioner::getRelativeTolerance() const
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::fetRelativeTolerance() not supported" << std::endl);
    return 0.0;
}// getRelativeTolerance

inline int
DiagonalScalingPreconditioner::getNumIterations() const
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::getNumIterations() not supported" << std::endl);
    return 0;
}// getNumIterations

inline double
DiagonalScalingPreconditioner::getResidualNorm() const
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::getResidualNorm() not supported" << std::endl);
    return 0.0;
}// getResidualNorm

inline void
DiagonalScalingPreconditioner::enableLogging(
    bool enabled)
{
    //TBOX_WARNING("DiagonalScalingPreconditioner::enableLogging() not supported" << std::endl);
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
